#include <splinesketch/paper_splinesketch.hpp>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <iostream>
#include <random>

namespace splinesketch {
struct PaperSplineSketchInspector {
  using Sketch = PaperSplineSketch;
  using Node = Sketch::Node;
  static inline std::uint64_t comparisons = 0;

  static void query(const std::vector<Node>& nodes, double x) {
    const auto result = Sketch::spline_rank(nodes, x);
    assert(std::isfinite(result));
    if (x < nodes.front().x) assert(result == 0);
    else if (x >= nodes.back().x) assert(result == nodes.back().prefix);
    else {
      const auto right = std::lower_bound(nodes.begin(), nodes.end(), x,
          [](const Node& node, double value) { return node.x < value; });
      if (right->x == x) assert(result == right->prefix);
      else {
        const auto& left = *std::prev(right);
        const auto fraction = Sketch::fraction(left.x, right->x, x);
        assert(std::isfinite(fraction) && fraction >= 0 && fraction <= 1);
        assert(result >= left.prefix && result <= right->prefix);
      }
    }
    ++comparisons;
  }

  static void interval(double a, double b, long double prefix, long double mass) {
    assert(std::isfinite(a) && std::isfinite(b) && a < b);
    std::vector<Node> nodes{Node{a, prefix}, Node{b, mass}};
    Sketch::rebuild(nodes);
    assert(nodes.back().prefix <= 0x1p80L);
    const double midpoint = Sketch::midpoint(a, b);
    const std::vector<double> queries{a, b, midpoint, std::nextafter(a, b),
                                    std::nextafter(b, a), -INFINITY, INFINITY};
    for (const auto& node : nodes)
      assert(std::isnan(node.slope) || node.slope >= 0);
    for (double x : queries) query(nodes, x);

    const std::vector<long double> slopes{0, std::numeric_limits<long double>::denorm_min(),
      1, std::numeric_limits<long double>::max(), INFINITY,
      std::numeric_limits<long double>::quiet_NaN()};
    // Exercise the query premise directly, including fallback after slope
    // overflow/NaN. These synthetic slopes are not reachability claims.
    for (auto left : slopes) for (auto right : slopes) {
      nodes[0].slope = left; nodes[1].slope = right;
      for (double x : queries) query(nodes, x);
    }
  }

  static void run() {
    const double tiny = std::numeric_limits<double>::denorm_min();
    const double maximum = std::numeric_limits<double>::max();
    const std::vector<std::pair<double, double>> spans{{-maximum, maximum},
      {0, 4 * tiny}, {-4 * tiny, 4 * tiny}, {0x1p-538, 0x1p-537},
      {0x1p-1022, 0x1p-1021}, {1, std::nextafter(std::nextafter(1.0, INFINITY), INFINITY)},
      {0x1p1022, maximum}, {-maximum, -0x1p1022}};
    const std::vector<std::pair<long double, long double>> ranks{{0, 0}, {0, 1},
      {1, 0}, {1, 0x1p-1022L}, {0x1p79L, 0x1p79L}, {0x1p53L, 1},
      {0, std::numeric_limits<long double>::denorm_min()}};
    for (auto span : spans) for (auto rank : ranks)
      interval(span.first, span.second, rank.first, rank.second);
    std::mt19937_64 random(0x4e554d45524943ULL);
    for (unsigned i = 0; i < 1000; ++i) {
      const auto first = random(), second = random();
      double a, b;
      std::memcpy(&a, &first, sizeof a); std::memcpy(&b, &second, sizeof b);
      if (!std::isfinite(a) || !std::isfinite(b) || a == b) continue;
      if (a > b) std::swap(a, b);
      interval(a, b, 0x1p53L, 1);
    }

    // Subtraction followed by re-addition need not recover an odd prefix.
    // Force both roundings, including when testing binary64 long double.
    long double unit = 1;
    for (int bit = 0; bit < std::numeric_limits<long double>::digits; ++bit) unit /= 2;
    volatile long double odd = 1 + 2 * unit;
    volatile long double small = unit;
    volatile long double difference = odd - small;
    volatile long double recovered = difference + small;
    assert(recovered != odd);
    assert(odd - recovered == 2 * unit);
  }
};
}

int main() {
  splinesketch::PaperSplineSketchInspector::run();
  std::cout << "Finite bracketed paper interpolation: "
            << splinesketch::PaperSplineSketchInspector::comparisons
            << " checks; odd-prefix reconstruction drift reproduced\n";
}
