#include <splinesketch/paper_splinesketch.hpp>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <iostream>
#include <random>
#include <string>

namespace splinesketch {
struct PaperSplineSketchInspector {
  using Sketch = PaperSplineSketch;
  using Node = Sketch::Node;
  using Item = std::pair<double, std::uint64_t>;

  static std::uint64_t bits(double x) {
    std::uint64_t result;
    std::memcpy(&result, &x, sizeof result);
    return result;
  }
  static double from_bits(std::uint64_t value) {
    double result;
    std::memcpy(&result, &value, sizeof result);
    return result;
  }
  // Exact dyadic serialization, including long-double bits that would be
  // lost by conversion to double. Avoid the platform frexpl ABI in the
  // -mlong-double-64 test. All trace values are nonnegative and below 2^64.
  static void dyadic(long double value) {
    assert(std::isfinite(value) && value >= 0 && value < 1000000);
    const auto whole = static_cast<std::uint64_t>(value);
    auto remainder = value - static_cast<long double>(whole);
    std::string digits;
    while (remainder != 0) {
      assert(digits.size() < 2048);
      remainder *= 2;
      const bool one = remainder >= 1;
      digits.push_back(one ? '1' : '0');
      if (one) remainder -= 1;
    }
    std::cout << '"' << whole << ':' << digits << '"';
  }
  static void items(const std::vector<Item>& values) {
    std::cout << '[';
    bool first = true;
    for (const auto& value : values) {
      if (!first) std::cout << ',';
      first = false;
      std::cout << '[' << bits(value.first) << ',' << value.second << ']';
    }
    std::cout << ']';
  }
  static void state(const char* kind, const std::vector<Node>& nodes,
                    const std::vector<Item>& incoming = {}) {
    std::cout << "{\"kind\":\"" << kind << "\",\"incoming\":";
    items(incoming);
    std::cout << ",\"nodes\":[";
    bool first = true;
    for (const auto& node : nodes) {
      if (!first) std::cout << ',';
      first = false;
      std::cout << '[' << bits(node.x) << ',';
      dyadic(node.mass);
      std::cout << ',';
      dyadic(node.prefix);
      std::cout << ',' << (node.protected_threshold ? "true" : "false") << ']';
    }
    std::cout << "]}\n";
  }
  static void geometry() {
    const auto emit = [](double a, double b) {
      if (a > b) std::swap(a, b);
      if (!(std::isfinite(a) && std::isfinite(b) && a < b)) return;
      std::cout << bits(a) << ' ' << bits(b) << ' ' << bits(Sketch::midpoint(a, b)) << '\n';
    };
    const double tiny = std::numeric_limits<double>::denorm_min();
    for (int a = -20; a <= 20; ++a)
      for (int b = a + 1; b <= 20; ++b) emit(a * tiny, b * tiny);
    for (int exponent = -1022; exponent <= 1023; ++exponent) {
      const double boundary = std::ldexp(1.0, exponent);
      std::vector<double> near{boundary};
      double left = boundary, right = boundary;
      for (int i = 0; i < 6; ++i) {
        left = std::nextafter(left, 0);
        right = std::nextafter(right, INFINITY);
        near.push_back(left); near.push_back(right);
      }
      for (std::size_t a = 0; a < near.size(); ++a)
        for (std::size_t b = a + 1; b < near.size(); ++b) {
          emit(near[a], near[b]); emit(-near[a], -near[b]);
        }
    }
    std::mt19937_64 random(0x50524f4f46ULL);
    for (int i = 0; i < 20000; ++i) {
      const double a = from_bits(random()), b = from_bits(random());
      emit(a, b);
    }
  }
  static void trace(bool heaps) {
    Sketch sketch(12);
    sketch.factor_ = 6;
    const double tiny = std::numeric_limits<double>::denorm_min();
    const double boundary = std::ldexp(1.0, -1021);
    // The initial exact ranks are roots. Include the double-rounding gap
    // beside the binade boundary and several ordinary interpolation gaps.
    std::vector<Item> initial{{-8, 1}, {-4, 2}, {-2, 3}, {-1, 4}, {0, 1},
      {boundary - tiny, 2}, {boundary + 4 * tiny, 3}, {1, 2}, {2, 1},
      {4, 3}, {8, 2}, {16, 1}};
    for (auto item : initial) sketch.count_ += item.second;
    Sketch::Snapshot exact{{}, initial, {}};
    exact.prepare(); sketch.initialize(exact);
    state("init", sketch.nodes_, initial);
    const std::vector<Item> incoming{{-16, 1}, {boundary, 2}, {3, 1}, {32, 1}};
    for (auto item : incoming) sketch.count_ += item.second;
    Sketch::Snapshot before{sketch.nodes_, incoming, {}};
    before.prepare();
    if (!heaps) {
      sketch.nodes_.insert(sketch.nodes_.begin(), Node{-16});
      sketch.nodes_.emplace_back(32);
      sketch.reestimate(before);
      state("update", sketch.nodes_, incoming);
      sketch.join_at(1); state("join", sketch.nodes_);
      sketch.join_at(sketch.nodes_.size() - 2); state("join", sketch.nodes_);
      // Split the binade-boundary interval repeatedly against one snapshot.
      auto it = std::lower_bound(sketch.nodes_.begin(), sketch.nodes_.end(), boundary + 4 * tiny,
          [](const Node& node, double x) { return node.x < x; });
      const auto at = static_cast<std::size_t>(it - sketch.nodes_.begin());
      sketch.join_at(2); state("join", sketch.nodes_);
      sketch.split_at(at - 1, before); state("split", sketch.nodes_);
      sketch.join_at(sketch.nodes_.size() - 2); state("join", sketch.nodes_);
      sketch.split_at(at, before); state("split", sketch.nodes_);
      Sketch::rebuild(sketch.nodes_);
    } else {
      Sketch::HeapRebalance heap(sketch, before, false);
      heap.reestimate(before); state("update", heap.snapshot_nodes(), incoming);
      heap.join(heap.links[heap.head].next); state("join", heap.snapshot_nodes());
      heap.join(heap.links[heap.tail].prev); state("join", heap.snapshot_nodes());
      auto right = heap.head;
      while (sketch.nodes_[right].x != boundary + 4 * tiny) right = heap.links[right].next;
      heap.join(heap.links[heap.links[heap.head].next].next);
      state("join", heap.snapshot_nodes());
      heap.split(right); state("split", heap.snapshot_nodes());
      heap.join(heap.links[heap.tail].prev); state("join", heap.snapshot_nodes());
      heap.split(right); state("split", heap.snapshot_nodes());
      heap.finish();
    }
    // Shift fractional prefixes across a significand binade. This introduces
    // actual summation/conversion defects, beyond interpolation uncertainty.
    const std::vector<Item> shifted{{-16, 64}};
    sketch.count_ += 64;
    Sketch::Snapshot next{sketch.nodes_, shifted, {}};
    next.prepare();
    if (heaps) {
      Sketch::HeapRebalance heap(sketch, next, false);
      heap.reestimate(next); state("update", heap.snapshot_nodes(), shifted);
      heap.finish();
    } else {
      sketch.reestimate(next); state("update", sketch.nodes_, shifted);
      Sketch::rebuild(sketch.nodes_);
    }
    state("finish", sketch.nodes_);
    std::vector<double> queries{-INFINITY, INFINITY, -3, 0.25, 1.5, 3, 6, 12};
    for (auto node : sketch.nodes_) {
      queries.push_back(node.x);
      queries.push_back(std::nextafter(node.x, -INFINITY));
      queries.push_back(std::nextafter(node.x, INFINITY));
    }
    for (double x : queries) {
      std::cout << "{\"kind\":\"query\",\"x\":" << bits(x) << ",\"rank\":";
      dyadic(sketch.rank(x));
      std::cout << "}\n";
    }
  }
  static void regressions() {
    const double tiny = std::numeric_limits<double>::denorm_min();
    const double boundary = std::ldexp(1.0, -1021);
    const double a = boundary - tiny, b = boundary + 4 * tiny;
    assert(Sketch::midpoint(a, b) == boundary);
    assert(Sketch::midpoint(-tiny, 3 * tiny) == 2 * tiny);
    assert(Sketch::midpoint(-3 * tiny, tiny) == -2 * tiny);
    assert(Sketch::midpoint(-std::numeric_limits<double>::max(),
                             std::numeric_limits<double>::max()) == 0);
    const auto adjacent = std::nextafter(1.0, INFINITY);
    const auto mid = Sketch::midpoint(1.0, adjacent);
    assert(mid == 1.0 || mid == adjacent);
  }
};
} // namespace splinesketch

int main(int argc, char** argv) {
  using Inspector = splinesketch::PaperSplineSketchInspector;
  Inspector::regressions();
  if (argc == 2 && std::string(argv[1]) == "--geometry") Inspector::geometry();
  else if (argc == 2 && std::string(argv[1]) == "--trace") {
    std::cout << "{\"kind\":\"format\",\"precision\":"
              << std::numeric_limits<long double>::digits << "}\n";
    Inspector::trace(false);
    std::cout << "{\"kind\":\"format\",\"precision\":"
              << std::numeric_limits<long double>::digits << "}\n";
    Inspector::trace(true);
  } else std::cout << "Paper proof midpoint regressions passed\n";
}
