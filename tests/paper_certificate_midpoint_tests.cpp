#include <splinesketch/paper_splinesketch.hpp>
#include "../benchmarks/comparison/experiments/paper_certificate_midpoint.hpp"
#include "paper_resize_history_fixture.hpp"
#include <cfenv>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using Sketch = splinesketch::CertifiedPaperSplineSketch;
using namespace splinesketch::experimental;

static void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

static std::uint64_t bits(double value) {
  std::uint64_t result;
  std::memcpy(&result, &value, sizeof result);
  return result;
}

namespace splinesketch {
struct PaperCertificateInspector {
  static void append_queries(const Sketch& sketch, std::vector<double>& queries) {
    for (const auto& node : sketch.nodes_) {
      queries.push_back(node.x);
      queries.push_back(std::nextafter(node.x, -INFINITY));
      queries.push_back(std::nextafter(node.x, INFINITY));
    }
  }
};
} // namespace splinesketch

static void arithmetic_checks() {
  for (std::uint64_t lower = 0; lower <= 64; ++lower)
    for (std::uint64_t upper = lower; upper <= 64; ++upper) {
      const auto r = certificate_midpoint(lower, upper);
      require(r.estimate == (lower + upper) / 2.0 && r.max_error == (upper - lower) / 2.0,
              "small exact midpoint failed");
    }
  const auto odd = certificate_midpoint((std::uint64_t{1} << 53) + 1,
                                       (std::uint64_t{1} << 53) + 2);
  require(odd.estimate == 0x1.0000000000001p53 && odd.max_error == 1,
          "half bit lost by double rounding");
  const auto top = std::numeric_limits<std::uint64_t>::max();
  const auto singleton = certificate_midpoint(top, top);
  require(singleton.estimate == 0x1p64 && singleton.max_error == 1, "2^64 boundary failed");
  const auto full = certificate_midpoint(0, top);
  require(full.estimate == 0x1p63 && full.max_error == 0x1p63, "full integer range failed");
  bool rejected = false;
  try { certificate_midpoint(2, 1); } catch (const std::invalid_argument&) { rejected = true; }
  require(rejected, "reversed bounds accepted");
  rejected = false;
  try { certificate_midpoint_uniform_error(2, 1); }
  catch (const std::invalid_argument&) { rejected = true; }
  require(rejected, "impossible uniform width accepted");
}

struct Comparison {
  unsigned better = 0, equal = 0, worse = 0;
  std::uint64_t queries = 0;
};

static void check_state(const Sketch& sketch, std::vector<double> data, Comparison& comparison) {
  std::sort(data.begin(), data.end());
  require(sketch.count() == data.size(), "count changed by query adapter");
  std::vector<double> queries{-INFINITY, INFINITY};
  for (std::size_t i = 0; i < data.size(); ++i) {
    queries.push_back(data[i]);
    queries.push_back(std::nextafter(data[i], -INFINITY));
    queries.push_back(std::nextafter(data[i], INFINITY));
    if (i) queries.push_back(data[i - 1] / 2 + data[i] / 2);
  }
  splinesketch::PaperCertificateInspector::append_queries(sketch, queries);
  std::sort(queries.begin(), queries.end());
  queries.erase(std::unique(queries.begin(), queries.end()), queries.end());
  const auto uniform = certificate_midpoint_max_error(sketch);
  double projected_worst = 0, midpoint_worst = 0;
  for (double query : queries) {
    const auto truth = static_cast<std::uint64_t>(std::upper_bound(data.begin(), data.end(), query) - data.begin());
    const auto projected = sketch.rank_with_error(query);
    const auto r = certificate_midpoint_rank(sketch, query);
    require(r.lower_rank == projected.lower_rank && r.upper_rank == projected.upper_rank,
            "adapter changed integer certificate");
    require(r.lower_rank <= truth && truth <= r.upper_rank, "certificate excludes truth");
    const double error = std::fabs(r.estimate - static_cast<double>(truth));
    require(error <= r.max_error && error <= uniform, "midpoint allowance too small");
    const auto projected_radius = std::max(std::fabs(projected.estimate - r.lower_rank),
                                         std::fabs(projected.estimate - r.upper_rank));
    require(r.max_error <= projected_radius, "midpoint lost certificate minimax comparison");
    require(r.estimate == (r.lower_rank + r.upper_rank) / 2.0 &&
            r.max_error == (r.upper_rank - r.lower_rank) / 2.0,
            "small reachable certificate midpoint not exact");
    projected_worst = std::max(projected_worst, std::fabs(projected.estimate - truth));
    midpoint_worst = std::max(midpoint_worst, error);
  }
  if (midpoint_worst < projected_worst) ++comparison.better;
  else if (midpoint_worst > projected_worst) ++comparison.worse;
  else ++comparison.equal;
  comparison.queries += queries.size();
}

static void history_checks(bool theoretical, Comparison& comparison) {
  const auto policy = theoretical ? Sketch::BoundPolicy::theoretical : Sketch::BoundPolicy::practical;
  std::mt19937_64 random(223);
  for (unsigned shape = 0; shape < 4; ++shape) {
    Sketch sketch(16, policy);
    std::vector<double> data;
    check_state(sketch, data, comparison);
    for (unsigned i = 0; i < 149; ++i) {
      const double value = shape == 0 ? static_cast<double>(random() % 17) - 8 :
          shape == 1 ? std::ldexp(static_cast<double>(random() % 31), -20) :
          shape == 2 ? static_cast<double>(i * i) : -static_cast<double>(i * i);
      sketch.add(value); data.push_back(value);
      if (i % 29 == 0) check_state(sketch, data, comparison); // Include nonempty exact buffers.
    }
    sketch.consolidate();
    check_state(sketch, data, comparison);
    for (unsigned step = 0; step < 8; ++step) {
      sketch.resize(step % 2 ? 16 : 8);
      check_state(sketch, data, comparison);
    }
    Sketch other(8, policy);
    for (unsigned i = 0; i < 47; ++i) { other.add(i / 4.0); data.push_back(i / 4.0); }
    sketch.merge(other);
    check_state(sketch, data, comparison);
    sketch.finalize();
    check_state(sketch, data, comparison);
  }
  Sketch cluster(16, policy);
  const std::vector<double> cluster_data(paper_test::resize_history_input.begin(),
                                       paper_test::resize_history_input.end());
  for (double x : cluster_data) cluster.add(x);
  cluster.consolidate();
  check_state(cluster, cluster_data, comparison);
  for (unsigned step = 0; step < 24; ++step) {
    cluster.resize(step % 2 ? 16 : 32);
    check_state(cluster, cluster_data, comparison);
  }

  // Both endpoints are reachable with identical active state (the separate
  // resize information test checks that premise). This makes radius 2.5 tight
  // at q=18, while also showing that a midpoint can hurt one actual stream.
  for (bool high : {false, true}) {
    Sketch witness(6, policy);
    std::vector<double> data;
    for (unsigned i = 0; i < 30; ++i) data.push_back(i * i);
    for (unsigned i = 1; i <= 5; ++i) data[i] = (high ? 0 : 18) + 3 * i;
    for (double value : data) witness.add(value);
    witness.consolidate();
    const auto r = certificate_midpoint_rank(witness, 18);
    const auto truth = high ? 6 : 1;
    require(r.lower_rank == 1 && r.upper_rank == 6 && r.estimate == 3.5 && r.max_error == 2.5 &&
            std::fabs(r.estimate - truth) == r.max_error, "tight reachable witness failed");
    if (high) require(std::fabs(witness.rank(18) - truth) < r.max_error,
                      "expected pointwise accuracy tradeoff missing");
    check_state(witness, data, comparison);
  }
  bool rejected = false;
  try { certificate_midpoint_rank(cluster, std::numeric_limits<double>::quiet_NaN()); }
  catch (const std::invalid_argument&) { rejected = true; }
  require(rejected, "NaN query accepted");
}

static void intervals() {
  std::uint64_t lower, upper;
  while (std::cin >> lower >> upper) {
    std::cout << "{\"lower\":" << lower << ",\"upper\":" << upper << ",\"modes\":[";
    bool first = true;
    for (int mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO}) {
      require(std::fesetround(mode) == 0, "rounding mode unavailable");
      const auto r = certificate_midpoint(lower, upper);
      if (!first) std::cout << ',';
      first = false;
      std::cout << '[' << bits(r.estimate) << ',' << bits(r.max_error) << ','
                << bits(certificate_midpoint_uniform_error(upper - lower, upper)) << ']';
    }
    std::cout << "]}\n";
  }
  require(std::cin.eof(), "invalid interval input");
  require(std::fesetround(FE_TONEAREST) == 0, "rounding mode restore failed");
}

static void large_counts(bool emit) {
  for (bool theoretical : {false, true}) {
    Sketch sketch(8, theoretical ? Sketch::BoundPolicy::theoretical : Sketch::BoundPolicy::practical);
    sketch.add(0);
    std::uint64_t zeros = 1, count = 1;
    for (unsigned step = 1; step <= 63; ++step) {
      // Public operations reach count 2^64-1 without materializing the stream.
      sketch.merge(sketch);
      const int value = step % 2;
      sketch.add(value);
      sketch.consolidate();
      zeros = 2 * zeros + (value == 0);
      count = 2 * count + 1;
      require(sketch.count() == count, "large public history count differs");
      if (step < 51) continue;
      for (int query : {-1, 0, 1}) {
        const auto truth = query < 0 ? 0 : query == 0 ? zeros : count;
        const auto r = certificate_midpoint_rank(sketch, query);
        require(r.lower_rank <= truth && truth <= r.upper_rank, "large certificate excludes truth");
        require(std::isfinite(r.estimate) && std::isfinite(r.max_error), "nonfinite large query");
        if (step == 63 && query == 1)
          require(r.estimate == 0x1p64 && r.max_error == 1, "public UINT64_MAX return allowance failed");
        if (emit)
          std::cout << "{\"step\":" << step << ",\"theoretical\":" << theoretical
                    << ",\"query\":" << query << ",\"count\":" << count << ",\"truth\":" << truth
                    << ",\"lower\":" << r.lower_rank << ",\"upper\":" << r.upper_rank
                    << ",\"estimate_bits\":" << bits(r.estimate) << ",\"radius_bits\":" << bits(r.max_error)
                    << ",\"uniform_bits\":" << bits(certificate_midpoint_max_error(sketch)) << "}\n";
      }
    }
  }
}

int main(int argc, char** argv) {
  if (argc == 2 && std::string(argv[1]) == "--intervals") { intervals(); return 0; }
  if (argc == 2 && std::string(argv[1]) == "--large-counts") { large_counts(true); return 0; }
  require(argc == 1, "usage: midpoint-tests [--intervals|--large-counts]");
  arithmetic_checks();
  large_counts(false);
  Comparison comparison;
  for (bool theoretical : {false, true}) history_checks(theoretical, comparison);
  require(comparison.queries > 10000, "missing history coverage");
  std::cout << "{\"queries\":" << comparison.queries << ",\"checkpoints_better\":" << comparison.better
            << ",\"checkpoints_equal\":" << comparison.equal << ",\"checkpoints_worse\":"
            << comparison.worse << "}\n";
}
