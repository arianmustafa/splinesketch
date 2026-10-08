#include <splinesketch/splinesketch.hpp>

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <iostream>
#include <random>

namespace splinesketch {
struct SplineSketchInvariantInspector {
  struct Interval {
    std::uint64_t lower, upper, absorbed;
  };

  static Interval interval(const SplineSketch& sketch, double x) {
    std::uint64_t stored = 0, below = 0;
    for (const auto& entry : sketch.heavy_) {
      stored += entry.second.exact;
      if (entry.first <= x) below += entry.second.exact;
    }
    for (const auto& entry : sketch.pending_) {
      stored += entry.second;
      if (entry.first <= x) below += entry.second;
    }
    assert(stored <= sketch.count_);
    const auto absorbed = sketch.count_ - stored;
    return {below, below + absorbed, absorbed};
  }

  static std::uint64_t held(const SplineSketch& sketch, double x) {
    const auto it = sketch.heavy_.find(x);
    return it == sketch.heavy_.end() ? 0 : it->second.exact;
  }
};
}  // namespace splinesketch

using splinesketch::SplineSketch;
using Inspector = splinesketch::SplineSketchInvariantInspector;

static void check(const SplineSketch& sketch, std::vector<double> input) {
  std::sort(input.begin(), input.end());
  std::vector<double> queries{-INFINITY, -0.5, 0.5, INFINITY};
  for (double x : input) {
    queries.push_back(x);
    queries.push_back(std::nextafter(x, -INFINITY));
    queries.push_back(std::nextafter(x, INFINITY));
  }
  std::sort(queries.begin(), queries.end());
  queries.erase(std::unique(queries.begin(), queries.end()), queries.end());
  for (double x : queries) {
    const auto exact = static_cast<std::uint64_t>(
        std::upper_bound(input.begin(), input.end(), x) - input.begin());
    const auto interval = Inspector::interval(sketch, x);
    assert(interval.lower <= exact && exact <= interval.upper);
    const double estimate = sketch.rank(x);
    if (std::isfinite(estimate)) {
      const double certificate = std::max(
          std::fabs(estimate - static_cast<double>(interval.lower)),
          std::fabs(estimate - static_cast<double>(interval.upper)));
      assert(std::fabs(estimate - static_cast<double>(exact)) <= certificate);
    }
    if (interval.absorbed == 0) assert(estimate == exact);
  }
  for (auto it = input.begin(); it != input.end();) {
    const auto end = std::upper_bound(it, input.end(), *it);
    const auto frequency = static_cast<std::uint64_t>(end - it);
    const auto held = Inspector::held(sketch, *it);
    assert(held <= frequency);
    assert(frequency - held <= sketch.count() / sketch.bucket_capacity());
    it = end;
  }
}

static void low_distinct_is_exact() {
  const double largest = std::numeric_limits<double>::max();
  const double values[] = {-largest, -0.0, 1.0,
                          std::nextafter(1.0, 2.0), largest};
  SplineSketch sketch(6);
  std::vector<double> input;
  std::mt19937_64 rng(173);
  check(sketch, input);
  for (int i = 0; i < 12000; ++i) {
    const double x = values[rng() % 5];
    sketch.add(x);
    input.push_back(x);
    if (i % 131 == 0) {
      sketch.consolidate();
      assert(sketch.bucket_count() == 0);
      assert(Inspector::interval(sketch, 0).absorbed == 0);
      check(sketch, input);
    }
  }
  check(sketch, input);
}

static void evictions_and_reinsertions() {
  std::mt19937_64 rng(771);
  for (std::size_t k : {6U, 16U, 32U}) {
    SplineSketch sketch(k);
    std::vector<double> input;
    for (int i = 0; i < 1800; ++i) {
      // Alternate heavy hitters and many distinct values so previously held
      // values are evicted, absorbed, and subsequently reinserted.
      const double x = (i / 100) % 3 == 0 ? i % 3 :
          static_cast<double>(rng() % (3 * k));
      sketch.add(x);
      input.push_back(x);
      if (i % 31 == 0) check(sketch, input);
      if (i % 47 == 0) {
        sketch.consolidate();
        check(sketch, input);
      }
    }
    sketch.consolidate();
    assert(Inspector::interval(sketch, 0).absorbed > 0);
    check(sketch, input);
  }
}

static void interval_survives_numeric_failure() {
  const double largest = std::numeric_limits<double>::max();
  const std::vector<double> input{
      -largest, -1.5e308, -1e308, 1e308, 1.5e308, largest};
  SplineSketch sketch(6);
  for (double x : input) sketch.add(x);
  const auto interval = Inspector::interval(sketch, 0);
  assert(interval.lower == 0 && interval.upper == 6);
  assert(interval.absorbed == 6);
  if (std::numeric_limits<long double>::max_exponent == 1024)
    assert(std::isnan(sketch.rank(0)));
  check(sketch, input);
}

int main() {
  low_distinct_is_exact();
  evictions_and_reinsertions();
  interval_survives_numeric_failure();
  std::cout << "Implementation proof checks passed\n";
}
