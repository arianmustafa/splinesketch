#include <splinesketch/splinesketch.hpp>

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

namespace splinesketch {

// Test-only access to completed states. This does not certify the paper's
// bucket-size or split-history lemmas, which are still open in the audit.
struct SplineSketchInvariantInspector {
  static void check(const SplineSketch& sketch) {
    assert(sketch.nodes_.size() <= sketch.capacity_);
    assert(sketch.heavy_.size() < sketch.capacity_);
    long double bucket_sum = 0;
    for (std::size_t i = 0; i < sketch.nodes_.size(); ++i) {
      const auto& node = sketch.nodes_[i];
      assert(std::isfinite(node.x));
      assert(std::isfinite(node.mass) && node.mass >= 0);
      assert(std::isfinite(node.prefix) && node.prefix >= 0);
      assert(std::isfinite(node.slope) && node.slope >= 0);
      if (i) assert(sketch.nodes_[i - 1].x < node.x);
      bucket_sum += node.mass;
      const long double tolerance = 64 * std::numeric_limits<long double>::epsilon() *
                                    static_cast<long double>(sketch.count_) + 1e-12L;
      assert(std::fabs(bucket_sum - node.prefix) <= tolerance);
    }
    long double total = bucket_sum;
    for (const auto& item : sketch.pending_) {
      assert(std::isfinite(item.first) && item.second > 0);
      total += item.second;
    }
    for (const auto& item : sketch.heavy_) {
      assert(std::isfinite(item.first));
      assert(item.second.residual > 0 && item.second.residual <= item.second.exact);
      total += item.second.exact;
    }
    const long double tolerance = 64 * std::numeric_limits<long double>::epsilon() *
                                  static_cast<long double>(sketch.count_) + 1e-12L;
    assert(std::fabs(total - static_cast<long double>(sketch.count_)) <= tolerance);
  }
};
}  // namespace splinesketch

using splinesketch::SplineSketch;
using splinesketch::SplineSketchInvariantInspector;

static void check_queries(const SplineSketch& sketch, const std::vector<double>& input) {
  SplineSketchInvariantInspector::check(sketch);
  std::vector<double> sorted_input = input;
  std::sort(sorted_input.begin(), sorted_input.end());
  std::vector<double> queries = input;
  queries.push_back(-std::numeric_limits<double>::infinity());
  queries.push_back(std::numeric_limits<double>::infinity());
  for (double x : input) {
    queries.push_back(std::nextafter(x, -INFINITY));
    queries.push_back(std::nextafter(x, INFINITY));
  }
  std::sort(queries.begin(), queries.end());
  queries.erase(std::unique(queries.begin(), queries.end()), queries.end());
  double previous = 0;
  for (double x : queries) {
    const double estimate = sketch.rank(x);
    const double exact = static_cast<double>(
        std::upper_bound(sorted_input.begin(), sorted_input.end(), x) - sorted_input.begin());
    assert(std::isfinite(estimate));
    assert(estimate >= 0 && estimate <= sketch.count());
    assert(estimate + 1e-9 >= previous);
    if (input.size() < sketch.bucket_capacity()) assert(estimate == exact);
    assert(std::fabs(estimate - exact) <= sketch.count());
    previous = estimate;
  }
  assert(sketch.rank(INFINITY) == input.size());
}

static void exhaustive_streams() {
  const std::vector<double> alphabet{-2, 0, 1, 3};
  for (std::size_t capacity : {6U, 8U, 12U}) {
    for (unsigned length = 0; length <= 7; ++length) {
      std::uint32_t cases = 1;
      for (unsigned i = 0; i < length; ++i) cases *= alphabet.size();
      for (std::uint32_t code = 0; code < cases; ++code) {
        SplineSketch sketch(capacity);
        std::vector<double> input;
        std::uint32_t digits = code;
        for (unsigned i = 0; i < length; ++i) {
          const double x = alphabet[digits % alphabet.size()];
          digits /= alphabet.size();
          sketch.add(x);
          input.push_back(x);
          check_queries(sketch, input);  // includes any automatic consolidation
        }
        sketch.consolidate();
        check_queries(sketch, input);
      }
    }
  }
}

static void edge_and_operation_histories() {
  const double near_one = std::nextafter(1.0, 2.0);
  const std::vector<double> alphabet{-1e308, -1.0, 0.0, 1.0, near_one, 1e308};
  std::mt19937_64 rng(0x51504f);
  for (int trial = 0; trial < 50; ++trial) {
    SplineSketch left(6), right(6);
    std::vector<double> left_data, right_data;
    for (int i = 0; i < 80; ++i) {
      const double x = alphabet[rng() % alphabet.size()];
      auto& sketch = i & 1 ? left : right;
      auto& data = i & 1 ? left_data : right_data;
      sketch.add(x);
      data.push_back(x);
      check_queries(sketch, data);
    }
    left.merge(right);  // equal-capacity, balanced inputs
    left_data.insert(left_data.end(), right_data.begin(), right_data.end());
    check_queries(left, left_data);
    left.resize(8);  // first and only resize on this merged sketch
    check_queries(left, left_data);
  }
}

static void longer_adversarial_streams() {
  for (std::size_t capacity : {6U, 16U, 32U}) {
    for (int shape = 0; shape < 3; ++shape) {
      SplineSketch sketch(capacity);
      std::vector<double> data;
      for (int i = 0; i < 3000; ++i) {
        const double x = shape == 0 ? (i % 19 == 0 ? 1e100 : i % 7 - 3.0) :
                         shape == 1 ? static_cast<double>(i) :
                         (i % 13 == 0 ? -1e100 : std::nextafter(1.0, 2.0));
        sketch.add(x);
        data.push_back(x);
        if (i % 53 == 0) check_queries(sketch, data);
      }
      sketch.consolidate();
      check_queries(sketch, data);
    }
  }
}

int main() {
  exhaustive_streams();
  edge_and_operation_histories();
  longer_adversarial_streams();
  std::cout << "Completed-state invariants and short-stream enumeration passed\n";
}
