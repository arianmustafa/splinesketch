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
#include <stdexcept>
#include <vector>

using splinesketch::SplineSketch;

static void check_exact_small() {
  SplineSketch sketch(16);
  assert(sketch.rank(0) == 0);
  bool threw = false;
  try { sketch.quantile(0.5); } catch (const std::logic_error&) { threw = true; }
  assert(threw);
  for (double x : {-4.0, -4.0, 0.0, 2.0, 2.0, 2.0, 9.0}) sketch.add(x);
  assert(sketch.count() == 7);
  assert(sketch.rank(-5) == 0);
  assert(sketch.rank(-4) == 2);
  assert(sketch.rank(1) == 3);
  assert(sketch.rank(2) == 6);
  assert(sketch.rank(9) == 7);
  assert(sketch.quantile(0) == -4);
  assert(sketch.quantile(0.5) == 2);
  assert(sketch.quantile(1) == 9);
  threw = false;
  try { sketch.add(std::numeric_limits<double>::quiet_NaN()); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);
}

static void check_stream_and_queries() {
  SplineSketch sketch(64);
  std::mt19937_64 rng(17);
  std::normal_distribution<double> normal(0, 1);
  std::vector<double> data;
  for (int i = 0; i < 30000; ++i) {
    const double x = i % 7 == 0 ? 3.0 : normal(rng);
    sketch.add(x);
    data.push_back(x);
  }
  sketch.consolidate();
  std::sort(data.begin(), data.end());
  assert(sketch.bucket_count() <= sketch.bucket_capacity());
  assert(sketch.rank(-std::numeric_limits<double>::infinity()) == 0);
  assert(sketch.rank(std::numeric_limits<double>::infinity()) == data.size());
  double previous = 0;
  double worst = 0;
  for (int i = -50; i <= 50; ++i) {
    const double x = i / 10.0;
    const double actual = std::upper_bound(data.begin(), data.end(), x) - data.begin();
    const double estimated = sketch.rank(x);
    assert(estimated >= previous - 1e-7);
    previous = estimated;
    worst = std::max(worst, std::abs(actual - estimated));
  }
  assert(worst < 0.025 * data.size());
  assert(sketch.rank(3.0) - sketch.rank(std::nextafter(3.0, -INFINITY)) > 4000);
  for (double q : {0.01, 0.1, 0.5, 0.9, 0.99}) {
    const double x = sketch.quantile(q);
    assert(std::isfinite(x));
    assert(x >= data.front() && x <= data.back());
    assert(sketch.rank(x) + 1e-6 >= std::ceil(q * data.size()));
  }
}

static void check_merge_and_resize() {
  SplineSketch left(48), right(48);
  std::mt19937_64 rng(29);
  std::uniform_real_distribution<double> uniform(-100, 100);
  for (int i = 0; i < 20000; ++i) {
    (i & 1 ? left : right).add(i % 5 == 0 ? 42 : uniform(rng));
  }
  left.merge(right);
  assert(left.count() == 20000);
  assert(left.rank(std::numeric_limits<double>::infinity()) == 20000);
  assert(left.bucket_count() <= 48);
  const double heavy_rank = left.rank(42) - left.rank(std::nextafter(42.0, -INFINITY));
  assert(heavy_rank >= 3900);
  left.resize(24);
  assert(left.bucket_count() <= 24);
  assert(left.heavy_hitter_count() < 24);
  left.resize(96);
  assert(left.bucket_count() <= 96);
  assert(left.rank(std::numeric_limits<double>::infinity()) == 20000);
  SplineSketch copy = left;
  copy.merge(copy);
  assert(copy.count() == 40000);
  assert(copy.rank(std::numeric_limits<double>::infinity()) == 40000);
}

static void check_extremes() {
  SplineSketch sketch(8);
  for (int i = 0; i < 100; ++i) sketch.add(i % 2 ? 1e308 : -1e308);
  assert(sketch.rank(-1e308) == 50);
  assert(sketch.rank(1e308) == 100);
  assert(sketch.quantile(0.5) == -1e308);
  assert(sketch.quantile(0.51) == 1e308);
  sketch.resize(6);
  assert(sketch.rank(1e308) == 100);
}

static void check_resize_with_full_heavy_table() {
  SplineSketch sketch(48);
  for (int i = 0; i < 30; ++i) sketch.add(i);
  assert(sketch.heavy_hitter_count() == 30);
  sketch.resize(6);
  assert(sketch.heavy_hitter_count() < 6);
  assert(sketch.rank(29) == 30);
  double previous = 0;
  for (int i = -1; i <= 30; ++i) {
    const double current = sketch.rank(i);
    assert(current >= previous);
    previous = current;
  }
  for (int i = 30; i < 1000; ++i) sketch.add(i);
  assert(sketch.rank(999) == 1000);
}

static void check_merge_tree() {
  std::mt19937_64 rng(101);
  std::normal_distribution<double> normal;
  std::vector<SplineSketch> sketches;
  std::vector<double> values;
  for (int part = 0; part < 8; ++part) {
    sketches.emplace_back(64);
    for (int i = 0; i < 4000; ++i) {
      const double value = i % 11 == 0 ? -2.0 : normal(rng) + part * 0.25;
      sketches.back().add(value);
      values.push_back(value);
    }
  }
  while (sketches.size() > 1) {
    std::vector<SplineSketch> next;
    for (std::size_t i = 0; i < sketches.size(); i += 2) {
      sketches[i].merge(sketches[i + 1]);
      next.push_back(std::move(sketches[i]));
    }
    sketches = std::move(next);
  }
  auto& sketch = sketches.front();
  std::sort(values.begin(), values.end());
  assert(sketch.count() == values.size());
  double previous = 0, worst = 0;
  for (int i = -40; i <= 60; ++i) {
    const double query = i / 10.0;
    const double actual = std::upper_bound(values.begin(), values.end(), query) - values.begin();
    const double estimated = sketch.rank(query);
    assert(estimated >= previous - 1e-7);
    previous = estimated;
    worst = std::max(worst, std::abs(actual - estimated));
  }
  assert(worst < 0.05 * values.size());
}

static void check_inverse(const SplineSketch& sketch) {
  for (double q : {std::nextafter(0.0, 1.0), 0.001, 0.1, 0.333, 0.5,
                   0.9, 0.999, std::nextafter(1.0, 0.0)}) {
    const long double target = std::ceil(q * static_cast<long double>(sketch.count()));
    const double x = sketch.quantile(q);
    assert(std::isfinite(x));
    assert(sketch.rank(x) >= target);
    assert(sketch.rank(std::nextafter(x, -INFINITY)) < target);
  }
}

static void check_prepared_quantiles() {
  // Exercise both sides of the small-table cutoff, including discontinuities,
  // adjacent doubles, signed zero, and values spanning the finite double range.
  for (int distinct : {1, 16, 17, 31}) {
    SplineSketch sketch(64);
    for (int i = 0; i < distinct; ++i) {
      double x = i - 15.0;
      if (i == 0) x = -1e308;
      if (i == 1) x = 1e308;
      if (i == 2) x = std::nextafter(1.0, 2.0);
      sketch.add(x);
      sketch.add(x);
    }
    sketch.add(-0.0);
    sketch.add(0.0);
    check_inverse(sketch);
  }

  SplineSketch sketch(32);
  for (int i = 0; i < 100; ++i) sketch.add(0);
  for (int i = 1; i <= 31; ++i) sketch.add(i);
  // One surviving heavy hitter and a not-yet-full pending buffer.
  assert(sketch.bucket_count() == 0);
  assert(sketch.heavy_hitter_count() == 1);
  check_inverse(sketch);
  sketch.consolidate();
  check_inverse(sketch);

  std::mt19937_64 rng(883);
  std::normal_distribution<double> normal;
  for (int i = 0; i < 8192; ++i) sketch.add(i % 5 ? normal(rng) : 0);
  check_inverse(sketch);
  sketch.consolidate();
  check_inverse(sketch);
  // Weighted merges expose rounding differences in prepared prefix counts.
  for (int i = 0; i < 39; ++i) sketch.merge(sketch);
  assert(sketch.count() < (std::uint64_t{1} << 53));
  check_inverse(sketch);
}

static void check_rebuilds_after_structural_edits() {
  SplineSketch sketch(64), other(96);
  std::mt19937_64 rng(922);
  std::normal_distribution<double> normal;
  for (int i = 0; i < 4096; ++i) {
    sketch.add(normal(rng));
    other.add(normal(rng) + 2);
  }
  sketch.consolidate();
  other.consolidate();
  // Both heavy tables and buffers are empty: merge must restore slopes even
  // when the final consolidate() has nothing to do.
  assert(sketch.heavy_hitter_count() == 0);
  // Finish the partial 96-item cycle before merging.
  for (int i = 4096; i < 4128; ++i) other.add(normal(rng) + 2);
  other.consolidate();
  assert(other.heavy_hitter_count() == 0);
  sketch.merge(other);
  check_inverse(sketch);
  for (std::size_t capacity : {6, 96, 12, 64}) {
    sketch.resize(capacity);
    assert(sketch.bucket_count() <= capacity);
    check_inverse(sketch);
    for (int i = 0; i < 131; ++i) sketch.add(normal(rng) - 2);
    sketch.consolidate();
    check_inverse(sketch);
  }
}

int main() {
  check_exact_small();
  check_stream_and_queries();
  check_merge_and_resize();
  check_extremes();
  check_resize_with_full_heavy_table();
  check_merge_tree();
  check_prepared_quantiles();
  check_rebuilds_after_structural_edits();
  std::cout << "SplineSketch tests passed\n";
}
