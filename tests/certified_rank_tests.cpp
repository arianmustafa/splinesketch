#include <splinesketch/certified_splinesketch.hpp>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <iostream>
#include <random>

using splinesketch::CertifiedSplineSketch;
using splinesketch::SplineSketch;
using splinesketch::detail::RankEnvelope;

static void inspect(const RankEnvelope& envelope, const std::vector<double>& data) {
  assert(envelope.count() == data.size());
  std::uint64_t prefix = 0;
  for (std::size_t i = 0; i < envelope.size(); ++i) {
    const auto entry = envelope.entry(i);
    assert(entry.prefix > prefix);
    const auto gap = entry.prefix - prefix;
    prefix = entry.prefix;
    assert(gap + entry.delta <= envelope.threshold());
    const auto strict = std::count_if(data.begin(), data.end(),
                                     [&](double x) { return x < entry.value; });
    const auto inclusive = std::count_if(data.begin(), data.end(),
                                        [&](double x) { return x <= entry.value; });
    // A retained occurrence's possible positions intersect the exact tie block.
    assert(entry.prefix <= static_cast<std::uint64_t>(inclusive));
    assert(entry.prefix + entry.delta > static_cast<std::uint64_t>(strict));
  }
  assert(prefix == envelope.count());
  assert(envelope.uncertainty() < envelope.threshold());
  for (double query : {-10., -1., -.5, 0., .5, 1., 1.5, 2., 2.5, 3., 4., 10.}) {
    const auto truth = std::count_if(data.begin(), data.end(),
                                     [&](double x) { return x <= query; });
    const auto bounds = envelope.bounds(query);
    assert(bounds.lower <= static_cast<std::uint64_t>(truth));
    assert(bounds.upper >= static_cast<std::uint64_t>(truth));
    assert(bounds.upper - bounds.lower <= envelope.uncertainty());
  }
}

static void same_envelope(const RankEnvelope& a, const RankEnvelope& b) {
  assert(a.count() == b.count() && a.size() == b.size());
  assert(a.uncertainty() == b.uncertainty());
  for (std::size_t i = 0; i < a.size(); ++i) {
    const auto x = a.entry(i), y = b.entry(i);
    assert(x.value == y.value && x.prefix == y.prefix && x.delta == y.delta);
  }
}

static void local_compaction_equivalence() {
  std::mt19937_64 random(72839461);
  for (std::size_t capacity : {6U, 8U, 16U, 32U, 64U, 128U, 257U}) {
    RankEnvelope local(capacity), full(capacity);
    for (unsigned i = 0; i < 12000; ++i) {
      double value = static_cast<double>(random() % 1000001) / 1000;
      if (i % 7 == 0) value = static_cast<int>(random() % 3) - 1;
      if (i % 53 == 0) value = -static_cast<double>(i) - 10000; // New minimum.
      if (i % 59 == 0) value = static_cast<double>(i) + 10000;  // New maximum.
      local.commit_add(local.prepare_add(value));
      full.commit_add(full.prepare_add(value));
      full.resize(capacity); // Force the original full reverse sweep.
      same_envelope(local, full);
    }
  }
}

static void exhaustive_envelopes() {
  // Capacity 2 exercises compression in short histories. Public sketches
  // retain their minimum capacity of 6.
  for (unsigned code = 0; code < 6561; ++code) {
    RankEnvelope direct(2), a(2), b(3);
    std::vector<double> data, adata, bdata;
    unsigned digits = code;
    for (unsigned i = 0; i < 8; ++i) {
      const double value = digits % 3;
      digits /= 3;
      direct.commit_add(direct.prepare_add(value));
      data.push_back(value);
      inspect(direct, data);
      auto& part = i % 2 ? a : b;
      auto& values = i % 2 ? adata : bdata;
      part.commit_add(part.prepare_add(value));
      values.push_back(value);
    }
    inspect(a.merged(b), data);
    inspect(b.merged(a), data);
    auto merged = a.merged(b);
    auto doubled = data;
    doubled.insert(doubled.end(), data.begin(), data.end());
    inspect(merged.merged(merged), doubled);
    merged.resize(1);
    inspect(merged, data);
    merged.resize(20);
    assert(merged.capacity() == 1);
    inspect(merged, data);
  }
}

static void verify(const CertifiedSplineSketch& certified, const SplineSketch& raw,
                   const std::vector<double>& data) {
  assert(certified.count() == data.size());
  assert(certified.max_rank_uncertainty() <
         std::max<std::uint64_t>(1, certified.count() / certified.certificate_capacity()));
  std::vector<double> queries = data;
  queries.insert(queries.end(), {-INFINITY, -1e308, -10, 0, .5, 10, 1e308, INFINITY});
  std::sort(queries.begin(), queries.end());
  const auto values = queries;
  for (double value : values) {
    queries.push_back(std::nextafter(value, -INFINITY));
    queries.push_back(std::nextafter(value, INFINITY));
  }
  for (double query : queries) {
    const auto truth = static_cast<std::uint64_t>(std::count_if(data.begin(), data.end(),
                                          [&](double x) { return x <= query; }));
    const double exact_rank = static_cast<double>(truth);
    const auto result = certified.rank_with_error(query);
    assert(result.lower_rank <= truth && truth <= result.upper_rank);
    assert(result.estimate == certified.rank(query));
    assert(std::isfinite(result.estimate) && std::isfinite(result.max_error));
    assert(std::fabs(result.estimate - exact_rank) <= result.max_error);
    assert(std::fabs(result.estimate - exact_rank) <= certified.max_rank_error());
    const double unguarded = raw.rank(query);
    if (std::isfinite(unguarded))
      assert(std::fabs(result.estimate - exact_rank) <= std::fabs(unguarded - exact_rank));
  }
}

static void operations() {
  std::mt19937_64 random(453879);
  for (unsigned trial = 0; trial < 60; ++trial) {
    const std::size_t capacity = 6 + random() % 27;
    CertifiedSplineSketch sketch(capacity);
    SplineSketch raw(capacity);
    std::vector<double> data;
    for (unsigned batch = 0; batch < 8; ++batch) {
      for (unsigned i = 0; i < 50; ++i) {
        const double value = i % 11 ? static_cast<int>(random() % 41) - 20
                                   : (random() % 2 ? 1e6 : -1e6);
        sketch.add(value); raw.add(value); data.push_back(value);
      }
      verify(sketch, raw, data);
      if (batch % 3 == 0) {
        const std::size_t next_capacity = 6 + random() % 27;
        sketch.resize(next_capacity); raw.resize(next_capacity);
      } else if (batch % 3 == 1) {
        sketch.consolidate(); raw.consolidate();
      } else {
        CertifiedSplineSketch other(6 + random() % 27);
        SplineSketch other_raw(other.bucket_capacity());
        for (unsigned i = 0; i < 45; ++i) {
          const double value = static_cast<int>(random() % 41) - 20;
          other.add(value); other_raw.add(value); data.push_back(value);
        }
        sketch.merge(other); raw.merge(other_raw);
      }
      verify(sketch, raw, data);
    }
    sketch.merge(sketch); raw.merge(raw);
    const auto copy = data;
    data.insert(data.end(), copy.begin(), copy.end());
    verify(sketch, raw, data);
    assert(sketch.quantile(0) == *std::min_element(data.begin(), data.end()));
    assert(sketch.quantile(1) == *std::max_element(data.begin(), data.end()));
    for (double q : {.1, .5, .9}) {
      const double value = sketch.quantile(q);
      assert(std::isfinite(value));
      assert(sketch.rank(value) >= std::ceil(q * static_cast<long double>(sketch.count())));
    }
  }
}

static void packed_boundaries() {
  CertifiedSplineSketch maximum(6), power(6);
  power.add(0);
  for (unsigned bit = 0; bit < 32; ++bit) {
    maximum.merge(power);
    if (bit != 31) power.merge(power);
  }
  assert(maximum.count() == std::numeric_limits<std::uint32_t>::max());
  maximum.add(0);  // Insertion promotes to 64-bit entries.
  CertifiedSplineSketch small(6);
  small.add(-1); small.add(1);
  auto a = maximum, b = small;
  a.merge(small); b.merge(maximum);  // Both directions of mixed storage merge.
  for (double query : {-2., -1., -.5, 0., .5, 1., 2.}) {
    const std::uint64_t truth = query < -1 ? 0 : query < 0 ? 1
        : query < 1 ? maximum.count() + 1 : maximum.count() + 2;
    for (const auto* sketch : {&a, &b}) {
      const auto result = sketch->rank_with_error(query);
      assert(result.lower_rank <= truth && truth <= result.upper_rank);
      assert(result.estimate == static_cast<double>(truth));
      assert(sketch->max_rank_uncertainty() == 0);
    }
  }
}

static void extremes_and_large_counts() {
  CertifiedSplineSketch sketch(6);
  SplineSketch raw(6);
  std::vector<double> data{-std::numeric_limits<double>::max(), -1.5e308, -1e308,
                          1e308, 1.5e308, std::numeric_limits<double>::max()};
  for (double x : data) { sketch.add(x); raw.add(x); }
  sketch.consolidate(); raw.consolidate();
  verify(sketch, raw, data);
  // Three exact runs, self-merged to a count above 2^53. Check with integer
  // truth: using double subtraction for this assertion would hide rounding.
  CertifiedSplineSketch large(6);
  for (double value : {-1., 0., 1.}) large.add(value);
  std::uint64_t repetitions = 1;
  for (unsigned i = 0; i < 61; ++i) {
    large.merge(large);
    repetitions *= 2;
    assert(large.max_rank_uncertainty() == 0);
    for (int query = -2; query <= 2; ++query) {
      const std::uint64_t truth = query < -1 ? 0 : query < 0 ? repetitions
                                         : query < 1 ? 2 * repetitions : 3 * repetitions;
      auto result = large.rank_with_error(query);
      assert(result.lower_rank <= truth && truth <= result.upper_rank);
      assert(std::fabs(static_cast<long double>(result.estimate) - truth) <= result.max_error);
      assert(std::fabs(static_cast<long double>(result.estimate) - truth) <= large.max_rank_error());
    }
  }
  bool overflow = false;
  try { large.merge(large); large.merge(large); }
  catch (const std::overflow_error&) { overflow = true; }
  assert(overflow);
  // Reach UINT64_MAX without allocating a stream of that size.
  CertifiedSplineSketch maximum(6), power(6);
  power.add(0);
  for (unsigned bit = 0; bit < 64; ++bit) {
    maximum.merge(power);
    if (bit != 63) power.merge(power);
  }
  assert(maximum.count() == std::numeric_limits<std::uint64_t>::max());
  const auto endpoint = maximum.rank_with_error(INFINITY);
  assert(endpoint.lower_rank == maximum.count() && endpoint.upper_rank == maximum.count());
  assert(std::fabs(static_cast<long double>(endpoint.estimate) - maximum.count()) <= endpoint.max_error);
  overflow = false;
  try { maximum.add(0); } catch (const std::overflow_error&) { overflow = true; }
  assert(overflow && maximum.count() == std::numeric_limits<std::uint64_t>::max());
  assert(maximum.certificate_size() < 100);
  assert(maximum.max_rank_uncertainty() == 0);
  auto odd = power;
  odd.add(-1); odd.add(0); odd.add(1);
  const std::uint64_t exact_middle = power.count() + 2;
  const auto middle = odd.rank_with_error(0);
  assert(middle.lower_rank <= exact_middle && middle.upper_rank >= exact_middle);
  assert(std::fabs(static_cast<long double>(middle.estimate) - exact_middle) <= middle.max_error);
  assert(std::fabs(static_cast<long double>(middle.estimate) - exact_middle) <= odd.max_rank_error());
  CertifiedSplineSketch empty;
  assert(empty.rank_with_error(0).max_error == 0);
  bool invalid = false;
  try { empty.rank(NAN); } catch (const std::invalid_argument&) { invalid = true; }
  assert(invalid);
}

int main() {
  exhaustive_envelopes();
  local_compaction_equivalence();
  operations();
  packed_boundaries();
  extremes_and_large_counts();
  std::cout << "Certified integer rank bounds and noninferiority passed\n";
}
