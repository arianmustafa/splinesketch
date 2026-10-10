#include <splinesketch/paper_splinesketch.hpp>
#include <splinesketch/certified_splinesketch.hpp>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cfenv>
#include <iostream>
#include <map>
#include <random>

namespace splinesketch {
struct PaperSplineSketchInspector {
  static void check(const PaperSplineSketch& s) {
    assert(s.nodes_.size() <= s.capacity_);
    assert(s.buffer_.size() < 5 * s.capacity_);
    assert(s.heavy_.size() < s.capacity_ || s.finalized_);
    long double sum = 0;
    for (std::size_t i = 0; i < s.nodes_.size(); ++i) {
      const auto& n = s.nodes_[i];
      assert(std::isfinite(n.x) && std::isfinite(n.mass) && n.mass >= 0);
      if (i) assert(s.nodes_[i - 1].x < n.x);
      sum += n.mass;
      assert(std::fabs(sum - n.prefix) <= 1e-9L);
    }
    for (const auto& h : s.heavy_) {
      assert(h.second.residual > 0 && h.second.residual <= h.second.exact);
      sum += h.second.exact;
    }
    sum += s.buffer_.size();
    assert(std::fabs(sum - static_cast<long double>(s.count_)) <= 1e-8L);
  }
  static void forwarded_bound(const PaperSplineSketch& s, const std::map<double, std::uint64_t>& truth,
                              std::size_t historical_min) {
    auto forwarded = truth;
    for (auto item : s.heavy_) forwarded[item.first] -= item.second.exact;
    for (double x : s.buffer_) --forwarded[x];
    for (auto item : forwarded) assert(item.second <= s.count_ / historical_min);
  }
  static void reject_illegal_join() {
    PaperSplineSketch s(6);
    s.count_ = 100;
    s.nodes_ = {{0, 1}, {1, 1}, {2, 1}};
    s.nodes_[1].protected_threshold = true;
    bool threw = false;
    try { s.join_at(1); } catch (const std::logic_error&) { threw = true; }
    assert(threw && s.nodes_.size() == 3);
    s.nodes_[1].protected_threshold = false;
    s.factor_ = 0.01L;
    threw = false;
    try { s.join_at(1); } catch (const std::logic_error&) { threw = true; }
    assert(threw && s.nodes_.size() == 3);
    // The heap backend must reject illegal joins even if heap membership is
    // stale: protection and mass guards are checked at the actual transition.
    s.factor_ = 3;
    s.nodes_[1].protected_threshold = true;
    PaperSplineSketch::rebuild(s.nodes_);
    const PaperSplineSketch::Snapshot before{s.nodes_, {}, {}};
    PaperSplineSketch::HeapRebalance heap(s, before, false);
    threw = false;
    try { heap.join(1); } catch (const std::logic_error&) { threw = true; }
    assert(threw && heap.size == 3 && heap.links[1].prev == 0 && heap.links[1].next == 2);
    s.nodes_[1].protected_threshold = false;
    s.factor_ = 0.01L;
    threw = false;
    try { heap.join(1); } catch (const std::logic_error&) { threw = true; }
    assert(threw && heap.size == 3 && heap.links[1].prev == 0 && heap.links[1].next == 2);
  }
  static void sparse_initialization() {
    PaperSplineSketch s(128);
    s.count_ = 20;
    PaperSplineSketch::Snapshot initial{{}, {{0, 10}, {1, 10}}, {}};
    initial.prepare();
    s.initialize(initial);
    assert(s.nodes_.size() == 128);
    for (auto n : s.nodes_) assert(!n.protected_threshold);
    // A new endpoint must be accommodated through legal joins, with no reset.
    s.count_ = 30;
    std::vector<PaperSplineSketch::Item> released{{2, 10}};
    s.incorporate(released);
    check(s);
    assert(s.rank(INFINITY) == 30);
  }
  static void storage_after_shrink(const PaperSplineSketch& s) {
    assert(s.nodes_.capacity() <= 2 * s.capacity_ + 2);
    assert(s.buffer_.capacity() <= 10 * s.capacity_);
  }
  static void exact_batch_limits() {
#ifdef __SIZEOF_INT128__
    PaperSplineSketch s(6);
    for (std::uint64_t n : {std::uint64_t{6000}, std::uint64_t{1} << 53,
                           (std::uint64_t{1} << 63) + 113, std::numeric_limits<std::uint64_t>::max()})
      for (std::size_t k : {std::size_t{6}, std::size_t{128}, std::size_t{100000000003ULL},
                           std::numeric_limits<std::size_t>::max() / 10})
        for (std::uint64_t factor : {std::uint64_t{3}, std::uint64_t{6}, std::uint64_t{1024},
                                    std::uint64_t{3} << 60}) {
          s.count_ = n; s.capacity_ = k; s.factor_ = static_cast<long double>(factor);
          const auto expected = static_cast<__uint128_t>(n) * factor / (2 * static_cast<__uint128_t>(k));
          const auto limit = std::max<std::uint64_t>(1, static_cast<std::uint64_t>(
              std::min(expected, static_cast<__uint128_t>(n))));
          assert(s.occurrence_batch_limit() == limit);
        }
#endif
  }
};
} // namespace splinesketch

using Sketch = splinesketch::PaperSplineSketch;
using Inspect = splinesketch::PaperSplineSketchInspector;

static std::uint64_t reference_target(double q, std::uint64_t count) {
#ifdef __SIZEOF_INT128__
  // Independent oracle: decode the binary64 bits and multiply directly in
  // 128 bits, rather than using frexp and the production word products.
  std::uint64_t bits;
  std::memcpy(&bits, &q, sizeof bits);
  const auto exponent = static_cast<unsigned>(bits >> 52);
  const auto significand = (bits & ((std::uint64_t{1} << 52) - 1)) |
      (exponent ? std::uint64_t{1} << 52 : 0);
  const auto product = static_cast<__uint128_t>(count) * significand;
  const auto shift = exponent ? 1075 - exponent : 1074;
  if (shift >= 128) return product != 0;
  return static_cast<std::uint64_t>((product >> shift) +
      ((product & ((static_cast<__uint128_t>(1) << shift) - 1)) != 0));
#else
  return splinesketch::detail::quantile_rank_target(q, count);
#endif
}

static void exact_quantile_arithmetic() {
  using splinesketch::detail::quantile_rank_target;
  using splinesketch::detail::rank_reaches_target;
  const auto maximum = std::numeric_limits<std::uint64_t>::max();
  std::vector<std::pair<double, std::uint64_t>> cases;
  for (std::uint64_t count : {std::uint64_t{0}, std::uint64_t{1}, std::uint64_t{5},
                            std::uint64_t{10}, maximum}) {
    for (double q : {0., 0.1, 0.5, std::nextafter(1., 0.), 1.,
                    std::numeric_limits<double>::denorm_min()}) cases.emplace_back(q, count);
  }
  for (unsigned bit = 0; bit < 64; ++bit) {
    const auto power = std::uint64_t{1} << bit;
    for (auto count : {power - 1, power, power + 1})
      for (double q : {0.1, 0.2, std::nextafter(0.5, 0.), 0.5, std::nextafter(0.5, 1.),
                      std::nextafter(1., 0.)}) cases.emplace_back(q, count);
  }
  // Every normal denominator exponent and nearby subnormal boundaries.
  for (unsigned exponent = 0; exponent < 1023; ++exponent) {
    for (std::uint64_t fraction : {std::uint64_t{0}, std::uint64_t{1},
                                  (std::uint64_t{1} << 52) - 1}) {
      const auto bits = (static_cast<std::uint64_t>(exponent) << 52) | fraction;
      double q;
      std::memcpy(&q, &bits, sizeof q);
      cases.emplace_back(q, maximum);
    }
  }
  std::mt19937_64 random(0x5155414e54494c45ULL);
  for (unsigned i = 0; i < 10000; ++i) {
    const auto bits = random() % 0x3ff0000000000001ULL;
    double q;
    std::memcpy(&q, &bits, sizeof q);
    cases.emplace_back(q, random());
  }
  const auto old_mode = std::fegetround();
  for (auto mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO}) {
    assert(std::fesetround(mode) == 0);
    for (auto input : cases)
      assert(quantile_rank_target(input.first, input.second) == reference_target(input.first, input.second));
    assert(quantile_rank_target(0.1, 10) == 2); // Exact binary64 0.1 exceeds 1/10.
    assert(quantile_rank_target(0.5, (std::uint64_t{1} << 53) + 1) == (std::uint64_t{1} << 52) + 1);
    assert(quantile_rank_target(std::numeric_limits<double>::denorm_min(), maximum) == 1);
    assert(quantile_rank_target(1, maximum) == maximum);
    assert(!rank_reaches_target(0x1p53, (std::uint64_t{1} << 53) + 1));
    assert(rank_reaches_target(0x1p53 + 2, (std::uint64_t{1} << 53) + 1));
    assert(rank_reaches_target(0x1p64, maximum));
    assert(!rank_reaches_target(std::nextafter(0x1p64, 0.), maximum));
    assert(rank_reaches_target(0.5, 0) && !rank_reaches_target(0.5, 1));
  }
  assert(std::fesetround(old_mode) == 0);
}

template<class S>
static void large_count_median(S zeros, S ones, unsigned exponent,
                               unsigned extra_zeros, unsigned extra_ones, double expected,
                               bool prepare_exact = false) {
  zeros.add(0); ones.add(1);
  for (unsigned i = 0; i < exponent; ++i) {
    zeros.merge(zeros); ones.merge(ones);
  }
  for (unsigned i = 0; i < extra_zeros; ++i) zeros.add(0);
  for (unsigned i = 0; i < extra_ones; ++i) ones.add(1);
  zeros.merge(ones);
  assert(zeros.count() == (std::uint64_t{1} << (exponent + 1)) + extra_zeros + extra_ones);
  assert(zeros.quantile(0.5) == expected);
  if (prepare_exact) {
    // More than 16 distinct exact values exercises the prepared-prefix
    // search as well. Equal additions on either side preserve the median.
    for (unsigned i = 0; i < 16; ++i) {
      zeros.add(-1. - i); zeros.add(2. + i);
    }
    assert(zeros.quantile(0.5) == expected);
  }
  zeros.consolidate();
  assert(zeros.quantile(0.5) == expected);
  if constexpr (std::is_same_v<S, splinesketch::PaperSplineSketch> ||
                std::is_same_v<S, splinesketch::CertifiedPaperSplineSketch>) {
    zeros.finalize();
    assert(zeros.quantile(0.5) == expected);
  }
}

static void large_count_quantiles() {
  using Core = splinesketch::SplineSketch;
  using Certified = splinesketch::CertifiedSplineSketch;
  using Paper = splinesketch::CertifiedPaperSplineSketch;
  // ceil(n/2) = 2^52+1. Converting the odd count to binary64 first
  // incorrectly selects the last zero instead of the first one.
  large_count_median(Core(8), Core(8), 52, 0, 1, 1);
  large_count_median(Certified(8), Certified(8), 52, 0, 1, 1);
  large_count_median(Core(64), Core(64), 52, 0, 1, 1, true);
  large_count_median(Certified(64), Certified(64), 52, 0, 1, 1, true);
  for (auto policy : {Sketch::BoundPolicy::practical, Sketch::BoundPolicy::theoretical}) {
    large_count_median(Sketch(8, policy), Sketch(8, policy), 52, 0, 1, 1);
    large_count_median(Sketch(64, policy), Sketch(64, policy), 52, 0, 1, 1, true);
    // The target itself is now unrepresentable as binary64. Its comparison
    // must retain the low integer bit after computing the correct target.
    large_count_median(Sketch(8, policy), Sketch(8, policy), 53, 0, 1, 1);
  }
  for (auto policy : {Paper::BoundPolicy::practical, Paper::BoundPolicy::theoretical}) {
    large_count_median(Paper(8, policy), Paper(8, policy), 52, 0, 1, 1);
    large_count_median(Paper(64, policy), Paper(64, policy), 52, 0, 1, 1, true);
    large_count_median(Paper(8, policy), Paper(8, policy), 53, 0, 1, 1);
    // A singleton integer certificate settles the crossing even when its
    // displayed rank rounds down, or rounds up to a target not yet reached.
    large_count_median(Paper(8, policy), Paper(8, policy), 53, 1, 0, 0);
    large_count_median(Paper(8, policy), Paper(8, policy), 53, 3, 5, 1);
  }
}

static void query_checks(const Sketch& s, const std::vector<double>& data) {
  Inspect::check(s);
  std::vector<double> queries = data;
  for (double x : data) {
    queries.push_back(std::nextafter(x, -INFINITY));
    queries.push_back(std::nextafter(x, INFINITY));
  }
  queries.push_back(-INFINITY);
  queries.push_back(INFINITY);
  std::sort(queries.begin(), queries.end());
  queries.erase(std::unique(queries.begin(), queries.end()), queries.end());
  double previous = 0;
  for (double x : queries) {
    const auto rank = s.rank(x);
    assert(std::isfinite(rank) && rank >= 0 && rank <= static_cast<double>(s.count()));
    assert(rank + 1e-9 >= previous);
    previous = rank;
  }
  for (double q : {0.0, 0.001, 0.2, 0.5, 0.9, 0.999, 1.0}) {
    const auto x = s.quantile(q);
    const auto target = reference_target(q, s.count());
    assert(std::isfinite(x));
    if (q > 0) assert(s.rank(x) >= target);
    if (q > 0 && q < 1 && x > *std::min_element(data.begin(), data.end()))
      assert(s.rank(std::nextafter(x, -INFINITY)) < target);
  }
}

static void streaming_and_merges() {
  std::mt19937_64 random(0x5041504552ULL);
  for (auto policy : {Sketch::BoundPolicy::practical, Sketch::BoundPolicy::theoretical})
  for (std::size_t k : {6U, 8U, 16U, 64U, 128U, 1024U}) {
    for (int shape = 0; shape < 4; ++shape) {
      Sketch left(k, policy), right(k, policy);
      std::vector<double> data;
      std::map<double, std::uint64_t> left_truth, right_truth;
      for (int i = 0; i < 12000; ++i) {
        const double x = shape == 0 ? static_cast<double>(random() % 41) :
            shape == 1 ? static_cast<double>(i) :
            shape == 2 ? (i % 29 ? static_cast<double>(random() % 1000) / 1000 : 1e100) :
            (i % 2 ? std::nextafter(1.0, 2.0) : 1.0);
        auto& s = i % 2 ? left : right;
        auto& truth = i % 2 ? left_truth : right_truth;
        s.add(x);
        ++truth[x];
        data.push_back(x);
        if (i % 113 == 0) {
          Inspect::check(s);
          Inspect::forwarded_bound(s, truth, k);
        }
      }
      // Merge unflushed buffers, then verify exact mass and MG filtering.
      left.merge(right);
      for (auto item : right_truth) left_truth[item.first] += item.second;
      Inspect::forwarded_bound(left, left_truth, k);
      query_checks(left, data);
      left.resize(std::max<std::size_t>(6, k / 2));
      Inspect::forwarded_bound(left, left_truth, std::max<std::size_t>(6, k / 2));
      query_checks(left, data);
      left.resize(k + 3);
      query_checks(left, data);
      left.finalize();
      query_checks(left, data);
      assert(left.finalized());
      bool threw = false;
      try { left.add(0); } catch (const std::logic_error&) { threw = true; }
      assert(threw);
    }
  }
}

static void edge_cases() {
  Sketch s(8);
  bool threw = false;
  try { s.quantile(0.5); } catch (const std::logic_error&) { threw = true; }
  assert(threw);
  std::vector<double> data{-std::numeric_limits<double>::max(), -1.5e308, -1e308,
                          0, std::numeric_limits<double>::denorm_min(), 1e308,
                          1.5e308, std::numeric_limits<double>::max()};
  for (int i = 0; i < 100; ++i) for (double x : data) s.add(x);
  s.consolidate();
  query_checks(s, data);
  s.merge(s);
  assert(s.count() == 1600);
  query_checks(s, data);
  Sketch a(16), b(8);
  for (int i = 0; i < 10; ++i) a.add(i);
  for (int i = 0; i < 20; ++i) b.add(i);
  a.merge(b);
  assert(a.bucket_capacity() == 8 && a.count() == 30);
  threw = false;
  try { a.add(INFINITY); } catch (const std::invalid_argument&) { threw = true; }
  assert(threw && a.count() == 30);
  Sketch growing(6);
  for (int i = 0; i < 1200; ++i) growing.add(i);
  growing.consolidate();
  growing.resize(128);
  growing.add(1e100);
  growing.consolidate();
  Inspect::check(growing);
  growing.resize(6);
  Inspect::storage_after_shrink(growing);
  Sketch maximum(6), power(6);
  power.add(0);
  for (unsigned bit = 0; bit < 64; ++bit) {
    maximum.merge(power);
    if (bit != 63) power.merge(power);
  }
  assert(maximum.count() == std::numeric_limits<std::uint64_t>::max());
  threw = false;
  try { maximum.add(0); } catch (const std::overflow_error&) { threw = true; }
  assert(threw && maximum.count() == std::numeric_limits<std::uint64_t>::max());
  threw = false;
  try { maximum.merge(power); } catch (const std::overflow_error&) { threw = true; }
  assert(threw && maximum.count() == std::numeric_limits<std::uint64_t>::max());
}

int main() {
  exact_quantile_arithmetic();
  large_count_quantiles();
  Inspect::reject_illegal_join();
  Inspect::sparse_initialization();
  Inspect::exact_batch_limits();
  edge_cases();
  streaming_and_merges();
  std::cout << "Paper baseline mass, MG, monotonicity, inverse, merge, resize and finalization checks passed\n";
}
