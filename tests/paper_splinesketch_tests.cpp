#include <splinesketch/paper_splinesketch.hpp>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
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
    assert(std::isfinite(x));
    if (q > 0) assert(s.rank(x) >= std::ceil(q * static_cast<long double>(s.count())));
    if (q > 0 && q < 1 && x > *std::min_element(data.begin(), data.end()))
      assert(s.rank(std::nextafter(x, -INFINITY)) < std::ceil(q * static_cast<long double>(s.count())));
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
  Inspect::reject_illegal_join();
  Inspect::sparse_initialization();
  Inspect::exact_batch_limits();
  edge_cases();
  streaming_and_merges();
  std::cout << "Paper baseline mass, MG, monotonicity, inverse, merge, resize and finalization checks passed\n";
}
