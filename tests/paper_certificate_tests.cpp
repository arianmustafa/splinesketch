#include <splinesketch/paper_splinesketch.hpp>
#include "paper_resize_history_fixture.hpp"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <iostream>
#include <map>
#include <random>
#include <type_traits>

using Certified = splinesketch::CertifiedPaperSplineSketch;
using Paper = splinesketch::PaperSplineSketch;

namespace splinesketch {
struct PaperCertificateInspector {
  template<class Sketch>
  static void protected_retention(const Sketch& before, const Sketch& after) {
    assert(before.epoch_end_ == after.epoch_end_);
    for (const auto& old : before.nodes_) if (old.protected_threshold) {
      const auto found = std::find_if(after.nodes_.begin(), after.nodes_.end(),
          [&](const auto& node) { return node.x == old.x; });
      assert(found != after.nodes_.end() && found->protected_threshold);
    }
  }
  template<class Sketch>
  static void extrema_reserve(const Sketch& sketch) {
    if (sketch.policy_ != Sketch::BoundPolicy::practical) return;
    std::size_t protected_count = 0;
    for (const auto& node : sketch.nodes_)
      protected_count += static_cast<std::size_t>(node.protected_threshold);
    assert(protected_count <= sketch.capacity_ - 2);
  }
  static void check(const Certified& sketch, const std::map<double, std::uint64_t>& truth) {
    extrema_reserve(sketch);
    auto absorbed = truth;
    for (const auto& item : sketch.heavy_) absorbed[item.first] -= item.second.exact;
    for (double value : sketch.buffer_) --absorbed[value];
    std::uint64_t total = 0;
    for (auto item : absorbed) total += item.second;
    for (std::size_t i = 0; i < sketch.nodes_.size(); ++i) {
      const auto& node = sketch.nodes_[i];
      if (i) {
        assert(sketch.nodes_[i - 1].lower + node.atom <= node.lower);
        assert(sketch.nodes_[i - 1].upper + node.atom <= node.upper);
      }
      std::uint64_t exact = 0;
      for (auto item : absorbed) if (item.first <= node.x) exact += item.second;
      assert(node.lower <= exact && exact <= node.upper && node.upper <= total);
      const auto found = absorbed.find(node.x);
      const auto frequency = found == absorbed.end() ? 0 : found->second;
      assert(node.atom <= frequency);
    }
    if (!sketch.nodes_.empty()) {
      assert(sketch.nodes_.back().lower == total && sketch.nodes_.back().upper == total);
      assert(sketch.nodes_.front().lower == absorbed[sketch.nodes_.front().x]);
      assert(sketch.nodes_.front().upper == sketch.nodes_.front().lower);
    }
  }
  static void numeric_fallback() {
    Certified sketch(6);
    sketch.count_ = 6;
    sketch.nodes_ = {{0, 3}, {1, 3}};
    sketch.nodes_[0].lower = sketch.nodes_[0].upper = 3;
    sketch.nodes_[0].atom = 3;
    sketch.nodes_[1].lower = sketch.nodes_[1].upper = 6;
    sketch.nodes_[1].atom = 3;
    Certified::rebuild(sketch.nodes_);
    sketch.nodes_[0].prefix = std::numeric_limits<long double>::quiet_NaN();
    const auto certificate = sketch.rank_with_error(0.5);
    assert(certificate.estimate == 3 && certificate.max_error == 0);
  }
  static void huge_initialization() {
    Certified sketch(6);
    sketch.count_ = std::numeric_limits<std::uint64_t>::max();
    const auto weight = sketch.count_ / 7;
    Certified::Snapshot initial{{}, {}, {}};
    for (unsigned i = 0; i < 7; ++i) initial.items.emplace_back(i, weight);
    initial.items.back().second += sketch.count_ % 7;
    initial.prepare(); sketch.initialize(initial);
    std::map<double, std::uint64_t> truth;
    for (auto item : initial.items) truth[item.first] = item.second;
    check(sketch, truth);
    assert(sketch.rank_bounds(INFINITY).lower == sketch.count_);
  }
};
} // namespace splinesketch

static std::uint64_t checks = 0;

static void queries(const Certified& sketch, const std::map<double, std::uint64_t>& truth,
                    const Paper* raw = nullptr) {
  splinesketch::PaperCertificateInspector::check(sketch, truth);
  std::vector<double> xs{-INFINITY, INFINITY, -100, -0.5, 0, 0.5, 1, 3.125, 100};
  for (auto item : truth) {
    xs.push_back(item.first);
    xs.push_back(std::nextafter(item.first, -INFINITY));
    xs.push_back(std::nextafter(item.first, INFINITY));
  }
  const auto maximum = sketch.max_rank_error();
  for (double x : xs) {
    std::uint64_t exact = 0;
    for (auto item : truth) if (item.first <= x) exact += item.second;
    const auto interval = sketch.rank_bounds(x);
    const auto certificate = sketch.rank_with_error(x);
    assert(interval.lower <= exact && exact <= interval.upper);
    assert(interval.lower == certificate.lower_rank && interval.upper == certificate.upper_rank);
    assert(std::isfinite(certificate.estimate) && certificate.estimate == sketch.rank(x));
    const auto error = std::fabs(static_cast<long double>(certificate.estimate) - exact);
    assert(error <= certificate.max_error && error <= maximum);
    if (raw) {
      const auto original = raw->rank(x);
      assert(std::isfinite(original));
      assert(error <= std::fabs(static_cast<long double>(original) - exact));
    }
    ++checks;
  }
  if (sketch.count())
    for (double q : {0.0, 0.01, 0.1, 0.5, 0.9, 0.99, 1.0}) {
      const auto x = sketch.quantile(q);
      assert(std::isfinite(x));
      const auto target = std::ceil(q * static_cast<long double>(sketch.count()));
      if (q > 0) assert(sketch.rank(x) >= target);
      if (q > 0 && q < 1 && x > truth.begin()->first)
        assert(sketch.rank(std::nextafter(x, -INFINITY)) < target);
    }
}

static void exhaustive() {
  // Enumerate suffixes after six distinct observations, with an explicit
  // consolidation schedule that exercises MG cancellation and initialization.
  for (std::uint64_t code = 0; code < 4096; ++code) {
    Certified sketch(6);
    Paper raw(6);
    std::map<double, std::uint64_t> truth;
    auto rest = code;
    for (int i = 0; i < 12; ++i) {
      const double value = static_cast<double>(i < 6 ? static_cast<std::uint64_t>(i) : rest % 8);
      if (i >= 6) rest /= 8;
      sketch.add(value); raw.add(value); ++truth[value];
      if (i == 4 || i == 8) { sketch.consolidate(); raw.consolidate(); }
    }
    queries(sketch, truth, &raw);
    sketch.consolidate(); raw.consolidate(); queries(sketch, truth, &raw);
  }
}

static void operations() {
  std::mt19937_64 random(0x424f554e4453ULL);
  const double maximum = std::numeric_limits<double>::max();
  const double tiny = std::numeric_limits<double>::denorm_min();
  const double boundary = std::ldexp(1.0, -1021);
  const std::vector<double> edge{-maximum, -1e308, -1, -tiny, -0.0,
    tiny, boundary - tiny, boundary, boundary + 4 * tiny, 1, 1e308, maximum};
  for (auto policy : {Certified::BoundPolicy::practical, Certified::BoundPolicy::theoretical})
    for (std::size_t k : {6U, 32U, 128U, 512U, 1024U})
      for (int shape = 0; shape < 3; ++shape) {
        Certified sketch(k, policy), other(k + 3, policy);
        Paper raw(k, policy == Certified::BoundPolicy::practical
          ? Paper::BoundPolicy::practical : Paper::BoundPolicy::theoretical);
        std::map<double, std::uint64_t> truth, other_truth;
        for (int i = 0; i < 14000; ++i) {
          const double value = shape == 0 ? static_cast<double>(random() % 151) / 7 :
              shape == 1 ? edge[static_cast<std::size_t>(random() % edge.size())] :
              (i % 19 ? static_cast<double>(i % 271) / 271 : 1e100);
          sketch.add(value); raw.add(value); ++truth[value];
          if (i % 907 == 0) queries(sketch, truth, &raw);
          if (i % 3 == 0) { other.add(value); ++other_truth[value]; }
        }
        queries(sketch, truth, &raw);
        sketch.merge(other);
        for (auto item : other_truth) truth[item.first] += item.second;
        queries(sketch, truth);
        sketch.resize(std::max<std::size_t>(6, k / 2)); queries(sketch, truth);
        sketch.resize(k + 7); queries(sketch, truth);
        sketch.add(-maximum); ++truth[-maximum];
        sketch.consolidate(); queries(sketch, truth);
        sketch.finalize(); queries(sketch, truth);
        sketch.finalize(); queries(sketch, truth);
      }
}

static void reserve_transitions() {
  using Inspect = splinesketch::PaperCertificateInspector;
  for (unsigned seed : {1U, 8U, 19U}) {
    std::mt19937_64 random(seed);
    Certified sketch(7); Paper raw(7);
    std::map<double, std::uint64_t> truth;
    for (unsigned i = 0; i < 1024; ++i) {
      const double value = i % 13 ? std::ldexp(static_cast<double>(random() % 31), -20) :
          (i % 2 ? -1.0 : 1.0) * (i + 1);
      sketch.add(value); raw.add(value); ++truth[value];
      if (i % 7 == 0) { sketch.consolidate(); raw.consolidate(); }
      if (i % 11 == 0) {
        // Walk adjacent capacities in both directions, including small shrinks
        // that keep their barriers and shrinks that must reset protection.
        const auto at = (i / 11) % 20;
        const auto k = 6 + (at <= 10 ? at : 20 - at);
        sketch.resize(k); raw.resize(k);
      }
      if (i % 47 == 0) {
        // Different capacities and pending buffers exercise union reduction.
        Certified other(6 + (i % 17)); Paper plain(6 + (i % 17));
        for (double x : {static_cast<double>(i + 10000), -static_cast<double>(i + 10000), 0.0}) {
          other.add(x); plain.add(x); ++truth[x];
        }
        sketch.merge(other); raw.merge(plain);
      }
      Inspect::extrema_reserve(sketch); Inspect::extrema_reserve(raw);
      if (i % 127 == 0) queries(sketch, truth, &raw);
    }
    sketch.consolidate(); raw.consolidate(); queries(sketch, truth, &raw);
    sketch.finalize(); raw.finalize();
    Inspect::extrema_reserve(raw); queries(sketch, truth, &raw);
  }
}

static void clustered_extreme_queries() {
  const double tiny = std::numeric_limits<double>::denorm_min();
  for (std::size_t k : {16U, 128U, 1024U})
    for (auto policy : {Certified::BoundPolicy::practical, Certified::BoundPolicy::theoretical}) {
      Certified sketch(k, policy);
      Paper raw(k, policy == Certified::BoundPolicy::practical
          ? Paper::BoundPolicy::practical : Paper::BoundPolicy::theoretical);
      const auto add = [&](double value) {
        sketch.add(value); raw.add(value); sketch.consolidate(); raw.consolidate();
      };
      for (std::size_t i = 0; i < k; ++i) add(static_cast<double>(i));
      for (std::size_t i = 1; i <= 4 * k; ++i) add(static_cast<double>(i) * tiny);
      const auto truth = [&](double x) -> std::uint64_t {
        if (x < 0) return 0;
        const auto clustered = x >= 4 * k * tiny ? 4 * k : static_cast<std::size_t>(x / tiny);
        const auto roots = x >= k - 1 ? k : static_cast<std::size_t>(x) + 1;
        return clustered + roots;
      };
      std::vector<double> probes{-INFINITY, 0, tiny, 2 * k * tiny, 4 * k * tiny,
          (4 * k + 1) * tiny, std::nextafter(1.0, -INFINITY), 1, static_cast<double>(k - 1), INFINITY};
      for (double q : {0.0, 0.1, 0.5, 0.9, 1.0}) probes.push_back(sketch.quantile(q));
      for (double x : probes) {
        const auto answer = sketch.rank_with_error(x);
        const auto exact = truth(x);
        assert(answer.lower_rank <= exact && exact <= answer.upper_rank);
        assert(std::isfinite(answer.estimate) && std::isfinite(raw.rank(x)));
        const auto error = std::fabs(static_cast<long double>(answer.estimate) - exact);
        assert(error <= answer.max_error && error <= sketch.max_rank_error());
        ++checks;
      }
      assert(sketch.count() == 5 * k && sketch.quantile(0) == 0 && sketch.quantile(1) == k - 1);
      splinesketch::PaperCertificateInspector::extrema_reserve(sketch);
      splinesketch::PaperCertificateInspector::extrema_reserve(raw);
    }
}

static void capacity_growth_preserves_certificates() {
  // Growing a grid cannot reconstruct observations already forgotten by the
  // smaller sketch. Validate its inherited uncertainty against input truth.
  for (auto policy : {Certified::BoundPolicy::practical, Certified::BoundPolicy::theoretical}) {
    Certified source(6, policy);
    Paper raw(6, policy == Certified::BoundPolicy::practical
        ? Paper::BoundPolicy::practical : Paper::BoundPolicy::theoretical);
    std::map<double, std::uint64_t> truth;
    const auto add = [&](double value) {
      source.add(value); raw.add(value); ++truth[value];
    };
    for (unsigned i = 0; i < 6; ++i) add(i);
    source.consolidate(); raw.consolidate();
    for (unsigned i = 1; i <= 24; ++i) {
      add(i * std::numeric_limits<double>::denorm_min());
      source.consolidate(); raw.consolidate();
    }
    auto grown = source;
    auto plain_grown = raw;
    grown.resize(8192); plain_grown.resize(8192);
    assert(grown.count() == 30 && grown.bucket_capacity() == 8192);
    queries(grown, truth, &plain_grown);

    Certified large(8192, policy);
    Paper plain_large(8192, policy == Certified::BoundPolicy::practical
        ? Paper::BoundPolicy::practical : Paper::BoundPolicy::theoretical);
    for (unsigned i = 0; i < 31; ++i) { large.add(100); plain_large.add(100); }
    large.consolidate(); plain_large.consolidate();
    auto reverse = source;
    auto plain_reverse = raw;
    reverse.merge(large); plain_reverse.merge(plain_large);
    large.merge(source); plain_large.merge(raw);
    truth[100] += 31;
    assert(large.bucket_capacity() == 8192 && reverse.bucket_capacity() == 8192);
    assert(large.count() == 61 && reverse.count() == 61);
    queries(large, truth, &plain_large);
    queries(reverse, truth, &plain_reverse);
    splinesketch::PaperCertificateInspector::extrema_reserve(grown);
    splinesketch::PaperCertificateInspector::extrema_reserve(plain_grown);
    splinesketch::PaperCertificateInspector::extrema_reserve(large);
    splinesketch::PaperCertificateInspector::extrema_reserve(plain_large);
  }
}

static void repeated_capacity_changes() {
  for (auto policy : {Certified::BoundPolicy::practical, Certified::BoundPolicy::theoretical}) {
    Certified sketch(16, policy);
    Paper raw(16, policy == Certified::BoundPolicy::practical
        ? Paper::BoundPolicy::practical : Paper::BoundPolicy::theoretical);
    std::map<double, std::uint64_t> truth;
    for (double x : paper_test::resize_history_input) {
      sketch.add(x); raw.add(x); ++truth[x];
    }
    sketch.consolidate(); raw.consolidate();
    const auto held = sketch.heavy_hitter_count();
    const auto verify = [&] {
      assert(sketch.count() == paper_test::resize_history_input.size());
      assert(raw.count() == sketch.count() && sketch.heavy_hitter_count() == held);
      const auto answer = sketch.rank_with_error(paper_test::resize_history_query);
      std::uint64_t exact = 0;
      for (const auto& item : truth) if (item.first <= paper_test::resize_history_query) exact += item.second;
      assert(answer.lower_rank <= exact && exact <= answer.upper_rank);
      assert(std::fabs(static_cast<long double>(answer.estimate) - exact) <= answer.max_error);
      queries(sketch, truth, &raw);
    };
    verify();
    for (unsigned i = 0; i < 24; ++i) {
      const auto capacity = i % 2 ? 16 : 32;
      sketch.resize(capacity); raw.resize(capacity);
      assert(sketch.bucket_capacity() == static_cast<std::size_t>(capacity));
      verify();
    }
  }
}

static void exact_discrete_capacity_changes() {
  // Separate seeded data exposed a resize join score that improved a clustered
  // witness while discarding exact atoms here. Preserve the existing exact
  // certificates through repeated shrink/grow cycles on these 17 keys.
  for (auto policy : {Certified::BoundPolicy::practical, Certified::BoundPolicy::theoretical}) {
    Certified sketch(32, policy);
    Paper raw(32, policy == Certified::BoundPolicy::practical
        ? Paper::BoundPolicy::practical : Paper::BoundPolicy::theoretical);
    std::mt19937_64 random(97);
    std::map<double, std::uint64_t> truth;
    for (unsigned i = 0; i < 269; ++i) {
      const double value = static_cast<double>(random() % 17) - 8;
      sketch.add(value); raw.add(value); ++truth[value];
    }
    sketch.consolidate(); raw.consolidate();
    const auto verify = [&] {
      assert(sketch.count() == 269 && raw.count() == 269);
      assert(sketch.max_rank_uncertainty() == 0 && sketch.max_rank_error() == 0);
      queries(sketch, truth, &raw);
    };
    verify();
    for (unsigned i = 0; i < 24; ++i) {
      sketch.resize(i % 2 ? 32 : 16); raw.resize(i % 2 ? 32 : 16);
      verify();
    }
  }
}

static void large_counts() {
  Certified sketch(6), power(6);
  power.add(0);
  for (unsigned bit = 0; bit < 64; ++bit) {
    sketch.merge(power);
    if (bit != 63) power.merge(power);
  }
  assert(sketch.count() == std::numeric_limits<std::uint64_t>::max());
  assert(sketch.rank_bounds(0).lower == sketch.count());
  assert(sketch.rank_bounds(-1).upper == 0);
  const auto result = sketch.rank_with_error(0);
  // UINT64_MAX rounds to 2^64, exactly one rank unit above its true rank.
  assert(result.estimate == std::ldexp(1.0, 64));
  assert(result.max_error == 1 && sketch.max_rank_error() >= 1);
}

template<class Sketch>
static Sketch repeated_value(double value, std::uint64_t count,
                             typename Sketch::BoundPolicy policy) {
  assert(count > 0);
  Sketch sketch(6, policy);
  sketch.add(value); sketch.consolidate();
  unsigned bit = 0;
  for (auto rest = count; rest > 1; rest >>= 1) ++bit;
  while (bit) {
    --bit;
    sketch.merge(sketch);
    if ((count >> bit) & 1) { sketch.add(value); sketch.consolidate(); }
  }
  assert(sketch.count() == count);
  return sketch;
}

template<class Sketch>
static void initialization_extrema(std::uint64_t total,
                                   typename Sketch::BoundPolicy policy) {
  // Seven equally weighted keys and a singleton maximum. Build the counts
  // entirely through public operations without iterating over observations.
  // Merging cancels all eight held keys and initializes the spline. Rounded
  // last-quantile positions formerly omitted the singleton maximum.
  assert((total - 1) % 7 == 0);
  const auto weight = (total - 1) / 7;
  Sketch left(6, policy), right(6, policy);
  for (unsigned key = 0; key < 7; ++key) {
    auto item = repeated_value<Sketch>(key, weight, policy);
    (key < 4 ? left : right).merge(item);
  }
  right.add(7); right.consolidate(); left.merge(right);
  const auto check = [&] {
    assert(left.count() == total);
    assert(left.quantile(0) == 0);
    assert(left.quantile(1) == 7);
    if constexpr (std::is_same_v<Sketch, Certified>) {
      const auto at_six = left.rank_bounds(6);
      assert(at_six.lower <= total - 1 && total - 1 <= at_six.upper);
      for (double x : {7.0, static_cast<double>(INFINITY)}) {
        const auto bounds = left.rank_bounds(x);
        assert(bounds.lower == total && bounds.upper == total);
      }
      assert(left.rank_bounds(-1).upper == 0);
    }
  };
  check(); left.consolidate(); check(); left.finalize(); check();
}

static void huge_initialization_extrema() {
  // Counterexamples for native extended precision and binary64, respectively.
  for (auto total : {18446744073709551608ULL, 18446744073709550586ULL}) {
    for (auto policy : {Paper::BoundPolicy::practical, Paper::BoundPolicy::theoretical})
      initialization_extrema<Paper>(total, policy);
    for (auto policy : {Certified::BoundPolicy::practical, Certified::BoundPolicy::theoretical})
      initialization_extrema<Certified>(total, policy);
  }
}

template<class Sketch>
static void shrink_extrema_reserve(typename Sketch::BoundPolicy policy) {
  using Inspect = splinesketch::PaperCertificateInspector;
  Sketch sketch(7, policy);
  std::map<double, std::uint64_t> truth;
  const auto add = [&](double value) {
    sketch.add(value); ++truth[value]; sketch.consolidate();
    Inspect::extrema_reserve(sketch);
  };
  for (int i = 0; i < 7; ++i) add(i);
  for (int i = 0; i < 14; ++i) {
    add(std::ldexp(1.0, -10) + std::ldexp(i, -20));
    if (i == 6) {
      // At count 14, the protected cuts still fit capacity six. That small
      // shrink must retain their protection and their coordinates.
      auto smaller = sketch;
      smaller.resize(6); Inspect::protected_retention(sketch, smaller);
      Inspect::extrema_reserve(smaller);
      if constexpr (std::is_same_v<Sketch, Certified>) queries(smaller, truth);
    }
  }
  auto growing = sketch;
  growing.resize(8); Inspect::protected_retention(sketch, growing);
  Inspect::extrema_reserve(growing);
  if constexpr (std::is_same_v<Sketch, Certified>) queries(growing, truth);
  // This small shrink retained five protected nodes at capacity six. The
  // sixth fresh key then released both extrema before the epoch ending at 35,
  // leaving no legal second deletion. Resize must restore the reserve.
  sketch.resize(6); Inspect::extrema_reserve(sketch);
  for (double value : {10000, -10001, 10002, -10003, 10004, -10005}) add(value);
  assert(sketch.count() == 27 && sketch.bucket_capacity() == 6 && sketch.bucket_count() <= 6);
  assert(sketch.quantile(0) == -10005 && sketch.quantile(1) == 10004);
  if constexpr (std::is_same_v<Sketch, Certified>) queries(sketch, truth);
  else assert(sketch.rank(-INFINITY) == 0 && sketch.rank(INFINITY) == 27);
}

template<class Sketch>
static void protected_grid_updates(typename Sketch::BoundPolicy policy) {
  using Inspect = splinesketch::PaperCertificateInspector;
  Sketch sketch(6, policy);
  std::map<double, std::uint64_t> truth;
  const auto verify = [&] {
    Inspect::extrema_reserve(sketch);
    if constexpr (std::is_same_v<Sketch, Certified>) queries(sketch, truth);
    else {
      assert(sketch.rank(-INFINITY) == 0 && sketch.rank(INFINITY) == sketch.count());
      assert(sketch.quantile(0) == truth.begin()->first);
      assert(sketch.quantile(1) == truth.rbegin()->first);
      for (const auto& item : truth) {
        const auto rank = sketch.rank(item.first);
        assert(std::isfinite(rank) && rank >= 0 && rank <= sketch.count());
      }
    }
  };
  // These public updates formerly failed at observation 24, when a mandatory
  // split had exhausted removable cuts. Raising the practical bound avoids
  // that split while retaining every protected cut.
  for (double value : {1, 2, 0, 3, 4, 10000, 0, 0, 0,
                       5, 6, 8, 10, 16, 5, 6, 8, 10, 16, 5, 6, 8, 10, 16}) {
    sketch.add(value); ++truth[value];
    const auto before = sketch;
    sketch.consolidate(); Inspect::protected_retention(before, sketch); verify();
  }
  // Both new extrema need space before the current epoch ends at 30.
  // This formerly failed because earlier splits left only one removable cut.
  for (double value : {-5, -4, -3, -2, -1}) {
    sketch.add(value); ++truth[value];
    const auto before = sketch;
    sketch.consolidate(); Inspect::protected_retention(before, sketch); verify();
  }
  sketch.add(10001); ++truth[10001];
  const auto before = sketch;
  sketch.consolidate(); Inspect::protected_retention(before, sketch);
  assert(sketch.count() == 30 && sketch.bucket_count() <= 6);
  verify();
  // A later observation advances the epoch and releases the protection.
  // Consolidation must keep all 30 earlier observations.
  sketch.add(0); ++truth[0]; sketch.consolidate(); verify();
  assert(sketch.count() == 31);
}

int main() {
  splinesketch::PaperCertificateInspector::numeric_fallback();
  splinesketch::PaperCertificateInspector::huge_initialization();
  large_counts(); huge_initialization_extrema(); exhaustive(); operations(); reserve_transitions();
  clustered_extreme_queries();
  capacity_growth_preserves_certificates();
  repeated_capacity_changes();
  exact_discrete_capacity_changes();
  for (auto policy : {Paper::BoundPolicy::practical, Paper::BoundPolicy::theoretical}) {
    protected_grid_updates<Paper>(policy);
    shrink_extrema_reserve<Paper>(policy);
  }
  for (auto policy : {Certified::BoundPolicy::practical, Certified::BoundPolicy::theoretical}) {
    protected_grid_updates<Certified>(policy);
    shrink_extrema_reserve<Certified>(policy);
  }
  std::cout << "Certified paper: " << checks << " exact rank/certificate checks passed\n";
}
