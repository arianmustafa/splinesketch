#include <splinesketch/paper_splinesketch.hpp>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <iostream>
#include <map>
#include <random>

using Certified = splinesketch::CertifiedPaperSplineSketch;
using Paper = splinesketch::PaperSplineSketch;

namespace splinesketch {
struct PaperCertificateInspector {
  static void check(const Certified& sketch, const std::map<double, std::uint64_t>& truth) {
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

int main() {
  splinesketch::PaperCertificateInspector::numeric_fallback();
  splinesketch::PaperCertificateInspector::huge_initialization();
  large_counts(); exhaustive(); operations();
  std::cout << "Certified paper: " << checks << " exact rank/certificate checks passed\n";
}
