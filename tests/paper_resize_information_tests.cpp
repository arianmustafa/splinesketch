#include <splinesketch/paper_splinesketch.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <vector>

using Certified = splinesketch::CertifiedPaperSplineSketch;
using Paper = splinesketch::PaperSplineSketch;

static void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

namespace splinesketch {
struct PaperCertificateInspector {
  template<class S> static void same_state(const S& a, const S& b) {
    // Compare every active value used by subsequent queries and transitions,
    // rather than a hash or just the returned ranks. Allocation addresses,
    // padding and inactive capacity bytes are not part of this logical state;
    // container capacities are also equal for this witness.
#ifdef SPLINESKETCH_RESIZE_REGRID_ONLY
    require(!a.rank_resize_ && !b.rank_resize_, "resize flag escaped the operation");
#endif
    require(a.capacity_ == b.capacity_ && a.count_ == b.count_ && a.epoch_end_ == b.epoch_end_ &&
            a.policy_ == b.policy_ && a.factor_ == b.factor_ && a.finalized_ == b.finalized_,
            "scalar state differs");
    require(a.buffer_.empty() && b.buffer_.empty() && a.heavy_.empty() && b.heavy_.empty(),
            "witness contains retained observations outside the grid");
    require(a.buffer_.capacity() == b.buffer_.capacity() && a.nodes_.capacity() == b.nodes_.capacity() &&
            a.heavy_.bucket_count() == b.heavy_.bucket_count() &&
            a.heavy_.max_load_factor() == b.heavy_.max_load_factor(), "container state differs");
    require(a.nodes_.size() == b.nodes_.size(), "grid size differs");
    for (std::size_t i = 0; i < a.nodes_.size(); ++i) {
      const auto& x = a.nodes_[i]; const auto& y = b.nodes_[i];
      require(x.x == y.x && x.mass == y.mass && x.prefix == y.prefix && x.slope == y.slope &&
              x.protected_threshold == y.protected_threshold, "curve state differs");
      if constexpr (std::is_same_v<S, Certified>)
        require(x.lower == y.lower && x.upper == y.upper && x.atom == y.atom, "integer state differs");
    }
  }
};
} // namespace splinesketch

static std::vector<double> input(bool high_rank) {
  std::vector<double> data;
  for (unsigned i = 0; i < 30; ++i) data.push_back(i * i);
  for (unsigned i = 1; i <= 5; ++i) data[i] = (high_rank ? 0 : 18) + 3 * i;
  return data;
}

static std::uint64_t exact_rank(const std::vector<double>& data, double query) {
  return static_cast<std::uint64_t>(std::count_if(data.begin(), data.end(),
      [query](double x) { return x <= query; }));
}

static std::uint64_t bits(double value) {
  std::uint64_t result;
  std::memcpy(&result, &value, sizeof result);
  return result;
}

template<class S> static S replay(const std::vector<double>& data, bool theoretical) {
  S result(6, theoretical ? S::BoundPolicy::theoretical : S::BoundPolicy::practical);
  for (double x : data) result.add(x);
  result.consolidate();
  return result;
}

template<class S> static void emit(bool theoretical) {
  constexpr double query = 18;
  const auto high_data = input(true), low_data = input(false);
  const auto high_truth = exact_rank(high_data, query), low_truth = exact_rank(low_data, query);
  require(high_truth == 6 && low_truth == 1, "unexpected independent ranks");
  auto high = replay<S>(high_data, theoretical), low = replay<S>(low_data, theoretical);
  require(high.count() == 30 && low.count() == 30, "incorrect initial count");
  splinesketch::PaperCertificateInspector::same_state(high, low);
  const auto initial = high;
  const auto old = high.rank(query);
  require(old == low.rank(query) && low_truth < old && old < high_truth, "invalid common estimate");

  // The family reaches both endpoints of the old certificate with identical
  // processing state. Any deterministic data-free replacement b must be the
  // same for both streams. If b > old it hurts rank 1; if b < old it hurts
  // rank 6. Keeping b == old is allowed: do not require resize distortion.
  for (unsigned capacity : {12, 6}) {
    high.resize(capacity); low.resize(capacity);
    require(high.bucket_capacity() == capacity && low.bucket_capacity() == capacity, "resize failed");
    require(high.count() == 30 && low.count() == 30, "data-free resize changed count");
    splinesketch::PaperCertificateInspector::same_state(high, low);
    require(high.rank(query) == low.rank(query), "identical states returned different ranks");
  }
  const auto changed = high.rank(query);
  require(std::isfinite(changed), "nonfinite changed estimate");
  if (changed > old)
    require(std::fabs(changed - low_truth) > std::fabs(old - low_truth), "upward change did not hurt low rank");
  if (changed < old)
    require(std::fabs(changed - high_truth) > std::fabs(old - high_truth), "downward change did not hurt high rank");

  std::cout << "{\"certified\":" << std::is_same_v<S, Certified> << ",\"theoretical\":" << theoretical
            << ",\"query_bits\":" << bits(query) << ",\"old_bits\":" << bits(old)
            << ",\"new_bits\":" << bits(changed) << ",\"low_truth\":" << low_truth
            << ",\"high_truth\":" << high_truth;
  if constexpr (std::is_same_v<S, Certified>) {
    const auto before = initial.rank_with_error(query), after = high.rank_with_error(query);
    require(before.lower_rank == low_truth && before.upper_rank == high_truth, "endpoints not attained");
    require(after.lower_rank <= low_truth && high_truth <= after.upper_rank, "certificate excludes a stream");
    std::cout << ",\"old_lower\":" << before.lower_rank << ",\"old_upper\":" << before.upper_rank
              << ",\"new_lower\":" << after.lower_rank << ",\"new_upper\":" << after.upper_rank
              << ",\"old_error_bits\":" << bits(before.max_error)
              << ",\"new_error_bits\":" << bits(after.max_error)
              << ",\"old_uniform_bits\":" << bits(initial.max_rank_error())
              << ",\"new_uniform_bits\":" << bits(high.max_rank_error());
  }
  std::cout << "}\n";
}

int main() {
  for (bool theoretical : {false, true}) { emit<Paper>(theoretical); emit<Certified>(theoretical); }
}
