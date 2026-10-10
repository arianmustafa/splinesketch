#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace splinesketch::experimental {

// Query-only experiment: consume the integer certificate without changing the
// grid, update rules, or production spline query policy.
struct CertificateMidpointEstimate {
  double estimate;
  std::uint64_t lower_rank, upper_rank;
  double max_error;
};

namespace midpoint_detail {
static_assert(std::numeric_limits<double>::is_iec559 &&
              std::numeric_limits<double>::radix == 2 &&
              std::numeric_limits<double>::digits == 53,
              "certificate midpoint requires IEEE binary64");
static_assert(sizeof(double) == sizeof(std::uint64_t), "binary64 storage required");

inline unsigned integer_shift(std::uint64_t value) noexcept {
  // At most 11 discarded bits for uint64_t. The remaining significand is
  // below 2^53, so incrementing it once still converts to double exactly.
  unsigned shift = 0;
  while ((value >> shift) >= (std::uint64_t{1} << 53)) ++shift;
  return shift;
}

// Least binary64 value >= the integer. The significand and scaling operations
// are exact, including the case rounding UINT64_MAX upward produces 2^64.
inline double integer_ceiling(std::uint64_t value) noexcept {
  const auto shift = integer_shift(value);
  const auto mask = (std::uint64_t{1} << shift) - 1;
  const auto significand = (value >> shift) + ((value & mask) != 0);
  return std::ldexp(static_cast<double>(significand), static_cast<int>(shift));
}

inline double distance(double estimate, std::uint64_t endpoint) noexcept {
  if (estimate < 0x1p52) {
    // For our midpoint estimates, both endpoints and both endpoint distances
    // are exact here: the midpoint and distances are half integers <= 2^52.
    return std::fabs(estimate - static_cast<double>(endpoint));
  }
  if (estimate == 0x1p64) {
    if (!endpoint) return estimate;
    return integer_ceiling(std::numeric_limits<std::uint64_t>::max() - endpoint + 1);
  }
  // Every binary64 value in [2^52,2^64) is an integer; the exceptional upper
  // endpoint above must be handled before the uint64_t conversion.
  const auto integer = static_cast<std::uint64_t>(estimate);
  return integer_ceiling(integer >= endpoint ? integer - endpoint : endpoint - integer);
}

inline std::uint64_t ordered_bits(double value) noexcept {
  std::uint64_t bits;
  std::memcpy(&bits, &value, sizeof bits);
  constexpr auto sign = std::uint64_t{1} << 63;
  return bits & sign ? ~bits : bits ^ sign;
}

inline double from_ordered_bits(std::uint64_t ordered) noexcept {
  constexpr auto sign = std::uint64_t{1} << 63;
  const auto bits = ordered & sign ? ordered ^ sign : ~ordered;
  double value;
  std::memcpy(&value, &bits, sizeof value);
  return value;
}

inline bool reaches(std::uint64_t lower, std::uint64_t upper, std::uint64_t target) noexcept {
  // Equivalent to (lower+upper)/2 >= target, without overflowing either sum
  // or converting integer ranks to floating point. Bounds must be ordered.
  if (lower >= target) return true;
  return upper >= target && target - lower <= upper - target;
}
} // namespace midpoint_detail

inline CertificateMidpointEstimate certificate_midpoint(std::uint64_t lower, std::uint64_t upper) {
  if (lower > upper) throw std::invalid_argument("reversed rank certificate");
  const auto width = upper - lower;
  const auto floor_midpoint = lower + width / 2; // No overflowing L+U.
  const auto half = width & 1;
  double estimate;
  if (floor_midpoint < (std::uint64_t{1} << 52)) {
    estimate = static_cast<double>(floor_midpoint) + (half ? 0.5 : 0.0);
  } else {
    const auto shift = midpoint_detail::integer_shift(floor_midpoint);
    const auto spacing = std::uint64_t{1} << shift;
    auto significand = floor_midpoint >> shift;
    // Compare the exact discarded remainder with half a spacing. Keep the
    // half-integer bit: converting the floor first can double-round odd widths.
    const auto twice_remainder = 2 * (floor_midpoint & (spacing - 1)) + half;
    if (twice_remainder > spacing || (twice_remainder == spacing && (significand & 1)))
      ++significand;
    estimate = std::ldexp(static_cast<double>(significand), static_cast<int>(shift));
  }
  const auto radius = std::max(midpoint_detail::distance(estimate, lower),
                               midpoint_detail::distance(estimate, upper));
  return {estimate, lower, upper, radius};
}

// Proof for any real rank t in [L,U]:
//   max_t |r-t| = max(|r-L|,|r-U|) = (U-L)/2 + |r-(L+U)/2|.
// Therefore a nearest binary64 midpoint minimizes this radius among all
// binary64 estimates; either nearest value is optimal in a tie. The integer
// rounding above chooses ties to even without relying on long double or the
// floating-point rounding mode. The returned radius is rounded upward exactly.
// This is minimax over the certificate endpoints, not over a smaller feasible
// rank set inferred from other sketch state, and not a pointwise accuracy claim.
//
// Monotonicity for CertifiedPaperSplineSketch: its ordered grid satisfies
// L_i<=U_i, a_i>=0, L_i-L_(i-1)>=a_i, U_i-U_(i-1)>=a_i. At each interior
// cut the envelope passes through gap [L_(i-1),U_i-a_i], knot [L_i,U_i], and
// next gap [L_i,U_(i+1)-a_(i+1)] with neither endpoint decreasing. Below the
// first cut it is [0,0]; at and above the last cut it is [U_last,U_last]. The
// first knot also respects this order. Pending exact counts add the same
// nondecreasing step function to both endpoints. Initialization from an
// exact cumulative count, sampling with the source's exact atoms, adding
// cumulative counts, summing merge sources, and deleting cuts preserve these
// inequalities. Splits sample the frozen envelope, so induction covers them
// too. Consequently M(x)=(L(x)+U(x))/2 is nondecreasing. Fixed nearest-even
// binary64 rounding is nondecreasing (its ordered rounding cells cannot
// reverse order); query-dependent choices between equally close tied values
// would not suffice. Hence the returned midpoint is nondecreasing for all
// ordered non-NaN queries, including adjacent keys and counts above 2^53.
template<class CertifiedSketch>
CertificateMidpointEstimate certificate_midpoint_rank(const CertifiedSketch& sketch, double query) {
  const auto bounds = sketch.rank_bounds(query);
  return certificate_midpoint(bounds.lower, bounds.upper);
}

inline double certificate_midpoint_uniform_error(std::uint64_t width, std::uint64_t count) {
  if (width > count) throw std::invalid_argument("certificate width exceeds count");
  // Every midpoint is between 0 and count. Up to 2^52 all half integers are
  // representable. Above it, half the largest spacing reached by count bounds
  // nearest rounding error, even if rounding crosses a binade boundary.
  // The bound is W/2 plus that allowance, rounded outward for binary64 returns.
  if (count <= (std::uint64_t{1} << 52)) return static_cast<double>(width) / 2;
  const double allowance = std::ldexp(1.0,
      static_cast<int>(midpoint_detail::integer_shift(count)) - 1);
  // Integer ceiling and division by two are outward/exact. One step upward
  // covers rounding of the addition under any IEEE rounding mode.
  return std::nextafter(midpoint_detail::integer_ceiling(width) / 2 + allowance, INFINITY);
}

template<class CertifiedSketch>
double certificate_midpoint_max_error(const CertifiedSketch& sketch) {
  return certificate_midpoint_uniform_error(sketch.max_rank_uncertainty(), sketch.count());
}

// Integer-rank quantile experiment. Invert the exact midpoint M=(L+U)/2,
// rather than its rounded display value: for some counts above 2^53 the
// displayed terminal rank rounds below count, making target=count unreachable.
// Targets are 1..count, with no floating-point q*count conversion. Return the
// first finite binary64 value in numerical order with M(x)>=target; zero is
// canonicalized to +0. This need not be an observed value.
//
// Proof for the certified paper envelopes: reaches(x) is monotone. At DBL_MAX the
// bounds equal count, so the high endpoint is feasible. Binary search preserves
// the first feasible point in [low,high] and terminates in at most 64 steps.
// Signed zeros have the same bounds. For p=nextafter(x,-infinity), M(p)<target
// and M(x)>=target, including x=-DBL_MAX, whose predecessor is -infinity.
// If W is max_rank_uncertainty and R is the true inclusive rank, the integer
// certificates give R(x)>=target-W/2 and R(p)<target+W/2. Thus the distance of
// target from the closed rank bracket [R(p),R(x)] is at most floor(W/2), since
// this distance is an integer. Duplicate atoms can make |R(x)-target| as
// large as count-target, even when W=0. No binary64 rank-return allowance is
// needed for this integer predicate.
template<class CertifiedSketch>
double certificate_midpoint_select(const CertifiedSketch& sketch, std::uint64_t target) {
  if (!sketch.count()) throw std::logic_error("selection from empty sketch");
  if (!target || target > sketch.count()) throw std::invalid_argument("rank target must be in [1,count]");
  auto low = midpoint_detail::ordered_bits(-std::numeric_limits<double>::max());
  auto high = midpoint_detail::ordered_bits(std::numeric_limits<double>::max());
  while (low < high) {
    const auto middle = low + (high - low) / 2;
    const auto bounds = sketch.rank_bounds(midpoint_detail::from_ordered_bits(middle));
    if (midpoint_detail::reaches(bounds.lower, bounds.upper, target)) high = middle;
    else low = middle + 1;
  }
  const auto value = midpoint_detail::from_ordered_bits(low);
  return value == 0 ? 0.0 : value;
}

} // namespace splinesketch::experimental
