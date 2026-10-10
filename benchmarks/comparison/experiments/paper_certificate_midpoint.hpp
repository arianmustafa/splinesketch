#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
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

} // namespace splinesketch::experimental
