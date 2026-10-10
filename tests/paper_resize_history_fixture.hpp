#pragma once
#include <array>

namespace paper_test {
// A reduced public stream for capacity changes with no new observations.
inline constexpr double unit = 0x1p-20;
inline constexpr double resize_history_query = 0x1.fffffffffffffp-21;
inline constexpr std::array<double, 124> resize_history_input{{
  29 * unit, 12 * unit, 5 * unit, 1, 9 * unit, 26 * unit, 24 * unit, 1 * unit,
  1 * unit, 7 * unit, 22 * unit, 19 * unit, 26 * unit, 21 * unit, 2 * unit, 18 * unit,
  17 * unit, 3 * unit, 8 * unit, -456, 0, 20 * unit, 23 * unit, 29 * unit,
  6 * unit, 12 * unit, -196, 17 * unit, 2 * unit, 15 * unit, 14 * unit, 2 * unit,
  18 * unit, 25 * unit, 3 * unit, 13 * unit, 12 * unit, 2 * unit, 495, 3 * unit,
  25 * unit, 30 * unit, -352, 30 * unit, 30 * unit, 19 * unit, 4 * unit, 1 * unit,
  12 * unit, 10 * unit, 24 * unit, -326, 8 * unit, 12 * unit, 10 * unit, 21 * unit,
  20 * unit, 23 * unit, 28 * unit, 0, 30 * unit, 5 * unit, 365, 26 * unit,
  2 * unit, 18 * unit, 3 * unit, 17 * unit, 6 * unit, 24 * unit, 5 * unit, 11 * unit,
  53, 235, 1 * unit, 12 * unit, 28 * unit, 18 * unit, 9 * unit, 8 * unit,
  17 * unit, 29 * unit, 25 * unit, 18 * unit, 17 * unit, 26 * unit, 23 * unit, 12 * unit,
  16 * unit, 26 * unit, 24 * unit, 27 * unit, 14 * unit, 10 * unit, 1 * unit, 14 * unit,
  22 * unit, 19 * unit, 3 * unit, 27 * unit, 17 * unit, 20 * unit, -118, 22 * unit,
  16 * unit, 17 * unit, 11 * unit, 21 * unit, 21 * unit, 30 * unit, 25 * unit, -170,
  24 * unit, 28 * unit, 29 * unit, 2 * unit, 28 * unit, 8 * unit, 0, 5 * unit,
  26 * unit, 12 * unit, 23 * unit, 1 * unit,
}};
} // namespace paper_test
