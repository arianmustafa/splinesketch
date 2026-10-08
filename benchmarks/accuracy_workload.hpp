#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace accuracy_workload {

// SplitMix64 and integer arithmetic keep the generated streams independent of
// the standard library's distribution implementations.
struct Random {
  std::uint64_t state;
  std::uint64_t next() {
    state += 0x9e3779b97f4a7c15ULL;
    std::uint64_t z = state;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }
};

constexpr std::size_t observations = 6000;
constexpr std::size_t seeds = 3;
const std::size_t capacities[] = {8, 16, 32, 64, 128};
const char* shapes[] = {"uniform", "clustered", "duplicates", "outliers", "ascending"};
const char* workflows[] = {"direct", "merge4", "resize"};

double generate(std::size_t shape, std::size_t index, Random& rng) {
  switch (shape) {
    case 0: return static_cast<double>(rng.next() % 1000001) / 1000.0 - 500.0;
    case 1: {
      std::int64_t sum = 0;
      for (int i = 0; i < 4; ++i) sum += static_cast<std::int64_t>(rng.next() % 2001);
      return static_cast<double>(sum - 4000) / 1000.0;
    }
    case 2: {
      const auto choice = rng.next() % 10;
      if (choice < 7) return static_cast<double>(static_cast<int>(choice % 3) - 1);
      return static_cast<double>(static_cast<int>(rng.next() % 129) - 64) / 8.0;
    }
    case 3: {
      if (rng.next() % 50 == 0) {
        // Sequence the random draws explicitly; operand evaluation order of
        // multiplication otherwise permits different streams across builds.
        const double sign = rng.next() & 1 ? 1.0 : -1.0;
        const double magnitude = static_cast<double>(10000 + rng.next() % 1000);
        return sign * magnitude;
      }
      return static_cast<double>(static_cast<int>(rng.next() % 2001) - 1000) / 1000.0;
    }
    default: return static_cast<double>(index) + static_cast<double>(rng.next() % 3) / 10.0;
  }
}

std::vector<double> make_queries(const std::vector<double>& sorted) {
  std::vector<double> queries;
  const double lo = sorted.front(), hi = sorted.back();
  queries.push_back(std::nextafter(lo, -std::numeric_limits<double>::infinity()));
  queries.push_back(std::nextafter(hi, std::numeric_limits<double>::infinity()));
  for (std::size_t i = 0; i <= 256; ++i) {
    const std::size_t at = i * (sorted.size() - 1) / 256;
    const double value = sorted[at];
    queries.push_back(value);
    queries.push_back(std::nextafter(value, -std::numeric_limits<double>::infinity()));
    if (at + 1 < sorted.size() && value < sorted[at + 1])
      queries.push_back(value + (sorted[at + 1] - value) / 2.0);
    queries.push_back(lo + (hi - lo) * static_cast<double>(i) / 256.0);
  }
  std::sort(queries.begin(), queries.end());
  queries.erase(std::unique(queries.begin(), queries.end()), queries.end());
  return queries;
}

}  // namespace accuracy_workload
