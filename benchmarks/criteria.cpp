#include <splinesketch/splinesketch.hpp>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <random>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  const std::size_t buckets = argc > 1 ? std::stoul(argv[1]) : 128;
  constexpr std::size_t parts = 10'000;
  constexpr std::size_t criteria = 30;

  std::mt19937_64 rng(1984);
  std::normal_distribution<double> normal;
  std::vector<double> input(parts * criteria);
  for (std::size_t part = 0; part < parts; ++part)
    for (std::size_t criterion = 0; criterion < criteria; ++criterion)
      input[part * criteria + criterion] = normal(rng) + criterion * 0.1;

  std::vector<splinesketch::SplineSketch> sketches;
  sketches.reserve(criteria);
  for (std::size_t i = 0; i < criteria; ++i) sketches.emplace_back(buckets);

  std::vector<double> milliseconds;
  milliseconds.reserve(parts);
  const auto total_start = std::chrono::steady_clock::now();
  for (std::size_t part = 0; part < parts; ++part) {
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t criterion = 0; criterion < criteria; ++criterion)
      sketches[criterion].add(input[part * criteria + criterion]);
    const auto end = std::chrono::steady_clock::now();
    milliseconds.push_back(std::chrono::duration<double, std::milli>(end - start).count());
  }
  const auto total_end = std::chrono::steady_clock::now();
  std::sort(milliseconds.begin(), milliseconds.end());
  const auto percentile = [&](double q) {
    return milliseconds[static_cast<std::size_t>((parts - 1) * q)];
  };
  const double mean = std::chrono::duration<double, std::milli>(total_end - total_start).count() / parts;
  std::cout << "30 criteria, " << buckets << " buckets, milliseconds per part:\n"
            << "mean " << mean << ", p99 " << percentile(0.99)
            << ", p99.9 " << percentile(0.999)
            << ", max " << milliseconds.back() << '\n';
}
