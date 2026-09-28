#include <splinesketch/splinesketch.hpp>

#include <chrono>
#include <cstddef>
#include <iostream>
#include <random>
#include <vector>

int main() {
  constexpr std::size_t observations = 1'000'000;
  constexpr std::size_t queries = 100'000;
  constexpr std::size_t quantiles = 10'000;
  std::mt19937_64 rng(12345);
  std::normal_distribution<double> normal;
  std::vector<double> input(observations);
  for (auto& value : input) value = normal(rng);

  splinesketch::SplineSketch sketch(128);
  const auto start = std::chrono::steady_clock::now();
  for (const double value : input) sketch.add(value);
  sketch.consolidate();
  const auto after_updates = std::chrono::steady_clock::now();
  volatile double sink = 0;
  for (std::size_t i = 0; i < queries; ++i)
    sink = sink + sketch.rank(input[i]);
  const auto after_queries = std::chrono::steady_clock::now();
  for (std::size_t i = 0; i < quantiles; ++i)
    sink = sink + sketch.quantile((i % 999 + 1) / 1000.0);
  const auto after_quantiles = std::chrono::steady_clock::now();

  const auto updates_us = std::chrono::duration<double, std::micro>(after_updates - start).count();
  const auto queries_us = std::chrono::duration<double, std::micro>(after_queries - after_updates).count();
  const auto quantiles_us = std::chrono::duration<double, std::micro>(after_quantiles - after_queries).count();
  std::cout << "update: " << updates_us / observations << " us/item\n"
            << "rank:   " << queries_us / queries << " us/query\n"
            << "quantile: " << quantiles_us / quantiles << " us/query ("
            << sketch.heavy_hitter_count() << " heavy hitters)\n"
            << "sink: " << sink << '\n';
}
