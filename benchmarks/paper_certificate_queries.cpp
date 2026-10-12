#include <splinesketch/certificate_midpoint.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using Sketch = splinesketch::CertifiedPaperSplineSketch;
static volatile double sink = 0;

template<class Query>
static double nanoseconds_per_query(std::size_t count, Query query) {
  volatile double sum = 0;
  const auto start = std::chrono::steady_clock::now();
  for (std::size_t i = 0; i < count; ++i) sum = sum + query(i);
  const auto end = std::chrono::steady_clock::now();
  sink = sink + sum;
  return std::chrono::duration<double, std::nano>(end - start).count() / count;
}

static double median(std::array<double, 5> samples) {
  std::sort(samples.begin(), samples.end());
  return samples[samples.size() / 2];
}

static void benchmark(std::size_t capacity, bool discrete, bool buffered) {
  constexpr std::size_t ranks = 8192, selections = 128;
  std::mt19937_64 random(239);
  std::normal_distribution<double> normal;
  std::vector<double> values;
  const auto next = [&] { return discrete ? static_cast<double>(random() % 17) : normal(random); };
  Sketch sketch(capacity);
  for (unsigned i = 0; i < 20000; ++i) {
    values.push_back(next()); sketch.add(values.back());
  }
  sketch.consolidate();
  // Fewer than 5*k pending additions keeps the buffered state unflushed.
  for (std::size_t i = 0; i < 2 * capacity - 1; ++i) {
    values.push_back(next()); sketch.add(values.back());
  }
  if (!buffered) sketch.consolidate();
  std::array<double, ranks> queries;
  for (auto& query : queries) query = values[random() % values.size()];
  std::array<double, selections> quantiles;
  std::array<std::uint64_t, selections> targets;
  for (std::size_t i = 0; i < selections; ++i) {
    quantiles[i] = static_cast<double>(i + 1) / (selections + 1);
    targets[i] = splinesketch::detail::quantile_rank_target(quantiles[i], sketch.count());
  }
  const auto spline_rank = [&](std::size_t i) {
    const auto result = sketch.rank_with_error(queries[i]);
    return result.estimate + result.max_error;
  };
  const auto midpoint_rank = [&](std::size_t i) {
    const auto result = splinesketch::certificate_midpoint_rank(sketch, queries[i]);
    return result.estimate + result.max_error;
  };
  const auto spline_quantile = [&](std::size_t i) { return sketch.quantile(quantiles[i]); };
  const auto midpoint_select = [&](std::size_t i) {
    return splinesketch::certificate_midpoint_select(sketch, targets[i]);
  };
  // Warm each path, then alternate paired measurement order to reduce bias.
  nanoseconds_per_query(ranks, spline_rank); nanoseconds_per_query(ranks, midpoint_rank);
  nanoseconds_per_query(selections, spline_quantile); nanoseconds_per_query(selections, midpoint_select);
  std::array<double, 5> rank_times, midpoint_times, quantile_times, selection_times;
  for (unsigned round = 0; round < rank_times.size(); ++round) {
    if (round % 2 == 0) {
      rank_times[round] = nanoseconds_per_query(ranks, spline_rank);
      midpoint_times[round] = nanoseconds_per_query(ranks, midpoint_rank);
      quantile_times[round] = nanoseconds_per_query(selections, spline_quantile);
      selection_times[round] = nanoseconds_per_query(selections, midpoint_select);
    } else {
      midpoint_times[round] = nanoseconds_per_query(ranks, midpoint_rank);
      rank_times[round] = nanoseconds_per_query(ranks, spline_rank);
      selection_times[round] = nanoseconds_per_query(selections, midpoint_select);
      quantile_times[round] = nanoseconds_per_query(selections, spline_quantile);
    }
  }
  std::cout << capacity << ',' << (discrete ? "discrete" : "normal") << ','
            << (buffered ? "buffered" : "consolidated") << ',' << sketch.count() << ','
            << sketch.max_rank_uncertainty() << ',' << median(rank_times) << ',' << median(midpoint_times) << ','
            << median(quantile_times) << ',' << median(selection_times) << '\n';
}

int main(int argc, char** argv) {
  if (argc > 2) throw std::invalid_argument("usage: paper-certificate-queries [capacity]");
  const std::vector<std::size_t> capacities = argc == 2 ?
      std::vector<std::size_t>{std::stoul(argv[1])} : std::vector<std::size_t>{16, 128, 1024};
  std::cout << std::setprecision(9)
            << "capacity,shape,state,count,width,spline_rank_with_error_ns,midpoint_rank_ns,spline_quantile_ns,midpoint_select_ns\n";
  for (auto capacity : capacities)
    for (bool discrete : {false, true})
      for (bool buffered : {false, true}) benchmark(capacity, discrete, buffered);
  return !std::isfinite(sink);
}
