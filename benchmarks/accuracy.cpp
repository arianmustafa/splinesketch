#include <splinesketch/splinesketch.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

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
      if (rng.next() % 50 == 0)
        return (rng.next() & 1 ? 1.0 : -1.0) *
               static_cast<double>(10000 + rng.next() % 1000);
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

struct Sample {
  double error;
  std::size_t capacity, shape, workflow, seed;
  double query, estimate;
  std::size_t exact;
};

struct Summary {
  std::vector<double> errors;
  Sample worst{};
  void add(const Sample& sample) {
    errors.push_back(sample.error);
    if (errors.size() == 1 || sample.error > worst.error) worst = sample;
  }
  void print(const std::string& group, const std::string& name) {
    std::sort(errors.begin(), errors.end());
    const auto percentile = [&](std::size_t numerator, std::size_t denominator) {
      return errors[(errors.size() - 1) * numerator / denominator];
    };
    // Percent of the stream, not percent of the estimated or exact rank.
    std::cout << group << ',' << name << ',' << errors.size() << ','
              << percentile(1, 2) * 100 << ',' << percentile(95, 100) * 100
              << ',' << errors.back() * 100 << '\n';
  }
};

}  // namespace

int main(int argc, char** argv) {
  const bool details = argc == 2 && std::string(argv[1]) == "--details";
  if (argc > 1 && !details) {
    std::cerr << "usage: " << argv[0] << " [--details]\n";
    return 2;
  }

  Summary overall, by_capacity[5], by_shape[5], by_workflow[3];
  std::cout << std::fixed << std::setprecision(6);
  if (details) std::cout << "group,name,queries,median_percent,p95_percent,max_percent\n";
  for (std::size_t capacity = 0; capacity < 5; ++capacity) {
    for (std::size_t shape = 0; shape < 5; ++shape) {
      for (std::size_t workflow = 0; workflow < 3; ++workflow) {
        Summary case_summary;
        for (std::size_t seed = 0; seed < seeds; ++seed) {
          Random rng{0x6a09e667f3bcc909ULL + seed * 0x100000001b3ULL + shape};
          std::vector<double> data;
          data.reserve(observations);
          for (std::size_t i = 0; i < observations; ++i)
            data.push_back(generate(shape, i, rng));

          splinesketch::SplineSketch sketch(capacities[capacity]);
          if (workflow == 0) {
            for (double x : data) sketch.add(x);
          } else if (workflow == 1) {
            std::vector<splinesketch::SplineSketch> parts;
            for (int i = 0; i < 4; ++i) parts.emplace_back(capacities[capacity]);
            for (std::size_t i = 0; i < data.size(); ++i) parts[i % 4].add(data[i]);
            parts[0].merge(parts[1]);
            parts[2].merge(parts[3]);
            parts[0].merge(parts[2]);
            sketch = std::move(parts[0]);
          } else {
            sketch.resize(capacities[capacity] * 2);
            for (std::size_t i = 0; i < data.size() / 2; ++i) sketch.add(data[i]);
            sketch.resize(std::max<std::size_t>(6, capacities[capacity] / 2));
            for (std::size_t i = data.size() / 2; i < data.size(); ++i) sketch.add(data[i]);
            sketch.resize(capacities[capacity]);
          }
          sketch.consolidate();
          if (sketch.count() != data.size()) throw std::runtime_error("count mismatch");

          std::sort(data.begin(), data.end());
          for (double x : make_queries(data)) {
            const std::size_t exact = std::upper_bound(data.begin(), data.end(), x) - data.begin();
            const double estimate = sketch.rank(x);
            if (!std::isfinite(estimate) || estimate < 0 || estimate > data.size())
              throw std::runtime_error("rank outside valid range");
            Sample sample{std::abs(estimate - static_cast<double>(exact)) / data.size(),
                          capacity, shape, workflow, seed, x, estimate, exact};
            case_summary.add(sample);
            by_capacity[capacity].add(sample);
            by_shape[shape].add(sample);
            by_workflow[workflow].add(sample);
            overall.add(sample);
          }
        }
        if (details) case_summary.print("case", std::to_string(capacities[capacity]) + "/" +
                                        shapes[shape] + "/" + workflows[workflow]);
      }
    }
  }

  if (!details) std::cout << "group,name,queries,median_percent,p95_percent,max_percent\n";
  for (std::size_t i = 0; i < 5; ++i)
    by_capacity[i].print("capacity", std::to_string(capacities[i]));
  for (std::size_t i = 0; i < 5; ++i) by_shape[i].print("shape", shapes[i]);
  for (std::size_t i = 0; i < 3; ++i) by_workflow[i].print("workflow", workflows[i]);
  overall.print("overall", "all");
  const auto& worst = overall.worst;
  std::cerr << std::setprecision(17) << "worst case: capacity=" << capacities[worst.capacity]
            << " shape=" << shapes[worst.shape] << " workflow=" << workflows[worst.workflow]
            << " seed=" << worst.seed << " query=" << worst.query << " exact=" << worst.exact
            << " estimate=" << worst.estimate << '\n';
}
