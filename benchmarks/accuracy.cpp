#include "sketch_variant.hpp"
#include "accuracy_workload.hpp"

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

using namespace accuracy_workload;

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

          Sketch sketch(capacities[capacity]);
          if (workflow == 0) {
            for (double x : data) sketch.add(x);
          } else if (workflow == 1) {
            std::vector<Sketch> parts;
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
#if defined(SPLINESKETCH_CERTIFIED) || defined(SPLINESKETCH_PAPER_CERTIFIED)
            const auto certificate = sketch.rank_with_error(x);
            if (certificate.lower_rank > exact || certificate.upper_rank < exact ||
                std::abs(estimate - static_cast<double>(exact)) > certificate.max_error ||
                std::abs(estimate - static_cast<double>(exact)) > sketch.max_rank_error())
              throw std::runtime_error("rank certificate violated");
#endif
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
