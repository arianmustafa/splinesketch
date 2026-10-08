#include <splinesketch/paper_splinesketch.hpp>
#include "../accuracy_workload.hpp"
#include <iostream>
#include <iomanip>
#include <string>

using Certified = splinesketch::CertifiedPaperSplineSketch;
using Paper = splinesketch::PaperSplineSketch;
using namespace accuracy_workload;

template<class Sketch>
void fill(Sketch& sketch, const std::vector<double>& data, unsigned workflow, std::size_t k) {
  if (workflow == 0 || workflow == 4) {
    for (double value : data) sketch.add(value);
  } else if (workflow == 1 || workflow == 2) {
    std::vector<Sketch> parts;
    for (unsigned i = 0; i < 4; ++i) parts.emplace_back(workflow == 2 && i % 2 ? k + 3 : k);
    for (std::size_t i = 0; i < data.size(); ++i) parts[i % 4].add(data[i]);
    parts[0].merge(parts[1]); parts[2].merge(parts[3]); parts[0].merge(parts[2]);
    sketch = std::move(parts[0]);
  } else {
    sketch.resize(k * 2);
    for (std::size_t i = 0; i < data.size() / 2; ++i) sketch.add(data[i]);
    sketch.resize(std::max<std::size_t>(6, k / 2));
    for (std::size_t i = data.size() / 2; i < data.size(); ++i) sketch.add(data[i]);
    sketch.resize(k);
  }
  sketch.consolidate();
  if (workflow == 4) sketch.finalize();
}

int main(int argc, char** argv) {
  const bool fresh = argc == 2 && std::string(argv[1]) == "--fresh";
  if (argc != 1 && !fresh) throw std::invalid_argument("usage: quality [--fresh]");
  const char* names[]{"stream", "balanced_merge", "mixed_merge", "resize", "finalize"};
  std::cout << "k,shape,seed,workflow,n,queries,error_median_percent,error_p95_percent,error_max_percent,"
               "raw_error_max_percent,point_bound_median_percent,point_bound_p95_percent,point_bound_max_percent,"
               "width_median_percent,width_p95_percent,width_max_percent,uniform_bound_percent\n"
            << std::setprecision(17);
  for (auto k : {128U, 256U, 512U, 1024U})
    for (unsigned shape = 0; shape < 5; ++shape)
      for (unsigned seed = fresh ? 71 : 0; seed < (fresh ? 74U : 3U); ++seed) {
        Random rng{0x6a09e667f3bcc909ULL + seed * 0x100000001b3ULL + shape};
        std::vector<double> data;
        for (unsigned i = 0; i < 60000; ++i) data.push_back(generate(shape, i, rng));
        auto sorted = data;
        std::sort(sorted.begin(), sorted.end());
        const auto queries = make_queries(sorted);
        for (unsigned workflow = 0; workflow < 5; ++workflow) {
          Certified sketch(k); Paper raw(k);
          fill(sketch, data, workflow, k); fill(raw, data, workflow, k);
          std::vector<double> errors, point, widths;
          double raw_max = 0;
          const auto uniform = sketch.max_rank_error();
          for (double x : queries) {
            const auto exact = static_cast<std::uint64_t>(
                std::upper_bound(sorted.begin(), sorted.end(), x) - sorted.begin());
            const auto result = sketch.rank_with_error(x);
            const auto error = std::abs(result.estimate - static_cast<double>(exact));
            const auto original_error = std::abs(raw.rank(x) - static_cast<double>(exact));
#ifdef SPLINESKETCH_ALLOW_ALGORITHM_ACCURACY_CHANGES
            constexpr bool require_original_accuracy = false;
#else
            constexpr bool require_original_accuracy = true;
#endif
            if (result.lower_rank > exact || result.upper_rank < exact ||
                error > result.max_error || error > uniform ||
                (require_original_accuracy && error > original_error))
              throw std::runtime_error("invalid paper certificate or accuracy regression");
            errors.push_back(100 * error / data.size());
            point.push_back(100 * result.max_error / data.size());
            widths.push_back(100.0 * static_cast<double>(result.upper_rank - result.lower_rank) / data.size());
            raw_max = std::max(raw_max, 100 * original_error / data.size());
          }
          const auto triple = [](std::vector<double>& values) {
            std::sort(values.begin(), values.end());
            std::cout << ',' << values[(values.size() - 1) / 2] << ','
                      << values[(values.size() - 1) * 95 / 100] << ',' << values.back();
          };
          std::cout << k << ',' << shapes[shape] << ',' << seed << ',' << names[workflow] << ','
                    << sketch.count() << ',' << queries.size();
          triple(errors); std::cout << ',' << raw_max;
          triple(point); triple(widths);
          std::cout << ',' << 100 * uniform / data.size() << '\n';
        }
      }
}
