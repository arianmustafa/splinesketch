#include <splinesketch/splinesketch.hpp>
#define splinesketch previous_splinesketch
#include SPLINESKETCH_PREVIOUS_CORE
#include SPLINESKETCH_PREVIOUS_PAPER
#undef splinesketch
#include <splinesketch/paper_splinesketch.hpp>
#include "../accuracy_workload.hpp"
#include <iostream>
#include <string>

using namespace accuracy_workload;

template<class Sketch>
static void fill(Sketch& s, const std::vector<double>& data, unsigned workflow, std::size_t k,
                 typename Sketch::BoundPolicy policy = Sketch::BoundPolicy::practical) {
  if (workflow == 3) {
    fill(s, data, 0, k, policy);
    s.finalize();
    return;
  }
  if (workflow == 0) for (double x : data) s.add(x);
  else if (workflow == 1) {
    std::vector<Sketch> parts;
    for (int i = 0; i < 4; ++i) parts.emplace_back(k, policy);
    for (std::size_t i = 0; i < data.size(); ++i) parts[i % 4].add(data[i]);
    parts[0].merge(parts[1]); parts[2].merge(parts[3]); parts[0].merge(parts[2]);
    s = std::move(parts[0]);
  } else {
    s.resize(2 * k);
    for (std::size_t i = 0; i < data.size() / 2; ++i) s.add(data[i]);
    s.resize(std::max<std::size_t>(6, k / 2));
    for (std::size_t i = data.size() / 2; i < data.size(); ++i) s.add(data[i]);
    s.resize(k);
  }
  s.consolidate();
}

int main(int argc, char** argv) {
  const bool extended = argc == 2 && std::string(argv[1]) == "--extended";
  if (argc != 1 && !extended) throw std::invalid_argument("usage: equivalence [--extended]");
  const std::vector<std::size_t> tested_capacities = extended ? std::vector<std::size_t>{32,128,256,512,1024}
      : std::vector<std::size_t>(std::begin(capacities), std::end(capacities));
  const std::size_t count = extended ? 24000 : observations;
  std::size_t cases = 0, queries = 0, inverse_queries = 0;
  for (unsigned policy = 0; policy < (extended ? 2U : 1U); ++policy)
  for (auto k : tested_capacities) for (unsigned shape = 0; shape < 5; ++shape)
    for (unsigned seed = extended ? 63 : 0; seed < (extended ? 67 : 63); ++seed) {
      if (seed >= 3 && seed < 23) continue;
      Random rng{0x6a09e667f3bcc909ULL + seed * 0x100000001b3ULL + shape};
      std::vector<double> data;
      for (std::size_t i = 0; i < count; ++i) data.push_back(generate(shape, i, rng));
      for (unsigned w = extended || seed < 3 ? 0 : 2; w < (extended ? 4U : 3U); ++w) {
        using Before = previous_splinesketch::PaperSplineSketch;
        using After = splinesketch::PaperSplineSketch;
        Before before(k, policy ? Before::BoundPolicy::theoretical : Before::BoundPolicy::practical);
        After after(k, policy ? After::BoundPolicy::theoretical : After::BoundPolicy::practical);
        fill(before, data, w, k, policy ? Before::BoundPolicy::theoretical : Before::BoundPolicy::practical);
        fill(after, data, w, k, policy ? After::BoundPolicy::theoretical : After::BoundPolicy::practical);
        auto sorted = data; std::sort(sorted.begin(), sorted.end());
        for (double x : make_queries(sorted)) {
          if (before.rank(x) != after.rank(x)) throw std::runtime_error("paper optimization rank mismatch");
          ++queries;
        }
        for (double q : {0., .001, .1, .5, .9, .999, 1.}) {
          if (before.quantile(q) != after.quantile(q)) throw std::runtime_error("paper optimization quantile mismatch");
          ++inverse_queries;
        }
        ++cases;
      }
    }
  // Include unflushed raw buffers and large exact MG tables in inverse checks.
  for (std::size_t k : {8U, 32U, 128U}) {
    previous_splinesketch::PaperSplineSketch before(k);
    splinesketch::PaperSplineSketch after(k);
    for (unsigned i = 0; i < 12 * k; ++i) {
      double x = static_cast<double>((i * 97) % 349) / 7;
      before.add(x); after.add(x);
      if (i % 37 == 0)
        for (double q : {.001, .1, .5, .9, .999}) {
          if (before.quantile(q) != after.quantile(q)) throw std::runtime_error("unflushed inverse mismatch");
          ++inverse_queries;
        }
    }
  }
  std::cout << "Exact ranks and quantiles match: " << cases << " cases, " << queries
            << " rank queries, " << inverse_queries << " inverse queries\n";
}
