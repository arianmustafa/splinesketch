// Compile with macros naming the reconstructed previous headers. Loading the
// current core first lets both certified headers include it harmlessly under
// pragma once; the previous core is explicitly loaded in a distinct namespace.
#include <splinesketch/splinesketch.hpp>
#define splinesketch previous_splinesketch
#include SPLINESKETCH_PREVIOUS_CORE
#include SPLINESKETCH_PREVIOUS_CERTIFIED
#undef splinesketch
#include <splinesketch/certified_splinesketch.hpp>
#include "../accuracy_workload.hpp"
#include <iostream>

using namespace accuracy_workload;

template<class Sketch>
static void fill(Sketch& sketch, const std::vector<double>& data,
                 unsigned workflow, std::size_t capacity) {
  if (workflow == 0) for (double x : data) sketch.add(x);
  else if (workflow == 1) {
    std::vector<Sketch> parts;
    for (unsigned i = 0; i < 4; ++i) parts.emplace_back(capacity);
    for (std::size_t i = 0; i < data.size(); ++i) parts[i % 4].add(data[i]);
    parts[0].merge(parts[1]); parts[2].merge(parts[3]); parts[0].merge(parts[2]);
    sketch = std::move(parts[0]);
  } else {
    sketch.resize(2 * capacity);
    for (std::size_t i = 0; i < data.size()/2; ++i) sketch.add(data[i]);
    sketch.resize(std::max<std::size_t>(6, capacity/2));
    for (std::size_t i = data.size()/2; i < data.size(); ++i) sketch.add(data[i]);
    sketch.resize(capacity);
  }
  sketch.consolidate();
}

int main() {
  std::size_t cases = 0, queries = 0;
  for (auto capacity : capacities) for (unsigned shape = 0; shape < 5; ++shape)
    for (unsigned seed = 0; seed < 43; ++seed) {
      if (seed >= 3 && seed < 23) continue;
      Random random{0x6a09e667f3bcc909ULL + seed*0x100000001b3ULL + shape};
      std::vector<double> data;
      for (std::size_t i = 0; i < observations; ++i) data.push_back(generate(shape,i,random));
      for (unsigned workflow = seed < 3 ? 0 : 2; workflow < 3; ++workflow) {
        previous_splinesketch::CertifiedSplineSketch before(capacity);
        splinesketch::CertifiedSplineSketch after(capacity);
        fill(before,data,workflow,capacity); fill(after,data,workflow,capacity);
        auto sorted = data; std::sort(sorted.begin(), sorted.end());
        for (double x : make_queries(sorted)) {
          const auto a = before.rank_with_error(x);
          const auto b = after.rank_with_error(x);
          if (a.estimate != b.estimate || a.lower_rank != b.lower_rank ||
              a.upper_rank != b.upper_rank || a.max_error != b.max_error ||
              before.max_rank_uncertainty() != after.max_rank_uncertainty())
            throw std::runtime_error("previous/candidate certificate mismatch");
          ++queries;
        }
        for (double q : {0., .1, .5, .9, 1.})
          if (before.quantile(q) != after.quantile(q))
            throw std::runtime_error("previous/candidate quantile mismatch");
        ++cases;
      }
    }
  std::cout << "Exact estimates, certificates and quantiles match: "
            << cases << " cases, " << queries << " queries\n";
}
