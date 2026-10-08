// Isolate the cost of the optional uniform-certificate query.
#include <splinesketch/paper_splinesketch.hpp>
#include "../accuracy_workload.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
using Sketch = splinesketch::CertifiedPaperSplineSketch;
static double query(const Sketch& sketch) { return sketch.max_rank_error(); }
int main() {
  std::cout << "k,n,nodes,width,uniform,ns_per_call\n" << std::setprecision(17);
  double (*volatile call)(const Sketch&) = query;
  for (auto k : {128U, 256U, 512U, 1024U}) {
    Sketch sketch(k); accuracy_workload::Random random{0x6a09e667f3bcc909ULL};
    for (unsigned i = 0; i < 60000; ++i) sketch.add(accuracy_workload::generate(0, i, random));
    sketch.consolidate(); volatile double sink = 0;
    constexpr unsigned repetitions = 10000;
    const auto start = std::chrono::steady_clock::now();
    for (unsigned i = 0; i < repetitions; ++i) sink = sink + call(sketch);
    const auto stop = std::chrono::steady_clock::now();
    const auto elapsed = std::chrono::duration<double, std::nano>(stop - start).count() / repetitions;
    std::cout << k << ',' << sketch.count() << ',' << sketch.bucket_count() << ','
              << sketch.max_rank_uncertainty() << ',' << sketch.max_rank_error() << ',' << elapsed << '\n';
  }
}
