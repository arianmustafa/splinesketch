// Linux/glibc allocation accounting, scoped to sketch operations.
#include <splinesketch/paper_splinesketch.hpp>
#include "../accuracy_workload.hpp"
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <malloc.h>
#include <new>
namespace memory { bool active = false; std::size_t live = 0, peak = 0; }
[[gnu::noinline]] void* operator new(std::size_t n) {
  void* p = std::malloc(n ? n : 1);
  if (!p) throw std::bad_alloc();
  if (memory::active) { memory::live += malloc_usable_size(p); memory::peak = std::max(memory::peak, memory::live); }
  return p;
}
[[gnu::noinline]] void operator delete(void* p) noexcept {
  if (p && memory::active) memory::live -= malloc_usable_size(p);
  std::free(p);
}
void operator delete(void* p, std::size_t) noexcept { ::operator delete(p); }
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete[](void* p) noexcept { ::operator delete(p); }
void operator delete[](void* p, std::size_t) noexcept { ::operator delete(p); }
int main() {
  using Sketch = splinesketch::CertifiedPaperSplineSketch;
  using Clock = std::chrono::steady_clock;
  std::cout << std::setprecision(17) << "k,shape,seed,n,resident_peak,update_peak,final_bytes,update_ns,rank_ns\n";
  constexpr unsigned count = 60000;
  for (unsigned k : {128U, 256U, 512U, 1024U}) for (unsigned shape = 0; shape < 5; ++shape)
    for (unsigned seed = 0; seed < 3; ++seed) {
      accuracy_workload::Random rng{0x6a09e667f3bcc909ULL + seed * 0x100000001b3ULL + shape};
      std::vector<double> data;
      for (unsigned i = 0; i < count; ++i) data.push_back(accuracy_workload::generate(shape, i, rng));
      auto sorted = data; std::sort(sorted.begin(), sorted.end());
      const auto queries = accuracy_workload::make_queries(sorted);
      std::size_t resident = 0, peak = 0, final = 0;
      memory::live = memory::peak = 0; memory::active = true;
      {
        Sketch sketch(k); resident = memory::live;
        for (double x : data) { sketch.add(x); resident = std::max(resident, memory::live); }
        sketch.consolidate(); resident = std::max(resident, memory::live);
        peak = memory::peak; final = memory::live;
      }
      memory::active = false;
      if (memory::live) throw std::runtime_error("allocation tracking imbalance");
      Sketch sketch(k);
      const auto start = Clock::now();
      for (double x : data) sketch.add(x);
      sketch.consolidate(); const auto updated = Clock::now();
      volatile double sink = 0;
      for (unsigned repeat = 0; repeat < 10; ++repeat) for (double x : queries) sink = sink + sketch.rank(x);
      const auto ranked = Clock::now();
      const auto ns = [](auto a, auto b) { return std::chrono::duration<double, std::nano>(b - a).count(); };
      std::cout << k << ',' << accuracy_workload::shapes[shape] << ',' << seed << ',' << count << ','
                << resident + sizeof(Sketch) << ',' << peak + sizeof(Sketch) << ',' << final + sizeof(Sketch) << ','
                << ns(start, updated) / count << ',' << ns(updated, ranked) / (10 * queries.size()) << '\n';
    }
}
