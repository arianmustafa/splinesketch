// Linux/glibc measurement adapter. Allocation accounting is enabled only for
// sketch operations; dataset storage, queries, and CSV output are excluded.
#include "../sketch_variant.hpp"
#include "../accuracy_workload.hpp"
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <malloc.h>
#include <new>
#include <string>

namespace memory {
bool active = false;
std::size_t live = 0, peak = 0;
}
void* operator new(std::size_t n) {
  void* p = std::malloc(n ? n : 1);
  if (!p) throw std::bad_alloc();
  if (memory::active) {
    memory::live += malloc_usable_size(p);
    memory::peak = std::max(memory::peak, memory::live);
  }
  return p;
}
void operator delete(void* p) noexcept {
  if (p && memory::active) memory::live -= malloc_usable_size(p);
  std::free(p);
}
void operator delete(void* p, std::size_t) noexcept { ::operator delete(p); }
void* operator new[](std::size_t n) { return ::operator new(n); }
void operator delete[](void* p) noexcept { ::operator delete(p); }
void operator delete[](void* p, std::size_t) noexcept { ::operator delete(p); }

int main(int argc, char** argv) {
  const bool grid = argc == 3 && std::string(argv[2]) == "--grid";
  const bool fresh = argc == 3 && std::string(argv[2]) == "--certificate-heldout";
  const bool heldout = fresh || (argc == 3 && std::string(argv[2]) == "--large-heldout");
  const bool large = heldout || (argc == 3 && std::string(argv[2]) == "--large");
  if (argc != 2 && !grid && !large) { std::cerr << "usage: local OUTPUT_DIRECTORY [--grid|--large|--large-heldout|--certificate-heldout]\n"; return 2; }
  using namespace accuracy_workload;
  using Clock = std::chrono::steady_clock;
  const std::size_t count = large ? 60000 : observations;
  std::cout << "variant,k,shape,seed,resident_peak,update_peak,final_bytes,update_ns,rank_ns,quantile_ns,median_percent,p95_percent,max_percent\n" << std::setprecision(17);
  for (std::size_t shape=0; shape<5; ++shape) for (std::size_t seed=fresh ? 71 : heldout ? 67 : 0;seed<(fresh ? 74U : heldout ? 70U : 3U);++seed) {
    Random rng{0x6a09e667f3bcc909ULL + seed * 0x100000001b3ULL + shape};
    std::vector<double> data;
    for (std::size_t i=0;i<count;++i) data.push_back(generate(shape,i,rng));
    auto sorted=data; std::sort(sorted.begin(),sorted.end());
    auto queries=make_queries(sorted);
    const std::string base=std::string(argv[1])+"/"+shapes[shape]+"-"+std::to_string(seed);
    { std::ofstream f(base+".data"); f<<std::setprecision(17); for(double x:data) f<<x<<'\n'; }
    { std::ofstream f(base+".queries"); f<<std::setprecision(17); for(double x:queries) f<<x<<'\n'; }
    const std::vector<std::size_t> tested_capacities = large
        ? std::vector<std::size_t>{128,256,512,1024}
        : grid
        ? std::vector<std::size_t>{6,8,12,16,24,32,48,64,80,96,128}
        : std::vector<std::size_t>(std::begin(capacities), std::end(capacities));
    for (auto k:tested_capacities) {
      std::size_t resident=0,peak=0,final=0;
      memory::live=memory::peak=0; memory::active=true;
      {
        Sketch sketch(k);
        resident=memory::live;
        for(double x:data) { sketch.add(x); resident=std::max(resident,memory::live); }
        sketch.consolidate(); resident=std::max(resident,memory::live);
        peak=memory::peak; final=memory::live;
      }
      memory::active=false;
      if(memory::live) throw std::runtime_error("allocation tracking imbalance");
      const auto object=sizeof(Sketch);
      std::vector<double> errors;
      Sketch sketch(k);
      auto start=Clock::now();
      for(double x:data) sketch.add(x);
      sketch.consolidate(); auto updated=Clock::now();
      volatile double sink=0;
      for(int rep=0;rep<10;++rep) for(double x:queries) sink=sink+sketch.rank(x);
      auto ranked=Clock::now();
      for(int i=1;i<1000;++i) sink=sink+sketch.quantile(i/1000.0);
      auto inverted=Clock::now();
      for(double x:queries) {
        double z=sketch.rank(x);
        if(!std::isfinite(z)||z<0||z>count) throw std::runtime_error("invalid rank");
        errors.push_back(100*std::abs(z-(std::upper_bound(sorted.begin(),sorted.end(),x)-sorted.begin()))/count);
      }
      std::sort(errors.begin(),errors.end());
      auto ns=[](auto a,auto b){return std::chrono::duration<double,std::nano>(b-a).count();};
#if defined(SPLINESKETCH_PAPER_CERTIFIED)
      constexpr const char* variant = "cpp_paper_certified,";
#elif defined(SPLINESKETCH_PAPER_THEORETICAL)
      constexpr const char* variant = "cpp_paper_theoretical,";
#elif defined(SPLINESKETCH_PAPER)
      constexpr const char* variant = "cpp_paper,";
#elif defined(SPLINESKETCH_CERTIFIED)
      constexpr const char* variant = "cpp_certified,";
#else
      constexpr const char* variant = "cpp,";
#endif
      std::cout<<variant<<k<<','<<shapes[shape]<<','<<seed<<','<<resident+object<<','<<peak+object<<','<<final+object<<','
        <<ns(start,updated)/count<<','<<ns(updated,ranked)/(10*queries.size())<<','<<ns(ranked,inverted)/999<<','
        <<errors[(errors.size()-1)/2]<<','<<errors[(errors.size()-1)*95/100]<<','<<errors.back()<<'\n';
    }
  }
}
