// Unseen-seed accuracy and shrink-cost comparison. Compile this same source
// against the baseline and candidate headers, then compare the emitted CSVs.
#include "../sketch_variant.hpp"
#include "../accuracy_workload.hpp"
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using SplineSketch = Sketch;
using namespace accuracy_workload;
using Clock = std::chrono::steady_clock;

double elapsed_us(Clock::time_point start, Clock::time_point end) {
  return std::chrono::duration<double, std::micro>(end-start).count();
}

int main(int argc,char** argv) {
  std::cout<<std::setprecision(17);
  if(argc==2 && std::string(argv[1])=="--timing") {
    std::cout<<"from,to,trial,resize_us\n";
    for(std::size_t k:{16U,32U,64U,128U,256U,512U}) {
      Random rng{0xa54ff53a5f1d36f1ULL};SplineSketch baseline(k);
      for(std::size_t i=0;i<6000;++i) baseline.add(generate(3,i,rng));
      baseline.consolidate();
      for(int trial=0;trial<5;++trial) {
        auto s=baseline;auto start=Clock::now();s.resize(std::max<std::size_t>(6,k/2));auto end=Clock::now();
        if(s.count()!=6000||s.bucket_count()>s.bucket_capacity())throw std::runtime_error("invalid resize");
        std::cout<<k<<','<<s.bucket_capacity()<<','<<trial<<','<<elapsed_us(start,end)<<'\n';
      }
    }
    return 0;
  }
  const bool fresh = argc == 2 && std::string(argv[1]) == "--fresh";
  const bool heldout = argc == 2 && std::string(argv[1]) == "--heldout";
  if(argc!=1 && !fresh && !heldout){std::cerr<<"usage: resize [--timing|--fresh|--heldout]\n";return 2;}
  const std::size_t first_seed = heldout ? 43 : fresh ? 23 : 3;
  std::cout<<"k,shape,seed,input_fingerprint,queries,median_percent,p95_percent,max_percent,shrink_us,grow_us\n";
  for(auto k:capacities)for(std::size_t shape=0;shape<5;++shape)for(std::size_t seed=first_seed;seed<first_seed+20;++seed) {
    Random rng{0x6a09e667f3bcc909ULL+seed*0x100000001b3ULL+shape};std::vector<double> data;
    for(std::size_t i=0;i<6000;++i)data.push_back(generate(shape,i,rng));
    std::uint64_t fingerprint=14695981039346656037ULL;
    for(double x:data) {
      std::uint64_t bits;std::memcpy(&bits,&x,sizeof bits);
      fingerprint=(fingerprint^bits)*1099511628211ULL;
    }
    SplineSketch s(2*k);
    for(std::size_t i=0;i<3000;++i)s.add(data[i]);
    auto start=Clock::now();s.resize(std::max<std::size_t>(6,k/2));auto shrunk=Clock::now();
    for(std::size_t i=3000;i<6000;++i)s.add(data[i]);
    auto grow_start=Clock::now();s.resize(k);auto grown=Clock::now();s.consolidate();
    if(s.count()!=6000||s.bucket_count()>k)throw std::runtime_error("invalid resize");
    std::sort(data.begin(),data.end());std::vector<double> errors;
    for(double x:make_queries(data)) {
      double rank=s.rank(x);
      if(!std::isfinite(rank)||rank<0||rank>6000)throw std::runtime_error("invalid rank");
      const auto exact=std::upper_bound(data.begin(),data.end(),x)-data.begin();
      errors.push_back(100*std::abs(rank-exact)/6000);
    }
    std::sort(errors.begin(),errors.end());
    std::cout<<k<<','<<shapes[shape]<<','<<seed<<','<<fingerprint<<','<<errors.size()<<','<<errors[(errors.size()-1)/2]<<','
      <<errors[(errors.size()-1)*95/100]<<','<<errors.back()<<','<<elapsed_us(start,shrunk)<<','<<elapsed_us(grow_start,grown)<<'\n';
  }
}
