// Read-only diagnostics. SPLINESKETCH_TESTING exposes the existing inspector
// friendship; this executable does not alter the production algorithm.
#define SPLINESKETCH_TESTING
#include <splinesketch/splinesketch.hpp>
#include "../accuracy_workload.hpp"
#include <iomanip>
#include <iostream>
#include <string>

namespace splinesketch {
struct SplineSketchInvariantInspector {
  static void dump(const SplineSketch& s) {
    std::cout << "n=" << s.count_ << " epoch=" << s.epoch_end_
              << " factor=" << double(s.bucket_bound_factor_)
              << " pending=" << s.pending_.size() << " held=" << s.heavy_.size() << '\n';
    for (const auto& n:s.nodes_)
      std::cout << n.x << ':' << double(n.mass) << ':' << n.protected_threshold << ' ';
    std::cout << '\n';
  }
};
}

using splinesketch::SplineSketch;
using splinesketch::SplineSketchInvariantInspector;
using namespace accuracy_workload;

bool fails(const std::vector<double>& data) {
  if(data.size()<8) return false;
  SplineSketch s(8); for(double x:data) s.add(x); s.consolidate();
  const auto exact=std::count_if(data.begin(),data.end(),[](double x){return x<=7.5;});
  return std::abs(s.rank(7.5)-exact)>0.20*data.size();
}

int main(int argc,char** argv) {
  std::cout << std::setprecision(17);
  if(argc==2 && std::string(argv[1])=="--minimize") {
    Random rng{0x6a09e667f3bcc909ULL+0x100000001b3ULL+3};
    std::vector<double> data;
    for(std::size_t i=0;i<observations;++i) data.push_back(generate(3,i,rng));
    if(!fails(data)) throw std::runtime_error("original witness no longer fails");
    // Delta debugging to deletion-minimality for this predicate, not a claim
    // of globally shortest input. Preserve order, capacity, and query.
    std::size_t granularity=2;
    while(data.size()>1) {
      bool removed=false; const auto chunk=(data.size()+granularity-1)/granularity;
      for(std::size_t start=0;start<data.size();start+=chunk) {
        auto candidate=data;
        candidate.erase(candidate.begin()+start,candidate.begin()+std::min(start+chunk,candidate.size()));
        if(fails(candidate)) {data=std::move(candidate);granularity=std::max<std::size_t>(2,granularity-1);removed=true;break;}
      }
      if(!removed) {if(granularity>=data.size())break;granularity=std::min(data.size(),granularity*2);}
    }
    for(double x:data) std::cout << x << '\n';
  } else if(argc==2 && std::string(argv[1])=="--stdin") {
    SplineSketch sketch(8); double x; std::size_t exact=0;
    while(std::cin>>x) {sketch.add(x);exact+=x<=7.5;if(sketch.count()%8==0)SplineSketchInvariantInspector::dump(sketch);}
    sketch.consolidate();
    std::cout << "query=7.5 exact=" << exact << " estimate=" << sketch.rank(7.5) << '\n';
  } else if(argc==1) {
    for(int workflow:{0,2}) for(std::size_t seed=0;seed<3;++seed) {
      Random rng{0x6a09e667f3bcc909ULL+seed*0x100000001b3ULL+3};
      SplineSketch sketch(workflow==2?16:8);
      std::vector<double> data; double previous_error=0;
      std::cout << "WORKFLOW " << workflow << " SEED " << seed << '\n';
      for(std::size_t i=0;i<observations;++i) {
        auto before=sketch;
        if(workflow==2 && i==3000) sketch.resize(6);
        const double x=generate(3,i,rng); data.push_back(x); sketch.add(x);
        const auto exact=std::count_if(data.begin(),data.end(),[](double v){return v<=0.999;});
        const double error=sketch.rank(0.999)-exact;
        if(std::abs(error-previous_error)>100) {
          std::cout << "JUMP " << i+1 << " value=" << x << " err=" << error << " delta=" << error-previous_error << '\n';
          SplineSketchInvariantInspector::dump(before); SplineSketchInvariantInspector::dump(sketch);
        }
        previous_error=error;
      }
      if(workflow==2) sketch.resize(8);
      sketch.consolidate();std::sort(data.begin(),data.end());
      double worst=0,query=0;
      for(double x:make_queries(data)) {
        const double error=std::abs(sketch.rank(x)-(std::upper_bound(data.begin(),data.end(),x)-data.begin()));
        if(error>worst){worst=error;query=x;}
      }
      std::cout << "FINAL max_percent=" << worst*100/data.size() << " query=" << query << '\n';
      SplineSketchInvariantInspector::dump(sketch);
    }
  } else {std::cerr << "usage: outliers [--minimize|--stdin]\n";return 2;}
}
