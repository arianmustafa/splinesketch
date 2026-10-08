// Check the experimental cubic-distance calculation against a dense rank
// oracle. Compile only against the resize_exact experimental header.
#define SPLINESKETCH_TESTING
#include <splinesketch/splinesketch.hpp>
#include "../accuracy_workload.hpp"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <iostream>

namespace splinesketch {
struct SplineSketchInvariantInspector {
  static void check_distance() {
    accuracy_workload::Random rng{123456};
    for (int trial=0;trial<100;++trial) {
      std::vector<SplineSketch::Node> reference;
      double x=0;
      for(int i=0;i<20;++i) {
        x+=std::pow(10.0,static_cast<int>(rng.next()%11)-5);
        reference.push_back(SplineSketch::Node{x,static_cast<long double>(rng.next()%100)});
      }
      SplineSketch::rebuild(reference);
      auto candidate=reference;
      for(int i=0;i<10;++i) {
        const auto at=1+rng.next()%(candidate.size()-2);
        candidate[at+1].mass+=candidate[at].mass;
        candidate.erase(candidate.begin()+static_cast<std::ptrdiff_t>(at));
      }
      SplineSketch::rebuild(candidate);
      const auto distance=SplineSketch::spline_distance(reference,candidate);
      long double sampled=0;
      assert(std::isfinite(distance));
      for(std::size_t j=1;j<reference.size();++j)for(int i=0;i<=1024;++i) {
        const double q=reference[j-1].x+(reference[j].x-reference[j-1].x)*i/1024.0;
        sampled=std::max(sampled,std::fabs(SplineSketch::spline_rank(reference,q)-
                                        SplineSketch::spline_rank(candidate,q)));
      }
      // Sampling approaches the maximum from below. Allow a small numerical
      // tolerance; this is an independent check, not a certified rounding bound.
      assert(distance+1e-8L>=sampled);
      assert(distance<=sampled+1e-3L);
      assert(SplineSketch::spline_distance(reference,reference)<1e-8L);
    }
  }
};
}

int main() {
  splinesketch::SplineSketchInvariantInspector::check_distance();
  std::cout<<"100 cubic-distance comparisons passed\n";
}
