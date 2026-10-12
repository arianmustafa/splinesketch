#include <splinesketch/certificate_midpoint.hpp>
#include <splinesketch/certified_splinesketch.hpp>

#ifdef SPLINESKETCH_TESTING
#error Public-header checks must compile without testing macros
#endif

double midpoint_from_other_translation_unit(const splinesketch::CertifiedPaperSplineSketch&);

int main() {
  splinesketch::SplineSketch core(64);
  splinesketch::CertifiedSplineSketch certified(64);
  splinesketch::PaperSplineSketch paper(64);
  splinesketch::CertifiedPaperSplineSketch sketch(64);
  for (int i = 0; i < 33; ++i) {
    core.add(i); certified.add(i); paper.add(i); sketch.add(i);
  }
  const auto rank = splinesketch::certificate_midpoint_rank(sketch, 16);
  const auto selected = splinesketch::certificate_midpoint_select(sketch, 17);
  return !(core.count() == 33 && certified.count() == 33 && paper.count() == 33 &&
           core.quantile(0.5) == 16 && certified.quantile(0.5) == 16 &&
           paper.quantile(0.5) == 16 && sketch.quantile(0.5) == 16 &&
           rank.estimate == 17 && rank.max_error == 0 && selected == 16 &&
           midpoint_from_other_translation_unit(sketch) == selected);
}
