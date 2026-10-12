#include <splinesketch/certified_splinesketch.hpp>
#include <splinesketch/paper_splinesketch.hpp>
#include <splinesketch/certificate_midpoint.hpp>

#ifdef SPLINESKETCH_TESTING
#error Public-header checks must compile without testing macros
#endif

double midpoint_from_other_translation_unit(const splinesketch::CertifiedPaperSplineSketch& sketch) {
  return splinesketch::certificate_midpoint_select(sketch, 17);
}
