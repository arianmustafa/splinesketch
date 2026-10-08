#include <splinesketch/splinesketch.hpp>

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

// Compiled only with SPLINESKETCH_VERIFY_GUARANTEES. This stream previously
// chose a pair rejected by the paper's Definition 1, despite eligible pairs.
// It now checks that capacity reduction chooses an eligible pair.
int main() {
  const int values[] = {27, 36, 61, 24, 21, 59, 3, 91, 56, 79, 51, 4,
                        50, 86, 82, 50, 35, 55, 34, 67, 46, 52, 7, 64,
                        60, 17, 7, 63, 30, 54, 98, 54, 51, 12, 91, 50};
  splinesketch::SplineSketch sketch(6);
  for (int i = 0; i < 35; ++i) sketch.add(values[i]);
  assert(sketch.count() == 35);
  sketch.add(values[35]);  // automatic consolidation
  assert(sketch.count() == 36);
}
