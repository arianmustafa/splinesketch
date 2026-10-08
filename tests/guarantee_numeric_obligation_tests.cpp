#include <splinesketch/splinesketch.hpp>

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <limits>
#include <stdexcept>
#include <string>

// Built with -mlong-double-64 on supported GCC/x86-64 targets. All six inputs
// are finite; the verification guard must catch an overflowing adjacent span.
int main() {
  static_assert(std::numeric_limits<long double>::max_exponent ==
                    std::numeric_limits<double>::max_exponent,
                "numeric obligation test requires binary64 long double range");
  splinesketch::SplineSketch sketch(6);
  const double values[] = {-std::numeric_limits<double>::max(), -1.5e308,
                           -1e308, 1e308, 1.5e308,
                           std::numeric_limits<double>::max()};
  for (int i = 0; i < 5; ++i) sketch.add(values[i]);
  bool caught = false;
  try {
    sketch.add(values[5]);
  } catch (const std::logic_error& error) {
    caught = std::string(error.what()) ==
             "finite endpoints overflow long double span";
  }
  assert(caught);
  assert(sketch.count() == 5);
}
