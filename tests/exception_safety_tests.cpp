#include <splinesketch/splinesketch.hpp>

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>

namespace {
int allocations_before_failure = -1;
}

void* operator new(std::size_t size) {
  if (allocations_before_failure == 0) throw std::bad_alloc();
  if (allocations_before_failure > 0) --allocations_before_failure;
  if (void* memory = std::malloc(size ? size : 1)) return memory;
  throw std::bad_alloc();
}

void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void* memory) noexcept { std::free(memory); }
void operator delete[](void* memory) noexcept { std::free(memory); }
void operator delete(void* memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void* memory, std::size_t) noexcept { std::free(memory); }

using splinesketch::SplineSketch;

static void check_same_state(const SplineSketch& actual, const SplineSketch& expected) {
  assert(actual.count() == expected.count());
  assert(actual.bucket_capacity() == expected.bucket_capacity());
  assert(actual.bucket_count() == expected.bucket_count());
  assert(actual.heavy_hitter_count() == expected.heavy_hitter_count());
  for (double x : {-100.0, -1.0, 0.0, 1.0, 5.0, 25.0, 100.0})
    assert(actual.rank(x) == expected.rank(x));
  if (expected.count()) {
    for (double q : {0.0, 0.1, 0.5, 0.9, 1.0})
      assert(actual.quantile(q) == expected.quantile(q));
  }
}

template <class Operation>
static void check_allocation_failures(const SplineSketch& baseline, Operation operation) {
  bool saw_failure = false;
  bool saw_success = false;
  for (int fail_at = 0; fail_at < 256; ++fail_at) {
    SplineSketch sketch = baseline;
    allocations_before_failure = fail_at;
    try {
      operation(sketch);
      allocations_before_failure = -1;
      saw_success = true;
      break;
    } catch (const std::bad_alloc&) {
      allocations_before_failure = -1;
      saw_failure = true;
      check_same_state(sketch, baseline);
      SplineSketch reference = baseline;
      sketch.add(1000.0);
      reference.add(1000.0);
      check_same_state(sketch, reference);
    }
  }
  assert(saw_failure && saw_success);
}

int main() {
  SplineSketch empty(6);
  check_allocation_failures(empty, [](SplineSketch& sketch) { sketch.add(1.0); });

  SplineSketch full(6);
  for (int i = 0; i < 5; ++i) full.add(i);
  check_allocation_failures(full, [](SplineSketch& sketch) { sketch.add(5.0); });

  SplineSketch buffered(6);
  for (int i = 0; i < 5; ++i) buffered.add(i);
  buffered.add(0.0);
  buffered.add(5.0);
  check_allocation_failures(buffered, [](SplineSketch& sketch) { sketch.consolidate(); });

  SplineSketch populated(32);
  for (int i = 0; i < 100; ++i) populated.add(i);
  check_allocation_failures(populated, [](SplineSketch& sketch) { sketch.resize(6); });
  check_allocation_failures(populated, [](SplineSketch& sketch) { sketch.resize(64); });

  std::cout << "Exception safety tests passed\n";
}
