#if defined(SPLINESKETCH_PAPER) || defined(SPLINESKETCH_PAPER_CERTIFIED)
#include <splinesketch/paper_splinesketch.hpp>
#elif defined(SPLINESKETCH_CERTIFIED)
#include <splinesketch/certified_splinesketch.hpp>
#else
#include <splinesketch/splinesketch.hpp>
#endif

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

#ifdef SPLINESKETCH_PAPER_CERTIFIED
using SplineSketch = splinesketch::CertifiedPaperSplineSketch;
#elif defined(SPLINESKETCH_PAPER)
using SplineSketch = splinesketch::PaperSplineSketch;
#elif defined(SPLINESKETCH_CERTIFIED)
using SplineSketch = splinesketch::CertifiedSplineSketch;
#else
using splinesketch::SplineSketch;
#endif


static void check_same_state(const SplineSketch& actual, const SplineSketch& expected) {
  assert(actual.count() == expected.count());
  assert(actual.bucket_capacity() == expected.bucket_capacity());
  assert(actual.bucket_count() == expected.bucket_count());
  assert(actual.heavy_hitter_count() == expected.heavy_hitter_count());
#if defined(SPLINESKETCH_PAPER) || defined(SPLINESKETCH_PAPER_CERTIFIED)
  assert(actual.finalized() == expected.finalized());
#endif
#ifdef SPLINESKETCH_CERTIFIED
  assert(actual.certificate_capacity() == expected.certificate_capacity());
  assert(actual.certificate_size() == expected.certificate_size());
  assert(actual.max_rank_uncertainty() == expected.max_rank_uncertainty());
#endif
#if defined(SPLINESKETCH_CERTIFIED) || defined(SPLINESKETCH_PAPER_CERTIFIED)
  assert(actual.max_rank_uncertainty() == expected.max_rank_uncertainty());
  assert(actual.max_rank_error() == expected.max_rank_error());
  for (double x : {-100.0, -1.0, 0.0, 1.0, 5.0, 25.0, 100.0}) {
    const auto a = actual.rank_with_error(x), b = expected.rank_with_error(x);
    assert(a.lower_rank == b.lower_rank && a.upper_rank == b.upper_rank);
    assert(a.max_error == b.max_error);
  }
#endif
  for (double x : {-100.0, -1.0, 0.0, 1.0, 5.0, 25.0, 100.0})
    assert(actual.rank(x) == expected.rank(x));
  if (expected.count()) {
    for (double q : {0.0, 0.1, 0.5, 0.9, 1.0})
      assert(actual.quantile(q) == expected.quantile(q));
  }
}

template <class Operation>
static void check_allocation_failures(const SplineSketch& baseline, Operation operation,
                                      const SplineSketch* expected_success = nullptr) {
  bool saw_failure = false;
  bool saw_success = false;
#if defined(SPLINESKETCH_PAPER) || defined(SPLINESKETCH_PAPER_CERTIFIED)
  constexpr int failure_limit = 1024; // weighted MG plus all legal resize joins
#else
  constexpr int failure_limit = 256;
#endif
  for (int fail_at = 0; fail_at < failure_limit; ++fail_at) {
    SplineSketch sketch = baseline;
    allocations_before_failure = fail_at;
    try {
      operation(sketch);
      allocations_before_failure = -1;
      if (expected_success) check_same_state(sketch, *expected_success);
      saw_success = true;
      break;
    } catch (const std::bad_alloc&) {
      allocations_before_failure = -1;
      saw_failure = true;
      check_same_state(sketch, baseline);
#if defined(SPLINESKETCH_PAPER) || defined(SPLINESKETCH_PAPER_CERTIFIED)
      if (baseline.finalized()) continue;
#endif
      SplineSketch reference = baseline;
      sketch.add(1000.0);
      reference.add(1000.0);
      check_same_state(sketch, reference);
    }
  }
  assert(saw_failure && saw_success);
}

static void copy_assignment() {
  SplineSketch source(128), destination(6);
  for (int i = 0; i < 1000; ++i) source.add(i);
  source.consolidate();
  source.add(25.5); // Include exact data alongside consolidated buckets.
  destination.add(-1);
  const SplineSketch original_source(source);
  check_allocation_failures(destination, [&](SplineSketch& sketch) {
    auto& assigned = (sketch = source);
    assert(&assigned == &sketch);
  }, &source);
  check_allocation_failures(source, [&](SplineSketch& sketch) { sketch = destination; },
                            &destination);
  check_same_state(source, original_source);

  SplineSketch self(source);
  const SplineSketch& alias = self;
  allocations_before_failure = 0;
  self = alias; // Self-assignment must preserve state without allocating.
  allocations_before_failure = -1;
  check_same_state(self, source);

  SplineSketch empty(32);
  destination = empty;
  check_same_state(destination, empty);
#if defined(SPLINESKETCH_PAPER) || defined(SPLINESKETCH_PAPER_CERTIFIED)
  SplineSketch finalized(source);
  finalized.finalize();
  check_allocation_failures(finalized, [&](SplineSketch& sketch) { sketch = source; },
                            &source);
  check_allocation_failures(source, [&](SplineSketch& sketch) { sketch = finalized; },
                            &finalized);
#endif
}

int main() {
  copy_assignment();
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

  check_allocation_failures(populated, [&](SplineSketch& sketch) { sketch.merge(populated); });
#if defined(SPLINESKETCH_PAPER) || defined(SPLINESKETCH_PAPER_CERTIFIED)
  check_allocation_failures(populated, [](SplineSketch& sketch) { sketch.finalize(); });
#endif
#ifdef SPLINESKETCH_CERTIFIED
  check_allocation_failures(populated, [](SplineSketch& sketch) { sketch.add(500.0); });
  // Cross the packed counter limit without inserting billions of items.
  SplineSketch narrow_maximum(6), power(6);
  power.add(0);
  for (unsigned bit = 0; bit < 32; ++bit) {
    narrow_maximum.merge(power);
    if (bit != 31) power.merge(power);
  }
  assert(narrow_maximum.count() == std::numeric_limits<std::uint32_t>::max());
  check_allocation_failures(narrow_maximum, [](SplineSketch& sketch) { sketch.add(0); });
  check_allocation_failures(narrow_maximum, [&](SplineSketch& sketch) { sketch.merge(power); });
#endif
  std::cout << "Exception safety tests passed\n";
}
