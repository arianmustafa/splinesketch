# SplineSketch

[![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C?logo=cplusplus&logoColor=white)](#installation)
[![Header only](https://img.shields.io/badge/header--only-yes-brightgreen)](include/splinesketch/splinesketch.hpp)
[![Dependencies: standard library](https://img.shields.io/badge/dependencies-standard%20library-blue)](#installation)
[![arXiv: 2504.01206](https://img.shields.io/badge/arXiv-2504.01206-b31b1b)](https://arxiv.org/abs/2504.01206v3)

A header-only C++17 library for estimating ranks and quantiles from streaming
data, without storing every observation.

SplineSketch combines adaptive buckets, cubic spline interpolation, and a
frequency table for repeated values. It is based on the
[SplineSketch paper](https://arxiv.org/abs/2504.01206v3) by Aleksander Łukasiewicz,
Jakub Tětek, and Pavel Veselý.

**Status:** Experimental. This implementation differs from the paper and does
not provide a certified worst-case error bound for a given capacity.

[Installation](#installation) · [Quick start](#quick-start) ·
[API](#api) · [Tests and benchmarks](#tests-and-benchmarks) ·
[Implementation notes](docs/implementation.md)

## Features

- Estimate ranks, medians, and other quantiles from finite `double` values.
- Merge sketches built from separate streams.
- Adjust bucket capacity as storage requirements change.
- Query at any time, including before buffered updates are consolidated.
- Use a single header with no dependencies beyond the C++ standard library.

## Installation

Requires a C++17 compiler and IEEE-754 binary64 `double`. CMake is optional;
the CMake build requires version 3.16 or newer.

### Single header

Copy [`include/splinesketch/splinesketch.hpp`](include/splinesketch/splinesketch.hpp)
into your project's include directory, preserving the `splinesketch/` folder:

```cpp
#include <splinesketch/splinesketch.hpp>
```

### CMake

Place this repository at `third_party/splinesketch`, then add the following
after defining your executable or library target:

```cmake
set(SPLINESKETCH_BUILD_TESTS OFF CACHE BOOL "" FORCE)
add_subdirectory(third_party/splinesketch)
target_link_libraries(your_target PRIVATE splinesketch)
```

## Quick start

```cpp
#include <splinesketch/splinesketch.hpp>
#include <iostream>

int main() {
    splinesketch::SplineSketch sketch(128);
    const double values[] = {1.0, 2.0, 2.0, 4.0, 8.0};
    for (double value : values) sketch.add(value);

    std::cout << sketch.count() << '\n';       // 5 observations
    std::cout << sketch.rank(2.0) << '\n';     // 3 observations <= 2
    std::cout << sketch.quantile(0.5) << '\n'; // median: 2
}
```

Save this as `example.cpp`. From the repository root, compile and run it with:

```sh
c++ -std=c++17 -O2 -I include example.cpp -o example
./example
```

The example is exact because all distinct values fit in the frequency table.
Once observations are incorporated into buckets, queries generally become
approximate. `rank()` returns a count; divide by `count()` for an estimated
cumulative probability when the sketch is nonempty.

### Merge and resize

```cpp
splinesketch::SplineSketch left(128), right(128);
left.add(1.0);
left.add(2.0);
right.add(3.0);
right.add(4.0);

left.merge(right); // left now represents all four observations
left.resize(64);   // reduce the bucket capacity
```

Merging retains the receiver's capacity and leaves a distinct source sketch
unchanged. Merge order can affect estimates. Increasing capacity later cannot
recover detail already lost to approximation.

## API

All methods belong to `splinesketch::SplineSketch`. The constructor's bucket
capacity defaults to `128` and must be between `6` and
`std::numeric_limits<std::size_t>::max() / 4`. Capacity controls storage and
bucket resolution; it is not an error tolerance.

| Method | Description |
| --- | --- |
| `add(double x)` | Insert one finite observation. Repeated values count separately. |
| `count() const` | Return the observation count as `std::uint64_t`. |
| `rank(double x) const` | Estimate the number of observations `<= x`, as `double` in `[0, count()]`. Return zero for an empty sketch. |
| `quantile(double q) const` | Estimate a quantile for `q` in `[0, 1]`. Endpoints return the minimum and maximum. Interior queries search for a value whose estimated rank reaches `ceil(q * count())`; the result need not be an observed value. |
| `consolidate()` | Incorporate pending entries into buckets, leaving the frequency table intact. Do nothing if there are no pending entries. |
| `merge(const SplineSketch& other)` | Combine observations into this sketch. Self-merge doubles the represented observations. |
| `resize(std::size_t k)` | Change capacity and adjust storage structures. Even resizing to the same capacity first consolidates pending entries. |
| `bucket_capacity() const` | Return the configured capacity. |
| `bucket_count() const` | Return the current number of bucket nodes, which can be less than the capacity or zero. |
| `heavy_hitter_count() const` | Return the number of entries in the frequency table. |

### Errors and thread safety

| Condition | Behavior |
| --- | --- |
| Non-finite input to `add()`, NaN input to `rank()`, invalid probability or capacity | Throw `std::invalid_argument`. |
| Valid quantile query on an empty sketch | Throw `std::logic_error`. |
| Observation count would exceed `UINT64_MAX` during addition or merge | Throw `std::overflow_error`. |
| `rank(-infinity)` / `rank(+infinity)` | Return zero / the total count converted to `double`. |

Allocating operations can also throw allocation exceptions. Failed mutating
operations leave the sketch unchanged. Inputs `-0.0` and
`+0.0` are treated as the same value. Const queries can run concurrently on an
otherwise unmodified sketch. Mutation requires external synchronization with
all other accesses to that sketch.

## Accuracy and limitations

- There is no certified worst-case error bound for this implementation.
  See the [differences from the paper](docs/implementation.md#differences-from-the-paper).
- Integer ranks above `2^53` may round to the same `double` value.
- Intermediate arithmetic uses `long double`. On platforms where its range
  equals `double`'s, extreme finite inputs can overflow intermediate calculations.
- Consolidation runs synchronously during some updates, so individual update
  latency can be much higher than the average.
- Resizing downward reduces logical capacity but does not shrink the allocated
  bucket or pending vectors. There is no serialization API.

See [implementation notes](docs/implementation.md) for complexity, memory costs,
numerical behavior, and the scope of the paper's guarantees.

## Tests and benchmarks

Build and run the tests:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
(cd build && ctest -C Release --output-on-failure)
```

Tests are enabled by default. The [test suite](tests/splinesketch_tests.cpp)
covers exact small inputs, monotonicity, quantile inversion, extreme values,
buffered updates, large weighted counts, merges, and resizing. Accuracy checks
use fixed synthetic streams and do not establish a general error bound.

To also build the benchmarks:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    -DSPLINESKETCH_BUILD_BENCHMARK=ON
cmake --build build --config Release

./build/splinesketch_benchmark
./build/splinesketch_criteria_benchmark 128
```

With a multi-configuration generator, executables are in `build/Release/`.

| Benchmark | Measures |
| --- | --- |
| [Throughput](benchmarks/throughput.cpp) | One million normal samples at capacity 128, followed by 100,000 rank and 10,000 quantile queries. Reports microseconds per operation; update timing includes a final consolidation. |
| [Criteria](benchmarks/criteria.cpp) | 10,000 parts, each updating 30 sketches. Reports mean, p99, p99.9, and maximum milliseconds per part. The optional argument sets capacity, defaulting to 128. Includes automatic consolidations, with no queries or final flush. |

Both benchmarks generate input before timing. Throughput queries use the
consolidated sketch with its frequency table retained. These are local synthetic
workloads, not reproductions of the paper's accuracy or storage comparisons.

## References

Aleksander Łukasiewicz, Jakub Tětek, and Pavel Veselý.
[*SplineSketch: Even More Accurate Quantiles with Error Guarantees*](https://arxiv.org/abs/2504.01206v3),
arXiv:2504.01206v3.

The [implementation notes](docs/implementation.md#differences-from-the-paper)
describe how this library differs from the algorithm and prototype in the paper.

## License

Licensed under the MIT License. See [LICENSE](LICENSE) for the full text.
