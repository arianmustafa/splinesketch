# Reproducible rank accuracy results

[Back to the README](../README.md)

The [accuracy benchmark](../benchmarks/accuracy.cpp) compares `rank(x)` with
the exact count of observations `<= x`, found with `std::upper_bound` on the
sorted input. These are measurements of this implementation, not a certified
error bound or a reproduction of the paper's experiments.

## Reproduce

From the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSPLINESKETCH_BUILD_BENCHMARK=ON
cmake --build build --config Release --target splinesketch_accuracy_benchmark
./build/splinesketch_accuracy_benchmark
./build/splinesketch_accuracy_benchmark --details
```

The accuracy runner is also part of CTest when tests are enabled. Its passing
status checks count preservation and valid estimated ranks; the numeric error
statistics are observational results rather than pass/fail guarantees.

With a multi-configuration generator, use
`./build/Release/splinesketch_accuracy_benchmark`. Default output is CSV
summarized by capacity, shape, and workflow. `--details` also prints every
capacity/shape/workflow combination. The worst-case query and its exact and
estimated ranks are printed to standard error.

The run below used GCC 11.4.0, C++17, Release (`-O3 -DNDEBUG`) on x86-64 Linux,
where `long double` has a 64-bit significand. The source uses SplitMix64 with
three fixed seeds (`0`, `1`, `2` in its seed formula), integer-generated input,
and no clock-dependent data. Floating-point results can still vary across
compilers and `long double` formats.

## Workload and metric

- **Capacities:** 8, 16, 32, 64, and 128 buckets.
- **Data:** 6,000 observations per case. Uniform values span approximately
  `[-500, 500]`; clustered values sum four discrete uniform draws; duplicates
  draw 70% of values from three repeated values; outliers put about 2% of
  observations around `±10,000` and the rest in `[-1, 1]`; ascending data
  arrive in increasing order. Generation formulas are in the benchmark source.
- **Workflows:** `direct` inserts into one sketch; `merge4` builds four
  interleaved sketches and merges them in a balanced tree; `resize` starts at
  twice the reported capacity, shrinks to `max(6, capacity/2)` halfway through,
  then grows to the reported capacity. Every final sketch is consolidated.
- **Queries:** For each sorted stream, 257 equally spaced order-statistic
  positions, the immediately preceding representable value for each, midpoints
  to the next distinct value, and 257 evenly spaced values across the observed
  range. The values immediately outside the observed range are included.
  Duplicate query values are removed. There were 187,020 evaluated queries
  across 225 cases (5 capacities × 5 shapes × 3 workflows × 3 seeds).

For each query, error is `abs(estimated rank - exact rank) / 6000`.
The tables express this as **percent of stream size**. “Median” and “p95”
are empirical query-error percentiles; “max” is the largest error at any
tested query. Repeated queries across capacities and workflows contribute
separately. These statistics describe the selected workload and query grid,
not all possible streams or query values.

## Observed results

| Capacity | Median error | p95 error | Max error |
| ---: | ---: | ---: | ---: |
| 8 | 0.182% | 4.250% | 27.923% |
| 16 | 0.091% | 1.208% | 20.072% |
| 32 | 0.062% | 0.395% | 6.160% |
| 64 | 0.042% | 0.239% | 2.328% |
| 128 | 0.030% | 0.154% | 0.975% |

| Data shape | Median error | p95 error | Max error |
| --- | ---: | ---: | ---: |
| Uniform | 0.081% | 0.413% | 0.863% |
| Clustered | 0.085% | 0.723% | 2.623% |
| Duplicates | 0.099% | 0.355% | 0.794% |
| Outliers | 0.148% | 4.753% | 27.923% |
| Ascending | 0.011% | 0.043% | 0.435% |

| Workflow | Median error | p95 error | Max error |
| --- | ---: | ---: | ---: |
| Direct | 0.044% | 0.580% | 27.923% |
| Four-way merge | 0.066% | 0.752% | 20.715% |
| Resize | 0.063% | 1.464% | 20.072% |

Across all queries, the median was **0.057%**, p95 was **0.818%**, and the
worst observed error was **27.923%**. That worst query was in the capacity-8,
outlier, direct case with seed 1: at `x = 7.5`, the exact rank was 5,945 and
the estimate was 4,269.638242813302. Large gaps between central values and
outliers are especially difficult at low capacity in this workload.

The [implementation notes](implementation.md#numerical-behavior-and-error-bounds)
explain why the paper's guarantee has not been established for this code.
