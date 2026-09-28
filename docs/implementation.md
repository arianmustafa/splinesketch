# Implementation notes

[Back to the README](../README.md)

Details of the data structures, computational costs, numerical behavior, and
differences from the SplineSketch paper. For installation and usage, see the
[README](../README.md).

## Representation and updates

Let `k` be the configured bucket capacity. The sketch partitions represented
observations among three structures:

| Structure | Stored information |
| --- | --- |
| Bucket vector | Up to `k` ordered thresholds, estimated bucket masses, cumulative masses, interpolation slopes, and protection flags. |
| Misra–Gries table | Up to `k - 1` distinct values. Each has a residual counter for eviction and an exact count of the observations currently held for that value. |
| Pending vector | Weighted `(value, count)` entries awaiting consolidation. During streaming updates its size is `O(k)`. |

`add()` first checks the frequency table. An existing entry is incremented;
an absent value is inserted if there is room. If the table is full and the
value is absent, all residual counters are decremented. Entries reaching zero
are evicted, and their held counts are forwarded to the pending vector together
with the new observation. A streaming update triggers consolidation once the
pending vector contains at least `k` entries.

An exact held count does not necessarily equal the value's total frequency:
earlier occurrences may have been evicted and incorporated into buckets.
Rank queries add the spline estimate to the pending and held counts at or
below the query value.

Consolidation sorts and combines pending entries, then updates bucket masses
by walking the sorted thresholds and entries. Existing thresholds reuse their
snapshot prefix sums. Oversized buckets are split at representable midpoints;
adjacent buckets are joined to satisfy capacity. Further splits and joins are
selected using the density-change heuristic. Split counts are estimated from
the snapshot taken before consolidation.

Joins and splits update cumulative masses immediately. Interpolation slopes
are rebuilt after structural edits, before the next snapshot or query uses
them. This avoids repeated slope calculations within a batch.

## Time and memory costs

Let `b` be the bucket count, `p` the pending-entry count, and `h` the frequency
table size; write `e = p + h`.

| Operation | Cost |
| --- | --- |
| Count and capacity accessors | `O(1)` |
| `add()` without a full-table event | Expected `O(1)` hash-table work; inserting a new value can allocate. |
| Full-table event | `O(k)` counter processing, plus consolidation if triggered. |
| `rank()` | `O(log(b + 1) + e)` time; no allocation. |
| `consolidate()` | Up to `O(k^2)` time at streaming sizes, with `O(k)` temporary storage. The sorted mass-update pass is linear after sorting; repeated bucket scans and vector edits account for the quadratic bound. |

Quantiles use at most 64 search steps over ordered IEEE-754 bit patterns.
With at most 16 exact entries, each step calls `rank()` and the query allocates
no memory. With more entries, the query first constructs sorted prefix counts
in `O(e log e)` time and `O(e)` temporary storage. A prepared step costs
`O(log(b + 1) + log e)`, except near a rounding boundary, where it calls
`rank()` to preserve that function's accumulation behavior. Preparation is
repeated for each quantile query.

Calling `consolidate()` before a series of queries removes the pending scan;
it does not remove the frequency-table scan. Consolidation runs synchronously
inside the triggering `add()`, so individual update latency can be much higher
than the average.

Memory use includes vector capacity, hash-table allocations, and temporary
snapshots. Shrinking with `resize()` reduces logical capacity but does not shrink the bucket
or pending vectors' allocated capacity. Merging temporarily holds copies of
both inputs and their combined thresholds. There is no serialization API.

## Numerical behavior and error bounds

Thresholds and results use IEEE-754 binary64 `double`. Observation counts and
exact weights use `std::uint64_t`. Bucket masses, prefixes, slopes, and
interpolation arithmetic use `long double`. Its precision and exponent range
depend on the compiler and platform; where its range equals `double`'s,
subtracting extreme finite endpoints can overflow intermediate calculations.

All integer ranks up to `2^53` are representable in the return type. Above
that count, adjacent integer ranks may round to the same `double`. This limit
concerns numeric representation; bucket approximation error also exists below
it. Conversion of input data to `double` can itself lose distinctions before
the sketch sees the values.

This implementation provides no certified worst-case error bound for a given
`k`. In particular, the capacity-reduction path ignores the bucket-size bound
and can remove protected thresholds if needed. Rebalancing is limited to
`2 * k` oversized-bucket passes and `k` heuristic passes, and buffer processing
does not enforce the batch-size condition used in the paper's proof.
The [theorems in §3.3 and Appendix B](https://arxiv.org/pdf/2504.01206v3#page=12)
therefore require a separate applicability analysis for this code.

## Differences from the paper

References below are to [version 3](https://arxiv.org/pdf/2504.01206v3).
The practical bound factor of 3, its doubling when a required split has no
eligible join, and its reset at an epoch transition follow §4.1.

| Area | Paper | This implementation |
| --- | --- | --- |
| Frequency filtering (§4.2) | Describes MG processing in batches and versions with and without MG. | Processes each observation through MG; filtering is always enabled. |
| Consolidation (§4.1) | Describes heap-based processing and an iteration-based prototype that batches non-overlapping edits. | Selects individual edits through sequential scans of a vector. |
| Quantile inversion (§4) | Locates a bucket, then numerically inverts its interpolation. | Searches the represented value range and includes pending and MG counts at each step. |
| Merge (§3.2) | Combines input buffers before consolidation; discusses capacity selection when inputs differ. | Consolidates copies first and always keeps the receiver's capacity. |
| Final storage (§4.2) | Compacts MG entries and adjusts the bucket budget. | Retains separate bucket and MG capacities; `consolidate()` performs no final storage compaction. |

The sorted buffer sweep follows the approach described in §4.1. Prefix reuse,
deferred slope rebuilds, and per-query preparation reduce repeated work without
changing this implementation's bucket-selection rules or interpolation formulas.
