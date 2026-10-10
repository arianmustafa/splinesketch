"""Exact checks of a reachable shared-state resize witness.

For two attainable true ranks L < U and an old estimate a in (L,U),
any replacement b != a worsens one attainable truth: b>a hurts L,
b<a hurts U. This is a pointwise claim at the demonstrated shared state,
not a theorem that every resize worsens error or every heuristic is useless.

More generally, squared error changes by
    (b-t)^2 - (a-t)^2 = (b-a)*(b+a-2*t).
This is affine in t, so nonpositive values at L and U are necessary and
sufficient for nonworsening on the whole interval [L,U]. If a is inside
that interval, the two inequalities force b=a. C++ checks that this
particular active sketch state can actually arise with either endpoint truth.
The indistinguishability claim concerns deterministic decisions using those
active fields, which are exactly what the current query/resize paths read.
"""
from fractions import Fraction
import json
import struct
import subprocess
import sys


def number(bits):
    return Fraction(struct.unpack('!d', struct.pack('!Q', bits))[0])


def check(row):
    q, a, b = (number(row[key]) for key in ('query_bits', 'old_bits', 'new_bits'))
    # Reconstruct public inputs independently; no truth count from C++ is used
    # to decide which endpoint stream the changed estimate harms.
    streams = []
    for offset in (18, 0):
        values = [Fraction(i * i) for i in range(30)]
        values[1:6] = [Fraction(offset + 3 * i) for i in range(1, 6)]
        streams.append(values)
    low, high = [sum(value <= q for value in values) for values in streams]
    assert q == 18 and (low, high) == (1, 6)
    assert (row['low_truth'], row['high_truth']) == (low, high)
    assert low < a < high
    changes = []
    for truth in (low, high):
        old_error, new_error = abs(a - truth), abs(b - truth)
        squared = (b - truth) ** 2 - (a - truth) ** 2
        assert squared == (b - a) * (b + a - 2 * truth)
        changes.append(new_error - old_error)
        if row['certified']:
            assert row['old_lower'] <= truth <= row['old_upper']
            assert row['new_lower'] <= truth <= row['new_upper']
            assert old_error <= number(row['old_error_bits']) and old_error <= number(row['old_uniform_bits'])
            assert new_error <= number(row['new_error_bits']) and new_error <= number(row['new_uniform_bits'])
    if b > a:
        assert changes[0] > 0
    elif b < a:
        assert changes[1] > 0
    else:
        assert changes == [0, 0]
    assert (changes[0] <= 0 and changes[1] <= 0) == (a == b)
    if row['certified']:
        assert (row['old_lower'], row['old_upper']) == (low, high)
    return max(changes)


def main():
    if len(sys.argv) not in (2, 3):
        raise SystemExit('usage: checks reference-executable [candidate-executable]')
    output = subprocess.check_output([sys.argv[1]], text=True)
    rows = [json.loads(line) for line in output.splitlines()]
    assert len(rows) == 4
    assert {(row['certified'], row['theoretical']) for row in rows} == {(a, b) for a in (0, 1) for b in (0, 1)}
    largest = max(check(row) for row in rows)
    print(f'Four shared-state witnesses checked exactly; largest pointwise increase {float(largest):.17g}')
    if len(sys.argv) == 3:
        candidates = [json.loads(line) for line in subprocess.check_output([sys.argv[2]], text=True).splitlines()]
        assert len(candidates) == 4
        reference = {(row['certified'], row['theoretical']): row for row in rows}
        assert {(row['certified'], row['theoretical']) for row in candidates} == set(reference)
        for candidate in candidates:
            check(candidate)
            key = candidate['certified'], candidate['theoretical']
            row = reference[key]
            assert row['query_bits'] == candidate['query_bits'] and row['old_bits'] == candidate['old_bits']
            a, b = number(row['new_bits']), number(candidate['new_bits'])
            low, high = row['low_truth'], row['high_truth']
            assert low < a < high
            changes = [abs(b - truth) - abs(a - truth) for truth in (low, high)]
            # Both algorithms return one common value for both histories.
            # Thus any different candidate estimate must hurt an attainable
            # endpoint, even though it may improve the other endpoint.
            assert (max(changes) <= 0) == (a == b)
            print(f'certified={key[0]} theoretical={key[1]} candidate error changes '
                  f'at ranks 1 and 6: {float(changes[0]):.17g}, {float(changes[1]):.17g}')


if __name__ == '__main__':
    main()
