#!/usr/bin/env python3
"""Exact rational oracle for the experimental certificate minimax query."""
from fractions import Fraction
import json
import math
import random
import struct
import subprocess
import sys


MAX_COUNT = (1 << 64) - 1


def unpack(bits):
    return struct.unpack('>d', struct.pack('>Q', bits))[0]


def fraction(bits):
    value = unpack(bits)
    assert math.isfinite(value), ('nonfinite result', bits)
    return Fraction.from_float(value)


def ceiling(value):
    result = float(value)
    if Fraction.from_float(result) < value:
        result = math.nextafter(result, math.inf)
    return Fraction.from_float(result)


def intervals():
    cases = {(lower, upper) for lower in range(65) for upper in range(lower, 65)}

    def add(a, b):
        if 0 <= a <= b <= MAX_COUNT:
            cases.add((a, b))

    for exponent in range(65):
        power = 1 << exponent
        spacing = 1 << max(0, exponent - 52)
        # Around every power of two, rounding ties, and both parities of the
        # retained significand; include half-integer midpoints and singletons.
        offsets = {-5, -2, -1, 0, 1, 2, 5, -spacing, spacing,
                   spacing // 2 - 1, spacing // 2, spacing // 2 + 1,
                   3 * spacing // 2 - 1, 3 * spacing // 2, 3 * spacing // 2 + 1}
        endpoints = sorted(power + delta for delta in offsets)
        for a in endpoints:
            for b in endpoints:
                add(a, b)
            add(0, a)
            add(a, MAX_COUNT)
    for offset in range(4097):
        add(MAX_COUNT - offset, MAX_COUNT)
        add(MAX_COUNT - offset, MAX_COUNT - offset)
    randomizer = random.Random(227)
    for _ in range(10000):
        a, b = sorted((randomizer.getrandbits(64), randomizer.getrandbits(64)))
        add(a, b)
    return sorted(cases)


def check(row, expected_interval):
    lower, upper = row['lower'], row['upper']
    assert (lower, upper) == expected_interval, ('changed input', row, expected_interval)
    midpoint = Fraction(lower + upper, 2)
    nearest = Fraction.from_float(float(midpoint))
    exact_radius = max(abs(nearest - lower), abs(nearest - upper))
    rounded_radius = ceiling(exact_radius)
    # Independent neighboring representable estimates verify the minimax
    # choice and recognize both optimal answers if the real midpoint is a tie.
    alternatives = [math.nextafter(float(midpoint), direction) for direction in (-math.inf, math.inf)]
    for alternative in alternatives:
        alternative = Fraction.from_float(alternative)
        assert exact_radius <= max(abs(alternative - lower), abs(alternative - upper))
    allowance = Fraction(0) if upper <= 1 << 52 else Fraction(2) ** (upper.bit_length() - 54)
    real_uniform = Fraction(upper - lower, 2) + allowance
    assert len(row['modes']) == 4
    for mode in row['modes']:
        estimate, radius, uniform = map(fraction, mode)
        assert estimate == nearest, ('nearest midpoint', lower, upper, mode, nearest)
        assert radius == rounded_radius, ('outward radius', lower, upper, mode, exact_radius)
        assert radius >= exact_radius and uniform >= real_uniform and uniform >= exact_radius, (
            'insufficient allowance', lower, upper, mode, real_uniform)
        assert exact_radius == Fraction(upper - lower, 2) + abs(estimate - midpoint), ('radius identity', row)
        # Returning an unnecessarily huge allowance would conceal errors.
        spacing = Fraction.from_float(math.ulp(float(real_uniform)))
        assert uniform <= ceiling(real_uniform) + 2 * spacing, ('loose allowance', row)
    assert len({tuple(mode[:2]) for mode in row['modes']}) == 1, ('rounding-mode dependence', row)


def main():
    assert len(sys.argv) == 2, 'usage: midpoint-checks executable'
    executable = sys.argv[1]
    history = json.loads(subprocess.check_output([executable], text=True))
    assert history['partition_queries'] > 10000, ('missing complete certificate partitions', history)
    cases = intervals()
    completed = subprocess.run([executable, '--intervals'], check=True, text=True, capture_output=True,
                               input=''.join(f'{lower} {upper}\n' for lower, upper in cases))
    rows = [json.loads(line) for line in completed.stdout.splitlines()]
    assert len(rows) == len(cases), ('missing output', len(rows), len(cases))
    for expected, row in zip(cases, rows):
        check(row, expected)
    ordered_pairs = 0
    for keys in (('lower', 'upper'), ('upper', 'lower')):
        previous = None
        for row in sorted(rows, key=lambda row: tuple(row[key] for key in keys)):
            if previous is not None and previous['lower'] <= row['lower'] and previous['upper'] <= row['upper']:
                for old_mode, new_mode in zip(previous['modes'], row['modes']):
                    assert fraction(old_mode[0]) <= fraction(new_mode[0]), ('rounded midpoint decreased', previous, row)
                ordered_pairs += 1
            previous = row
    assert ordered_pairs > 10000, ('missing ordered intervals', ordered_pairs)
    large_rows = [json.loads(line) for line in subprocess.check_output(
        [executable, '--large-counts'], text=True).splitlines()]
    assert len(large_rows) == 2 * 13 * 3
    assert {(row['theoretical'], row['step'], row['query']) for row in large_rows} == {
        (theory, step, query) for theory in (0, 1) for step in range(51, 64) for query in (-1, 0, 1)}
    previous = None
    for row in large_rows:
        count = (1 << (row['step'] + 1)) - 1
        zeros = 1
        for step in range(1, row['step'] + 1):
            zeros = 2 * zeros + (step % 2 == 0)
        truth = 0 if row['query'] < 0 else zeros if row['query'] == 0 else count
        assert (row['count'], row['truth']) == (count, truth), ('large history reconstruction', row)
        assert row['lower'] <= truth <= row['upper'], ('large certificate', row)
        estimate, radius, uniform = map(fraction, (row['estimate_bits'], row['radius_bits'], row['uniform_bits']))
        midpoint = Fraction(row['lower'] + row['upper'], 2)
        assert estimate == Fraction.from_float(float(midpoint)), ('large midpoint', row)
        exact_radius = max(abs(estimate - row['lower']), abs(estimate - row['upper']))
        assert radius == ceiling(exact_radius), ('large radius', row)
        assert abs(estimate - truth) <= radius and abs(estimate - truth) <= uniform, ('large error bound', row)
        if previous is not None and (row['theoretical'], row['step']) == (previous['theoretical'], previous['step']):
            assert row['query'] > previous['query'] and estimate >= fraction(previous['estimate_bits']), (
                'large public midpoint decreased', previous, row)
        previous = row
    # Reachable endpoint witness: the certificate minimax radius is attained
    # by both streams, rather than merely by arbitrary real endpoints.
    truths = []
    for high in (False, True):
        data = [i * i for i in range(30)]
        data[1:6] = [(0 if high else 18) + 3 * i for i in range(1, 6)]
        truths.append(sum(x <= 18 for x in data))
    assert truths == [1, 6]
    assert [abs(Fraction(7, 2) - truth) for truth in truths] == [Fraction(5, 2)] * 2
    print(json.dumps({'intervals': len(cases), 'rounding_modes': 4, 'history': history,
                      'ordered_interval_pairs': ordered_pairs,
                      'large_count_queries': len(large_rows),
                      'reachable_minimax_radius': '5/2'}, sort_keys=True))


if __name__ == '__main__':
    main()
