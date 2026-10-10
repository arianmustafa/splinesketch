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
    for target, reached in row['targets']:
        assert bool(reached) == (lower + upper >= 2 * target), ('integer midpoint comparison', row, target)


def check_selection(row):
    data = [(Fraction.from_float(unpack(value)), weight) for value, weight in row['data']]
    assert sum(weight for _, weight in data) == row['count'], ('selection count', row)
    if row['case'].startswith('large-'):
        step = int(row['case'].split('-')[1])
        zeros = 1
        for stage in range(1, step + 1):
            zeros = 2 * zeros + (stage % 2 == 0)
        assert row['count'] == (1 << (step + 1)) - 1
        assert data == [(0, zeros), (1, row['count'] - zeros)]
    if row['case'].startswith('rounded-terminal-'):
        assert data == [(7, (1 << 53) + 1)]

    def rank(value):
        if math.isinf(value):
            return row['count'] if value > 0 else 0
        value = Fraction.from_float(value)
        return sum(weight for x, weight in data if x <= value)

    profile = [(unpack(value), lower, upper) for value, lower, upper in row['profile']]
    assert profile[0] == (-math.inf, 0, 0) and profile[-1] == (math.inf, row['count'], row['count'])
    for value, lower, upper in profile:
        assert lower <= rank(value) <= upper, ('selection profile certificate', row['case'], value)
    previous_selection = -math.inf
    for selection in row['selections']:
        target = selection['target']
        # Linear cell walk with exact fractions, independent of the C++
        # binary search, ordered-bit mapping and integer predicate.
        expected = next(value for value, lower, upper in profile
                        if math.isfinite(value) and Fraction(lower + upper, 2) >= target)
        value = unpack(selection['value_bits'])
        assert value == expected and math.isfinite(value), ('first selection crossing', row['case'], selection, expected)
        assert value != 0 or selection['value_bits'] == 0, ('noncanonical zero', selection)
        assert previous_selection <= value, ('selection target order', row['case'], selection)
        previous_selection = value
        assert len(selection['modes']) == 4 and set(selection['modes']) == {selection['value_bits']}
        midpoint_at = Fraction(selection['lower'] + selection['upper'], 2)
        midpoint_before = Fraction(selection['before_lower'] + selection['before_upper'], 2)
        assert midpoint_before < target <= midpoint_at, ('exact target crossing', row['case'], selection)
        before = rank(math.nextafter(value, -math.inf))
        at = rank(value)
        assert selection['before_lower'] <= before <= selection['before_upper']
        assert selection['lower'] <= at <= selection['upper']
        allowance = Fraction(row['width'], 2)
        assert at >= target - allowance and before < target + allowance, ('rank bracket', row['case'], selection)
        distance = max(before - target, target - at, 0)
        assert distance <= row['width'] // 2
        if row['case'].startswith('tight-'):
            assert row['width'] == 5 and distance == row['width'] // 2 == 2
    if row['case'].startswith('duplicates-'):
        assert row['width'] == 0 and row['selections'][0]['target'] == 1 and rank(7) - 1 == 98
    return len(row['selections'])


def main():
    assert len(sys.argv) == 2, 'usage: midpoint-checks executable'
    executable = sys.argv[1]
    history = json.loads(subprocess.check_output([executable], text=True))
    assert history['partition_queries'] > 10000, ('missing complete certificate partitions', history)
    assert history['selections'] > 10000, ('missing history selections', history)
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
    selection_rows = [json.loads(line) for line in subprocess.check_output(
        [executable, '--selection-fixtures'], text=True).splitlines()]
    expected_cases = {f'{case}-{policy}' for policy in ('practical', 'theoretical')
                      for case in ('cluster', 'tight', 'duplicates', 'rounded-terminal', 'raw', 'finalized', 'large-53', 'large-63')}
    assert len(selection_rows) == len(expected_cases) and {row['case'] for row in selection_rows} == expected_cases
    selection_checks = sum(check_selection(row) for row in selection_rows)
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
                      'selection_fixtures': len(selection_rows), 'exact_selection_checks': selection_checks,
                      'reachable_minimax_radius': '5/2'}, sort_keys=True))


if __name__ == '__main__':
    main()
