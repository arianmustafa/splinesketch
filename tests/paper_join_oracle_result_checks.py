#!/usr/bin/env python3
"""Cross-check saved oracle results, exact objectives and source provenance."""
import argparse
import csv
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import re
import struct

ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def error(row):
    return Fraction.from_float(float(row['error_high'])) + Fraction.from_float(float(row['error_low']))


def binary64(bits):
    return Fraction.from_float(struct.unpack('>d', int(bits).to_bytes(8, 'big'))[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--supplemental', action='store_true', help='Require saved plain-alias and sanitizer checks too.')
    args = parser.parse_args()
    path = args.results
    metadata = json.loads((path / 'metadata.json').read_text())
    for name, saved in metadata['source_sha256'].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == saved, name
    for name, saved in metadata['input_sha256'].items():
        source = ROOT / 'tests/data/paper_width_counterexample.hpp' if name == 'witness_header' else path / name
        assert hashlib.sha256(source.read_bytes()).hexdigest() == saved, name
    report = {'source_hashes_verified': len(metadata['source_sha256']), 'oracle_arms': 0,
              'oracle_certificate_queries': 0, 'exact_rational_queries': 0,
              'original_choice_controls': 0, 'failed_arms': 0,
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    for key, summary in metadata['oracle_summary'].items():
        data = rows(path / f'{key}.csv')
        baseline, final = data[0], data[-1]
        assert baseline['phase'] == 'baseline' and final['phase'] == 'oracle'
        queries = int(baseline['queries'])
        events = rows(path / f'{key}-events.csv')
        singles = [r for r in data if r['phase'] == 'single']
        assert len(singles) == summary['single_arms']
        for event in events:
            arms = [r for r in singles if r['event'] == event['event']]
            assert len(arms) == int(event['legal_choices'])
            assert len({r['choice_bits'] for r in arms}) == len(arms)
            original = [r for r in arms if r['is_original'] == '1']
            assert len(original) == 1 and original[0]['success'] == '1'
            assert original[0]['query_hash'] == baseline['query_hash']
            assert error(original[0]) == error(baseline)
            report['original_choice_controls'] += 1
        # Reconstruct accepted lookahead objectives from every arm in order.
        current = error(baseline)
        lookahead = [r for r in data if r['phase'] == 'lookahead']
        for _, arms in itertools.groupby(lookahead, key=lambda r: r['event']):
            current = min([current] + [error(r) for r in arms if r['success'] == '1'])
        assert current == error(final)
        assert error(final) <= error(baseline)
        assert len((path / f'{key}-plan.txt').read_text().splitlines()) == int(final['event'])
        cdf = rows(path / f'cdf-{key}.csv')
        assert len(cdf) == 1 and int(cdf[0]['queries']) == queries
        for row in data + cdf:
            assert int(row['queries']) == queries
            if row['success'] != '1':
                report['failed_arms'] += 1
                continue
            exact = abs(binary64(row['estimate_bits']) - int(row['truth']))
            assert exact == error(row), (key, row)
            assert exact <= Fraction.from_float(float(row['uniform']))
        report['oracle_arms'] += len(singles) + len(lookahead)
        diagnostic = (path / f'{key}.stderr').read_text().strip()
        match = re.fullmatch(r'Oracle: (\d+) baseline joins; (\d+) changed choices; (\d+) exact-error certificate checks', diagnostic)
        assert match and int(match[1]) == len(events) and int(match[2]) == int(final['event'])
        # Final is a reused result; the unhooked control adds one replay.
        assert int(match[3]) == sum(queries for r in data if r['success'] == '1')
        report['oracle_certificate_queries'] += int(match[3])
    exact_files = []
    for key in metadata['oracle_summary']:
        exact_files.extend(path / f'exact-{key}-{label}.txt' for label in ('baseline', 'single', 'lookahead'))
        exact_files.append(path / f'exact-cdf-{key}.txt')
    for file in exact_files:
        match = re.fullmatch(r'Oracle: (\d+) exact rational error/subtraction checks passed\n', file.read_text())
        assert match, file
        report['exact_rational_queries'] += int(match[1])
    report['heap_exact_rational_queries'] = 0
    for file in path.glob('heap-exact-cdf-*.txt'):
        match = re.fullmatch(r'Paper width: (\d+) exact rational query/error checks passed\n', file.read_text())
        assert match, file
        report['heap_exact_rational_queries'] += int(match[1])
    before, after = rows(path / 'quality-fresh-reference-native.csv'), rows(path / 'quality-fresh-cdf-native.csv')
    assert hashlib.sha256((path / 'quality-fresh-reference-native.csv').read_bytes()).hexdigest() == metadata['reference_quality']['sha256']
    provenance = json.loads((path / 'quality-fresh-cdf-native.csv.provenance.json').read_text())
    assert hashlib.sha256((path / 'quality-fresh-cdf-native.csv').read_bytes()).hexdigest() == provenance['csv_sha256']
    assert provenance['inputs']['compiler'] == metadata['compiler']
    assert provenance['inputs']['source_sha256']['1'] == metadata['instrumented_header_sha256']['cdf']
    assert len(before) == len(after) == 300
    assert all(all(a[f] == b[f] for f in ('k', 'shape', 'seed', 'workflow', 'n', 'queries', 'raw_error_max_percent'))
               for a, b in zip(before, after))
    report['fresh_quality_states'] = len(after)
    report['fresh_quality_queries'] = sum(int(r['queries']) for r in after)
    for file in path.glob('*.stderr'):
        if file.name.removesuffix('.stderr') in metadata['oracle_summary']:
            continue
        assert not file.read_text(), file
    if args.supplemental:
        report.update(plain_cases=0, plain_rank_queries=0, plain_inverse_queries=0)
        for mode in ('native', 'binary64'):
            match = re.fullmatch(r'Exact ranks and quantiles match: (\d+) cases, (\d+) rank queries, (\d+) inverse queries\n',
                                 (path / f'plain-cdf-{mode}.txt').read_text())
            assert match, mode
            for field, value in zip(('plain_cases', 'plain_rank_queries', 'plain_inverse_queries'), match.groups()):
                report[field] += int(value)
        assert (path / 'sanitizer-oracle.txt').read_text() == (path / 'witness-native.stderr').read_text()
        assert re.fullmatch(r'Paper width: 77912 query checks; split/join/merge inequalities passed\n',
                            (path / 'sanitizer-cdf.txt').read_text())
        report['address_and_undefined_sanitizers'] = 'passed: full witness oracle and forced-heap CDF width suite'
    report['status'] = 'passed'
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
