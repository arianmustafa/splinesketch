#!/usr/bin/env python3
"""Reproduce the certified-envelope comparison on Linux/glibc.

Run from any directory. No third-party Python packages are required. The
compiler builds the same workloads for the compact and certified classes.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import statistics
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def rows(path):
    with path.open() as source:
        return list(csv.DictReader(source))


def summarize(output):
    direct = {name: rows(output / f'local-{name}.csv')
              for name in ('baseline', 'certified')}
    with (output / 'direct-summary.csv').open('w', newline='') as dest:
        writer = csv.writer(dest)
        writer.writerow(['variant', 'k', 'resident_peak', 'update_peak',
                         'max_percent', 'median_update_ns', 'median_rank_ns'])
        for name, measurements in direct.items():
            for capacity in (8, 16, 32, 64, 128):
                subset = [r for r in measurements if int(r['k']) == capacity]
                writer.writerow([name, capacity,
                    max(int(r['resident_peak']) for r in subset),
                    max(int(r['update_peak']) for r in subset),
                    max(float(r['max_percent']) for r in subset),
                    statistics.median(float(r['update_ns']) for r in subset),
                    statistics.median(float(r['rank_ns']) for r in subset)])
    baseline = rows(output / 'fresh-baseline.csv')
    certified = rows(output / 'fresh-certified.csv')
    assert len(baseline) == len(certified) == 500
    for a, b in zip(baseline, certified):
        for key in ('k', 'shape', 'seed', 'input_fingerprint', 'queries'):
            assert a[key] == b[key], (key, a, b)
        for key in ('median_percent', 'p95_percent', 'max_percent'):
            assert float(b[key]) <= float(a[key]) + 1e-13, (key, a, b)
    with (output / 'fresh-summary.csv').open('w', newline='') as dest:
        writer = csv.writer(dest)
        writer.writerow(['group', 'name', 'cases', 'baseline_max_percent',
                         'certified_max_percent'])
        for field, names in [('shape', ('uniform', 'clustered', 'duplicates',
                                       'outliers', 'ascending')),
                             ('k', ('8', '16', '32', '64', '128'))]:
            for name in names:
                a = [r for r in baseline if r[field] == name]
                b = [r for r in certified if r[field] == name]
                writer.writerow([field, name, len(a),
                    max(float(r['max_percent']) for r in a),
                    max(float(r['max_percent']) for r in b)])
    grid = rows(output / 'local-certified-grid.csv')
    with (output / 'matched-memory.csv').open('w', newline='') as dest:
        writer = csv.writer(dest)
        writer.writerow(['budget_bytes', 'baseline_k', 'baseline_max_percent',
                         'certified_k', 'certified_resident_peak',
                         'certified_max_percent'])
        for capacity in (8, 16, 32, 64, 128):
            a = [r for r in direct['baseline'] if int(r['k']) == capacity]
            budget = max(int(r['resident_peak']) for r in a)
            fits = []
            for k in sorted({int(r['k']) for r in grid}):
                b = [r for r in grid if int(r['k']) == k]
                resident = max(int(r['resident_peak']) for r in b)
                if resident <= budget:
                    fits.append((k, resident, max(float(r['max_percent']) for r in b)))
            chosen = max(fits) if fits else ('', '', '')
            writer.writerow([budget, capacity,
                             max(float(r['max_percent']) for r in a), *chosen])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=ROOT / 'docs/results/2026-09-30/certified')
    parser.add_argument('--compiler', default='g++')
    parser.add_argument('--summarize-only', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.summarize_only:
        summarize(output)
        return
    flags = ['-O3', '-std=c++17', '-I', str(ROOT / 'include')]
    with tempfile.TemporaryDirectory(prefix='splinesketch-certified-', dir='/tmp') as temp:
        build = Path(temp)
        datasets = build / 'datasets'
        datasets.mkdir()
        for variant in ('baseline', 'certified'):
            for target, source in [('accuracy', 'benchmarks/accuracy.cpp'),
                                   ('local', 'benchmarks/comparison/local.cpp'),
                                   ('resize', 'benchmarks/comparison/resize.cpp')]:
                executable = build / f'{target}-{variant}'
                definitions = ['-DSPLINESKETCH_CERTIFIED'] if variant == 'certified' else []
                subprocess.run([args.compiler, *flags, *definitions,
                                str(ROOT / source), '-o', str(executable)], check=True)
                arguments = {'accuracy': ['--details'],
                             'local': [str(datasets)], 'resize': ['--fresh']}[target]
                filename = 'fresh' if target == 'resize' else target
                with (output / f'{filename}-{variant}.csv').open('w') as dest:
                    subprocess.run([str(executable), *arguments], stdout=dest, check=True)
                if target == 'local' and variant == 'certified':
                    with (output / 'local-certified-grid.csv').open('w') as dest:
                        subprocess.run([str(executable), str(datasets), '--grid'],
                                       stdout=dest, check=True)
        sources = ['include/splinesketch/splinesketch.hpp',
                   'include/splinesketch/certified_splinesketch.hpp',
                   'benchmarks/accuracy_workload.hpp', 'benchmarks/accuracy.cpp',
                   'benchmarks/comparison/local.cpp', 'benchmarks/comparison/resize.cpp',
                   'benchmarks/comparison/certified_experiments.py']
        metadata = {'compiler': subprocess.check_output([args.compiler, '--version'], text=True),
                    'flags': flags, 'platform': platform.platform(),
                    'fresh_seeds': [23, 42], 'observations': 6000,
                    'accuracy_cases': 225, 'fresh_resize_cases': 500,
                    'source_sha256': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                                      for name in sources},
                    'memory': 'Linux glibc malloc_usable_size; includes root object',
                    'timing': 'Single local run; descriptive, not a regression threshold'}
        (output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    summarize(output)
    print(f'Results written to {output}')


if __name__ == '__main__':
    main()
