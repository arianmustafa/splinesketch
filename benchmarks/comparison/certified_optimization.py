#!/usr/bin/env python3
"""Compare optimized certificates against the reconstructed previous headers."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
PATCH = ROOT / 'benchmarks/comparison/experiments/certified_optimization.patch'


def read_rows(path):
    with path.open() as source:
        return list(csv.DictReader(source))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=ROOT / 'docs/results/2026-09-30/certified-optimization')
    parser.add_argument('--compiler', default='g++')
    parser.add_argument('--repeats', type=int, default=5)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    flags = ['-O3', '-std=c++17', '-DSPLINESKETCH_CERTIFIED']
    with tempfile.TemporaryDirectory(prefix='splinesketch-opt-', dir='/tmp') as temp:
        build = Path(temp)
        before = build / 'previous'
        headers = ['splinesketch.hpp', 'certified_splinesketch.hpp']
        (before / 'include/splinesketch').mkdir(parents=True)
        for header in headers:
            shutil.copyfile(ROOT / 'include/splinesketch' / header,
                            before / 'include/splinesketch' / header)
        subprocess.run(['patch', '-R', '-p1', '-i', str(PATCH)], cwd=before, check=True)
        datasets = build / 'datasets'
        datasets.mkdir()
        for variant, include in [('before', before/'include'), ('after', ROOT/'include')]:
            for target, source in [('local', 'benchmarks/comparison/local.cpp'),
                                   ('accuracy', 'benchmarks/accuracy.cpp'),
                                   ('fresh', 'benchmarks/comparison/resize.cpp')]:
                subprocess.run([args.compiler, *flags, '-I', str(include),
                                str(ROOT/source), '-o', str(build/f'{variant}-{target}')], check=True)
        definitions = [f'-DSPLINESKETCH_PREVIOUS_CORE={json.dumps(str(before / "include/splinesketch/splinesketch.hpp"))}',
                       f'-DSPLINESKETCH_PREVIOUS_CERTIFIED={json.dumps(str(before / "include/splinesketch/certified_splinesketch.hpp"))}']
        subprocess.run([args.compiler, *flags, '-I', str(ROOT/'include'), *definitions,
                        str(ROOT/'benchmarks/comparison/certified_equivalence.cpp'),
                        '-o', str(build/'equivalence')], check=True)
        with (output/'equivalence.txt').open('w') as dest:
            subprocess.run([str(build/'equivalence')], stdout=dest, check=True)
        for variant in ('before', 'after'):
            for target, option in [('accuracy', '--details'), ('fresh', '--fresh')]:
                with (output/f'{variant}-{target}.csv').open('w') as dest:
                    subprocess.run([str(build/f'{variant}-{target}'), option], stdout=dest, check=True)
        assert (output/'before-accuracy.csv').read_bytes() == (output/'after-accuracy.csv').read_bytes()
        old = read_rows(output/'before-fresh.csv')
        new = read_rows(output/'after-fresh.csv')
        assert len(old) == len(new) == 500
        for a, b in zip(old, new):
            for field in a:
                if field not in ('shrink_us', 'grow_us'):
                    assert a[field] == b[field], (field, a, b)
        measurements = {'before': [], 'after': []}
        for repeat in range(args.repeats):
            order = ('before', 'after') if repeat % 2 == 0 else ('after', 'before')
            for variant in order:
                path = output/f'{variant}-local-{repeat}.csv'
                with path.open('w') as dest:
                    subprocess.run([str(build/f'{variant}-local'), str(datasets)], stdout=dest, check=True)
                sample = read_rows(path)
                assert len(sample) == 75
                measurements[variant].append(sample)
        with (output/'summary.csv').open('w', newline='') as dest:
            writer = csv.writer(dest)
            writer.writerow(['variant', 'k', 'resident_peak', 'update_peak', 'max_percent',
                             'median_update_ns', 'median_rank_ns', 'median_quantile_ns'])
            for variant, runs in measurements.items():
                for capacity in (8, 16, 32, 64, 128):
                    cases = [(shape, str(seed)) for shape in ('uniform', 'clustered', 'duplicates',
                                                             'outliers', 'ascending') for seed in range(3)]
                    medians = {field: [] for field in ('update_ns', 'rank_ns', 'quantile_ns')}
                    subset = [r for run in runs for r in run if int(r['k']) == capacity]
                    for shape, seed in cases:
                        repeated = [r for r in subset if r['shape'] == shape and r['seed'] == seed]
                        assert len(repeated) == args.repeats
                        for field in medians:
                            medians[field].append(statistics.median(float(r[field]) for r in repeated))
                    writer.writerow([variant, capacity,
                                     max(int(r['resident_peak']) for r in subset),
                                     max(int(r['update_peak']) for r in subset),
                                     max(float(r['max_percent']) for r in subset),
                                     *[statistics.median(medians[field]) for field in medians]])
        files = ['include/splinesketch/splinesketch.hpp',
                 'include/splinesketch/certified_splinesketch.hpp',
                 'benchmarks/comparison/certified_optimization.py',
                 'benchmarks/comparison/certified_equivalence.cpp',
                 'benchmarks/comparison/experiments/certified_optimization.patch',
                 'benchmarks/comparison/local.cpp', 'benchmarks/comparison/resize.cpp',
                 'benchmarks/accuracy.cpp', 'benchmarks/accuracy_workload.hpp']
        metadata = {'compiler': subprocess.check_output([args.compiler, '--version'], text=True),
                    'flags': flags, 'platform': platform.platform(), 'repeats': args.repeats,
                    'order': 'alternating before/after; no concurrent benchmark jobs',
                    'timing_summary': 'median across cases of per-case medians across repetitions',
                    'memory': 'Linux/glibc malloc_usable_size plus root object; peaks across cases/runs',
                    'accuracy': 'identical original matrix summaries and fresh 500-case summaries',
                    'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files},
                    'previous_header_sha256': {name: hashlib.sha256((before/'include/splinesketch'/name).read_bytes()).hexdigest()
                                               for name in headers}}
        (output/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(f'Paired checks and repeated measurements saved to {output}')


if __name__ == '__main__':
    main()
