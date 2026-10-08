#!/usr/bin/env python3
"""Reproduce the width search, reduced witness, tighter bounds and API query cost."""
import argparse
import csv
import datetime
import hashlib
import json
from pathlib import Path
import platform
import shutil
import statistics
import tempfile
from paper_baseline import capture, output_run, rows, run

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'benchmarks/comparison'
PATCH = SOURCE / 'experiments/paper_certificate_refinement.patch'
MEDIAN_PATCH = SOURCE / 'experiments/paper_median_split_experiment.patch'


def compare(before, after, varying):
    a, b = rows(before), rows(after)
    assert len(a) == len(b), (before, after)
    assert all(all(x[key] == y[key] for key in x if key not in varying)
               for x, y in zip(a, b)), (before, after)
    for x, y in zip(a, b):
        for field in varying - {'ns_per_call'}:
            assert float(y[field]) <= float(x[field]), (field, x, y)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--trials', type=int, default=5)
    args = parser.parse_args()
    if args.trials < 1:
        parser.error('--trials must be positive')
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    compiler = shutil.which('c++')
    if not compiler:
        raise SystemExit('C++17 compiler required')
    flags = ['-std=c++17', '-O3', '-DNDEBUG']
    paths = [ROOT / 'include/splinesketch' / name for name in ('splinesketch.hpp', 'paper_splinesketch.hpp')]
    paths += [Path(__file__).resolve(), PATCH, MEDIAN_PATCH, SOURCE / 'paper_baseline.py',
              SOURCE / 'paper_width_search.cpp', SOURCE / 'paper_width_cost.cpp',
              SOURCE / 'paper_certificate_quality.cpp', ROOT / 'benchmarks/accuracy_workload.hpp',
              ROOT / 'tests/data/paper_width_counterexample.hpp']
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    metadata = {'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'compiler': capture([compiler, '--version']), 'flags': flags,
                'platform': platform.platform(), 'trials': args.trials,
                'source_sha256': hashes, 'development_seed': 71, 'fresh_seed': 97,
                'short_observations_per_case': 200000, 'long_observations_per_case': 2000000,
                'long_order': 'sequential before/after; exact query hashes and metrics compared',
                'timing_order': 'alternating before/after; sequential after all correctness jobs',
                'cost_scope': 'max_rank_error() only; estimate/update paths and storage unchanged'}
    with tempfile.TemporaryDirectory(dir='/tmp', prefix='splinesketch-width-') as directory:
        build = Path(directory); reference = build / 'reference'; headers = reference / 'include/splinesketch'
        headers.mkdir(parents=True)
        for name in ('splinesketch.hpp', 'paper_splinesketch.hpp'):
            shutil.copyfile(ROOT / 'include/splinesketch' / name, headers / name)
        run(['patch', '--batch', '-R', '-p1', '-i', PATCH], cwd=reference)
        metadata['reference_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in headers.iterdir()}
        modes = [('native', [])]
        if platform.machine() == 'x86_64' and 'gcc version' in capture([compiler, '-v']):
            modes += [('binary64', ['-mlong-double-64'])]
        for mode, extra in modes:
            for variant in ('before', 'after'):
                include = reference / 'include' if variant == 'before' else ROOT / 'include'
                binary = build / f'quality-{variant}-{mode}'
                run([compiler, *flags, *extra, '-I', include, SOURCE / 'paper_certificate_quality.cpp', '-o', binary])
                for suite, options in [('original', []), ('fresh', ['--fresh'])]:
                    output_run([binary, *options], output / f'quality-{suite}-{variant}-{mode}.csv')
            for suite in ('original', 'fresh'):
                compare(output / f'quality-{suite}-before-{mode}.csv', output / f'quality-{suite}-after-{mode}.csv',
                        {'uniform_bound_percent'})
        # Search includes prefixes at the largest normalized width and final states.
        for variant in ('before', 'after'):
            include = reference / 'include' if variant == 'before' else ROOT / 'include'
            binary = build / f'search-{variant}'
            run([compiler, *flags, '-I', include, SOURCE / 'paper_width_search.cpp', '-o', binary])
            for name, options in [('short', []), ('long', ['--long', '97']), ('histories', ['--histories'])]:
                output_run([binary, *options], output / f'{name}-{variant}.csv')
            output_run([binary, '--reduce'], output / f'reduction-{variant}.txt')
        for suite in ('short', 'long'):
            compare(output / f'{suite}-before.csv', output / f'{suite}-after.csv', {'peak_uniform', 'final_uniform'})
        compare(output / 'histories-before.csv', output / 'histories-after.csv', {'uniform'})
        assert (output / 'reduction-before.txt').read_bytes() == (output / 'reduction-after.txt').read_bytes()
        # Development experiment: compare changed cuts on the same fixed
        # streams. Truth containment/error checks still run; estimate accuracy
        # may regress and is measured rather than asserted away.
        experiment = build / 'median'
        shutil.copytree(reference, experiment)
        run(['patch', '--batch', '-p1', '-i', MEDIAN_PATCH], cwd=experiment)
        binary = build / 'median-quality'
        run([compiler, *flags, '-DSPLINESKETCH_ALLOW_ALGORITHM_ACCURACY_CHANGES',
             '-I', experiment / 'include', SOURCE / 'paper_certificate_quality.cpp', '-o', binary])
        output_run([binary], output / 'median-quality-native.csv')
        metadata['median_experiment'] = {'shipped': False, 'suite': 'original development seeds 0-2',
            'headers_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in (experiment / 'include/splinesketch').iterdir()}}
        for variant in ('before', 'after'):
            include = reference / 'include' if variant == 'before' else ROOT / 'include'
            run([compiler, *flags, '-I', include, SOURCE / 'paper_width_cost.cpp', '-o', build / f'cost-{variant}'])
        for trial in range(args.trials):
            order = ('before', 'after') if trial % 2 == 0 else ('after', 'before')
            for variant in order:
                output_run([build / f'cost-{variant}'], output / f'cost-{variant}-{trial}.csv')
            compare(output / f'cost-before-{trial}.csv', output / f'cost-after-{trial}.csv', {'uniform', 'ns_per_call'})
        with (output / 'cost-summary.csv').open('w', newline='') as stream:
            writer = csv.writer(stream); writer.writerow(['variant', 'k', 'median_ns_per_call'])
            for variant in ('before', 'after'):
                samples = [r for i in range(args.trials) for r in rows(output / f'cost-{variant}-{i}.csv')]
                for k in (128, 256, 512, 1024):
                    writer.writerow([variant, k, statistics.median(float(r['ns_per_call']) for r in samples if int(r['k']) == k)])
    assert hashes == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}, 'sources changed during comparison'
    (output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print('Results:', output)


if __name__ == '__main__':
    main()
