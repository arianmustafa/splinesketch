#!/usr/bin/env python3
"""Reproduce certificate-aware join experiments without editing production."""
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
VARIANTS = ('width', 'guarded', 'balanced-low', 'balanced-mid', 'balanced-high')


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
    sources = [ROOT / 'include/splinesketch' / name for name in ('splinesketch.hpp', 'paper_splinesketch.hpp')]
    sources += [Path(__file__).resolve(), SOURCE / 'paper_baseline.py', SOURCE / 'paper_certificate_quality.cpp',
                SOURCE / 'paper_join_paired.cpp', SOURCE / 'paper_join_cost.cpp', SOURCE / 'paper_equivalence.cpp',
                ROOT / 'benchmarks/accuracy_workload.hpp', ROOT / 'tests/data/paper_width_counterexample.hpp',
                ROOT / 'tests/paper_join_selection_checks.cpp', ROOT / 'tests/paper_width_tests.cpp',
                ROOT / 'tests/paper_width_exact_checks.py']
    sources += [SOURCE / 'experiments' / f'paper_join_{name}.patch' for name in VARIANTS]
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in sources}
    metadata = {'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'compiler': capture([compiler, '--version']), 'flags': flags, 'platform': platform.platform(),
                'source_sha256': hashes, 'trials': args.trials, 'production_modified': False,
                'development_suite': 'original seeds 0-2; all five candidates',
                'heldout_suite': 'seeds 71-73; width and guarded candidates; no retuning',
                'paired_stream_seed': 103, 'paired_long_seeds': [97, 103],
                'paired_queries': 'common workload queries plus union of both node query sets',
                'paired_inputs': 'reference-generated, including feedback attacks; identical for both algorithms',
                'timing': 'five rotated sequential trials after correctness; normal fixed-capacity streaming',
                'allocation_scope': 'glibc allocated bytes during sketch operations; excludes datasets and queries',
                'proof_review': 'written local proof and finite checks; no independent or formal review'}
    with tempfile.TemporaryDirectory(dir='/tmp', prefix='splinesketch-join-') as directory:
        build = Path(directory); reference = build / 'reference'
        includes = {}
        for name in ('reference', *VARIANTS):
            headers = build / name / 'include/splinesketch'; headers.mkdir(parents=True)
            for file in ('splinesketch.hpp', 'paper_splinesketch.hpp'):
                shutil.copyfile(ROOT / 'include/splinesketch' / file, headers / file)
            if name != 'reference':
                run(['patch', '--batch', '-p1', '-i', SOURCE / 'experiments' / f'paper_join_{name}.patch'], cwd=build / name)
            includes[name] = headers.parent
        # GCC can deduplicate byte-identical pragma-once headers at distinct
        # paths. The namespace comparison needs a second core definition.
        core = includes['reference'] / 'splinesketch/splinesketch.hpp'
        core.write_text(core.read_text() + '\n// Reference namespace include marker.\n')
        metadata['reference_core_with_include_marker_sha256'] = digest(core)
        metadata['candidate_header_sha256'] = {name: digest(include / 'splinesketch/paper_splinesketch.hpp')
                                                for name, include in includes.items()}
        previous = [f'-DSPLINESKETCH_PREVIOUS_CORE="{reference}/include/splinesketch/splinesketch.hpp"',
                    f'-DSPLINESKETCH_PREVIOUS_PAPER="{reference}/include/splinesketch/paper_splinesketch.hpp"']
        modes = [('native', [])]
        if platform.machine() == 'x86_64' and 'gcc version' in capture([compiler, '-v']):
            modes.append(('binary64', ['-mlong-double-64']))
        for mode, extra in modes:
            names = ('reference', *VARIANTS) if mode == 'native' else ('reference', 'width', 'guarded')
            for name in names:
                quality = build / f'{name}-quality-{mode}'
                run([compiler, *flags, *extra, '-DSPLINESKETCH_ALLOW_ALGORITHM_ACCURACY_CHANGES',
                     '-I', includes[name], SOURCE / 'paper_certificate_quality.cpp', '-o', quality])
                output_run([quality], output / f'quality-original-{name}-{mode}.csv')
                if name in ('reference', 'width', 'guarded'):
                    output_run([quality, '--fresh'], output / f'quality-fresh-{name}-{mode}.csv')
                print('Quality complete:', mode, name, flush=True)
            for name, number in (('width', 1), ('guarded', 2)):
                checks = build / f'{name}-checks-{mode}'
                run([compiler, *flags, *extra, '-DSPLINESKETCH_TESTING', '-DSPLINESKETCH_PAPER_FORCE_HEAPS',
                     '-DSPLINESKETCH_VERIFY_PAPER_HEAPS', f'-DSPLINESKETCH_JOIN_EXPERIMENT={number}',
                     '-I', includes[name], ROOT / 'tests/paper_join_selection_checks.cpp', '-o', checks])
                output_run([checks], output / f'selection-{name}-{mode}.txt')
                width_checks = build / f'{name}-width-{mode}'
                run([compiler, *flags, *extra, '-DSPLINESKETCH_TESTING', '-DSPLINESKETCH_PAPER_FORCE_HEAPS',
                     '-DSPLINESKETCH_VERIFY_PAPER_HEAPS', '-I', includes[name], ROOT / 'tests/paper_width_tests.cpp', '-o', width_checks])
                output_run(['python3', ROOT / 'tests/paper_width_exact_checks.py', width_checks],
                           output / f'exact-{name}-{mode}.txt')
                paired = build / f'{name}-paired-{mode}'
                run([compiler, *flags, *extra, *previous, '-I', includes[name], SOURCE / 'paper_join_paired.cpp', '-o', paired])
                output_run([paired], output / f'paired-{name}-{mode}.csv')
                if mode == 'native':
                    output_run([paired, '--long'], output / f'long-{name}.csv')
                equivalence = build / f'{name}-plain-{mode}'
                run([compiler, *flags, *extra, *previous, '-I', includes[name], SOURCE / 'paper_equivalence.cpp', '-o', equivalence])
                output_run([equivalence, '--extended'], output / f'plain-{name}-{mode}.txt')
                print('Checks and paired streams complete:', mode, name, flush=True)
        for name in ('reference', 'width', 'guarded'):
            run([compiler, *flags, '-I', includes[name], SOURCE / 'paper_join_cost.cpp', '-o', build / f'{name}-cost'])
        for trial in range(args.trials):
            names = ['reference', 'width', 'guarded']; order = names[trial % 3:] + names[:trial % 3]
            for name in order:
                output_run([build / f'{name}-cost'], output / f'cost-{name}-{trial}.csv')
            print('Timing trial complete:', trial + 1, flush=True)
        with (output / 'cost-summary.csv').open('w', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['variant', 'k', 'median_update_ns', 'median_rank_ns', 'median_resident_peak', 'median_update_peak', 'median_final_bytes'])
            for name in ('reference', 'width', 'guarded'):
                trials = [rows(output / f'cost-{name}-{i}.csv') for i in range(args.trials)]
                keys = [('k', 'shape', 'seed')]
                assert all(len(trial) == 60 for trial in trials)
                assert all(all(a[field] == b[field] for field in keys[0]) for trial in trials for a, b in zip(trials[0], trial))
                for k in (128, 256, 512, 1024):
                    values = []
                    for field in ('update_ns', 'rank_ns', 'resident_peak', 'update_peak', 'final_bytes'):
                        per_case = [statistics.median(float(trial[i][field]) for trial in trials)
                                    for i, row in enumerate(trials[0]) if int(row['k']) == k]
                        values.append(statistics.median(per_case))
                    writer.writerow([name, k, *values])
    assert hashes == {str(p.relative_to(ROOT)): digest(p) for p in sources}, 'sources changed during experiment'
    (output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print('Results:', output, flush=True)


if __name__ == '__main__':
    main()
