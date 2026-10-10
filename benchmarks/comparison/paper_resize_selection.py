#!/usr/bin/env python3
"""Compare a resize-only certificate loss prototype without modifying production."""
import argparse
import csv
import datetime
import hashlib
import json
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'benchmarks/comparison'
PATCHES = {
    'sampled': SOURCE / 'experiments/paper_resize_certificate_loss.patch',
    'regrid': SOURCE / 'experiments/paper_resize_regrid_loss.patch',
}


def run(command, **kwargs):
    subprocess.run([str(arg) for arg in command], check=True, **kwargs)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize(path, suite):
    with path.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    result = {}
    for group, expected in [(suite, 288), ('cluster_witness', 2), ('exact_discrete', 2),
                            ('released_discrete', 2), ('mixed', 6)]:
        cases = [row for row in rows if row['group'] == group]
        if len(cases) != expected:
            raise RuntimeError(f'{group}: expected {expected} cases, got {len(cases)}')
        counts = {'better': 0, 'equal': 0, 'worse': 0}
        for row in cases:
            before, after = float(row['before_error']), float(row['after_error'])
            counts['better' if after < before else 'worse' if after > before else 'equal'] += 1
        result[group] = {
            'cases': len(cases), 'final_measured_error': counts,
            'queries_checked': sum(int(row['queries']) for row in cases),
            'cases_with_any_worse_checkpoint': sum(int(row['worse_steps']) > 0 for row in cases),
            'worst_error_regression': max(cases, key=lambda row: float(row['after_error']) - float(row['before_error'])),
            'median_exploratory_resize_time_ratio': statistics.median(
                float(row['after_resize_ns']) / float(row['before_resize_ns']) for row in cases),
        }
    return result, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--compiler', default='c++')
    parser.add_argument('--variant', choices=PATCHES, default='regrid')
    parser.add_argument('--suite', choices=('development', 'fresh'), default='fresh')
    args = parser.parse_args()
    compiler = shutil.which(args.compiler)
    if not compiler:
        parser.error('C++17 compiler required')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    patch = PATCHES[args.variant]
    sources = [ROOT / 'include/splinesketch' / name for name in ('splinesketch.hpp', 'paper_splinesketch.hpp')]
    sources += [Path(__file__).resolve(), SOURCE / 'paper_resize_selection.cpp', patch,
                ROOT / 'tests/paper_resize_history_fixture.hpp']
    hashes = {str(path.relative_to(ROOT)): digest(path) for path in sources}
    version = subprocess.check_output([compiler, '--version'], text=True)
    verbose_version = subprocess.run([compiler, '-v'], capture_output=True, text=True, check=True).stderr
    flags = ['-std=c++17', '-O2']
    if args.variant == 'regrid':
        flags.append('-DSPLINESKETCH_RESIZE_REGRID_ONLY')
    modes = [('native-scanner', []), ('native-heap', ['-DSPLINESKETCH_PAPER_FORCE_HEAPS',
                                                  '-DSPLINESKETCH_VERIFY_PAPER_HEAPS'])]
    if platform.machine() == 'x86_64' and 'gcc version' in verbose_version:
        modes.append(('binary64-heap', ['-mlong-double-64', '-DSPLINESKETCH_PAPER_FORCE_HEAPS',
                                      '-DSPLINESKETCH_VERIFY_PAPER_HEAPS']))
    metadata = {
        'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'compiler': version, 'flags': flags, 'platform': platform.platform(), 'source_sha256': hashes,
        'variant': args.variant, 'suite': args.suite, 'development_seeds': [113, 127, 149],
        'fresh_seeds': [179, 191, 211], 'shapes': list(range(8)), 'capacities': [16, 32, 64],
        'resize_operations_per_case': 24, 'mixed_seed': 163,
        'historical_cases': ['124-input resize fixture', '269-input discrete stream, seed 97',
                             '141-input discrete stream with heavy-hitter releases, seed 113'],
        'queries': 'input keys and neighbours, quarterpoints, union of both grids and neighbours, infinities',
        'timing': 'exploratory single measurements; reference first; heap modes include invariant verification',
        'proof_scope': 'local integer width identity and finite certificate checks; no capacity-only rate',
        'production_modified': False, 'summaries': {},
    }
    with tempfile.TemporaryDirectory(dir='/tmp', prefix='splinesketch-resize-selection-') as directory:
        build = Path(directory)
        for name in ('reference', 'candidate'):
            headers = build / name / 'include/splinesketch'
            headers.mkdir(parents=True)
            for source in sources[:2]:
                shutil.copyfile(source, headers / source.name)
        run(['patch', '--batch', '-p1', '-i', patch], cwd=build / 'candidate')
        reference = build / 'reference/include/splinesketch'
        # GCC may deduplicate byte-identical pragma-once files. Both namespaces
        # must have their own core definition for this paired comparison.
        core = reference / 'splinesketch.hpp'
        core.write_text(core.read_text() + '\n// Reference namespace include marker.\n')
        previous = [f'-DSPLINESKETCH_PREVIOUS_CORE="{core}"',
                    f'-DSPLINESKETCH_PREVIOUS_PAPER="{reference / "paper_splinesketch.hpp"}"']
        metadata['candidate_header_sha256'] = digest(build / 'candidate/include/splinesketch/paper_splinesketch.hpp')
        paired_rows = {}
        for mode, extra in modes:
            executable = build / mode
            run([compiler, *flags, *extra, *previous, '-I', build / 'candidate/include',
                 SOURCE / 'paper_resize_selection.cpp', '-o', executable])
            path = output / f'{mode}.csv'
            with path.open('w') as stream:
                run([executable, *(['--fresh'] if args.suite == 'fresh' else [])], stdout=stream)
            metadata['summaries'][mode], paired_rows[mode] = summarize(path, args.suite)
            print(mode, metadata['summaries'][mode][args.suite]['final_measured_error'], flush=True)
        # Native scanner and verified heaps must agree on all deterministic
        # results, including common query sets. Timing naturally differs.
        scanner, heap = paired_rows['native-scanner'], paired_rows['native-heap']
        comparable = lambda rows: [{key: value for key, value in row.items() if not key.endswith('_ns')}
                                   for row in rows]
        if comparable(scanner) != comparable(heap):
            raise RuntimeError('native scanner and heap disagree')
    if hashes != {str(path.relative_to(ROOT)): digest(path) for path in sources}:
        raise RuntimeError('sources changed during experiment')
    (output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print('Results:', output, flush=True)


if __name__ == '__main__':
    main()
