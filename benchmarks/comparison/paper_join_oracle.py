#!/usr/bin/env python3
"""Full-data join lookahead and a data-free local CDF-distortion experiment."""
import argparse
import csv
import datetime
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile
from paper_baseline import capture, output_run, rows, run

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'benchmarks/comparison'
FLAGS = ['-std=c++17', '-O3', '-DNDEBUG']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def quality_provenance(include, compiler):
    paths = [include / 'splinesketch/splinesketch.hpp', include / 'splinesketch/paper_splinesketch.hpp',
             SOURCE / 'paper_certificate_quality.cpp', ROOT / 'benchmarks/accuracy_workload.hpp']
    return {'compiler': capture([compiler, '--version']), 'flags': FLAGS + ['-DSPLINESKETCH_ALLOW_ALGORITHM_ACCURACY_CHANGES'],
            'suite': 'fresh seeds 71-73, 60000 observations, native',
            'source_sha256': {str(i): digest(path) for i, path in enumerate(paths)}}


def exact_error(row):
    return Fraction.from_float(float(row['error_high'])) + Fraction.from_float(float(row['error_low']))


def validate_oracle(path, events_path, plan_path):
    data = rows(path); events = rows(events_path)
    baseline = data[0]; assert baseline['phase'] == 'baseline' and baseline['success'] == '1'
    assert int(baseline['joins']) == len(events)
    singles = [row for row in data if row['phase'] == 'single']
    for event in events:
        arms = [r for r in singles if r['event'] == event['event']]
        assert len(arms) == int(event['legal_choices'])
        original = [r for r in arms if r['is_original'] == '1']
        assert len(original) == 1 and original[0]['success'] == '1'
        assert original[0]['query_hash'] == baseline['query_hash']
    final = data[-1]; assert final['phase'] == 'oracle' and final['success'] == '1'
    assert exact_error(final) <= exact_error(baseline)
    successful = [r for r in singles if r['success'] == '1']
    best = min(successful, key=exact_error) if successful else baseline
    event_states = {r['event']: r['state_hash'] for r in events}
    single_plan = plan_path.with_name(plan_path.name.replace('-plan.txt', '-single-plan.txt'))
    if best['phase'] == 'single' and exact_error(best) < exact_error(baseline):
        single_plan.write_text(f"{best['event']} {event_states[best['event']]} {best['choice_bits']}\n")
    else:
        single_plan.write_text('')
    return {'baseline_joins': len(events), 'single_arms': len(singles),
            'improving_single_arms': sum(exact_error(r) < exact_error(baseline) for r in successful),
            'failed_single_arms': len(singles) - len(successful), 'baseline_error_percent': float(baseline['error_percent']),
            'best_single_error_percent': float(best['error_percent']), 'lookahead_error_percent': float(final['error_percent']),
            'lookahead_changed_choices': int(final['event']), 'lookahead_arms': sum(r['phase'] == 'lookahead' for r in data),
            'queries': int(baseline['queries']), 'single_plan': str(single_plan)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reuse-fresh-native', type=Path,
                        help='Reuse a completed CSV with a matching .provenance.json source manifest.')
    args = parser.parse_args()
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=True)
    compiler = shutil.which('c++')
    if not compiler:
        raise SystemExit('C++17 compiler required')
    paths = [ROOT / 'include/splinesketch' / file for file in ('splinesketch.hpp', 'paper_splinesketch.hpp')]
    paths += [Path(__file__).resolve(), SOURCE / 'paper_baseline.py', SOURCE / 'paper_join_oracle.cpp',
              SOURCE / 'paper_certificate_quality.cpp', SOURCE / 'paper_join_paired.cpp', SOURCE / 'paper_equivalence.cpp',
              ROOT / 'benchmarks/accuracy_workload.hpp', ROOT / 'tests/data/paper_width_counterexample.hpp',
              ROOT / 'tests/paper_join_oracle_exact_checks.py', ROOT / 'tests/paper_width_tests.cpp',
              ROOT / 'tests/paper_width_exact_checks.py']
    paths += [SOURCE / 'experiments' / file for file in ('paper_join_oracle.patch', 'paper_join_distortion.patch')]
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in paths}
    metadata = {'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'compiler': capture([compiler, '--version']), 'flags': FLAGS, 'platform': platform.platform(),
                'source_sha256': hashes, 'production_modified': False, 'capacity': 32,
                'query_objective': 'exact binary64 returned-rank error on a fixed full-data probe set; not an all-query supremum',
                'search_scope': 'all single alternatives on baseline path, then sequential default-tail lookahead; not exhaustive global search',
                'distortion_score': 'data-free, local translated CDF, 17 probes per affected interval; no worst-case numeric bound claimed',
                'fresh_quality_reused': args.reuse_fresh_native is not None, 'oracle_summary': {}}
    with tempfile.TemporaryDirectory(dir='/tmp', prefix='splinesketch-oracle-') as directory:
        build = Path(directory); includes = {}
        for variant in ('reference', 'cdf'):
            headers = build / variant / 'include/splinesketch'; headers.mkdir(parents=True)
            for file in ('splinesketch.hpp', 'paper_splinesketch.hpp'):
                shutil.copyfile(ROOT / 'include/splinesketch' / file, headers / file)
            if variant == 'cdf':
                run(['patch', '--batch', '-p1', '-i', SOURCE / 'experiments/paper_join_distortion.patch'], cwd=build / variant)
            run(['patch', '--batch', '-p1', '-i', SOURCE / 'experiments/paper_join_oracle.patch'], cwd=build / variant)
            includes[variant] = headers.parent
        metadata['instrumented_header_sha256'] = {name: digest(inc / 'splinesketch/paper_splinesketch.hpp') for name, inc in includes.items()}
        modes = [('native', [])]
        if platform.machine() == 'x86_64' and 'gcc version' in capture([compiler, '-v']):
            modes.append(('binary64', ['-mlong-double-64']))
        binaries = {}
        for mode, extra in modes:
            for variant in ('reference', 'cdf'):
                binary = build / f'oracle-{variant}-{mode}'
                run([compiler, *FLAGS, *extra, '-I', includes[variant], SOURCE / 'paper_join_oracle.cpp', '-o', binary])
                binaries[variant, mode] = binary
        # Freeze feedback inputs from the native reference once, for all modes.
        inputs = [('witness', [])]
        for seed in (71, 97, 103):
            path = output / f'mass-{seed}.bits'
            output_run([binaries['reference', 'native'], '--generate-mass', seed], path)
            inputs.append((f'mass-{seed}', ['--data', path]))
        metadata['input_sha256'] = {p.name: digest(p) for p in output.glob('*.bits')}
        metadata['input_sha256']['witness_header'] = digest(ROOT / 'tests/data/paper_width_counterexample.hpp')
        for mode, extra in modes:
            for name, options in inputs:
                prefix = output / f'{name}-{mode}'
                csv_path = output / f'{name}-{mode}.csv'
                output_run([binaries['reference', mode], *options, '--output-prefix', prefix], csv_path)
                summary = validate_oracle(csv_path, Path(str(prefix) + '-events.csv'), Path(str(prefix) + '-plan.txt'))
                metadata['oracle_summary'][f'{name}-{mode}'] = summary
                for label, plan in [('baseline', None), ('single', summary['single_plan']), ('lookahead', str(prefix) + '-plan.txt')]:
                    plan_args = [] if plan is None else ['--plan', plan]
                    output_run(['python3', ROOT / 'tests/paper_join_oracle_exact_checks.py', binaries['reference', mode],
                                *options, *plan_args], output / f'exact-{name}-{mode}-{label}.txt')
                output_run([binaries['cdf', mode], *options, '--evaluate-only'], output / f'cdf-{name}-{mode}.csv')
                output_run(['python3', ROOT / 'tests/paper_join_oracle_exact_checks.py', binaries['cdf', mode], *options],
                           output / f'exact-cdf-{name}-{mode}.txt')
                print('Oracle and proxy complete:', name, mode, flush=True)
            checks = build / f'width-cdf-{mode}'
            run([compiler, *FLAGS, *extra, '-DSPLINESKETCH_TESTING', '-DSPLINESKETCH_PAPER_FORCE_HEAPS',
                 '-DSPLINESKETCH_VERIFY_PAPER_HEAPS', '-I', includes['cdf'], ROOT / 'tests/paper_width_tests.cpp', '-o', checks])
            output_run(['python3', ROOT / 'tests/paper_width_exact_checks.py', checks], output / f'heap-exact-cdf-{mode}.txt')
        # The source-identical prior production phase supplies reference quality.
        reference_dir = ROOT / 'docs/results/2026-10-01/paper-width-analysis'
        prior = json.loads((reference_dir / 'metadata.json').read_text())
        for path in (ROOT / 'include/splinesketch/splinesketch.hpp', ROOT / 'include/splinesketch/paper_splinesketch.hpp',
                     SOURCE / 'paper_certificate_quality.cpp', ROOT / 'benchmarks/accuracy_workload.hpp'):
            assert digest(path) == prior['source_sha256'][str(path.relative_to(ROOT))]
        reference_quality = reference_dir / 'quality-fresh-after-native.csv'
        shutil.copyfile(reference_quality, output / 'quality-fresh-reference-native.csv')
        metadata['reference_quality'] = {'reused': True, 'path': str(reference_quality.relative_to(ROOT)), 'sha256': digest(reference_quality)}
        target_quality = output / 'quality-fresh-cdf-native.csv'
        provenance = quality_provenance(includes['cdf'], compiler)
        if args.reuse_fresh_native:
            source = args.reuse_fresh_native.resolve()
            saved = json.loads(Path(str(source) + '.provenance.json').read_text())
            assert saved['inputs'] == provenance and saved['csv_sha256'] == digest(source), 'quality artifact provenance mismatch'
            shutil.copyfile(source, target_quality)
        else:
            binary = build / 'quality-cdf-native'
            run([compiler, *FLAGS, '-DSPLINESKETCH_ALLOW_ALGORITHM_ACCURACY_CHANGES', '-I', includes['cdf'],
                 SOURCE / 'paper_certificate_quality.cpp', '-o', binary])
            output_run([binary, '--fresh'], target_quality)
        Path(str(target_quality) + '.provenance.json').write_text(json.dumps({'inputs': provenance, 'csv_sha256': digest(target_quality)}, indent=2) + '\n')
        before, after = rows(output / 'quality-fresh-reference-native.csv'), rows(target_quality)
        assert len(before) == len(after) == 300
        assert all(all(x[field] == y[field] for field in ('k', 'shape', 'seed', 'workflow', 'n', 'queries', 'raw_error_max_percent'))
                   for x, y in zip(before, after))
        metadata['fresh_quality_counts'] = {field: {'improved': sum(float(y[field]) < float(x[field]) - 1e-9 for x, y in zip(before, after)),
                                                   'worsened': sum(float(y[field]) > float(x[field]) + 1e-9 for x, y in zip(before, after))}
                                           for field in ('error_max_percent', 'uniform_bound_percent')}
    assert hashes == {str(p.relative_to(ROOT)): digest(p) for p in paths}, 'sources changed during oracle analysis'
    (output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print('Results:', output, flush=True)


if __name__ == '__main__':
    main()
