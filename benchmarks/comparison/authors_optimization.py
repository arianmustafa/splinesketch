#!/usr/bin/env python3
"""Apply and evaluate a direct patch to the pinned authors' Java sources."""
import argparse
import csv
import datetime
import hashlib
import json
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'benchmarks/comparison'
REVISION = '320019fe9d156a42ec6c1edba880d060de8a1b13'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def capture(command):
    return subprocess.check_output([str(x) for x in command], stderr=subprocess.STDOUT, text=True).strip()


def run(command, output, error):
    with output.open('w') as out, error.open('w') as err:
        subprocess.run([str(x) for x in command], check=True, stdout=out, stderr=err)


def summarize(output, trials):
    with (output / 'benchmark.csv').open() as stream:
        data = list(csv.DictReader(stream))
    assert len(data) == 2 * 2 * 5 * 15 * trials
    groups = {}
    for row in data:
        key = tuple(row[f] for f in ('variant', 'mg', 'k', 'shape', 'seed'))
        groups.setdefault(key, []).append(row)
    for key, group in groups.items():
        assert len(group) == trials and len({r['trial'] for r in group}) == trials
        peer = groups[('optimized' if key[0] == 'original' else 'original', *key[1:])]
        assert len({r['query_hash'] for r in group + peer}) == 1, key
        assert len({r['resident_bytes'] for r in group + peer}) == 1, key
    fields = ('update_ns', 'batch_ns', 'single_ns', 'update_allocated_bytes',
              'batch_allocated_bytes', 'single_allocated_bytes', 'resident_bytes')
    with (output / 'summary.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['mg', 'k', 'metric', 'original_median', 'optimized_median', 'median_original_over_optimized'])
        for mg in ('false', 'true'):
            for k in ('8', '32', '128', '512', '1024'):
                for field in fields:
                    before, after = [], []
                    for shape in range(5):
                        for seed in (71, 72, 73):
                            for variant, target in (('original', before), ('optimized', after)):
                                group = groups[(variant, mg, k, str(shape), str(seed))]
                                target.append(statistics.median(float(r[field]) for r in group))
                    writer.writerow([mg, k, field, statistics.median(before), statistics.median(after),
                                     statistics.median(a / b for a, b in zip(before, after))])
    return {'rows': len(data), 'paired_streams': len(groups) // 2,
            'identical_query_hashes': True, 'identical_retained_graph_sizes': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--authors', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--trials', type=int, default=5)
    args = parser.parse_args()
    if args.trials < 1:
        parser.error('--trials must be positive')
    authors, output = args.authors.resolve(), args.output.resolve()
    assert capture(['git', '-C', authors, 'rev-parse', 'HEAD']) == REVISION
    assert not capture(['git', '-C', authors, 'status', '--porcelain', '--untracked-files=no'])
    output.mkdir(parents=True, exist_ok=True)
    sources = [Path(__file__).resolve(), SOURCE / 'AuthorsOptimization.java',
               SOURCE / 'experiments/authors_query_optimization.patch',
               ROOT / 'include/splinesketch/splinesketch.hpp', ROOT / 'include/splinesketch/paper_splinesketch.hpp']
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in sources}
    reference_hashes = {name: digest(authors / name) for name in ('SplineSketch.java', 'SplineSketchMG.java', 'LICENSE')}
    folders = {name: output / name for name in ('reference', 'candidate', 'classes')}
    for folder in folders.values():
        folder.mkdir(exist_ok=True)
    for name in ('SplineSketch', 'SplineSketchMG'):
        text = (authors / f'{name}.java').read_text()
        # Only class identifiers change in the control so both coexist in one JVM.
        (folders['reference'] / f'Original{name}.java').write_text(re.sub(r'\b' + name + r'\b', 'Original' + name, text))
        (folders['candidate'] / f'{name}.java').write_text(text)
    shutil.copyfile(authors / 'LICENSE', folders['candidate'] / 'LICENSE')
    subprocess.run(['patch', '--batch', '-p1', '-i', str(sources[2])], cwd=folders['candidate'], check=True)
    metadata = {
        'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'authors_revision': REVISION, 'authors_source_sha256': reference_hashes, 'authors_license': 'MIT',
        'source_sha256': hashes, 'java': capture(['java', '-version']), 'javac': capture(['javac', '-version']),
        'platform': platform.platform(), 'trials': args.trials, 'production_modified': False,
        'reference_transformation': 'LF line endings and class identifier renaming only; same original methods',
        'candidate_source_sha256': {p.name: digest(p) for p in folders['candidate'].glob('*.java')},
        'suite': 'five existing distribution shapes, seeds 71-73, 6000 observations, capacities 8/32/128/512/1024',
        'timing': 'same JVM; per-capacity warmup; five alternating paired trials; no CPU pinning',
        'allocation': 'ThreadMXBean cumulative current-thread bytes; not peak live workspace or RSS',
        'resident': 'Instrumentation reachable object graph; excludes query results and runtime/class storage',
    }
    java_files = [*folders['reference'].glob('*.java'), *folders['candidate'].glob('*.java'), SOURCE / 'AuthorsOptimization.java']
    run(['javac', '-d', folders['classes'], *java_files], output / 'compile.txt', output / 'compile.stderr')
    flags = ['-XX:-UsePerfData', '-Xms128m', '-Xmx512m']
    run(['java', *flags, '-cp', folders['classes'], 'AuthorsOptimization', '--check'], output / 'checks.json', output / 'checks.stderr')
    metadata['checks'] = json.loads((output / 'checks.json').read_text())
    print('Exact state, rank and derivative checks complete:', metadata['checks'], flush=True)
    (output / 'manifest').write_text('Premain-Class: AuthorsOptimization\n')
    run(['jar', 'cfm', output / 'agent.jar', output / 'manifest', '-C', folders['classes'], '.'], output / 'jar.txt', output / 'jar.stderr')
    flags += ['--add-opens', 'java.base/java.util=ALL-UNNAMED', '--add-opens', 'java.base/java.lang=ALL-UNNAMED',
              '-javaagent:' + str(output / 'agent.jar')]
    metadata['java_flags'] = flags
    run(['java', *flags, '-cp', folders['classes'], 'AuthorsOptimization', '--benchmark', args.trials],
        output / 'benchmark.csv', output / 'benchmark.stderr')
    metadata['benchmark_checks'] = summarize(output, args.trials)
    assert hashes == {str(p.relative_to(ROOT)): digest(p) for p in sources}, 'sources changed during evaluation'
    assert reference_hashes == {name: digest(authors / name) for name in reference_hashes}
    assert not capture(['git', '-C', authors, 'status', '--porcelain', '--untracked-files=no'])
    (output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print('Results:', output, flush=True)


if __name__ == '__main__':
    main()
