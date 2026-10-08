#!/usr/bin/env python3
"""Measure threshold certificates against the preceding and current paper sketches."""
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
SOURCE = ROOT / "benchmarks/comparison"
PATCH = SOURCE / "experiments/paper_certificate.patch"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--trials", type=int, default=5)
    args = parser.parse_args()
    if args.trials < 1:
        parser.error("--trials must be positive")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    compiler = shutil.which("c++")
    if not compiler:
        raise SystemExit("C++17 compiler required")
    flags = ["-std=c++17", "-O3", "-DNDEBUG"]
    sources = [ROOT / "include/splinesketch" / name for name in
               ("splinesketch.hpp", "paper_splinesketch.hpp")]
    sources += [SOURCE / "experiments/paper_certificate_refinement.patch", Path(__file__).resolve(), PATCH, SOURCE / "paper_baseline.py",
                SOURCE / "paper_equivalence.cpp", SOURCE / "paper_certificate_quality.cpp",
                SOURCE / "local.cpp", ROOT / "benchmarks/accuracy_workload.hpp",
                ROOT / "benchmarks/sketch_variant.hpp"]
    hashes = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
    metadata = {
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "compiler": capture([compiler, "--version"]), "flags": flags,
        "platform": platform.platform(), "trials": args.trials,
        "timing_order": "alternating before/plain/certified; sequential, no competing benchmark jobs",
        "timing_summary": "median across workloads of per-workload trial medians",
        "suites": {"original": {"observations": 60000, "seeds": [0, 2]},
                   "fresh": {"observations": 60000, "seeds": [71, 73]}},
        "quality_workflows": ["stream", "balanced_merge", "mixed_merge", "resize", "finalize"],
        "source_sha256": hashes,
    }
    with tempfile.TemporaryDirectory(prefix="splinesketch-paper-certificate-", dir="/tmp") as tmp:
        build = Path(tmp)
        reference = build / "reference"
        headers = reference / "include/splinesketch"
        headers.mkdir(parents=True)
        for name in ("splinesketch.hpp", "paper_splinesketch.hpp"):
            shutil.copyfile(ROOT / "include/splinesketch" / name, headers / name)
        run(["patch", "--batch", "-R", "-p1", "-i", SOURCE / "experiments/paper_certificate_refinement.patch"], cwd=reference)
        run(["patch", "--batch", "-R", "-p1", "-i", PATCH], cwd=reference)
        metadata["reference_sha256"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                         for path in headers.iterdir()}
        definitions = ["-DSPLINESKETCH_PREVIOUS_CORE=" + json.dumps(str(headers / "splinesketch.hpp")),
                       "-DSPLINESKETCH_PREVIOUS_PAPER=" + json.dumps(str(headers / "paper_splinesketch.hpp"))]
        modes = [("native", [])]
        if platform.machine() == "x86_64" and "gcc version" in capture([compiler, "-v"]):
            modes.append(("binary64", ["-mlong-double-64"]))
        modes += [(name + "-forced", [*extra, "-DSPLINESKETCH_PAPER_FORCE_HEAPS"])
                  for name, extra in modes[:]]
        for mode, extra in modes:
            binary = build / ("equivalence-" + mode)
            run([compiler, *flags, *extra, "-I", ROOT / "include", *definitions,
                 SOURCE / "paper_equivalence.cpp", "-o", binary])
            output_run([binary], output / ("equivalence-" + mode + ".txt"))
            output_run([binary, "--extended"], output / ("equivalence-extended-" + mode + ".txt"))
        # Quality measurements are untimed and complete before timing begins.
        for mode, extra in modes[:2]:
            if mode.endswith("forced"):
                continue
            binary = build / ("quality-" + mode)
            run([compiler, *flags, *extra, "-I", ROOT / "include",
                 SOURCE / "paper_certificate_quality.cpp", "-o", binary])
            output_run([binary], output / ("quality-original-" + mode + ".csv"))
            output_run([binary, "--fresh"], output / ("quality-fresh-" + mode + ".csv"))
        variants = ("before", "plain", "certified")
        for variant in variants:
            include = reference / "include" if variant == "before" else ROOT / "include"
            definition = "-DSPLINESKETCH_PAPER_CERTIFIED" if variant == "certified" else "-DSPLINESKETCH_PAPER"
            run([compiler, *flags, definition, "-I", include, SOURCE / "local.cpp", "-o", build / variant])
        inputs = {}
        with (output / "timing-summary.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["suite", "variant", "k", "resident_peak", "update_peak", "max_percent",
                             "update_ns", "rank_ns", "quantile_ns"])
            for suite, option in [("original", "--large"), ("fresh", "--certificate-heldout")]:
                data = build / suite
                data.mkdir()
                for trial in range(args.trials):
                    order = variants if trial % 2 == 0 else tuple(reversed(variants))
                    for variant in order:
                        output_run([build / variant, data, option],
                                   output / f"{suite}-{variant}-{trial}.csv")
                    before = rows(output / f"{suite}-before-{trial}.csv")
                    plain = rows(output / f"{suite}-plain-{trial}.csv")
                    certified = rows(output / f"{suite}-certified-{trial}.csv")
                    fields = ("k", "shape", "seed", "median_percent", "p95_percent", "max_percent")
                    assert all(all(a[field] == b[field] for field in fields) for a, b in zip(before, plain))
                    assert len(before) == len(plain) == len(certified)
                    assert all(float(c["max_percent"]) <= float(p["max_percent"])
                               for c, p in zip(certified, plain))
                for variant in variants:
                    samples = [row for trial in range(args.trials)
                               for row in rows(output / f"{suite}-{variant}-{trial}.csv")]
                    for k in sorted({int(row["k"]) for row in samples}):
                        subset = [row for row in samples if int(row["k"]) == k]
                        keys = sorted({(row["shape"], row["seed"]) for row in subset})
                        metrics = [statistics.median(statistics.median(float(row[field]) for row in subset
                                   if (row["shape"], row["seed"]) == key) for key in keys)
                                   for field in ("update_ns", "rank_ns", "quantile_ns")]
                        writer.writerow([suite, variant, k,
                                         max(int(row["resident_peak"]) for row in subset),
                                         max(int(row["update_peak"]) for row in subset),
                                         max(float(row["max_percent"]) for row in subset), *metrics])
                inputs[suite] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                 for path in sorted(data.iterdir())}
        metadata["input_sha256"] = inputs
    assert hashes == {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                      for path in sources}, "sources changed during benchmarking"
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print("Results:", output)


if __name__ == "__main__":
    main()
