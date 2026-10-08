#!/usr/bin/env python3
"""Reproduce the paper-based update optimization with exact output comparisons.

Requires C++17 and Linux/glibc. Reconstructs the preceding implementation by
reversing the checked-in patch; benchmarks execute sequentially.
"""
import argparse
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
PATCH = SOURCE / "experiments/paper_update_optimization.patch"


def main(phase="update"):
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
        raise SystemExit("a C++17 compiler is required")
    flags = ["-std=c++17", "-O3", "-DNDEBUG"]
    patch = SOURCE / ("experiments/paper_" + ("heap" if phase == "heap" else "update") + "_optimization.patch")
    suites = [("normal", "--grid"), ("large", "--large")]
    if phase == "heap":
        suites.append(("large-heldout", "--large-heldout"))
    sources = [ROOT / "include/splinesketch" / name for name in ("splinesketch.hpp", "paper_splinesketch.hpp")]
    sources += [Path(__file__).resolve(), SOURCE / "paper_baseline.py", PATCH,
                SOURCE / "experiments/paper_heap_optimization.patch",
                SOURCE / "experiments/paper_certificate.patch", SOURCE / "experiments/paper_certificate_refinement.patch",
                SOURCE / "paper_equivalence.cpp", SOURCE / "local.cpp",
                ROOT / "benchmarks/accuracy_workload.hpp", ROOT / "benchmarks/sketch_variant.hpp"]
    if phase == "heap":
        sources.append(SOURCE / "paper_heap_optimization.py")
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    metadata = {
        "phase": phase,
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "compiler": capture([compiler, "--version"]), "flags": flags,
        "platform": platform.platform(), "trials": args.trials,
        "timing_order": "alternating before/after/compact; no simultaneous benchmark jobs",
        "timing_summary": "median across workloads of per-workload trial medians",
        "suites": {"normal": {"observations": 6000, "seeds": [0, 2]},
                   "large": {"observations": 60000, "seeds": [0, 2]}},
        "extended_equivalence": {"observations": 24000, "seeds": [63, 66],
                                 "policies": ["practical", "theoretical"],
                                 "workflows": ["direct", "merge", "resize", "finalize"]},
        "source_sha256": hashes,
    }
    if phase == "heap":
        metadata["suites"]["large-heldout"] = {"observations": 60000, "seeds": [67, 69]}
    with tempfile.TemporaryDirectory(prefix="splinesketch-paper-update-", dir="/tmp") as tmp:
        build = Path(tmp)
        reference = build / "reference"
        headers = reference / "include/splinesketch"
        headers.mkdir(parents=True)
        for name in ("splinesketch.hpp", "paper_splinesketch.hpp"):
            shutil.copyfile(ROOT / "include/splinesketch" / name, headers / name)
        run(["patch", "--batch", "-R", "-p1", "-i", SOURCE / "experiments/paper_certificate_refinement.patch"], cwd=reference)
        run(["patch", "--batch", "-R", "-p1", "-i", SOURCE / "experiments/paper_certificate.patch"], cwd=reference)
        if phase == "update":
            run(["patch", "--batch", "-R", "-p1", "-i", SOURCE / "experiments/paper_heap_optimization.patch"], cwd=reference)
        run(["patch", "--batch", "-R", "-p1", "-i", patch], cwd=reference)
        metadata["reference_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in headers.iterdir()}
        definitions = ["-DSPLINESKETCH_PREVIOUS_CORE=" + json.dumps(str(headers / "splinesketch.hpp")),
                       "-DSPLINESKETCH_PREVIOUS_PAPER=" + json.dumps(str(headers / "paper_splinesketch.hpp"))]
        modes = [("native", [])]
        # This GCC/x86 ABI mode exercises the separate binary64 numeric path.
        if platform.machine() == "x86_64" and "gcc version" in capture([compiler, "-v"]):
            modes.append(("binary64", ["-mlong-double-64"]))
        if phase == "heap":
            modes += [(name + "-forced", [*extra, "-DSPLINESKETCH_PAPER_FORCE_HEAPS"])
                      for name, extra in modes[:]]
        for mode, extra in modes:
            binary = build / ("equivalence-" + mode)
            run([compiler, *flags, *extra, "-I", ROOT / "include", *definitions,
                 SOURCE / "paper_equivalence.cpp", "-o", binary])
            output_run([binary], output / ("equivalence-" + mode + ".txt"))
            output_run([binary, "--extended"], output / ("equivalence-extended-" + mode + ".txt"))
        for variant in ("before", "after", "compact"):
            include = reference / "include" if variant == "before" else ROOT / "include"
            definitions = [] if variant == "compact" else ["-DSPLINESKETCH_PAPER"]
            run([compiler, *flags, *definitions, "-I", include, SOURCE / "local.cpp",
                 "-o", build / variant])
        inputs = {}
        with (output / "timing-summary.csv").open("w", newline="") as stream:
            import csv
            writer = csv.writer(stream)
            writer.writerow(["suite", "variant", "k", "resident_peak", "update_peak", "max_percent",
                             "update_ns", "rank_ns", "quantile_ns"])
            for suite, option in suites:
                data = build / suite
                data.mkdir()
                variants = ("before", "after", "compact")
                for trial in range(args.trials):
                    order = variants if trial % 2 == 0 else tuple(reversed(variants))
                    for variant in order:
                        output_run([build / variant, data, option],
                                   output / (suite + "-" + variant + "-" + str(trial) + ".csv"))
                    before = rows(output / (suite + "-before-" + str(trial) + ".csv"))
                    after = rows(output / (suite + "-after-" + str(trial) + ".csv"))
                    fields = ("k", "shape", "seed", "median_percent", "p95_percent", "max_percent")
                    assert len(before) == len(after)
                    assert all(all(a[k] == b[k] for k in fields) for a, b in zip(before, after))
                for variant in variants:
                    samples = [r for trial in range(args.trials)
                               for r in rows(output / (suite + "-" + variant + "-" + str(trial) + ".csv"))]
                    for k in sorted({int(r["k"]) for r in samples}):
                        subset = [r for r in samples if int(r["k"]) == k]
                        keys = sorted({(r["shape"], r["seed"]) for r in subset})
                        metrics = [statistics.median(statistics.median(float(r[field]) for r in subset
                            if (r["shape"], r["seed"]) == key) for key in keys)
                            for field in ("update_ns", "rank_ns", "quantile_ns")]
                        writer.writerow([suite, variant, k,
                            max(int(r["resident_peak"]) for r in subset),
                            max(int(r["update_peak"]) for r in subset),
                            max(float(r["max_percent"]) for r in subset), *metrics])
                inputs[suite] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(data.iterdir())}
        metadata["input_sha256"] = inputs
    assert hashes == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}, \
        "sources changed while the benchmark was running"
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print("Results:", output)


if __name__ == "__main__":
    main()
