#!/usr/bin/env python3
"""Compare compact code, the paper baseline, its optimizations, and authors' Java.

Requires C++17 and Linux/glibc; Java comparisons additionally require JDK 21
and the unmodified authors' checkout at the revision below.
"""
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
SOURCE = ROOT / "benchmarks/comparison"
REVISION = "320019fe9d156a42ec6c1edba880d060de8a1b13"
JAVA_GRID = "6,8,12,16,19,24,32,40,48,64,80,96,128,144,152,160,168,192,256"


def run(args, **kwargs):
    subprocess.run([str(a) for a in args], check=True, **kwargs)


def capture(args):
    return subprocess.check_output([str(a) for a in args], text=True, stderr=subprocess.STDOUT).strip()


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def output_run(args, path):
    with path.open("w") as dest, path.with_suffix(".stderr").open("w") as err:
        run(args, stdout=dest, stderr=err)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--authors", type=Path)
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
    metadata = {
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "compiler": capture([compiler, "--version"]), "flags": flags,
        "platform": platform.platform(), "trials": args.trials,
        "timing_order": "alternating variants; no simultaneous benchmark jobs",
        "timing_summary": "median across workloads of per-workload trial medians",
        "authors_revision": REVISION if args.authors else None,
        "seed_ranges": {"original": [0, 2], "fresh": [23, 42], "heldout": [43, 62]},
    }
    with tempfile.TemporaryDirectory(prefix="splinesketch-paper-repro-", dir="/tmp") as tmp:
        build = Path(tmp)
        reference = build / "reference"
        headers = reference / "include/splinesketch"
        headers.mkdir(parents=True)
        for name in ("splinesketch.hpp", "paper_splinesketch.hpp"):
            shutil.copyfile(ROOT / "include/splinesketch" / name, headers / name)
        # Reconstruct the original baseline across all optimization phases.
        run(["patch", "--batch", "-R", "-p1", "-i", SOURCE / "experiments/paper_certificate_refinement.patch"], cwd=reference)
        run(["patch", "--batch", "-R", "-p1", "-i", SOURCE / "experiments/paper_certificate.patch"], cwd=reference)
        run(["patch", "--batch", "-R", "-p1", "-i", SOURCE / "experiments/paper_heap_optimization.patch"], cwd=reference)
        run(["patch", "--batch", "-R", "-p1", "-i", SOURCE / "experiments/paper_update_optimization.patch"], cwd=reference)
        run(["patch", "--batch", "-R", "-p1", "-i", SOURCE / "experiments/paper_optimization.patch"], cwd=reference)
        metadata["reference_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in headers.iterdir()}
        definitions = ["-DSPLINESKETCH_PREVIOUS_CORE=" + json.dumps(str(headers / "splinesketch.hpp")),
                       "-DSPLINESKETCH_PREVIOUS_PAPER=" + json.dumps(str(headers / "paper_splinesketch.hpp"))]
        run([compiler, *flags, "-I", ROOT / "include", *definitions,
             SOURCE / "paper_equivalence.cpp", "-o", build / "equivalence"])
        output_run([build / "equivalence"], output / "equivalence.txt")
        variants = ("compact", "reference", "paper", "theoretical")
        for variant in variants:
            include = reference / "include" if variant == "reference" else ROOT / "include"
            definitions = [] if variant == "compact" else ["-DSPLINESKETCH_PAPER"]
            if variant == "theoretical":
                definitions.append("-DSPLINESKETCH_PAPER_THEORETICAL")
            for target, path in (("accuracy", ROOT / "benchmarks/accuracy.cpp"),
                                 ("resize", SOURCE / "resize.cpp"), ("local", SOURCE / "local.cpp")):
                if variant == "theoretical" and target == "local":
                    continue
                run([compiler, *flags, *definitions, "-I", include, path, "-o", build / (variant + "-" + target)])
            output_run([build / (variant + "-accuracy"), "--details"], output / (variant + "-accuracy.csv"))
            for group in ("fresh", "heldout"):
                output_run([build / (variant + "-resize"), "--" + group], output / (variant + "-" + group + ".csv"))
        assert (output / "reference-accuracy.csv").read_bytes() == (output / "paper-accuracy.csv").read_bytes()
        for group in ("fresh", "heldout"):
            old, new = rows(output / ("reference-" + group + ".csv")), rows(output / ("paper-" + group + ".csv"))
            assert len(old) == len(new) == 500
            for a, b in zip(old, new):
                assert all(a[k] == b[k] for k in a if k not in ("shrink_us", "grow_us"))
        data = build / "data"
        data.mkdir()
        variants = variants[:3]
        for trial in range(args.trials):
            order = variants if trial % 2 == 0 else tuple(reversed(variants))
            for variant in order:
                output_run([build / (variant + "-local"), data, "--grid"],
                           output / (variant + "-local-" + str(trial) + ".csv"))
        with (output / "timing-summary.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["variant", "k", "resident_peak", "update_peak", "max_percent",
                             "update_ns", "rank_ns", "quantile_ns"])
            for variant in variants:
                samples = [r for trial in range(args.trials)
                           for r in rows(output / (variant + "-local-" + str(trial) + ".csv"))]
                for k in sorted({int(r["k"]) for r in samples}):
                    subset = [r for r in samples if int(r["k"]) == k]
                    keys = sorted({(r["shape"], r["seed"]) for r in subset})
                    metrics = [statistics.median(statistics.median(float(r[field]) for r in subset
                        if (r["shape"], r["seed"]) == key) for key in keys)
                        for field in ("update_ns", "rank_ns", "quantile_ns")]
                    writer.writerow([variant, k, max(int(r["resident_peak"]) for r in subset),
                        max(int(r["update_peak"]) for r in subset),
                        max(float(r["max_percent"]) for r in subset), *metrics])
        if args.authors:
            authors = args.authors.resolve()
            if capture(["git", "-C", authors, "rev-parse", "HEAD"]) != REVISION:
                raise SystemExit("authors' checkout must be pinned to " + REVISION)
            if capture(["git", "-C", authors, "status", "--porcelain", "--untracked-files=no"]):
                raise SystemExit("authors' checkout must be unmodified")
            classes = build / "classes"
            classes.mkdir()
            run(["javac", "-d", classes, authors / "SplineSketch.java", authors / "SplineSketchMG.java",
                 SOURCE / "AuthorsComparison.java"])
            (build / "manifest").write_text("Premain-Class: AuthorsComparison\n")
            run(["jar", "cfm", build / "agent.jar", build / "manifest", "-C", classes, "."])
            java_flags = ["-Xms128m", "-Xmx512m", "--add-opens", "java.base/java.util=ALL-UNNAMED",
                          "--add-opens", "java.base/java.lang=ALL-UNNAMED"]
            output_run(["java", *java_flags, "-javaagent:" + str(build / "agent.jar"), "-cp", classes,
                        "AuthorsComparison", data, JAVA_GRID], output / "authors.csv")
            metadata["java"] = capture(["java", "-version"])
            metadata["java_flags"] = java_flags
            metadata["java_capacity_grid"] = JAVA_GRID
            summary = rows(output / "timing-summary.csv")
            cpp = rows(output / "paper-local-0.csv")
            java = rows(output / "authors.csv")
            with (output / "matched-budgets.csv").open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(["budget", "variant", "k", "resident_peak", "max_percent"])
                for k in (8, 16, 32, 64, 128):
                    baseline = next(r for r in summary if r["variant"] == "compact" and int(r["k"]) == k)
                    budget = int(baseline["resident_peak"])
                    writer.writerow([budget, "compact", k, budget, baseline["max_percent"]])
                    for name, source in (("paper", cpp), ("java_mg", [r for r in java if r["variant"] == "java_mg"])):
                        groups = [[r for r in source if int(r["k"]) == cap]
                                  for cap in sorted({int(r["k"]) for r in source})]
                        fits = [group for group in groups if max(int(r["resident_peak"]) for r in group) <= budget]
                        if not fits:
                            writer.writerow([budget, name, "none", "", ""])
                            continue
                        group = fits[-1]
                        writer.writerow([budget, name, group[0]["k"], max(int(r["resident_peak"]) for r in group),
                                         max(float(r["max_percent"]) for r in group)])
        metadata["input_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(data.iterdir())}
    sources = [ROOT / "include/splinesketch/splinesketch.hpp", ROOT / "include/splinesketch/paper_splinesketch.hpp",
               Path(__file__).resolve(), SOURCE / "paper_equivalence.cpp", SOURCE / "experiments/paper_optimization.patch",
               SOURCE / "experiments/paper_update_optimization.patch", SOURCE / "experiments/paper_heap_optimization.patch",
               SOURCE / "experiments/paper_certificate.patch", SOURCE / "experiments/paper_certificate_refinement.patch",
               SOURCE / "local.cpp", SOURCE / "resize.cpp", SOURCE / "AuthorsComparison.java",
               ROOT / "benchmarks/accuracy.cpp", ROOT / "benchmarks/accuracy_workload.hpp",
               ROOT / "benchmarks/sketch_variant.hpp"]
    metadata["source_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print("Results:", output)


if __name__ == "__main__":
    main()
