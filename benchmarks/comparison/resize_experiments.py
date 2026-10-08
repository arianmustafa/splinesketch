#!/usr/bin/env python3
"""Reproduce rejected resize experiments against the recorded baseline header.

All experimental headers live in OUTPUT; the repository header is untouched.
Requires C++17, Python 3, and GNU patch. CSV inputs use unseen seeds 3..22.
"""
import argparse
import csv
import datetime
import hashlib
import json
from pathlib import Path
import platform
import statistics
import subprocess

ROOT = Path(__file__).resolve().parents[2]
BASELINE_SHA256 = "15b9d19b6425b7262417de367635a66aaf732dd496aa1261be31ba944263e1cd"
FLAGS = ["-std=c++17", "-O3", "-DNDEBUG"]
VARIANTS = ("baseline", "resize_shape", "resize_exact", "resize_chord")


def run(command, **kwargs):
    subprocess.run([str(x) for x in command], check=True, **kwargs)


def load(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def summarize(output):
    baseline = load(output / "baseline/holdout.csv")
    if len(baseline) != 500:
        raise RuntimeError("incomplete baseline grid")
    key = lambda r: (r["k"], r["shape"], r["seed"])
    original = {key(r): r for r in baseline}
    if len(original) != 500:
        raise RuntimeError("duplicate baseline cases")
    summary = []
    for variant in VARIANTS:
        rows = load(output / variant / "holdout.csv")
        if len(rows) != 500 or {key(r) for r in rows} != set(original):
            raise RuntimeError("incomplete variant grid")
        for r in rows:
            before = original[key(r)]
            if (r["input_fingerprint"], r["queries"]) != (before["input_fingerprint"], before["queries"]):
                raise RuntimeError("input/query mismatch")
        for dimension, values in (("shape", ("uniform", "clustered", "duplicates", "outliers", "ascending")),
                                  ("k", ("8", "16", "32", "64", "128"))):
            for value in values:
                group = [r for r in rows if r[dimension] == value]
                errors = [float(r["max_percent"]) for r in group]
                deltas = [float(r["max_percent"]) - float(original[key(r)]["max_percent"]) for r in group]
                summary.append([variant, dimension, value, len(group), max(errors), statistics.median(errors),
                                sum(d < -1e-9 for d in deltas), sum(d > 1e-9 for d in deltas)])
    with (output / "holdout-summary.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(["variant", "dimension", "value", "cases", "worst_percent", "median_case_max_percent", "improved", "worsened"])
        writer.writerows(summary)
    with (output / "timing-summary.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(["variant", "from", "to", "median_resize_us"])
        for variant in VARIANTS:
            rows = load(output / variant / "resize-timing.csv")
            for capacity in (16, 32, 64, 128, 256, 512):
                group = [r for r in rows if int(r["from"]) == capacity]
                writer.writerow([variant, capacity, group[0]["to"], statistics.median(float(r["resize_us"]) for r in group)])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    header = ROOT / "include/splinesketch/splinesketch.hpp"
    digest = hashlib.sha256(header.read_bytes()).hexdigest()
    if digest != BASELINE_SHA256:
        raise SystemExit("baseline header changed; review and rebase the experiment patches first")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = ROOT / "benchmarks/comparison"
    metadata = {
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "baseline_header_sha256": digest,
        "workload_sha256": hashlib.sha256((ROOT / "benchmarks/accuracy_workload.hpp").read_bytes()).hexdigest(),
        "resize_runner_sha256": hashlib.sha256((source / "resize.cpp").read_bytes()).hexdigest(),
        "compiler": subprocess.check_output(["c++", "--version"], text=True).strip(),
        "flags": FLAGS, "platform": platform.platform(),
        "holdout_seeds": list(range(3, 23)), "observations": 6000,
        "cpu": next((line.split(":", 1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines()
                     if line.startswith("model name")), "unknown"),
        "patch_sha256": {},
    }
    for variant in VARIANTS:
        target = output / variant
        include = target / "include/splinesketch"
        include.mkdir(parents=True, exist_ok=True)
        (include / "splinesketch.hpp").write_bytes(header.read_bytes())
        if variant != "baseline":
            patch = source / "experiments" / (variant + ".patch")
            metadata["patch_sha256"][variant] = hashlib.sha256(patch.read_bytes()).hexdigest()
            run(["patch", "--batch", "--forward", "--fuzz=0", "-p1", "-d", target, "-i", patch])
        for filename in ("resize", "accuracy"):
            cpp = source / "resize.cpp" if filename == "resize" else ROOT / "benchmarks/accuracy.cpp"
            run(["c++", *FLAGS, "-I", target / "include", cpp, "-o", target / filename])
        for command, filename in (([target / "resize"], "holdout.csv"),
                                  ([target / "resize", "--timing"], "resize-timing.csv")):
            with (target / filename).open("w") as stream:
                run(command, stdout=stream)
        with (target / "accuracy.csv").open("w") as stream, (target / "worst.txt").open("w") as error:
            run([target / "accuracy", "--details"], stdout=stream, stderr=error)
        print("Completed", variant, flush=True)
    for variant in ("resize_exact", "resize_chord"):
        for name in ("splinesketch_tests", "exception_safety_tests"):
            target = output / variant
            run(["c++", "-std=c++17", "-O2", "-I", target / "include",
                 ROOT / "tests" / (name + ".cpp"), "-o", target / name])
            with (target / (name + ".txt")).open("w") as stream:
                run([target / name], stdout=stream)
    target = output / "resize_exact"
    run(["c++", "-std=c++17", "-O2", "-I", target / "include", source / "resize_metric.cpp",
         "-o", target / "resize_metric"])
    with (target / "resize_metric.txt").open("w") as stream:
        run([target / "resize_metric"], stdout=stream)
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    summarize(output)
    print("Summary:", output / "holdout-summary.csv")


if __name__ == "__main__":
    main()
