#!/usr/bin/env python3
"""Test Section 4.1's batched edits against the current compact implementation.

Builds temporary headers; never changes the installed implementation. Requires
Linux/glibc for the existing allocation benchmark and a C++17 compiler.
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

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "benchmarks/comparison"


def run(args, **kwargs):
    subprocess.run([str(a) for a in args], check=True, **kwargs)


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--trials", type=int, default=5)
    args = parser.parse_args()
    if args.trials < 1:
        parser.error("--trials must be positive")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    data = output / "data"
    data.mkdir(exist_ok=True)
    header = (ROOT / "include/splinesketch/splinesketch.hpp").read_text()
    method = (SOURCE / "experiments/paper_rebalance.inc").read_text()
    start = header.index("  void rebalance(const Snapshot& before) {")
    end = header.index("  void insert_summary(", start)
    compiler = shutil.which("c++")
    if compiler is None:
        raise SystemExit("a C++17 compiler is required")
    flags = ["-std=c++17", "-O3", "-DNDEBUG"]
    hashes = {}
    for variant in ("baseline", "candidate"):
        include = output / variant / "include/splinesketch"
        include.mkdir(parents=True, exist_ok=True)
        content = header if variant == "baseline" else header[:start] + method + header[end:]
        (include / "splinesketch.hpp").write_text(content)
        hashes[variant] = hashlib.sha256(content.encode()).hexdigest()
        for name, source in (("accuracy", ROOT / "benchmarks/accuracy.cpp"),
                             ("fresh", SOURCE / "resize.cpp"),
                             ("local", SOURCE / "local.cpp"),
                             ("invariants", ROOT / "tests/guarantee_invariant_tests.cpp")):
            executable = output / variant / name
            definitions = ["-DSPLINESKETCH_TESTING"] if name == "invariants" else []
            run([compiler, *flags, *definitions, "-I", include.parent, source, "-o", executable])
            if name == "invariants":
                with (output / (variant + "-invariants.txt")).open("w") as stream:
                    run([executable], stdout=stream)
            elif name != "local":
                option = "--details" if name == "accuracy" else "--fresh"
                with (output / (variant + "-" + name + ".csv")).open("w") as stream, \
                     (output / (variant + "-" + name + ".stderr")).open("w") as error:
                    run([executable, option], stdout=stream, stderr=error)
    for trial in range(args.trials):
        order = ("baseline", "candidate") if trial % 2 == 0 else ("candidate", "baseline")
        for variant in order:
            with (output / (variant + "-local-" + str(trial) + ".csv")).open("w") as stream:
                run([output / variant / "local", data], stdout=stream)
    summary = {}
    for variant in ("baseline", "candidate"):
        accuracy = rows(output / (variant + "-accuracy.csv"))
        fresh = rows(output / (variant + "-fresh.csv"))
        timing = [rows(output / (variant + "-local-" + str(t) + ".csv")) for t in range(args.trials)]
        summary[variant] = {
            "workflow_max_percent": {w: max(float(r["max_percent"]) for r in accuracy
                if r["group"] == "case" and r["name"].endswith("/" + w))
                for w in ("direct", "merge4", "resize")},
            "fresh_max_percent": {s: max(float(r["max_percent"]) for r in fresh if r["shape"] == s)
                for s in ("uniform", "clustered", "duplicates", "outliers", "ascending")},
            "cost": {k: {
                "update_ns": statistics.median(statistics.median(float(r["update_ns"])
                    for r in trial if int(r["k"]) == k) for trial in timing),
                "resident_peak_bytes": max(int(r["resident_peak"]) for trial in timing for r in trial if int(r["k"]) == k),
                "update_peak_bytes": max(int(r["update_peak"]) for trial in timing for r in trial if int(r["k"]) == k)}
                for k in (8, 16, 32, 64, 128)}}
    metadata = {
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "authors_revision": "320019fe9d156a42ec6c1edba880d060de8a1b13",
        "header_sha256": hashes,
        "method_sha256": hashlib.sha256(method.encode()).hexdigest(),
        "compiler": subprocess.check_output([compiler, "--version"], text=True).strip(),
        "flags": flags, "trials": args.trials, "platform": platform.platform(),
        "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__).resolve(), ROOT / "benchmarks/accuracy.cpp", SOURCE / "resize.cpp",
                      SOURCE / "local.cpp", ROOT / "benchmarks/accuracy_workload.hpp",
                      ROOT / "tests/guarantee_invariant_tests.cpp")}}
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
