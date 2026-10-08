#!/usr/bin/env python3
"""Reproduce the fixed-capacity comparison; requires Linux/glibc, C++17, JDK 21+."""
import argparse
import csv
import datetime
import hashlib
import json
import platform
from pathlib import Path
import shutil
import statistics
import subprocess

REVISION = "320019fe9d156a42ec6c1edba880d060de8a1b13"
CAPACITIES = "6,8,9,12,16,19,20,24,32,36,38,40,42,44,45,48,64,72,76,80,84,88,94,96,128,144,152,160,168,176,190,192,256,383"
ROOT = Path(__file__).resolve().parents[2]


def run(args, **kwargs):
    return subprocess.run([str(a) for a in args], check=True, **kwargs)


def capture(args):
    return subprocess.check_output([str(a) for a in args], stderr=subprocess.STDOUT, text=True).strip()


def summarize(output):
    with (output / "comparison-cpp.csv").open() as f:
        cpp = list(csv.DictReader(f))
    with (output / "comparison-java.csv").open() as f:
        java = list(csv.DictReader(f))
    with (output / "matched-budgets.csv").open("w") as f:
        writer = csv.writer(f)
        writer.writerow(["budget_bytes", "variant", "k", "resident_peak", "max_percent",
                         "outlier_max_percent", "median_update_ns", "median_rank_ns"])
        for k in (8, 16, 32, 64, 128):
            local = [r for r in cpp if int(r["k"]) == k]
            budget = max(int(r["resident_peak"]) for r in local)
            groups = [local]
            for variant in ("java", "java_mg"):
                candidates = []
                for jk in sorted({int(r["k"]) for r in java}):
                    group = [r for r in java if r["variant"] == variant and int(r["k"]) == jk]
                    if len(group) != 15:
                        raise RuntimeError("incomplete Java workload grid")
                    if max(int(r["resident_peak"]) for r in group) <= budget:
                        candidates.append(group)
                if not candidates:
                    raise RuntimeError("no reference capacity fits budget")
                groups.append(candidates[-1])
            for group in groups:
                writer.writerow([budget, group[0]["variant"], group[0]["k"],
                                 max(int(r["resident_peak"]) for r in group),
                                 max(float(r["max_percent"]) for r in group),
                                 max(float(r["max_percent"]) for r in group if r["shape"] == "outliers"),
                                 statistics.median(float(r["update_ns"]) for r in group),
                                 statistics.median(float(r.get("rank_ns", r.get("rank_single_ns"))) for r in group)])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("authors", type=Path, help="checkout of the pinned authors' repository")
    parser.add_argument("output", type=Path, help="build, input, and result directory")
    parser.add_argument("--ablations", action="store_true", help="also test three temporary rule changes")
    args = parser.parse_args()
    authors, output = args.authors.resolve(), args.output.resolve()
    if capture(["git", "-C", authors, "rev-parse", "HEAD"]) != REVISION:
        raise SystemExit("authors' checkout must be pinned to " + REVISION)
    if capture(["git", "-C", authors, "status", "--porcelain", "--untracked-files=no"]):
        raise SystemExit("authors' tracked files must be unmodified")
    output.mkdir(parents=True, exist_ok=True)
    data, classes = output / "data", output / "classes"
    data.mkdir(exist_ok=True)
    classes.mkdir(exist_ok=True)
    compiler = shutil.which("c++")
    flags = ["-std=c++17", "-O3", "-DNDEBUG"]
    source = ROOT / "benchmarks/comparison"
    metadata = {
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "authors_revision": REVISION, "authors_license": "MIT",
        "local_head": capture(["git", "-C", ROOT, "rev-parse", "HEAD"]),
        "header_sha256": hashlib.sha256((ROOT / "include/splinesketch/splinesketch.hpp").read_bytes()).hexdigest(),
        "compiler": capture([compiler, "--version"]), "flags": flags,
        "java": capture(["java", "-version"]), "platform": platform.platform(),
        "machine": platform.machine(),
        "numeric_macros": [line for line in capture([compiler, "-dM", "-E", "-x", "c++", "/dev/null"]).splitlines() if any(name in line for name in ("__LDBL_MANT_DIG__", "__LDBL_MAX_EXP__", "__SIZEOF_LONG_DOUBLE__"))],
        "java_flags": ["-Xms128m", "-Xmx512m", "--add-opens", "java.base/java.util=ALL-UNNAMED",
                       "--add-opens", "java.base/java.lang=ALL-UNNAMED"],
        "capacity_grid": CAPACITIES,
    }
    for name in ("local", "outliers"):
        run([compiler, *flags, "-I", ROOT / "include", source / (name + ".cpp"), "-o", output / name])
    with (output / "comparison-cpp.csv").open("w") as f:
        run([output / "local", data], stdout=f)
    run(["javac", "-d", classes, authors / "SplineSketch.java", authors / "SplineSketchMG.java", source / "AuthorsComparison.java"])
    (output / "manifest").write_text("Premain-Class: AuthorsComparison\n")
    run(["jar", "cfm", output / "agent.jar", output / "manifest", "-C", classes, "."])
    with (output / "comparison-java.csv").open("w") as f, (output / "java.stderr").open("w") as err:
        run(["java", *metadata["java_flags"], "-javaagent:" + str(output / "agent.jar"),
             "-cp", classes, "AuthorsComparison", data, CAPACITIES], stdout=f, stderr=err)
    metadata["input_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(data.iterdir())}
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    for option, filename in (([], "outlier-trace.txt"), (["--minimize"], "outlier-witness.data")):
        with (output / filename).open("w") as f:
            run([output / "outliers", *option], stdout=f)
    with (output / "outlier-witness.data").open() as inp, (output / "outlier-witness.txt").open("w") as out:
        run([output / "outliers", "--stdin"], stdin=inp, stdout=out)
    if args.ablations:
        header = (ROOT / "include/splinesketch/splinesketch.hpp").read_text()
        threshold = "nodes_[i].mass <= 2 * static_cast<long double>(count_) / capacity_"
        reserve = "if (joinable < capacity_ / 3 + 2) break;"
        if threshold not in header or reserve not in header:
            raise RuntimeError("ablation source anchors changed")
        for variant in ("baseline", "split-threshold", "join-reserve", "both"):
            content = header
            if variant in ("split-threshold", "both"):
                content = content.replace(threshold, "nodes_[i].mass <= 0.01L * bound()")
            if variant in ("join-reserve", "both"):
                content = content.replace(reserve, "if (joinable == 0) break;")
            include = output / variant / "splinesketch"
            include.mkdir(parents=True, exist_ok=True)
            (include / "splinesketch.hpp").write_text(content)
            executable = output / variant / "accuracy"
            run([compiler, *flags, "-I", include.parent, ROOT / "benchmarks/accuracy.cpp", "-o", executable])
            with (output / (variant + "-accuracy.csv")).open("w") as out, (output / (variant + "-worst.txt")).open("w") as err:
                run([executable, "--details"], stdout=out, stderr=err)
    summarize(output)
    print("Results:", output / "matched-budgets.csv")


if __name__ == "__main__":
    main()
