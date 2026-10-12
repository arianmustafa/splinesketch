#!/usr/bin/env python3
"""Publish CI summaries and retain machine-readable test/benchmark results."""

import argparse
import csv
import html
import io
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import xml.etree.ElementTree as ET


TIMINGS = {
    "spline_rank_with_error_ns": "Spline rank with error",
    "midpoint_rank_ns": "Midpoint rank",
    "spline_quantile_ns": "Spline quantile",
    "midpoint_select_ns": "Midpoint selection",
}


def cell(value):
    return html.escape(str(value)).replace("|", "&#124;").replace("\n", " ")


def publish(output, title, metadata, body):
    output.mkdir(parents=True, exist_ok=True)
    metadata_text = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    (output / "metadata.json").write_text(metadata_text)
    text = f"## {title}\n\n{body}\n"
    (output / "summary.md").write_text(text)
    print(text)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as summary:
            summary.write(text)


def environment(compiler):
    return {"compiler": subprocess.check_output([compiler, "--version"], text=True).splitlines()[0],
            "compiler_flags": os.environ.get("CXXFLAGS", ""),
            "platform": platform.platform(), "runner_arch": platform.machine()}


def test_results(junit):
    cases = list(ET.parse(junit).getroot().iter("testcase"))
    if not cases:
        raise ValueError("JUnit report contains no tests")
    results = []
    for case in cases:
        failed = case.find("failure") is not None or case.find("error") is not None
        skipped = case.find("skipped") is not None or case.get("status") == "notrun"
        results.append({"name": case.get("name", "unnamed"),
                        "status": "failed" if failed else "skipped" if skipped else "passed",
                        "seconds": float(case.get("time", "0"))})
    return results


def tests(args):
    metadata = environment(args.compiler)
    metadata.update(revision=args.revision, job_outcome=args.outcome)
    body = (f"Commit: `{cell(args.revision)}`\n\nCompiler: {cell(metadata['compiler'])}\n\n"
            f"Job outcome before reporting: **{cell(args.outcome)}**\n\n")
    if not args.junit.exists():
        body += "**Tests did not produce a report.** See the configuration/build logs for the failed or interrupted step.\n"
        metadata["tests"] = None
    else:
        results = test_results(args.junit)
        metadata["tests"] = results
        counts = {status: sum(r["status"] == status for r in results)
                  for status in ("passed", "failed", "skipped")}
        seconds = sum(r["seconds"] for r in results)
        body += ("| Passed | Failed | Skipped | Total | Sum of test durations |\n"
                 "| ---: | ---: | ---: | ---: | ---: |\n"
                 f"| {counts['passed']} | {counts['failed']} | {counts['skipped']} | {len(results)} | {seconds:.2f} s |\n\n")
        failures = [r for r in results if r["status"] != "passed"]
        if failures:
            body += "| Test | Result | Duration |\n| --- | --- | ---: |\n"
            for result in failures:
                body += f"| {cell(result['name'])} | {result['status']} | {result['seconds']:.2f} s |\n"
        body += "\nDownload the job artifact for JUnit XML, complete test output, and build logs.\n"
    publish(args.output, "Test results", metadata, body)


def benchmark_rows(text):
    rows = {}
    for row in csv.DictReader(io.StringIO(text)):
        key = (int(row["capacity"]), row["shape"], row["state"])
        if key in rows:
            raise ValueError(f"duplicate benchmark workload: {key}")
        result = {"count": int(row["count"]), "width": int(row["width"])}
        for metric in TIMINGS:
            value = float(row[metric])
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"invalid benchmark timing: {key}, {metric}")
            result[metric] = value
        rows[key] = result
    expected = {(capacity, shape, state) for capacity in (16, 128, 1024)
                for shape in ("normal", "discrete") for state in ("buffered", "consolidated")}
    if set(rows) != expected:
        raise ValueError("benchmark did not report all twelve workloads")
    return rows


def medians(samples):
    first = samples[0]
    for sample in samples:
        if sample.keys() != first.keys():
            raise ValueError("benchmark workload set changed between runs")
        for key in first:
            for field in ("count", "width"):
                if sample[key][field] != first[key][field]:
                    raise ValueError("benchmark state changed between runs")
    return {key: {**first[key], **{metric: statistics.median(s[key][metric] for s in samples)
                                 for metric in TIMINGS}} for key in first}


def benchmark(args):
    args.output.mkdir(parents=True, exist_ok=True)
    metadata = environment(args.compiler)
    metadata.update(baseline=args.base_revision, revision=args.revision,
                    harness_revision=args.revision, repetitions=3,
                    flags=["-std=c++17", "-O3", "-DNDEBUG"])
    body = (f"Baseline: `{cell(args.base_revision)}` → changed: `{cell(args.revision)}`\n\n"
            f"Compiler: {cell(metadata['compiler'])}\n\n"
            "The changed commit's query harness is compiled against each revision's headers with "
            "C++17, `-O3 -DNDEBUG`. Both executables run serially on the same runner, in alternating "
            "order over three repetitions. Each repetition reports a median over five rounds; "
            "the table takes the median of those three results. Times are microseconds per query. "
            "Positive changes mean slower queries. Timing differences are informational.\n\n")
    try:
        executables = {}
        for label, source in (("head", args.head), ("base", args.base)):
            executable = args.output / f"{label}-benchmark"
            command = [args.compiler, *metadata["flags"], "-I", str(source / "include"),
                       str(args.head / "benchmarks/paper_certificate_queries.cpp"), "-o", str(executable)]
            completed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            (args.output / f"{label}-build.log").write_text(completed.stdout)
            if completed.returncode:
                if label == "head":
                    raise RuntimeError("changed benchmark failed to compile; see head-build.log")
                metadata["baseline_unavailable"] = "Current query harness could not compile against baseline headers."
            else:
                executables[label] = executable
        samples = {label: [] for label in executables}
        for repetition in range(3):
            for label in (("base", "head") if repetition % 2 == 0 else ("head", "base")):
                if label not in executables:
                    continue
                completed = subprocess.run([str(executables[label])], text=True,
                                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                (args.output / f"{label}-{repetition}.csv").write_text(completed.stdout)
                (args.output / f"{label}-{repetition}.log").write_text(completed.stderr)
                if completed.returncode:
                    raise RuntimeError(f"{label} benchmark failed with exit {completed.returncode}")
                samples[label].append(benchmark_rows(completed.stdout))
        head = medians(samples["head"])
        base = medians(samples["base"]) if "base" in samples else {}
        if not base:
            body += "**Baseline comparison unavailable:** current query harness is incompatible with baseline headers. See base-build.log.\n\n"
        comparisons = []
        for metric, name in TIMINGS.items():
            body += (f"### {name}\n\n| Capacity | Input | State | Before (µs) | After (µs) | Change | Speedup |\n"
                     "| ---: | --- | --- | ---: | ---: | ---: | ---: |\n")
            for key in sorted(head):
                after = head[key][metric]
                before = base[key][metric] if key in base else None
                change = (after / before - 1) * 100 if before else None
                comparisons.append({"capacity": key[0], "shape": key[1], "state": key[2],
                                    "metric": metric, "before_ns": before, "after_ns": after,
                                    "change_percent": change, "head_width": head[key]["width"],
                                    "base_width": base[key]["width"] if key in base else None})
                values = f"{before / 1000:.3f} | {after / 1000:.3f} | {change:+.1f}% | {before / after:.2f}×" if before else f"N/A | {after / 1000:.3f} | N/A | N/A"
                body += f"| {key[0]} | {key[1]} | {key[2]} | {values} |\n"
            body += "\n"
        (args.output / "comparison.json").write_text(json.dumps(comparisons, indent=2) + "\n")
        body += "Download the artifact for raw CSV samples, comparison JSON, build logs, and runner metadata.\n"
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        metadata["error"] = str(error)
        publish(args.output, "Query benchmarks", metadata, body + f"**Benchmark failed:** {cell(error)}\n")
        raise
    publish(args.output, "Query benchmarks", metadata, body)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    test_parser = commands.add_parser("tests")
    test_parser.add_argument("--junit", type=Path, required=True)
    test_parser.add_argument("--outcome", required=True)
    test_parser.set_defaults(run=tests)
    benchmark_parser = commands.add_parser("benchmark")
    benchmark_parser.add_argument("--head", type=Path, required=True)
    benchmark_parser.add_argument("--base", type=Path, required=True)
    benchmark_parser.add_argument("--base-revision", required=True)
    benchmark_parser.set_defaults(run=benchmark)
    for command in (test_parser, benchmark_parser):
        command.add_argument("--output", type=Path, required=True)
        command.add_argument("--compiler", default="g++")
        command.add_argument("--revision", required=True)
    args = parser.parse_args()
    args.output = args.output.resolve()
    if args.command == "benchmark":
        args.head, args.base = args.head.resolve(), args.base.resolve()
    args.run(args)


if __name__ == "__main__":
    main()
