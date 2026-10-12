#!/usr/bin/env python3
"""Exercise reporting failures, skipped checks and benchmark regressions."""

import importlib.util
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock

spec = importlib.util.spec_from_file_location("ci_report", Path(__file__).resolve().parents[1] / "tools/ci_report.py")
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def sample_csv(timing="100"):
    lines = ["capacity,shape,state,count,width," + ",".join(report.TIMINGS)]
    for capacity in (16, 128, 1024):
        for shape in ("normal", "discrete"):
            for state in ("buffered", "consolidated"):
                lines.append(f"{capacity},{shape},{state},20000,5," + ",".join([timing] * 4))
    return "\n".join(lines) + "\n"


class ReportTests(unittest.TestCase):
    def test_junit_failure_error_skipped_and_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            junit = Path(directory) / "tests.xml"
            junit.write_text('''<testsuites><testsuite>
              <testcase name="ok" time="1.25"/>
              <testcase name="bad" time="2"><failure>assertion</failure></testcase>
              <testcase name="crash"><error>signal</error></testcase>
              <testcase name="skip"><skipped/></testcase>
              <testcase name="disabled" status="notrun"/>
            </testsuite></testsuites>''')
            results = report.test_results(junit)
            self.assertEqual([r["status"] for r in results],
                             ["passed", "failed", "failed", "skipped", "skipped"])
            self.assertEqual(results[0]["seconds"], 1.25)

    def test_empty_junit_is_not_a_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            junit = Path(directory) / "tests.xml"
            junit.write_text("<testsuite/>")
            with self.assertRaises(ValueError):
                report.test_results(junit)

    def test_incomplete_and_duplicate_benchmarks_rejected(self):
        lines = sample_csv().splitlines()
        for invalid in ("\n".join(lines[:-1]), "\n".join(lines + [lines[1]])):
            with self.assertRaises(ValueError):
                report.benchmark_rows(invalid)

    def test_invalid_timings_rejected(self):
        for timing in ("nan", "inf", "0", "-1"):
            with self.assertRaises(ValueError):
                report.benchmark_rows(sample_csv(timing))

    def test_medians_keep_all_workloads_and_detect_state_changes(self):
        samples = [report.benchmark_rows(sample_csv(timing)) for timing in ("100", "300", "200")]
        result = report.medians(samples)
        self.assertEqual(len(result), 12)
        self.assertTrue(all(r["midpoint_select_ns"] == 200 for r in result.values()))
        samples[1][(16, "normal", "buffered")]["width"] = 6
        with self.assertRaises(ValueError):
            report.medians(samples)

    def test_markdown_cells_escape_output(self):
        self.assertEqual(report.cell("<value>|\n"), "&lt;value&gt;&#124; ")

    def test_missing_test_report_explains_failed_setup(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = argparse.Namespace(junit=root / "absent.xml", compiler="g++",
                                      output=root, revision="abc", outcome="failure")
            with mock.patch.object(report, "environment", return_value={"compiler": "test compiler"}), \
                    mock.patch.dict(os.environ, {"GITHUB_STEP_SUMMARY": ""}), contextlib.redirect_stdout(io.StringIO()):
                report.tests(args)
            text = (root / "summary.md").read_text()
            self.assertIn("Tests did not produce a report", text)
            self.assertIn("**failure**", text)
            self.assertIsNone(json.loads((root / "metadata.json").read_text())["tests"])

    def run_benchmark_report(self, root, fail_build=None):
        args = argparse.Namespace(output=root, compiler="g++", revision="head-sha",
                                  base_revision="base-sha", head=root / "head", base=root / "base")

        def process(command, **kwargs):
            if command[0] == "g++":
                label = Path(command[-1]).name.split("-")[0]
                return subprocess.CompletedProcess(command, int(label == fail_build), "build output")
            label = Path(command[0]).name.split("-")[0]
            return subprocess.CompletedProcess(command, 0, sample_csv("150" if label == "head" else "100"), "")

        with mock.patch.object(report, "environment", return_value={"compiler": "test compiler"}), \
                mock.patch.object(report.subprocess, "run", side_effect=process), \
                mock.patch.dict(os.environ, {"GITHUB_STEP_SUMMARY": ""}), contextlib.redirect_stdout(io.StringIO()):
            report.benchmark(args)

    def test_regressions_remain_visible_for_every_workload_and_metric(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.run_benchmark_report(root)
            comparisons = json.loads((root / "comparison.json").read_text())
            self.assertEqual(len(comparisons), 48)
            self.assertTrue(all(r["change_percent"] == 50 for r in comparisons))
            text = (root / "summary.md").read_text()
            self.assertEqual(text.count("+50.0% | 0.67×"), 48)
            self.assertEqual(len(list(root.glob("*.csv"))), 6)

    def test_unavailable_baseline_still_reports_changed_timings(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.run_benchmark_report(root, fail_build="base")
            comparisons = json.loads((root / "comparison.json").read_text())
            self.assertEqual(len(comparisons), 48)
            self.assertTrue(all(r["before_ns"] is None for r in comparisons))
            self.assertIn("Baseline comparison unavailable", (root / "summary.md").read_text())

    def test_changed_build_failure_publishes_error_and_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaises(RuntimeError):
                self.run_benchmark_report(root, fail_build="head")
            self.assertIn("Benchmark failed", (root / "summary.md").read_text())
            self.assertIn("error", json.loads((root / "metadata.json").read_text()))


if __name__ == "__main__":
    unittest.main()
