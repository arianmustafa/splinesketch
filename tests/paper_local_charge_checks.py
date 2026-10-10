"""Exact checks of frozen-interval charges, including a failed shrink lemma.

The finite trace checks support the separate mathematical argument; they do
not prove a universal error rate. Observer calls are inserted into a temporary
header copy and never change the repository's production headers.
"""
import argparse
import bisect
from fractions import Fraction as F
import json
import math
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile


def number(bits):
    return F(struct.unpack("!d", struct.pack("!Q", bits))[0])


def encoding(value):
    return struct.unpack("!Q", struct.pack("!d", float(value)))[0]


def dyadic(text):
    whole, digits = text.split(":")
    return F(int(whole)) + F(int(digits or "0", 2), 1 << len(digits))


def arithmetic_horizon(count, capacity, precision):
    """Integer-only structural upper bounds; roundoff comparisons stay exact."""
    epoch_end, epochs = 5 * capacity, 1
    while epoch_end < count:
        epoch_end = min((1 << 64) - 1, epoch_end + (epoch_end + 3) // 4)
        epochs += 1
    pass_bound = count // capacity + (4 * capacity * count.bit_length() + 2) // 3
    operation_bound = 3 * pass_bound + capacity * (2 * epochs + 3)
    theta = F(10 * (capacity + 3) * operation_bound, 1 << precision)
    return epochs, pass_bound, operation_bound, theta


def instrument(root, destination, force_scans=False):
    target = destination / "splinesketch"
    target.mkdir(parents=True)
    for file in (root / "include/splinesketch").glob("*.hpp"):
        shutil.copy2(file, target / file.name)
    path = target / "paper_splinesketch.hpp"
    text = path.read_text()
    edits = [
        ("  void initialize(const Snapshot& before) {",
         "  void initialize(const Snapshot& before) {\n    PaperCertificateInspector::begin(*this, before);"),
        ("        const auto at = std::upper_bound(before.sums.begin(), before.sums.end(), offset);",
         "        PaperCertificateInspector::selected(i, offset, total, capacity_);\n        const auto at = std::upper_bound(before.sums.begin(), before.sums.end(), offset);"),
        ("    const long double value = std::clamp(before.rank(mid), left, right);\n    nodes_.insert",
         "    const long double value = std::clamp(before.rank(mid), left, right);\n    PaperCertificateInspector::born(*this, before, mid, nodes_[i - 1].x, nodes_[i].x, value, left, right, nodes_[i].mass);\n    nodes_.insert"),
        ("      const auto value = std::clamp(before->rank(mid), low, high);\n      const auto affected",
         "      const auto value = std::clamp(before->rank(mid), low, high);\n      PaperCertificateInspector::born(sketch, *before, mid, sketch.nodes_[left].x, sketch.nodes_[right].x, value, low, high, sketch.nodes_[right].mass);\n      const auto affected"),
        ("    void reestimate(const Snapshot& snapshot) {\n      before =",
         "    void reestimate(const Snapshot& snapshot) {\n      PaperCertificateInspector::begin(sketch, snapshot);\n      before ="),
        ("    const Snapshot& before = supplied ? snapshot : local;\n    if (use_heap_backend())",
         "    const Snapshot& before = supplied ? snapshot : local;\n    PaperCertificateInspector::begin(*this, before);\n    if (use_heap_backend())"),
        ("      if (materialize) finish();",
         "      PaperCertificateInspector::end(sketch, *before, snapshot_nodes());\n      if (materialize) finish();"),
        ("  void join_at(std::size_t i) {",
         "  void join_at(std::size_t i) {\n    PaperCertificateInspector::joined(*this, nodes_[i].mass + nodes_[i + 1].mass, 0.75L * bound());"),
        ("    void join(Id id) {",
         "    void join(Id id) {\n      PaperCertificateInspector::joined(sketch, sketch.nodes_[id].mass + sketch.nodes_[links[id].next].mass, 0.75L * sketch.bound());"),
        ("    // Restore slopes before another batch snapshot or public query uses them.\n    rebuild(nodes_);",
         "    // Restore slopes before another batch snapshot or public query uses them.\n    rebuild(nodes_);\n    PaperCertificateInspector::end(*this, before, nodes_);"),
        ("    rebuild(result.nodes_);",
         "    rebuild(result.nodes_);\n    PaperCertificateInspector::merged(result, *this, other);"),
    ]
    for old, new in edits:
        assert text.count(old) == 1, "observer anchor changed: " + old
        text = text.replace(old, new)
    start = text.index("  void initialize(const Snapshot& before) {")
    end = text.index("  void reestimate(const Snapshot& before) {", start)
    region = text[start:end]
    anchor = "    clear_protection();\n  }"
    assert region.count(anchor) == 1
    region = region.replace(anchor, "    clear_protection();\n    PaperCertificateInspector::end(*this, before, nodes_);\n  }")
    text = text[:start] + region + text[end:]
    if force_scans:
        anchor = "  bool use_heap_backend() const {"
        assert text.count(anchor) == 1
        text = text.replace(anchor, anchor + "\n    return false; // Test the production scanner on the same stream.")
    path.write_text(text)


def nodes(rows):
    return [(number(x), dyadic(h), low, high, atom, identity)
            for x, h, low, high, atom, identity in rows]


def round_binary(value, precision, min_quantum):
    """Exact nearest/ties-even rounding, including the subnormal grid."""
    assert value >= 0
    if not value:
        return F(0)
    exponent = value.numerator.bit_length() - value.denominator.bit_length()
    if value < F(2) ** exponent:
        exponent -= 1
    quantum = F(2) ** max(exponent - precision + 1, min_quantum)
    scaled = value / quantum
    integer, remainder = divmod(scaled.numerator, scaled.denominator)
    twice = 2 * remainder
    if twice > scaled.denominator or (twice == scaled.denominator and integer % 2):
        integer += 1
    return integer * quantum


def check_operation_cases(binary, binary64=False):
    states = merges = samples = dips = rank_nodes = queries = useful_queries = 0
    tighter_passes = 0
    directional_useful = tighter_queries = 0
    resize_phase = None
    resize_histories = amplified_histories = 0
    pure_projections = exact_refinements = union_checks = interval_queries = 0
    projection_source = union_sources = None
    max_theta = F(0)
    credits, directional = {}, {}
    credit = over = under = F(0)
    # This second induction uses exact local transfer defects. It needs no
    # assumed per-edit rounding exponent or small arithmetic horizon.
    transfers = {}
    transfer_over = transfer_under = F(0)
    transfer_queries = transfer_useful = 0
    expected_transfer_nodes = None

    def grid(rows):
        return [(number(row[0]), dyadic(row[1])) for row in rows]

    def cell_span(source, x):
        if not source or x <= source[0][0] or x >= source[-1][0]:
            return F(0)
        at = bisect.bisect_left([point for point, _ in source], x)
        return F(0) if source[at][0] == x else source[at][1] - source[at - 1][1]

    def max_span(source):
        return max((b[1] - a[1] for a, b in zip(source, source[1:])
                    if math.nextafter(float(a[0]), math.inf) < float(b[0])), default=F(0))

    def offsets(source, x, sample):
        """One-sided cell charges, without using the true CDF at the sample."""
        if not source or x < source[0][0]:
            assert sample == 0
            return F(0), F(0)
        if x >= source[-1][0]:
            assert sample == source[-1][1]
            return F(0), F(0)
        at = bisect.bisect_left([point for point, _ in source], x)
        if source[at][0] == x:
            assert sample == source[at][1]
            return F(0), F(0)
        low, high = source[at - 1][1], source[at][1]
        assert low <= sample <= high, "compiled sample escaped its prefix bracket"
        return sample - low, high - sample

    def numeric_budget(scale, cost):
        # Twice the conservation exponent also covers sampled additions and
        # reconstruction at individual cuts, rather than only the last cut.
        rho = F(2 * cost, 1 << precision)
        assert rho <= F(1, 2)
        return scale * rho / (1 - rho)

    def envelope(rows):
        return [(number(row[0]), *row[-3:]) for row in rows]

    def interval(source, x):
        if not source or x < source[0][0]:
            return 0, 0
        if x >= source[-1][0]:
            return source[-1][2], source[-1][2]
        at = bisect.bisect_left([node[0] for node in source], x)
        if source[at][0] == x:
            return source[at][1], source[at][2]
        return source[at - 1][1], source[at][2] - source[at][3]

    def partition(*sources):
        points = {-math.inf, math.inf}
        for source in sources:
            for x, *_ in source:
                value = float(x)
                points.update((value, math.nextafter(value, -math.inf), math.nextafter(value, math.inf)))
        return [F(x) if math.isfinite(x) else x for x in sorted(points)]

    process = subprocess.Popen([str(binary), "--operation-cases"], stdout=subprocess.PIPE, text=True)
    try:
        for line in process.stdout:
            row = json.loads(line)
            if row["kind"] == "format":
                precision, min_quantum = row["precision"], row["min_quantum"]
            elif row["kind"] == "rank_context":
                credit = sum((credits.get(parent, F(0)) for parent in row["parents"]), F(0))
                parents = [directional.get(parent, (F(0), F(0))) for parent in row["parents"]]
                over = sum((pair[0] for pair in parents), F(0))
                under = sum((pair[1] for pair in parents), F(0))
                transfer_parents = [transfers.get(parent, (F(0), F(0))) for parent in row["parents"]]
                transfer_over = sum(pair[0] for pair in transfer_parents)
                transfer_under = sum(pair[1] for pair in transfer_parents)
            elif row["kind"] == "rank_pass":
                projection_source = envelope(row["source"]) if row["released"] == 0 else None
                if resize_phase is not None:
                    assert row["released"] == 0, "unchanged-count resize released observations"
                    resize_phase["passes"] += 1
                source = grid(row["source"])
                positive, negative = [F(0)], [F(0)]
                expected_transfer_nodes = []
                for x, stored, sample, incoming in row["transfers"]:
                    stored, sample = dyadic(stored), dyadic(sample)
                    pos, neg = offsets(source, number(x), sample)
                    defect = stored - sample - incoming
                    positive.append(pos + defect)
                    negative.append(neg - defect)
                    expected_transfer_nodes.append((x, stored))
                transfer_over += max(positive)
                transfer_under += max(negative)
                scale = (source[-1][1] if source else F(0)) + row["released"]
                assert len(row["births"]) == row["splits"]
                spatial = max((cell_span(source, number(x)) for x, _ in row["births"]), default=F(0))
                coarse = max_span(source) if row["splits"] else F(0)
                assert spatial <= coarse
                tighter_passes += spatial < coarse
                credit += spatial
                sampling = [offsets(source, number(x), dyadic(sample)) for x, sample in row["births"]]
                numeric = numeric_budget(scale, row["cost"])
                credit += numeric
                over += max((pair[0] for pair in sampling), default=F(0)) + numeric
                under += max((pair[1] for pair in sampling), default=F(0)) + numeric
            elif row["kind"] == "merge_samples":
                union_sources = envelope(row["left"]), envelope(row["right"])
                left, right = grid(row["left"]), grid(row["right"])
                positive, negative = [F(0)], [F(0)]
                expected_transfer_nodes = []
                for x, values in zip(row["cuts"], row["samples"]):
                    stored, a, b = dyadic(values[2]), dyadic(values[3]), dyadic(values[4])
                    ap, an = offsets(left, number(x), a)
                    bp, bn = offsets(right, number(x), b)
                    defect = stored - a - b
                    positive.append(ap + bp + defect)
                    negative.append(an + bn - defect)
                    expected_transfer_nodes.append((x, stored))
                transfer_over += max(positive)
                transfer_under += max(negative)
                # Every union coordinate is retained in at least one input.
                spatial = max((cell_span(left, number(x)) + cell_span(right, number(x))
                               for x in row["cuts"]), default=F(0))
                assert spatial <= max(max_span(left), max_span(right))
                scale = (left[-1][1] if left else F(0)) + (right[-1][1] if right else F(0))
                numeric = numeric_budget(scale, row["cost"])
                credit += spatial + numeric
                sampling = []
                assert len(row["cuts"]) == len(row["samples"])
                for x, values in zip(row["cuts"], row["samples"]):
                    a = offsets(left, number(x), dyadic(values[3]))
                    b = offsets(right, number(x), dyadic(values[4]))
                    sampling.append((a[0] + b[0], a[1] + b[1]))
                over += max((pair[0] for pair in sampling), default=F(0)) + numeric
                under += max((pair[1] for pair in sampling), default=F(0)) + numeric
                previous = prefix = F(0)
                last_sample = F(0)
                for sample, mass, stored_prefix, _, _ in row["samples"]:
                    sample, mass, stored_prefix = map(dyadic, (sample, mass, stored_prefix))
                    dips += sample < last_sample
                    target = max(previous, sample)
                    assert mass == round_binary(target - previous, precision, min_quantum)
                    prefix = round_binary(prefix + mass, precision, min_quantum)
                    assert stored_prefix == prefix
                    previous, last_sample = target, sample
                    samples += 1
                merges += 1
            elif row["kind"] == "ledger_state":
                theta = F(row["q"], 1 << precision)
                assert theta <= F(1, 2), "operation-aware arithmetic horizon exceeded"
                absorbed, prefix = row["absorbed"], dyadic(row["top"])
                assert absorbed <= row["count"]
                assert abs(prefix - absorbed) <= absorbed * theta / (1 - theta)
                assert 0 <= over <= credit and 0 <= under <= credit
                for _, stored, exact, *_ in row["nodes"]:
                    assert abs(dyadic(stored) - exact) <= credit, "retained-cut rank ledger failed"
                    assert -under <= dyadic(stored) - exact <= over, "directional retained-cut ledger failed"
                    assert -transfer_under <= dyadic(stored) - exact <= transfer_over, "exact transfer retained-cut proof failed"
                    rank_nodes += 1
                if row["stage"] in ("pass", "merge"):
                    assert expected_transfer_nodes == [(node[0], dyadic(node[1])) for node in row["nodes"]]
                    expected_transfer_nodes = None
                out_envelope = envelope(row["nodes"])
                if row["stage"] == "pass" and projection_source is not None:
                    refinement = {node[0] for node in projection_source} <= {node[0] for node in out_envelope}
                    for x in partition(projection_source, out_envelope):
                        old, new = interval(projection_source, x), interval(out_envelope, x)
                        assert new[0] <= old[0] and old[1] <= new[1], "pure pass narrowed its source envelope"
                        if refinement:
                            assert new == old, "retained-cut refinement changed integer uncertainty"
                        interval_queries += 1
                    pure_projections += 1
                    exact_refinements += refinement
                    projection_source = None
                elif row["stage"] == "merge":
                    assert union_sources is not None
                    left, right = union_sources
                    for x in partition(left, right, out_envelope):
                        a, b = interval(left, x), interval(right, x)
                        assert interval(out_envelope, x) == (a[0] + b[0], a[1] + b[1]), "union construction widened source intervals"
                        interval_queries += 1
                    union_checks += 1
                    union_sources = None
                curve_bound = credit + max_span(grid(row["nodes"]))
                last_nodes = row["nodes"]
                if row["stage"] == "complete":
                    credits[row["id"]] = credit
                    directional[row["id"]] = over, under
                    transfers[row["id"]] = transfer_over, transfer_under
                count = row["count"]
                last_q = row["q"]
                if resize_phase is not None:
                    assert count == resize_phase["count"] and absorbed == resize_phase["absorbed"]
                max_theta = max(max_theta, theta)
                states += 1
            elif row["kind"] == "rank_queries":
                # Held keys and buffer entries are added exactly in the
                # mathematical CDF; bound their compiled additions/conversions.
                rho = F(row["additions"] + 2, 1 << precision)
                assert rho <= F(1, 2)
                rounding = (prefix + count) * rho / (1 - rho)
                rounding += F(count, 1 << 53) * (1 + F(1, 1 << precision))
                certified_rounding = F(count, 1 << 53)
                source = grid(last_nodes)
                for raw, certified, exact, x, sample, pending, lower, upper in row["queries"]:
                    raw_error, certified_error = dyadic(raw) - exact, dyadic(certified) - exact
                    if resize_phase is not None and x == resize_phase["query"]:
                        assert exact == resize_phase["truth"]
                        resize_phase["error"] = abs(certified_error)
                    assert abs(raw_error) <= curve_bound + rounding, "raw rank ledger failed"
                    assert abs(certified_error) <= curve_bound + rounding + certified_rounding
                    coordinate = struct.unpack("!d", struct.pack("!Q", x))[0]
                    if math.isfinite(coordinate):
                        coordinate = number(x)
                    pos, neg = offsets(source, coordinate, dyadic(sample))
                    # These costs depend only on sketch state, not exact truth.
                    defect = dyadic(raw) - dyadic(sample) - pending
                    raw_positive = max(F(0), transfer_over + pos + defect)
                    raw_negative = max(F(0), transfer_under + neg - defect)
                    assert -raw_negative <= raw_error <= raw_positive, "exact transfer raw proof failed"
                    cert_positive = max(raw_positive, F(float(lower)) - lower)
                    cert_negative = max(raw_negative, upper - F(float(upper)))
                    assert lower <= exact <= upper
                    assert -cert_negative <= certified_error <= cert_positive, "exact transfer certified proof failed"
                    transfer_queries += 1
                    transfer_useful += max(cert_positive, cert_negative) < count
                    positive, negative = over + pos + rounding, under + neg + rounding
                    assert -negative <= raw_error <= positive, "directional raw rank ledger failed"
                    assert -negative - certified_rounding <= certified_error <= positive + certified_rounding
                    assert max(positive, negative) <= curve_bound + rounding
                    directional_useful += max(positive, negative) + certified_rounding < count
                    tighter_queries += max(positive, negative) < curve_bound + rounding
                    queries += 1
                    useful_queries += curve_bound + rounding + certified_rounding < count
            elif row["kind"] == "resize_history":
                if row["begin"]:
                    assert resize_phase is None
                    resize_phase = {"count": count, "absorbed": absorbed, "start_q": last_q,
                                    "passes": 0, "query": row["query"], "truth": row["truth"],
                                    "source_bound": dyadic(row["source_bound"])}
                    assert abs(dyadic(row["source_rank"]) - row["truth"]) <= resize_phase["source_bound"]
                else:
                    assert resize_phase is not None and resize_phase["passes"] >= 24
                    assert last_q > resize_phase["start_q"]
                    # Observe this known accuracy loss; do not require future
                    # implementations to retain it for the regression to pass.
                    numerical = count * max_theta / (1 - max_theta)
                    amplified_histories += resize_phase["error"] > resize_phase["source_bound"] + numerical
                    resize_histories += 1
                    resize_phase = None
        assert process.wait() == 0
    finally:
        if process.poll() is None:
            process.kill()
        process.wait()
    assert states > 1000 and merges > 50 and samples > 300
    assert rank_nodes > 10000 and queries > 50000 and useful_queries > 1000
    assert tighter_passes > 0 and directional_useful > useful_queries and tighter_queries > 0
    assert resize_phase is None and resize_histories == 2
    assert pure_projections > 100 and exact_refinements > 50 and union_checks == merges
    assert expected_transfer_nodes is None and transfer_queries == queries and transfer_useful > 1000
    if binary64:
        assert dips > 0, "public decreasing-sample regression was not exercised"
    print(f"Operation-aware conservation: {states} states, {merges} merges, {samples} exact samples, "
          f"{dips} decreasing samples; maximum theta {float(max_theta):.6g}")
    print(f"Compositional rank ledger: {rank_nodes} retained cuts, {queries} raw/certified query pairs; "
          f"{useful_queries} query allowances below the observation count")
    print(f"Split sampling: {tighter_passes} passes charge less than the largest source cell")
    print(f"Directional rank ledger: {directional_useful} query allowances below the observation count; "
          f"{tighter_queries} strictly tighter than the symmetric cell bound")
    print(f"Unchanged-count resize: {resize_histories} histories with no releases; "
          f"{amplified_histories} exceed their source uniform error allowance")
    print(f"Integer envelopes: {pure_projections} pure projections, {exact_refinements} exact refinements, "
          f"{union_checks} exact merge unions; {interval_queries} exhaustive partition queries")
    print(f"Exact transfer proof: {rank_nodes} retained cuts, {transfer_queries} query pairs; "
          f"{transfer_useful} allowances below count; no rounding-horizon premise")


def check_initialization_cases(binary):
    """Isolated compiled states test local numeric premises, not reachability."""
    cases = offsets = spans = 0
    process = subprocess.Popen([str(binary), "--initialization-cases"], stdout=subprocess.PIPE, text=True)
    try:
        for line in process.stdout:
            row = json.loads(line)
            if row["kind"] == "format":
                unit = F(1, 1 << row["precision"])
                beta8 = (1 + unit) ** 8 - 1
            elif row["kind"] == "initial_offset":
                k, total, index, offset = (row[key] for key in ("capacity", "total", "index", "offset"))
                target = F(index * (total - 1), k - 1)
                assert 0 <= offset < total
                assert abs(offset - target) <= F(1, 2) + beta8 * total
                if index == 0:
                    assert offset == 0
                if index == k - 1:
                    assert offset == total - 1
                offsets += 1
            else:
                assert row["kind"] == "initialization"
                k, count = row["capacity"], row["n"]
                weights = [(number(x), weight) for x, weight in row["items"]]
                assert max(weight for _, weight in weights) <= count // k
                assert sum(weight for _, weight in weights) <= count
                theta = 10 * (k + 3) * (k + 1) * unit
                assert theta <= F(1, 2)
                budget = count * theta / (1 - theta)
                roots = [(number(x), dyadic(h), low, high, atom)
                         for x, h, low, high, atom in row["nodes"]]
                assert roots[0][0] == weights[0][0] and roots[-1][0] == weights[-1][0]
                for x, h, low, high, atom in roots:
                    exact = sum(weight for value, weight in weights if value <= x)
                    assert low == exact == high
                    assert atom == sum(weight for value, weight in weights if value == x)
                    assert abs(h - exact) <= budget
                for a, b in zip(roots, roots[1:]):
                    assert 0 <= b[1] - a[1] <= F(3 * count, k) + 3 * budget
                    spans += 1
                cases += 1
        assert process.wait() == 0
    finally:
        if process.poll() is None:
            process.kill()
        process.wait()
        process.stdout.close()
    assert cases == 296 and offsets > 0 and spans > 0
    print(f"Compiled initialization: {cases} isolated cases, {offsets} rounded offsets, {spans} cell spans")


class LiveCharges:
    """Directional live charges and exact short-gap telescoping checks.

    Stored-prefix drift is measured explicitly. In the ideal specialization
    it vanishes; these compiled checks do not assume that specialization.
    """

    def __init__(self, capacity):
        self.capacity = capacity
        self.births, self.states = {}, {}
        self.barriers, self.epoch = set(), None
        self.comparisons = self.tail_checks = self.edges = self.tail_depth = 0
        self.final_bound = F(0)

    def support(self, bits):
        self.keys = [number(value) for value in bits]
        assert all(a < b for a, b in zip(self.keys, self.keys[1:]))
        self.spacing = min((b - a for a, b in zip(self.keys, self.keys[1:])), default=None)
        self.released = [0] * len(self.keys)
        self.tree = [0] * (len(self.keys) + 1)
        self.long_limit = 0
        if self.spacing:
            span = self.keys[-1] - self.keys[0]
            while span >= self.spacing and self.long_limit < self.capacity:
                self.long_limit += 1
                span *= F(7, 8)

    def prefix(self, x):
        at, total = bisect.bisect_right(self.keys, x), 0
        while at:
            total += self.tree[at]
            at -= at & -at
        return total

    def interval(self, a, b):
        return self.prefix(b) - self.prefix(a)

    def event(self, row):
        if row["kind"] == "live_birth":
            if self.epoch != row["epoch"]:
                self.barriers.clear()
                self.epoch = row["epoch"]
            a, x, b = (number(row[key]) for key in ("a", "x", "b"))
            assert a < x < b and F(1, 8) <= (x - a) / (b - a) <= F(7, 8)
            self.barriers.update((row["a"], row["x"], row["b"]))
            assert self.spacing is not None
            for side in ("lower", "upper"):
                parent = self.births[row[side]]
                assert row["id"] > row[side]
                assert parent["x"] == row["a" if side == "lower" else "b"]
                long_count = int(b - a >= self.spacing)
                if parent["kind"] == "live_birth" and parent["epoch"] == row["epoch"]:
                    pa, px, pb = (number(parent[key]) for key in ("a", "x", "b"))
                    left, right = (px, pb) if side == "lower" else (pa, px)
                    assert left <= a < x < b <= right
                    assert b - a <= F(7, 8) * (pb - pa)
                    long_count += parent["long_" + side]
                    self.edges += 1
                row["long_" + side] = long_count
                assert long_count <= self.long_limit
        self.births[row["id"]] = row

    def pass_(self, row, source, out, width):
        if self.epoch != row["epoch"]:
            self.barriers.clear()
            self.epoch = row["epoch"]
        assert self.barriers <= set(row["protected"])
        assert len(source) == len(row["source_live_ids"])
        assert len(out) == len(row["out_live_ids"]) == len(row["evaluations"])
        for bits, weight in row["items"]:
            key = number(bits)
            at = bisect.bisect_left(self.keys, key)
            assert at < len(self.keys) and self.keys[at] == key
            self.released[at] += weight
            assert self.released[at] <= row["n"] // self.capacity
            at += 1
            while at < len(self.tree):
                self.tree[at] += weight
                at += at & -at
        old_nodes = dict(zip(row["source_live_ids"], source))
        new_nodes = dict(zip(row["out_live_ids"], out))
        evaluations = dict(zip(row["out_live_ids"], row["evaluations"]))
        states = {}
        for identity in sorted(new_nodes):
            x, h, low, high, _, _ = new_nodes[identity]
            birth = self.births[identity]
            if identity in old_nodes:
                old = old_nodes[identity]
                added = evaluations[identity][1]
                assert low == old[2] + added and high == old[3] + added
                delta = h - old[1] - added
                previous = self.states[identity]
                tails = [dict(tail) if tail else None for tail in previous["tail"]]
                state = {"credit": [previous["credit"][0] + delta, previous["credit"][1] - delta],
                         "drift": previous["drift"] + delta, "tail": tails}
                for side, sign in enumerate((1, -1)):
                    if tails[side]:
                        tails[side]["value"] += sign * delta
            elif birth["kind"] == "live_root":
                state = {"credit": [max(F(0), h - low), max(F(0), high - h)],
                         "drift": F(0), "tail": [None, None]}
            else:
                a_id, b_id = birth["lower"], birth["upper"]
                a, b = new_nodes[a_id], new_nodes[b_id]
                assert a[0] < x < b[0]
                assert a[2] <= low and high <= b[3] - b[4]
                assert a[1] <= h <= b[1]
                state = {"credit": [], "drift": F(0), "tail": []}
                for side, parent in enumerate((a_id, b_id)):
                    charge = h - a[1] if side == 0 else b[1] - h
                    state["credit"].append(states[parent]["credit"][side] + charge)
                    tail = states[parent]["tail"][side]
                    short = b[0] - a[0] < self.spacing
                    if tail and tail["epoch"] == row["epoch"]:
                        assert short
                        tail = dict(tail)
                        tail["value"] += charge
                        tail["length"] += 1
                    elif short:
                        assert bisect.bisect_right(self.keys, b[0]) - bisect.bisect_right(self.keys, a[0]) <= 1
                        tail = {"epoch": row["epoch"], "a": a_id, "b": b_id,
                                "mass": b[1] - a[1], "released": self.interval(a[0], b[0]),
                                "drift_a": states[a_id]["drift"], "drift_b": states[b_id]["drift"],
                                "value": charge, "length": 1}
                    else:
                        tail = None
                    state["tail"].append(tail)
            assert h - low <= state["credit"][0] and high - h <= state["credit"][1]
            states[identity] = state
            self.comparisons += 1
        for state in states.values():
            for side, tail in enumerate(state["tail"]):
                if not tail or tail["epoch"] != row["epoch"]:
                    continue
                a_id, b_id = tail["a"], tail["b"]
                a, b = new_nodes[a_id], new_nodes[b_id]
                arrivals = self.interval(a[0], b[0]) - tail["released"]
                assert 0 <= arrivals <= row["n"] // self.capacity
                drift = (states[b_id]["drift"] - tail["drift_b"] if side == 0
                         else tail["drift_a"] - states[a_id]["drift"])
                assert tail["value"] <= tail["mass"] + arrivals + drift
                self.tail_checks += 1
                self.tail_depth = max(self.tail_depth, tail["length"])
        bound = max((sum(states[i]["credit"]) for i in new_nodes), default=F(0))
        for a_id, b_id in zip(row["out_live_ids"], row["out_live_ids"][1:]):
            a, b = new_nodes[a_id], new_nodes[b_id]
            if math.nextafter(float(a[0]), math.inf) < float(b[0]):
                bound = max(bound, states[a_id]["credit"][0] + states[b_id]["credit"][1]
                            + b[1] - a[1] - b[4])
        assert width <= bound
        self.final_bound = bound
        self.states = states


def check_gap_telescope_model():
    """Stress the local lemma with exact abstract transitions, not a claim
    that these synthetic initial prefixes are public-stream reachable.
    Both sides, incoming mass at one key, and signed prefix drift are covered.
    """
    for side in ("lower", "upper"):
        for rounding in (False, True):
            checker = LiveCharges(32)
            checker.support([encoding(x) for x in (0, F(1, 2), 1)])
            a, b = F(3, 8), F(5, 8)
            rows = [(a, F(100), 0, 1024, 0, 0), (b, F(200), 0, 1024, 0, 1)]
            for identity, x in enumerate((a, b)):
                checker.event({"kind": "live_root", "id": identity, "x": encoding(x)})
            checker.pass_({"epoch": 2000, "n": 1024, "protected": [], "items": [],
                           "source_live_ids": [], "out_live_ids": [0, 1],
                           "evaluations": [(None, 0, None)] * 2}, [], rows, 1024)
            left_id, right_id = 0, 1
            for step in range(12):
                identity = step + 2
                old = {node[5]: node for node in rows}
                left, right = old[left_id], old[right_id]
                x = (left[0] + right[0]) / 2
                checker.event({"kind": "live_birth", "id": identity, "x": encoding(x),
                               "a": encoding(left[0]), "b": encoding(right[0]),
                               "lower": left_id, "upper": right_id, "epoch": 2000})
                source = rows
                rows, evaluations = [], {}
                for node in source:
                    point, h, low, high, atom, old_id = node
                    added = int(point >= F(1, 2))
                    drift = F((-1) ** (step + old_id), 1 << 16) if rounding else F(0)
                    rows.append((point, h + added + drift, low + added, high + added,
                                 atom + int(point == F(1, 2)), old_id))
                    evaluations[old_id] = (None, added, None)
                added = int(x >= F(1, 2))
                h = (left[1] + right[1]) / 2 + added
                rows.append((x, h, left[2] + added, right[3] - right[4] + added,
                             int(x == F(1, 2)), identity))
                evaluations[identity] = (None, added, None)
                rows.sort()
                width = max(node[3] - node[2] for node in rows)
                width = max(width, max(v[3] - v[4] - u[2] for u, v in zip(rows, rows[1:])))
                checker.pass_({"epoch": 2000, "n": 1025 + step,
                               "protected": [encoding(node[0]) for node in rows],
                               "items": [(encoding(F(1, 2)), 1)],
                               "source_live_ids": [node[5] for node in source],
                               "out_live_ids": [node[5] for node in rows],
                               "evaluations": [evaluations[node[5]] for node in rows]}, source, rows, width)
                if side == "lower":
                    left_id = identity
                else:
                    right_id = identity
            assert checker.tail_depth == 12 and checker.tail_checks > 0
    print("Abstract short-gap telescopes: both sides, 12 nested births, exact and signed-drift modes")


def check(binary, require_witness, input_path=None, *, capacity=512, count=15360,
          seed=8, theoretical=False, pattern="contracting", consolidate_each=False):
    live_charges = LiveCharges(capacity)
    arithmetic_comparisons = 0
    arithmetic_budget_sum = F(0)
    largest_survivor_defect = F(0)
    largest_sample_displacement = F(0)
    largest_factor = F(1024 if theoretical else 3)
    split_mass_checks = join_mass_checks = completed_mass_checks = 0
    queries = []
    newborn_comparisons = 0
    total_splits = total_joins = 0
    births, credits = {}, {}
    passes = comparisons = failures = protected_edges = 0
    active_barriers = set()
    active_epoch = None
    final_width = 0
    final_bound = F(0)
    largest_ratio = F(0)
    witness = None
    deepest = 0
    command = [str(binary), "--capacity", str(capacity), "--count", str(count), "--seed", str(seed),
               "--pattern", pattern]
    if theoretical:
        command.append("--theoretical")
    if consolidate_each:
        command.append("--consolidate-each")
    expected_count = count
    if input_path:
        command.extend(["--input", str(input_path)])
        expected_count = len(input_path.read_text().splitlines())
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
    try:
        for line in process.stdout:
            row = json.loads(line)
            if row["kind"] == "format":
                precision = row["precision"]
                unit = F(1, 1 << precision)
                epochs, pass_bound, operation_bound, theta = arithmetic_horizon(expected_count, capacity, precision)
                assert theta <= F(1, 2), "structural arithmetic horizon exceeded"
                universal_budget = expected_count * theta / (1 - theta)
            elif row["kind"] == "support":
                live_charges.support(row["values"])
            elif row["kind"] in ("mass_split", "mass_join"):
                factor = dyadic(row["factor"])
                largest_factor = max(largest_factor, factor)
                base = F(row["n"], capacity) * largest_factor
                if row["kind"] == "mass_split":
                    span = dyadic(row["high"]) - dyadic(row["low"])
                    assert span >= 0
                    assert abs(span - dyadic(row["mass"])) <= universal_budget
                    assert span <= (F(3 * row["n"], capacity) + 3 * universal_budget
                                    if row["initial"] else F(3, 2) * base + 7 * universal_budget)
                    split_mass_checks += 1
                else:
                    combined, limit = dyadic(row["sum"]), dyadic(row["limit"])
                    assert combined <= limit <= F(3, 4) * base + universal_budget
                    join_mass_checks += 1
            elif row["kind"] in ("live_root", "live_birth"):
                live_charges.event(row)
            elif row["kind"] == "query_checks":
                truth = [(number(bits), weight) for bits, weight in row["truth"]]
                assert sum(weight for _, weight in truth) == expected_count
                coordinates, cumulative = [], [0]
                for x, weight in truth:
                    coordinates.append(x)
                    cumulative.append(cumulative[-1] + weight)
                for x_bits, estimate_bits in row["queries"]:
                    x = struct.unpack("!d", struct.pack("!Q", x_bits))[0]
                    exact = cumulative[bisect.bisect_right(coordinates, x)]
                    queries.append(abs(number(estimate_bits) - exact))
            elif row["kind"] == "root":
                births[row["id"]] = row
            elif row["kind"] == "birth":
                a, x, b = (number(row[key]) for key in ("live_a", "x", "live_b"))
                assert a < x < b
                assert F(1, 8) <= (x - a) / (b - a) <= F(7, 8)
                if active_epoch != row["epoch"]:
                    active_barriers.clear()
                    active_epoch = row["epoch"]
                active_barriers.update((row["live_a"], row["x"], row["live_b"]))
                for side in ("lower", "upper"):
                    parent = births[row[side]]
                    assert parent["x"] == row["a" if side == "lower" else "b"]
                    depth = 1
                    if parent["kind"] == "birth" and parent["epoch"] == row["epoch"]:
                        # Frozen source cells need not shrink. The LIVE split
                        # interval is trapped by the parent's protected half.
                        pa, px, pb = (number(parent[key]) for key in ("live_a", "x", "live_b"))
                        lower, upper = (px, pb) if side == "lower" else (pa, px)
                        assert lower <= a < x < b <= upper
                        assert b - a <= F(7, 8) * (pb - pa)
                        protected_edges += 1
                        depth += parent["depth_" + side]
                        ratio = ((number(row["b"]) - number(row["a"])) /
                                 (number(parent["b"]) - number(parent["a"])))
                        if ratio > F(7, 8):
                            failures += 1
                            if ratio > largest_ratio:
                                largest_ratio, witness = ratio, {"parent": parent, "child": row, "side": side}
                    row["depth_" + side] = depth
                    deepest = max(deepest, depth)
                    assert depth <= min(capacity, 10897)
                births[row["id"]] = row
            elif row["kind"] == "pass":
                source, out = nodes(row["source"]), nodes(row["out"])
                factor = dyadic(row["factor"])
                if not theoretical:
                    assert len(row["protected"]) <= capacity - 2
                largest_factor = max(largest_factor, factor)
                if not theoretical:
                    assert largest_factor < 8 * capacity
                masses = [dyadic(value) for value in row["masses"]]
                assert len(masses) == len(out)
                if theoretical and capacity <= 512 and theta <= F(1, 4):
                    # In this regime every unprotected pair is affordable and
                    # mandatory splits cannot consume protection capacity.
                    assert (1 + unit) * sum(masses) <= F(10 * row["n"], 7)
                    if len(out) == capacity:
                        assert capacity - len(row["protected"]) >= 2
                for index, (a, b) in enumerate(zip(out, out[1:]), 1):
                    span = b[1] - a[1]
                    assert abs(span - masses[index]) <= universal_budget
                    if math.nextafter(float(a[0]), math.inf) < float(b[0]):
                        assert span <= largest_factor * F(row["n"], capacity) + 3 * universal_budget
                        completed_mass_checks += 1
                exponent = 10 * (capacity + 3) * (1 + row["splits"] + row["joins"])
                total_splits += row["splits"]
                total_joins += row["joins"]
                relative = exponent * unit
                assert relative <= F(1, 2), "operation-count arithmetic horizon exceeded"
                scale = (source[-1][1] if source else F(0)) + sum(weight for _, weight in row["items"])
                assert scale <= 1 << 80
                arithmetic_budget = scale * relative / (1 - relative)
                arithmetic_budget_sum += arithmetic_budget
                assert not out or out[-1][1] <= scale + arithmetic_budget
                old_live = dict(zip(row["source_live_ids"], source))
                for identity, node, (value, added, _) in zip(row["out_live_ids"], out, row["evaluations"]):
                    if identity not in old_live:
                        birth = live_charges.births[identity]
                        if birth["kind"] == "live_birth":
                            live_nodes = dict(zip(row["out_live_ids"], out))
                            parent_span = live_nodes[birth["upper"]][1] - live_nodes[birth["lower"]][1]
                            assert parent_span <= (F(3, 2) * largest_factor * F(row["n"], capacity)
                                                   + 9 * universal_budget)
                            chosen = dyadic(birth["chosen"])
                            displacement = abs(chosen - dyadic(value) - added)
                            largest_sample_displacement = max(largest_sample_displacement, displacement)
                            assert abs(node[1] - chosen) <= arithmetic_budget
                            assert abs(node[1] - dyadic(value) - added) <= arithmetic_budget + displacement
                            newborn_comparisons += 1
                        else:
                            assert node[2] == node[3]
                            assert abs(node[1] - node[2]) <= universal_budget
                        continue
                    defect = abs(node[1] - old_live[identity][1] - added)
                    assert defect <= arithmetic_budget
                    largest_survivor_defect = max(largest_survivor_defect, defect)
                    arithmetic_comparisons += 1
                if not source:
                    for node, (_, added, _) in zip(out, row["evaluations"]):
                        assert abs(node[1] - added) <= arithmetic_budget
                if active_epoch != row["epoch"]:
                    active_barriers.clear()
                    active_epoch = row["epoch"]
                protected = set(row["protected"])
                assert active_barriers <= protected
                assert protected <= {bits for bits, *_ in row["out"]}
                assert all(a[0] < b[0] and a[1] <= b[1] for a, b in zip(source, source[1:]))
                assert all(a[0] < b[0] and a[1] <= b[1] for a, b in zip(out, out[1:]))
                coordinates = [node[0] for node in source]
                gamma = max((abs(node[1] - dyadic(value) - added)
                             for node, (value, added, _) in zip(out, row["evaluations"])), default=F(0))
                new_credits = {}
                for node, (interpolation, added, atom) in zip(out, row["evaluations"]):
                    x, h, low, high, known, identity = node
                    assert known == atom
                    at = bisect.bisect_left(coordinates, x)
                    if not source or at == 0 and x < coordinates[0]:
                        l = u = 0; p = q = F(0); span = F(0)
                        assert dyadic(interpolation) == 0
                    elif at == len(source):
                        old = source[-1]; l = u = old[3]
                        p, q = credits[old[5]]; span = F(0)
                        assert dyadic(interpolation) == old[1]
                    elif coordinates[at] == x:
                        old = source[at]; l, u = old[2:4]
                        p, q = credits[old[5]]; span = F(0)
                        assert dyadic(interpolation) == old[1]
                    else:
                        a, b = source[at - 1], source[at]
                        l, u = a[2], b[3] - b[4]
                        p, q = credits[a[5]][0], credits[b[5]][1]
                        span = b[1] - a[1]
                        assert a[1] <= dyadic(interpolation) <= b[1]
                        born = births[identity]
                        assert born["kind"] == "birth"
                        assert born["lower"] == a[5] and born["upper"] == b[5]
                        assert number(born["a"]) == a[0] and number(born["b"]) == b[0]
                    assert low == l + added and high == u + added
                    p, q = p + span + gamma, q + span + gamma
                    assert h - low <= p and high - h <= q
                    new_credits[identity] = (p, q)
                    comparisons += 1
                credits = new_credits
                width = max((node[3] - node[2] for node in out), default=0)
                bound = max((sum(credits[node[5]]) for node in out), default=F(0))
                for a, b in zip(out, out[1:]):
                    if math.nextafter(float(a[0]), math.inf) < float(b[0]):
                        width = max(width, b[3] - b[4] - a[2])
                        bound = max(bound, credits[a[5]][0] + credits[b[5]][1] + b[1] - a[1] - b[4])
                assert width <= bound
                forgotten = out[-1][3] - sum(node[4] for node in out) if out else 0
                assert width <= forgotten
                final_width, final_bound = width, min(bound, forgotten)
                live_charges.pass_(row, source, out, width)
                passes += 1
            elif row["kind"] == "summary":
                assert row["passes"] == passes
                assert row["count"] == expected_count
                assert row["violations"] == failures
                assert row["width"] == final_width
        assert process.wait() == 0
    finally:
        if process.poll() is None:
            process.kill()
        process.wait()
        process.stdout.close()
    assert passes > 0 and comparisons > 0
    epochs, pass_bound, operation_bound, theta = arithmetic_horizon(expected_count, capacity, precision)
    assert passes <= pass_bound
    assert total_splits <= capacity * (epochs + 1)
    assert total_joins <= total_splits + 2 * passes + capacity
    operations = passes + total_splits + total_joins
    assert operations <= operation_bound
    assert theta <= F(1, 2), "structural arithmetic horizon exceeded"
    universal_budget = expected_count * theta / (1 - theta)
    actual_theta = 10 * (capacity + 3) * operations * unit
    assert arithmetic_budget_sum <= expected_count * actual_theta / (1 - actual_theta)
    depth = live_charges.long_limit
    side_rate = F(25, 4) * (F(3, 2) * largest_factor * (depth + 1) + 1) * F(expected_count, capacity)
    side_numeric = (9 * depth * epochs + 10 * epochs + 2) * universal_budget
    for state in live_charges.states.values():
        assert max(state["credit"]) <= side_rate + side_numeric
    width_rate = largest_factor * F(expected_count, capacity) + 2 * side_rate
    width_numeric = (18 * depth * epochs + 20 * epochs + 7) * universal_budget
    assert live_charges.final_bound <= width_rate + width_numeric
    raw_rate = largest_factor * F(expected_count, capacity) + side_rate
    raw_numeric = (9 * depth * epochs + 10 * epochs + 7) * universal_budget
    output_allowance = F(3 * expected_count, 1 << 53)
    assert queries and max(queries) <= raw_rate + raw_numeric + output_allowance
    if require_witness:
        assert failures > 0, "saved stream no longer disproves 7/8 frozen-ancestry contraction"
    print(f"Capacity {capacity}, {'theoretical' if theoretical else 'practical'} policy, "
          f"{expected_count} observations, seed {seed}")
    print(f"Exact local charges: {passes} passes, {comparisons} node inequalities; "
          f"{failures} shrink-lemma counterexamples; largest ratio {largest_ratio}")
    print(f"Maximum observed one-sided depth within an epoch: {deepest}")
    print(f"Protected live-interval contraction: {protected_edges} exact ancestry edges")
    print(f"Directional live charges: {live_charges.comparisons} nodes, {live_charges.edges} ancestry edges; "
          f"{live_charges.tail_checks} short-gap telescope checks, maximum tail depth {live_charges.tail_depth}")
    print(f"Final directional live bound: {live_charges.final_bound}")
    print(f"Operation-count arithmetic: {arithmetic_comparisons} surviving cuts; "
          f"largest measured drift {float(largest_survivor_defect):.6g}; "
          f"sum of pass budgets {float(arithmetic_budget_sum):.6g}")
    print(f"Newborn reconstruction: {newborn_comparisons} cuts; "
          f"largest separate sampling displacement {float(largest_sample_displacement):.6g}")
    print(f"Structural arithmetic horizon: {operations}/{operation_bound} operations; "
          f"universal conservation budget {float(universal_budget):.6g}")
    print(f"Compiled mass ceilings: {split_mass_checks} splits, {join_mass_checks} joins, "
          f"{completed_mass_checks} completed cells; maximum factor {largest_factor}")
    print(f"Compiled input-aspect width bound: {float(width_rate):.6g} + "
          f"{float(width_numeric):.6g} numeric allowance; E={epochs}, D={depth}")
    print(f"Raw public queries: {len(queries)} exact-truth comparisons; "
          f"largest measured error {float(max(queries)):.6g}")
    print(f"Final width {final_width}; charged bound capped by forgotten mass {final_bound}")
    if witness:
        print(json.dumps(witness, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiler", required=True)
    parser.add_argument("--force-heaps", action="store_true")
    parser.add_argument("--force-scans", action="store_true")
    parser.add_argument("--binary64", action="store_true")
    parser.add_argument("--input", type=Path)
    parser.add_argument("--pattern", choices=("contracting", "grid", "staged-grid"))
    parser.add_argument("--capacity", type=int)
    parser.add_argument("--count", type=int)
    parser.add_argument("--seed", type=int)
    args = parser.parse_args()
    check_gap_telescope_model()
    root = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="splinesketch-local-charge-") as directory:
        temporary = Path(directory)
        instrument(root, temporary / "include", args.force_scans)
        binary = temporary / "check"
        command = [args.compiler, "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror",
                   "-I", str(temporary / "include"), str(root / "tests/paper_local_charge_tests.cpp"),
                   "-o", str(binary)]
        if args.force_heaps:
            command.extend(["-DSPLINESKETCH_PAPER_FORCE_HEAPS", "-DSPLINESKETCH_VERIFY_PAPER_HEAPS"])
        if args.binary64:
            command.append("-mlong-double-64")
        subprocess.run(command, check=True)
        check_operation_cases(binary, args.binary64)
        check_initialization_cases(binary)
        custom = args.input or any(value is not None for value in (args.pattern, args.capacity, args.count, args.seed))
        check(binary, require_witness=not args.binary64 and not custom, input_path=args.input,
              capacity=args.capacity if args.capacity is not None else 512,
              count=args.count if args.count is not None else 15360,
              seed=args.seed if args.seed is not None else 8,
              pattern=args.pattern if args.pattern is not None else "contracting")
        if not custom:
            for capacity, count, seed, theoretical in ((6, 2048, 3, False),
                                                       (32, 4096, 11, False),
                                                       (512, 8192, 1, True)):
                check(binary, require_witness=False, capacity=capacity, count=count,
                      seed=seed, theoretical=theoretical)
            check(binary, require_witness=False, capacity=128, count=16384, seed=17, pattern="grid")
            protected_grid = temporary / "protected-grid.bits"
            protected_grid.write_text("".join(f"{encoding(x)}\n" for x in
                (1, 2, 0, 3, 4, 10000, 0, 0, 0, 5, 6, 8, 10, 16,
                 5, 6, 8, 10, 16, 5, 6, 8, 10, 16,
                 -5, -4, -3, -2, -1, 10001)))
            check(binary, require_witness=False, input_path=protected_grid, capacity=6,
                  consolidate_each=True)


if __name__ == "__main__":
    main()
