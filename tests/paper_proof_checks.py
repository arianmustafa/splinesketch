"""Exact-rational checks of the code mapping; these are not the theorem's proof."""
import collections
from fractions import Fraction as F
import json
import math
import struct
import subprocess
import sys


def floating(bits):
    return struct.unpack("!d", struct.pack("!Q", bits))[0]


def number(bits):
    return F(floating(bits))


def dyadic(text):
    whole, digits = text.split(":")
    return F(int(whole)) + F(int(digits or "0", 2), 1 << len(digits))


def geometry(binary):
    rows = subprocess.check_output([binary, "--geometry"], text=True).splitlines()
    accepted = 0
    one_fifth = False
    for row in rows:
        av, bv, mv = map(int, row.split())
        a, b, mid = map(number, (av, bv, mv))
        assert a <= mid <= b, (a, b, mid)
        if a < mid < b:
            t = (mid - a) / (b - a)
            assert F(1, 8) <= t <= F(7, 8), t
            one_fifth |= t == F(1, 5)
            accepted += 1
        else:
            # No interior binary64 query exists exactly when the endpoints
            # are adjacent (canonical signed zeros denote the same value).
            assert math.nextafter(float(a), math.inf) == float(b), (a, b, mid)
    assert one_fifth and accepted > 100000
    return len(rows), accepted


def traces(binary):
    rows = [json.loads(line) for line in
            subprocess.check_output([binary, "--trace"], text=True).splitlines()]
    checked = 0
    ledgers = []
    recorded = []
    certificate_bounds = []
    for row in rows:
        kind = row["kind"]
        if kind == "format":
            precision = row["precision"]
            truth = collections.Counter()
            old = {}
            costs = {}
            envelopes = {}
            ledger = F(0)
            continue
        if kind == "query":
            x = floating(row["x"])
            estimate = dyadic(row["rank"])
            exact = sum(weight for value, weight in truth.items() if float(value) <= x)
            spans = [p2 - p1 for (_, p1), (_, p2) in zip(old.items(), list(old.items())[1:])]
            max_span = max(spans, default=F(0))
            n = sum(truth.values())
            # Here heavy and buffer are empty; only clamp/count conversion
            # and the binary64 return can add query arithmetic error.
            u, ud = F(1, 1 << precision), F(1, 1 << 53)
            query_error = 2 * u * n + ud * (1 + u) * n
            bound = ledger + max(costs.values(), default=F(0)) + max_span + query_error
            assert abs(estimate - exact) <= bound, (estimate, exact, bound)
            keys = list(old)
            if x < float(keys[0]):
                low = high = 0
            elif x >= float(keys[-1]):
                low = high = n
            elif F(x) in old:
                low, high = envelopes[F(x)]
            else:
                left = max(value for value in old if float(value) < x)
                right = min(value for value in old if float(value) > x)
                low, high = envelopes[left][0], envelopes[right][1]
            assert low <= exact <= high
            point_bound = max(abs(estimate - low), abs(estimate - high))
            assert abs(estimate - exact) <= point_bound
            assert abs(estimate - exact) <= certificate_bounds[-1] + query_error
            checked += 1
            continue
        nodes = {number(x): dyadic(prefix) for x, mass, prefix, protected in row["nodes"]}
        assert list(nodes) == sorted(nodes)
        assert all(a <= b for a, b in zip(nodes.values(), list(nodes.values())[1:]))
        incoming = {number(x): weight for x, weight in row["incoming"]}
        truth.update(incoming)
        exact_rank = lambda x: sum(weight for value, weight in truth.items() if value <= x)
        if kind == "init":
            ledger = max((abs(prefix - exact_rank(x)) for x, prefix in nodes.items()), default=F(0))
            costs = dict.fromkeys(nodes, F(0))
            envelopes = {x: (exact_rank(x), exact_rank(x)) for x in nodes}
        elif kind == "update":
            defect = F(0)
            for x, prefix in nodes.items():
                if x in old:
                    added = sum(weight for value, weight in incoming.items() if value <= x)
                    target = old[x] + added
                    low, high = envelopes[x]
                    envelopes[x] = (low + added, high + added)
                else:
                    assert x < min(old) or x > max(old)
                    target = exact_rank(x)
                    costs[x] = F(0)
                    envelopes[x] = (exact_rank(x), exact_rank(x))
                defect = max(defect, abs(prefix - target))
            ledger += defect
        elif kind in ("join", "split", "finish"):
            defect = max((abs(prefix - old[x]) for x, prefix in nodes.items() if x in old), default=F(0))
            if kind == "split":
                new = set(nodes) - set(old)
                assert len(new) == 1
                x = new.pop()
                left = max(value for value in old if value < x)
                right = min(value for value in old if value > x)
                low, high = old[left], old[right]
                defect = max(defect, low - nodes[x], nodes[x] - high)
                # Track the strongest of the two endpoint ancestry costs.
                costs[x] = max(costs[left], costs[right]) + high - low
                envelopes[x] = (envelopes[left][0], envelopes[right][1])
                by_x = {number(value): protected for value, mass, prefix, protected in row["nodes"]}
                assert by_x[left] and by_x[x] and by_x[right]
            elif kind == "join":
                assert len(set(old) - set(nodes)) == 1
            ledger += defect
        else:
            raise AssertionError(kind)
        costs = {x: costs[x] for x in nodes}
        envelopes = {x: envelopes[x] for x in nodes}
        for x, prefix in nodes.items():
            assert abs(prefix - exact_rank(x)) <= ledger + costs[x], (kind, x)
            assert envelopes[x][0] <= exact_rank(x) <= envelopes[x][1]
        old = nodes
        checked += len(nodes)
        if kind == "finish":
            ledgers.append(ledger)
            recorded.append(nodes)
            certificate = max((max(abs(prefix - envelopes[x][0]), abs(prefix - envelopes[x][1]))
                               for x, prefix in nodes.items()), default=F(0))
            for left, right in zip(nodes, list(nodes)[1:]):
                if math.nextafter(float(left), math.inf) < float(right):
                    certificate = max(certificate, abs(nodes[left] - envelopes[right][1]),
                                      abs(nodes[right] - envelopes[left][0]))
            certificate_bounds.append(certificate)
    assert recorded[0] == recorded[1], "scan and heap transition traces differ"
    assert max(ledgers) > 0, "fixture must exercise nonzero rounding defects"
    return checked, [str(value) for value in ledgers], [str(value) for value in certificate_bounds]


if __name__ == "__main__":
    total, accepted = geometry(sys.argv[1])
    checks, ledgers, certificates = traces(sys.argv[1])
    print(f"Exact checks: {total} midpoint pairs, {accepted} interior splits, "
          f"{checks} threshold/query inequalities; numeric ledgers {ledgers}; "
          f"uniform integer-envelope bounds {certificates}")
