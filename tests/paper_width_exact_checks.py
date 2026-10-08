#!/usr/bin/env python3
"""Check compiled rank/error pairs as exact rational values, avoiding LD rounding."""
import math
import struct
import subprocess
import sys
from fractions import Fraction


def rational(bits):
    value = struct.unpack(">d", int(bits).to_bytes(8, "big"))[0]
    assert math.isfinite(value), value
    return Fraction.from_float(value)


def main():
    result = subprocess.run([sys.argv[1], "--exact-queries"], check=True,
                            text=True, stdout=subprocess.PIPE)
    checks = 0
    for line in result.stdout.splitlines():
        if line.startswith("Paper width:"):
            continue
        estimate, truth, lower, upper, point, uniform = line.split(",")
        truth = int(truth)
        assert int(lower) <= truth <= int(upper)
        error = abs(rational(estimate) - truth)
        assert error <= rational(point), (line, error)
        assert error <= rational(uniform), (line, error)
        checks += 1
    assert checks >= 50000, checks
    print(f"Paper width: {checks} exact rational query/error checks passed")


if __name__ == "__main__":
    main()
