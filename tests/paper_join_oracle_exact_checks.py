#!/usr/bin/env python3
"""Verify oracle queries and error-free subtraction with exact fractions."""
import math
import struct
import subprocess
import sys
from fractions import Fraction


def rational(bits):
    value = struct.unpack('>d', int(bits).to_bytes(8, 'big'))[0]
    assert math.isfinite(value), value
    return Fraction.from_float(value)


def main():
    result = subprocess.run([sys.argv[1], '--exact-queries', *sys.argv[2:]], check=True,
                            text=True, stdout=subprocess.PIPE)
    checks = 0
    for line in result.stdout.splitlines():
        if line.startswith('Oracle exact:'):
            continue
        estimate, truth, lower, upper, point, uniform, high, low = line.split(',')
        truth = int(truth)
        assert int(lower) <= truth <= int(upper)
        error = abs(rational(estimate) - truth)
        assert error == rational(high) + rational(low), (line, error)
        assert error <= rational(point) and error <= rational(uniform), (line, error)
        checks += 1
    assert checks >= 5, checks
    print(f'Oracle: {checks} exact rational error/subtraction checks passed')


if __name__ == '__main__':
    main()
