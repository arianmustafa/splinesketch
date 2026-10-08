#!/usr/bin/env python3
"""Compare the hybrid heap consolidation against the preceding scanning code.

Runs native/binary64 and forced-heap exact-output checks, five alternating
timing trials, and an additional unseen-seed timing/accuracy suite.
Requires C++17 and Linux/glibc.
"""
from paper_update_optimization import main

if __name__ == "__main__":
    main(phase="heap")
