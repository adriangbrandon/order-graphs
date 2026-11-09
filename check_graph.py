#!/usr/bin/env python3
"""Check that order graphs given as edge lists are complete.

Usage: python3 check_graph.py [--uni] graph_file ...

A graph file has one edge per line, "order label" (e.g. "123456 3"); lines starting with # are
ignored. If a file name contains "uni", or --uni is given, the unidirectional model is used.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ordergraph_check import missing_requirements  # noqa: E402

ok = True
for path in [a for a in sys.argv[1:] if not a.startswith("--")]:
    E = set()
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        o, a = line.split()
        E.add((tuple(int(c) for c in o), int(a)))
    d = len(next(iter(E))[0])
    uni = "--uni" in sys.argv or "uni" in os.path.basename(path)
    miss = missing_requirements(E, d, uni)
    print(f"{path}: d={d}, {len(E)} edges, {'unidirectional' if uni else 'bidirectional'}: "
          f"{'complete' if not miss else f'NOT complete ({len(miss)} missing)'}")
    ok &= not miss
sys.exit(0 if ok else 1)
