#!/usr/bin/env python3
"""Decode a set of De Bruijn strings (output of greedy_ordergraph.cpp) into an order graph and verify it.

Usage: python3 decode_verify.py strings_file [--uni] [--out graph_file]
Needs greedy_ordergraph.py in the same directory (for decode and the completeness checker).
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from greedy_ordergraph import decode, missing_requirements, write_graph  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("strings")
ap.add_argument("--uni", action="store_true")
ap.add_argument("--out")
args = ap.parse_args()

sel = set()
for line in open(args.strings):
    line = line.strip()
    if line and not line.startswith("#"):
        sel.add(tuple(int(c) for c in line))
d = len(next(iter(sel)))
E = decode(sel, d)
miss = missing_requirements(E, d, args.uni)
model = "uni" if args.uni else "bi"
print(f"d={d} {model}: {len(sel)} strings -> {len(E)} edges, complete: {not miss}")
if miss:
    print("missing (first 5):", miss[:5])
    sys.exit(1)
if args.out:
    write_graph(args.out, E, d, args.uni)
    print("written", args.out)
