"""Decode a DistinctTail selection into De Bruijn strings, then into an order graph, and verify.

Usage: python3 tail_decode.py tails_file d [--uni] [--out graph_file]
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from greedy_ordergraph import decode, missing_requirements, write_graph  # noqa: E402
from tail_ilp import delta  # noqa: E402

f, d = sys.argv[1], int(sys.argv[2])
uni = "--uni" in sys.argv
out = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else None
sel = [tuple(int(c) for c in l.strip()) for l in open(f) if l.strip()]
S = set(sel)
pred = {}
for t in sel:
    c = t[-1]
    ps = [tp for tp in sel if tp != t and delta(tp, c) == t]
    assert ps, ("no predecessor", t)
    pred[t] = ps[0]
strings = set()
for t in sel:
    labels, v = [], t
    while len(labels) < d:
        labels.append(v[-1])
        v = pred[v]
    s = tuple(reversed(labels))
    # the distinct tail of s must be t
    seen, k = set(), 0
    for c in reversed(s):
        if c in seen:
            break
        seen.add(c); k += 1
    assert s[d - k:] == t, (s, t)
    strings.add(s)
E = decode(strings, d)
miss = missing_requirements(E, d, uni)
print(f"d={d} {'uni' if uni else 'bi'}: {len(sel)} tails -> {len(strings)} strings -> {len(E)} edges, "
      f"complete: {not miss}")
if out and not miss:
    write_graph(out, E, d, uni)
