#!/usr/bin/env python3
"""Small complete order graphs via a connection-aware greedy set cover on the De Bruijn graph.

Self-contained (standard library only; runs faster under PyPy).

Model (Arroyuelo et al., TODS 2024, Defs. 7.1-7.2): attributes 1..d, nodes are orders, an edge
(pi, a) goes from pi to a.pi (a moved to the front). A requirement (S, x) is served by a path with
|S|+1 edges labeled S + {x}, with x first or last (only last in the unidirectional model).

Method (see Appendix A of the paper for the De Bruijn view):
  * A node of the De Bruijn graph is a string s in [d]^d, the last d labels read when arriving
    at a node of the order graph. It covers the requirements whose witness is a suffix of s with
    distinct labels: (U - {last}, last) and, in the bidirectional model, (U - {first}, first).
  * Every selected string needs a selected predecessor. We start from the ring 12..d (its d
    rotations) and, to add a target string s, add the shortest De Bruijn path from the selection
    to s: if a selected string ends with s[1..o], only d-o new strings are needed (the
    intermediate ones also cover requirements).
  * Each step picks the target maximizing (newly covered requirements) / (new strings)**alpha,
    breaking ties at random (seed > 0) or by enumeration order (seed 0).
  * A reverse-delete pass removes strings that are no longer needed.
  * The selection is decoded into an order graph (proof of Lemma 13: keep one selected
    predecessor per string; the order of a node is given by the distinct labels of its infinite
    backward walk, most recent first, followed by the unseen attributes in increasing order),
    and its completeness is verified.

Usage:
  python3 greedy_ordergraph.py 6                         # bidirectional, seed 0
  python3 greedy_ordergraph.py 6 --uni --seeds 10         # unidirectional, seeds 0..9
  pypy3 greedy_ordergraph.py 8 --seeds 4 --out graphs/ -v
"""
import argparse
import itertools
import os
import random
import sys
import time
from math import comb


# ---------------------------------------------------------------------------------------------
# Order graphs and completeness (same checker as in Appendix C)

def move_front(pi, a):
    return (a,) + tuple(z for z in pi if z != a)


def longest_into(edges, d):
    into = {}
    nodes = set()
    for (u, y) in edges:
        v = move_front(u, y)
        into.setdefault(v, []).append((u, y))
        nodes |= {u, v}
    h = {v: 0 for v in nodes}
    changed = True
    while changed:
        changed = False
        for v, inc in into.items():
            best = max(min(d, 1 + min(h[u], u.index(y))) for (u, y) in inc)
            if best > h[v]:
                h[v] = best
                changed = True
    return h


def missing_requirements(edges, d, unidirectional=False):
    h = longest_into(edges, d)
    served = set()
    for (v, a) in edges:
        p = v.index(a) + 1
        for m in range(1, min(p, h[v] + 1) + 1):
            served.add((frozenset(v[:m - 1]), a))
            if m >= 2 and not unidirectional:
                served.add((frozenset(v[:m - 2]) | {a}, v[m - 2]))
    missing = []
    for k in range(d):
        for S in itertools.combinations(range(1, d + 1), k):
            for x in range(1, d + 1):
                if x not in S and (frozenset(S), x) not in served:
                    missing.append((set(S), x))
    return missing


# ---------------------------------------------------------------------------------------------
# Coverage of De Bruijn strings

def req_index(d):
    idx = {}
    for k in range(d):
        for S in itertools.combinations(range(1, d + 1), k):
            for x in range(1, d + 1):
                if x not in S:
                    idx[(frozenset(S), x)] = len(idx)
    return idx


def covers(s, d, idx, uni):
    mask = 0
    for m in range(1, d + 1):
        t = s[d - m:]
        if len(set(t)) < m:
            break
        U = frozenset(t)
        mask |= 1 << idx[(U - {t[-1]}, t[-1])]
        if not uni:
            mask |= 1 << idx[(U - {t[0]}, t[0])]
    return mask


# ---------------------------------------------------------------------------------------------
# Greedy with shortest-path connections

def greedy(d, cov, full, rng, alpha, verbose):
    cands = [s for s, c in cov.items() if c]
    base = tuple(range(1, d + 1))
    sel = {base[i:] + base[:i] for i in range(d)}
    covered = 0
    for s in sel:
        covered |= cov[s]
    witness = [dict() for _ in range(d)]          # witness[o][suffix of length o] = selected string

    def add(s):
        for o in range(d):
            witness[o].setdefault(s[d - o:], s)

    for s in sel:
        add(s)
    t0 = time.time()
    while covered != full:
        cands = [s for s in cands if s not in sel and (cov[s] & ~covered)]
        best, choices = None, []
        for s in cands:
            o = d - 1
            while s[:o] not in witness[o]:
                o -= 1
            src = witness[o][s[:o]]
            X = src + s[o:]
            new = [w for w in (X[j - d + 1:j + 1] for j in range(d, len(X))) if w not in sel]
            g = 0
            for w in new:
                g |= cov[w]
            g = bin(g & ~covered).count("1")
            score = g / len(new) ** alpha
            if best is None or score > best + 1e-12:
                best, choices = score, [new]
            elif rng is not None and score > best - 1e-12:
                choices.append(new)
        new = rng.choice(choices) if rng is not None else choices[0]
        for w in new:
            sel.add(w)
            covered |= cov[w]
            add(w)
        if verbose:
            left = bin(full & ~covered).count("1")
            print(f"    {len(sel)} strings, {left} requirements left, "
                  f"{len(cands)} candidates, {time.time() - t0:.0f}s", flush=True)
    return sel


def has_pred(w, sel, d):
    return any((c,) + w[:-1] in sel for c in range(1, d + 1))


def reverse_delete(sel, cov, full, d, rng):
    sel = set(sel)
    order = sorted(sel)
    rng.shuffle(order)
    changed = True
    while changed:
        changed = False
        for w in order:
            if w not in sel:
                continue
            rest = sel - {w}
            c = 0
            for s in rest:
                c |= cov[s]
            if c != full:
                continue
            succ = (w[1:] + (c2,) for c2 in range(1, d + 1))
            if all(has_pred(v, rest, d) for v in succ if v in rest):
                sel = rest
                changed = True
    return sel


def decode(sel, d):
    pred = {}
    for w in sel:
        for c in range(1, d + 1):
            v = (c,) + w[:-1]
            if v in sel:
                pred[w] = v
                break
        assert w in pred, "a selected string has no selected predecessor"

    def order(s):
        walk, v = list(reversed(s)), s
        for _ in range(len(sel) + 1):
            v = pred[v]
            walk.append(v[0])
        seen, out = set(), []
        for c in walk:
            if c not in seen:
                seen.add(c)
                out.append(c)
        return tuple(out + [c for c in range(1, d + 1) if c not in seen])

    ords = {s: order(s) for s in sel}
    E = set()
    for w in sel:
        pv, c = ords[pred[w]], w[-1]
        if pv[0] != c:
            assert move_front(pv, c) == ords[w]
            E.add((pv, c))
    return E


# ---------------------------------------------------------------------------------------------

def lb(d, uni):
    if uni:
        return d * comb(d - 1, (d - 1) // 2)
    return max(comb(d, m) * ((m + 1) // 2) for m in range(1, d + 1))


def write_graph(path, E, d, uni):
    with open(path, "w") as f:
        f.write(f"# complete order graph, d={d}, {'unidirectional' if uni else 'bidirectional'} model, "
                f"{len(E)} edges\n# one edge per line: source order, label (target = label moved to the front)\n")
        for (pi, a) in sorted(E):
            f.write("".join(map(str, pi)) + " " + str(a) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("d", type=int, nargs="+")
    ap.add_argument("--uni", action="store_true", help="unidirectional model")
    ap.add_argument("--seeds", type=int, default=1, help="runs with seeds 0..seeds-1 (default 1)")
    ap.add_argument("--alpha", type=float, default=1.0, help="score = gain / cost**alpha (default 1)")
    ap.add_argument("--out", help="directory for the best graph of each d")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()
    if args.out:
        os.makedirs(args.out, exist_ok=True)
    model = "uni" if args.uni else "bi"
    for d in args.d:
        t0 = time.time()
        idx = req_index(d)
        full = (1 << len(idx)) - 1
        cov = {s: covers(s, d, idx, args.uni) for s in itertools.product(range(1, d + 1), repeat=d)}
        if args.verbose:
            print(f"d={d}: {len(cov)} strings, {len(idx)} requirements, setup {time.time() - t0:.0f}s",
                  flush=True)
        best = None
        for seed in range(args.seeds):
            rng = random.Random(seed) if seed > 0 else None
            sel = greedy(d, cov, full, rng, args.alpha, args.verbose)
            sel2 = reverse_delete(sel, cov, full, d, random.Random(seed))
            E = decode(sel2, d)
            ok = not missing_requirements(E, d, args.uni)
            print(f"  d={d} {model} seed={seed}: greedy {len(sel)} -> pruned {len(sel2)} -> "
                  f"{len(E)} edges, complete: {ok}", flush=True)
            if ok and (best is None or len(E) < len(best)):
                best = E
        print(f"d={d} {model}: best {len(best)} edges (lower bound {lb(d, args.uni)}), "
              f"{time.time() - t0:.0f}s", flush=True)
        if args.out:
            write_graph(os.path.join(args.out, f"greedy_{model}_d{d}_{len(best)}.txt"), best, d, args.uni)


if __name__ == "__main__":
    main()
