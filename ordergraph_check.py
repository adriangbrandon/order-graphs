#!/usr/bin/env python3
"""Completeness checker for order graphs, and the SCD-based construction G_d.

Model (Arroyuelo et al., TODS 2024, Defs. 7.1-7.2):
  * Attributes are 1..d; an order is a permutation, written as a tuple.
  * An edge is (pi, a) with a != pi[0]; it goes from pi to a.pi (a moved to the front).
  * A requirement (S, x), with x not in S, is served if there is a path with |S|+1 edges
    whose labels are exactly S + {x}, with the edge labeled x first or last
    (only last in the unidirectional model).

Checking (Observation 3 and Lemma 4 of the paper):
  * h(v) = length of the longest path with distinct labels ending at v. Such a path has
    labels v[k-1], ..., v[0], so paths of every length <= h(v) end at v. It satisfies
        h(v) = max over edges (u, y) into v of 1 + min(h(u), pos_u(y) - 1),
    capped at d, and is computed as a least fixpoint.
  * Edge e = (v, a), with a at (1-based) position p of v, ends a path with m edges and
    distinct labels iff m <= p and m - 1 <= h(v). That path has labels v[m-2..0], a, and
    serves (v[0..m-2], a) (x last) and, for m >= 2, (v[0..m-3] + {a}, v[m-2]) (x first).

Usage:
    python3 ordergraph_check.py            # runs all checks reported in the paper
    python3 ordergraph_check.py 10         # also builds and checks G_d for d up to 10
"""
import itertools
import sys
from math import comb


def move_front(pi, a):
    """The order a.pi: move attribute a to the front of pi."""
    return (a,) + tuple(z for z in pi if z != a)


def nodes_of(edges):
    nodes = set()
    for (pi, a) in edges:
        nodes.add(pi)
        nodes.add(move_front(pi, a))
    return nodes


def longest_into(edges, d):
    """h(v) for every node, as defined in the module docstring."""
    into = {}
    for (u, y) in edges:
        into.setdefault(move_front(u, y), []).append((u, y))
    h = {v: 0 for v in nodes_of(edges)}
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
    """Requirements (S, x) not served by the order graph; empty list iff it is complete."""
    h = longest_into(edges, d)
    served = set()
    for (v, a) in edges:
        p = v.index(a) + 1                        # 1-based position of a in v
        for m in range(1, min(p, h[v] + 1) + 1):  # paths with m edges ending with (v, a)
            served.add((frozenset(v[:m - 1]), a))                     # x = a, last
            if m >= 2 and not unidirectional:
                served.add((frozenset(v[:m - 2]) | {a}, v[m - 2]))   # x = v(m-1), first
    missing = []
    for k in range(d):
        for S in itertools.combinations(range(1, d + 1), k):
            for x in range(1, d + 1):
                if x not in S and (frozenset(S), x) not in served:
                    missing.append((set(S), x))
    return missing


def ring(order):
    """The ring of the given order: re-sort repeatedly by the last attribute."""
    edges, node = set(), tuple(order)
    for _ in range(len(order)):
        edges.add((node, node[-1]))
        node = move_front(node, node[-1])
    assert node == tuple(order)
    return edges


def symmetric_chain_decomposition(d):
    """Greene-Kleitman SCD of the subsets of {1..d}, by parenthesis matching.
    Element i is ')' if in the set and '(' otherwise; matched pairs are fixed along a chain."""
    chains = {}
    for mask in range(1 << d):
        inset = [(mask >> i) & 1 for i in range(d)]
        stack, matched = [], []
        for i in range(d):
            if not inset[i]:
                stack.append(i)
            elif stack:
                matched += [stack.pop(), i]
        key = tuple(sorted((i, inset[i]) for i in matched))
        chains.setdefault(key, []).append(frozenset(i + 1 for i in range(d) if inset[i]))
    return [sorted(c, key=len) for c in chains.values()]


def construction(d):
    """The order graph G_d of Theorem 9: a history path and hairs for each chain."""
    edges = set()
    gadget_total = 0
    for chain in symmetric_chain_decomposition(d):
        beta = len(chain[0])
        order = sorted(chain[0])
        for s in chain[1:]:
            order += sorted(s - set(order))
        order += sorted(set(range(1, d + 1)) - set(order))
        pi = tuple(order)
        k_minus, k_plus = max(1, beta), min(d - 1, d - beta)
        # history: starts at pi[h+1..d] pi[1..h], labels pi(h), ..., pi(1), ends at pi
        node = pi[k_plus:] + pi[:k_plus]
        for i in range(k_plus - 1, -1, -1):
            assert node[0] != pi[i]
            edges.add((node, pi[i]))
            node = move_front(node, pi[i])
        assert node == pi
        # hairs: edges (pi, pi(p)) for p = k_minus+1..d
        for p in range(k_minus, d):
            edges.add((pi, pi[p]))
        gadget_total += k_plus + (d - k_minus)
    return edges, gadget_total


def parse(pairs):
    return {(tuple(int(c) for c in order), label) for order, label in pairs}


# The 40-edge complete order graph for d = 6 (Figure 1), as (source order, label).
GRAPH_D6 = parse([
    ("412563", 6), ("641253", 5), ("564123", 2), ("256413", 1), ("125643", 4),
    ("412563", 3), ("341256", 6), ("634125", 5), ("563412", 1), ("156342", 2),
    ("125643", 3), ("215634", 6), ("215634", 3), ("321564", 4), ("432156", 5),
    ("321564", 6), ("632154", 1), ("632154", 4), ("163254", 5), ("516324", 3),
    ("516324", 4), ("351624", 2), ("235164", 6), ("623514", 5), ("463215", 2),
    ("463215", 5), ("246315", 1), ("246315", 5), ("524631", 3), ("352461", 6),
    ("546321", 3), ("354621", 1), ("135462", 4), ("413562", 5), ("541362", 2),
    ("135462", 6), ("613542", 2), ("261354", 4), ("426135", 3), ("613542", 4),
])


def lb_bi(d):
    return max(comb(d, m) * ((m + 1) // 2) for m in range(1, d + 1))


def lb_uni(d):
    return d * comb(d - 1, (d - 1) // 2)


def report(name, edges, d, unidirectional=False):
    miss = missing_requirements(edges, d, unidirectional)
    status = "complete" if not miss else f"NOT complete ({len(miss)} missing, e.g. {miss[0]})"
    model = "unidirectional" if unidirectional else "bidirectional"
    print(f"{name}: d={d}, {len(edges)} edges, {model}: {status}")
    return not miss


def main():
    dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 9

    print("Sanity checks against Arroyuelo et al. (TODS 2024):")
    report("one ring 123", ring((1, 2, 3)), 3)
    report("one ring 1234", ring((1, 2, 3, 4)), 4)
    report("rings 1234 + 1324", ring((1, 2, 3, 4)) | ring((1, 3, 2, 4)), 4)

    print("\nFigure 1:")
    assert all(a != pi[0] for (pi, a) in GRAPH_D6)
    report("40-edge graph", GRAPH_D6, 6)
    print(f"  lower bound LB_bi(6) = {lb_bi(6)}")

    print("\nConstruction G_d (Theorem 9):")
    for d in range(3, dmax + 1):
        edges, gadgets = construction(d)
        bound = (d - 1) * comb(d, d // 2) + 2 ** d
        ok = not missing_requirements(edges, d, unidirectional=True)
        print(f"  d={d}: |E|={len(edges)}, gadget edges={gadgets}, Thm. 9 bound={bound}, "
              f"LB_bi={lb_bi(d)}, LB_uni={lb_uni(d)}, unidirectionally complete={ok}")

    print("\nFigure 2 without the four redundant history edges:")
    g4, _ = construction(4)
    redundant = parse([("1324", 4), ("4132", 2), ("2314", 4), ("4231", 1)])
    assert redundant <= g4
    report("G_4 minus redundant edges", g4 - redundant, 4)


if __name__ == "__main__":
    main()
