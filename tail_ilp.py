"""ILP on the DistinctTail graph (Ruben's reduction), solved with HiGHS via scipy.

Nodes: sequences of distinct symbols of [d] (lengths 1..d). Transition delta(t, c) = longest
distinct suffix of t.c. A node covers the requirement (S,x) of level m if |t| >= m and its last m
symbols form S+{x} with x last (or first, bidirectional). Constraints: cover all requirements;
every selected node has a selected predecessor (t' with delta(t', last(t)) = t).
"""
import itertools
import sys
import time

import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix


def tails(d):
    out = []
    for k in range(1, d + 1):
        out += list(itertools.permutations(range(1, d + 1), k))
    return out


def delta(t, c):
    if c in t:
        i = t.index(c)
        return t[i + 1:] + (c,)
    return t + (c,)


def build(d, uni):
    T = tails(d)
    idx = {t: i for i, t in enumerate(T)}
    reqs = {}
    rows, cols = [], []
    for j, t in enumerate(T):
        for m in range(1, len(t) + 1):
            suf = t[len(t) - m:]
            U = frozenset(suf)
            keys = {(U - {suf[-1]}, suf[-1])}
            if not uni:
                keys.add((U - {suf[0]}, suf[0]))
            for key in keys:
                r = reqs.setdefault(key, len(reqs))
                rows.append(r)
                cols.append(j)
    nreq = len(reqs)
    assert nreq == d * 2 ** (d - 1), nreq
    A_cov = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(nreq, len(T)))
    # in-degree: sum_{t' in pred(t)} y_t' - y_t >= 0
    prow, pcol, pval = [], [], []
    for j, t in enumerate(T):
        prow.append(j); pcol.append(j); pval.append(-1.0)
    for jp, tp in enumerate(T):
        for c in range(1, d + 1):
            t = delta(tp, c)
            j = idx[t]
            if j != jp:                      # a self-transition does not give a predecessor
                prow.append(j); pcol.append(jp); pval.append(1.0)
    A_in = coo_matrix((pval, (prow, pcol)), shape=(len(T), len(T)))
    return T, A_cov, A_in


def set_cuts(d, uni, T):
    """For each m and m-set U: the nodes whose last m symbols form U cover at most 2 (1 if uni) of
    the m requirements with set U, so at least ceil(m/2) (m) of them are selected."""
    rows, cols, rhs = [], [], []
    rid = {}
    for j, t in enumerate(T):
        for m in range(2, len(t) + 1):
            key = (m, frozenset(t[len(t) - m:]))
            r = rid.setdefault(key, len(rid))
            rows.append(r); cols.append(j)
    for (m, U), r in sorted(rid.items(), key=lambda kv: kv[1]):
        rhs.append(m if uni else (m + 1) // 2)
    A = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(rid), len(T)))
    return A, np.array(rhs, dtype=float)


def solve(d, uni, time_limit, fix=True):
    T, A_cov, A_in = build(d, uni)
    n = len(T)
    A_cut, rhs = set_cuts(d, uni, T)
    cons = [LinearConstraint(A_cov.tocsr(), lb=1, ub=np.inf),
            LinearConstraint(A_in.tocsr(), lb=0, ub=np.inf),
            LinearConstraint(A_cut.tocsr(), lb=rhs, ub=np.inf)]
    lb = np.zeros(n)
    if fix:   # symmetry breaking: some node has a full-length tail; relabel it to 12...d
        lb[T.index(tuple(range(1, d + 1)))] = 1
    t0 = time.time()
    res = milp(c=np.ones(n), constraints=cons, integrality=np.ones(n), bounds=Bounds(lb, 1),
               options={"time_limit": time_limit, "disp": False})
    sel = [T[i] for i in range(n) if res.x is not None and res.x[i] > 0.5]
    return res, sel, n, time.time() - t0


if __name__ == "__main__":
    d = int(sys.argv[1])
    uni = "--uni" in sys.argv
    tl = float(sys.argv[sys.argv.index("--time") + 1]) if "--time" in sys.argv else 600
    res, sel, n, dt = solve(d, uni, tl)
    print(f"d={d} {'uni' if uni else 'bi'}: {n} variables, status={res.status} ({res.message}), "
          f"objective={res.fun}, dual bound={getattr(res, 'mip_dual_bound', None)}, "
          f"gap={getattr(res, 'mip_gap', None)}, {dt:.1f}s", flush=True)
    if sel:
        with open(f"tail_ilp_{'uni' if uni else 'bi'}_d{d}.txt", "w") as f:
            for t in sel:
                f.write("".join(map(str, t)) + "\n")
