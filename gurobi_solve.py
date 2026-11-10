import itertools
import json
from collections import defaultdict, deque

import gurobipy as gp
from gurobipy import GRB

D = 7
BIDIRECTIONAL = True
TIME_LIMIT = 3600

A = [str(i) for i in range(1, D + 1)]
NAME = f"d{D}_{'bi' if BIDIRECTIONAL else 'uni'}"


def step(s, c):
    j = s.find(c)
    return s + c if j < 0 else s[j + 1:] + c


def covers(s):
    out = set()
    for k in range(2, len(s) + 1):
        t = s[-k:]
        out.add((t[-1], "".join(sorted(t[:-1]))))
        if BIDIRECTIONAL:
            out.add((t[0], "".join(sorted(t[1:]))))
    return out


def instance():
    states = ["".join(p) for k in range(1, D + 1) for p in itertools.permutations(A, k)]
    pred = defaultdict(set)
    for s in states:
        for c in A:
            pred[step(s, c)].add(s)
    cov = {s: covers(s) for s in states}
    by_req = defaultdict(list)
    for s in states:
        for r in cov[s]:
            by_req[r].append(s)
    return states, pred, cov, by_req


def model(states, pred, by_req):
    m = gp.Model()
    x = m.addVars(states, vtype=GRB.BINARY)
    m.setObjective(x.sum(), GRB.MINIMIZE)
    for group in by_req.values():
        m.addConstr(gp.quicksum(x[s] for s in group) >= 1)
    for s in states:
        m.addConstr(gp.quicksum(x[p] for p in pred[s]) >= x[s])
    if BIDIRECTIONAL:
        for k in range(3, D + 1, 2):
            ending = defaultdict(list)
            for s in states:
                if len(s) >= k:
                    ending[frozenset(s[-k:])].append(s)
            for group in ending.values():
                m.addConstr(gp.quicksum(x[s] for s in group) >= (k + 1) // 2)
    x["".join(A)].LB = 1
    return m, x


def save(X, pred):
    p = {s: min(pred[s] & X) for s in X}
    word = {}
    for s0 in sorted(X):
        path, s = [], s0
        while s not in path and s not in word:
            path.append(s)
            s = p[s]
        if s in word:
            continue
        cycle = path[path.index(s):][::-1]
        w = "".join(c[-1] for c in cycle) * (D // len(cycle) + 2)
        for i, c in enumerate(cycle):
            end = len(w) - len(cycle) + i + 1
            word[c] = w[end - D:end]
    children = defaultdict(list)
    for s, q in p.items():
        children[q].append(s)
    queue = deque(word)
    while queue:
        q = queue.popleft()
        for s in children[q]:
            if s not in word:
                word[s] = word[q][1:] + s[-1]
                queue.append(s)
    json.dump(sorted(X), open(f"solution_{NAME}.json", "w"))
    json.dump(sorted(word.values()), open(f"words_{NAME}.json", "w"))


if __name__ == "__main__":
    states, pred, cov, by_req = instance()
    m, x = model(states, pred, by_req)
    m.Params.TimeLimit = TIME_LIMIT
    m.optimize()
    save({s for s in states if x[s].X > 0.5}, pred)
