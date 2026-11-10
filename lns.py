import json
import math
import random
import time
from collections import deque

from gurobipy import GRB

from gurobi_solve import A, NAME, instance, model, save, step

MINUTES = 120
STEP_SECONDS = 300
SIZE = 1.4
SEED = 1

states, pred, cov, by_req = instance()
reqs = list(by_req)
near = {s: pred[s] | {step(s, c) for c in A} for s in states}
start = "".join(A)
rng = random.Random(SEED)


def ball(seeds, radius):
    out = set(seeds)
    for _ in range(radius):
        out |= {w for u in out for w in near[u]}
    return out


def destroy(X):
    size = SIZE * math.sqrt(len(X) / 77)
    kind = rng.randrange(3)
    if kind == 0:
        k = rng.randint(round(20 * size), round(44 * size))
        s0 = rng.choice(sorted(X))
        drop, queue = {s0}, deque([s0])
        while queue and len(drop) < k:
            for w in near[queue.popleft()]:
                if w in X and w not in drop:
                    drop.add(w)
                    queue.append(w)
        free = ball(drop, 2)
    elif kind == 1:
        R = set(rng.sample(reqs, min(len(reqs), rng.randint(round(40 * size),
                                                            round(109 * size)))))
        drop = {s for s in X if cov[s] & R}
        free = ball(drop, 1) | {s for r in R for s in by_req[r]}
        free |= {p for s in list(free) for p in pred[s]}
    else:
        drop = set(rng.sample(sorted(X), min(len(X), rng.randint(round(12 * size),
                                                                 round(29 * size)))))
        free = ball(drop, 2)
    drop.discard(start)
    return drop, free | X


X = set(json.load(open(f"solution_{NAME}.json")))
m, x = model(states, pred, by_req)
xs = [x[s] for s in states]
m.Params.OutputFlag = 0
m.Params.TimeLimit = STEP_SECONDS

t0, it = time.time(), 0
while time.time() - t0 < 60 * MINUTES:
    it += 1
    drop, free = destroy(X)
    m.setAttr("LB", xs, [int(s in X and s not in drop) for s in states])
    m.setAttr("UB", xs, [int(s in free) for s in states])
    m.setAttr("Start", xs, [int(s in X) for s in states])
    x[start].LB = 1
    m.optimize()
    if m.Status == GRB.INTERRUPTED:
        break
    if m.SolCount:
        Y = {s for s, v in zip(states, m.getAttr("X", xs)) if v > 0.5}
        if len(Y) < len(X):
            print(f"step {it}: {len(X)} -> {len(Y)} ({time.time() - t0:.0f} s)", flush=True)
            save(Y, pred)
        if len(Y) <= len(X):
            X = Y
