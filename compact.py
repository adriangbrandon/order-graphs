"""Compact text encoding of order graphs (for the appendix), with a parser to check round trips.

Encoding of a weakly connected component: a start order pi0, followed by a bracketed DFS of all
the edges reachable from pi0. An item "a[...]" is an edge labeled a from the current node to
a.(current node), followed by the bracketed items of that target; "a*" is an edge labeled a whose
target was already described (e.g., the edge closing the cycle). The target of an edge is always
obtained by moving its label to the front, so no orders other than pi0 need to be written.
"""
import sys


def move_front(pi, a):
    return (a,) + tuple(z for z in pi if z != a)


def read_graph(path):
    E = set()
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        o, a = line.split()
        E.add((tuple(int(c) for c in o), int(a)))
    return E


def encode(E):
    out = {}
    for (u, a) in E:
        out.setdefault(u, []).append(a)
    nodes = {u for u, _ in E} | {move_front(u, a) for u, a in E}
    indeg = {v: 0 for v in nodes}
    for (u, a) in E:
        indeg[move_front(u, a)] += 1
    seen = set()

    def dfs(u):
        seen.add(u)
        parts = []
        for a in sorted(out.get(u, [])):
            v = move_front(u, a)
            if v in seen:
                parts.append(f"{a}*")
            else:
                sub = dfs(v)
                parts.append(f"{a}[{sub}]" if sub else f"{a}")
        return "".join(parts)

    comps = []
    # roots first, then any unvisited node (cycles): pick a node on a cycle
    starts = sorted(v for v in nodes if indeg[v] == 0)
    for r in starts:
        if r not in seen:
            comps.append((r, dfs(r)))
    while len(seen) < len(nodes):
        v = min(nodes - seen)
        # walk back to a cycle node
        pred = {move_front(u, a): u for (u, a) in E}
        trail = set()
        while v not in trail and v in pred and pred[v] not in seen:
            trail.add(v)
            v = pred[v]
        comps.append((v, dfs(v)))
    return comps


def decode(comps):
    E = set()
    for pi0, s in comps:
        stack, cur, i = [], pi0, 0
        last = None
        while i < len(s):
            ch = s[i]
            if ch.isdigit():
                a = int(ch)
                E.add((cur, a))
                last = move_front(cur, a)
                if i + 1 < len(s) and s[i + 1] == "*":
                    i += 1
            elif ch == "[":
                stack.append(cur)
                cur = last
            elif ch == "]":
                cur = stack.pop()
            i += 1
    return E


def latex(comps, title):
    lines = [r"\paragraph{" + title + "}"]
    for pi0, s in comps:
        body = s.replace("]", "]\\allowbreak{}").replace("*", "$^*$\\allowbreak{}")
        lines.append(r"\noindent\texttt{" + "".join(map(str, pi0)) + r"}:~\texttt{" + body + r"}\par")
    return "\n".join(lines)


if __name__ == "__main__":
    for path in sys.argv[1:]:
        E = read_graph(path)
        comps = encode(E)
        assert decode(comps) == E, path
        n = sum(len(s) for _, s in comps)
        print(f"{path}: {len(E)} edges, {len(comps)} components, {n} characters")
        for pi0, s in comps:
            print("  ", "".join(map(str, pi0)), s)
