# Order graphs: code and graphs

Companion code for the paper on the size of complete order graphs (Becker, Gómez-Brandón,
Navarro, Prezza). Everything is Python 3 (standard library), except `tail_ilp.py` (needs
NumPy and SciPy ≥ 1.9 for HiGHS) and the C++17 version of the greedy heuristic.

Model: attributes are `1..d`, nodes are orders (permutations), and an edge labeled `a` goes from
`π` to `a·π` (`a` moved to the front). A requirement `(S,x)` is served by a path with `|S|+1`
edges labeled `S ∪ {x}`, with `x` first or last (only last in the unidirectional model).

## Files

| File | What it does |
|---|---|
| `ordergraph_check.py` | Completeness checker, lower bounds `LB_bi`, `LB_uni`, and the SCD-based construction `G_d`. Run without arguments to reproduce the checks in the paper; `python3 ordergraph_check.py 10` also builds and checks `G_d` up to `d = 10`. |
| `check_graph.py` | Checks graph files (one edge per line, `order label`). Uses the unidirectional model if the file name contains `uni` or `--uni` is given. |
| `greedy_ordergraph.cpp` | Connection-aware greedy set cover on the De Bruijn graph (C++, fast). Writes the selected strings. |
| `decode_verify.py` | Decodes the strings written by the C++ program into an order graph and checks it. |
| `greedy_ordergraph.py` | The same greedy heuristic in pure Python (slower; runs faster under PyPy). |
| `tail_ilp.py` | Exact ILP on distinct tails, with the per-set cuts and symmetry fixing, solved with HiGHS through SciPy. |
| `tail_decode.py` | Lifts an ILP solution (list of tails) to an order graph and checks it. |
| `compact.py` | Encoder/decoder of the compact bracket format used in the appendix. |
| `edges.py` | The short decoder printed in the appendix: `edges(start, s)` lists the edges of a component. |
| `graphs/` | The best order graphs of the paper, `{bi,uni}_d{d}_{edges}.txt`, and `compact.json` with their compact encodings. |

## Usage

```bash
# check the graphs of the paper
python3 check_graph.py graphs/*.txt

# greedy heuristic (C++)
g++ -O3 -march=native -std=c++17 -o greedy_ordergraph greedy_ordergraph.cpp
./greedy_ordergraph 8 --uni --seed 3 -v --out strings.txt
python3 decode_verify.py strings.txt --uni --out graph.txt

# greedy heuristic (Python)
python3 greedy_ordergraph.py 6 --uni --seeds 10

# exact ILP on distinct tails (writes tail_ilp_{bi|uni}_d{d}.txt)
python3 tail_ilp.py 6 --uni --time 600
python3 tail_decode.py tail_ilp_uni_d6.txt 6 --uni --out graph.txt

# compact format
python3 -c 'from edges import edges; print(edges("123456", "3[6[5]2]4"))'
```

The C++ program accepts `2 ≤ d ≤ 9`. Running times on a single core: seconds for `d = 7`,
minutes for `d = 8`, and several hours for `d = 9` (about 5 GB of memory).
For `d ≥ 7`, the optimal and best known graphs in `graphs/` were obtained by solving the
distinct-tail ILP with Gurobi (branch and bound, and its large neighborhood search warm-started
from branch-and-bound solutions); that Gurobi setup is not included here.
