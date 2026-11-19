# CliqueMITM — exact clique statistics with meet-in-the-middle

**Research Project, Laboratories 5–8 · Babeș–Bolyai University, Faculty of Mathematics and Computer Science**

This repository accompanies the four-page report *Meet-in-the-Middle Algorithms for Exact Clique Statistics*. The paper presents the model, principal experiments, and conclusions; this README provides the longer algorithmic explanation, implementation details, command-line instructions, validation procedure, and qualifications that could not all fit in the report.

**Scope.** The implemented programs count *all* cliques by cardinality, including the empty clique. This is different from listing maximal cliques or merely finding a single maximum clique. The balanced MITM decomposition is established rather than a newly invented algorithm. Our focus is its reusable implementation and the experimental support-block extension.

## 1. Problem and terminology

Let $`G=(V,E)`$ be an undirected simple graph with $`n=|V|`$. A subset $`C\subseteq V`$ is a clique when every two distinct vertices in $`C`$ are adjacent. Define

```math
 c_j(G)=|\{C\subseteq V : C\text{ is a clique},\ |C|=j\}|,
 \qquad P_G(x)=\sum_{j=0}^{n}c_j(G)x^j.
```

The goal is the **entire clique-size histogram** $`(c_0,c_1,\ldots,c_n)`$, with $`c_0=1`$ for the empty clique. The total number of cliques is $`P_G(1)`$ and the maximum clique *size* is the largest $`j`$ with $`c_j>0`$. These counting programs do not output a witness vertex set for a maximum clique.

A **maximal** clique cannot be extended by another vertex; a **maximum** clique has largest cardinality. A complete graph illustrates why maximal-clique enumeration is not an all-clique counting method: it has one maximal clique but $`2^n`$ total cliques (counting the empty one).

## 2. Repository guide

| Path | Purpose |
| --- | --- |
| `src/balanced_mitm.cpp` | Fixed left/right split, direct polynomial subset DP (`mitm-dp`), clique-indicator SOS/zeta transform (`mitm`), exhaustive counting (`enumerate`), and a **local** pivot-inspired baseline (`pivoter`). Graph limit: 63 vertices; DP side limit: 21 vertices. |
| `src/block_mitm.cpp` | Repeated support-block method, budgets 1–20, four selectors (`high`, `low`, `peel`, `clique3`). Graph limit: 120 vertices. |
| `src/block_mitm_compact.cpp` | Same block method with variable-width polynomial rows, supporting budgets up to 25. Graph limit: 120 vertices. |
| `scripts/test_balanced.py` | Randomized comparison of fixed MITM configurations against exhaustive enumeration. |
| `scripts/test_blocks.py` | Randomized brute-force comparisons for block counting and all four selectors. |
| `scripts/test_compact_blocks.py` | Randomized brute-force comparisons of the compact block implementation and several budgets. |
| `short_runs.py` | Rebuilds executables and runs a short benchmark campaign; also expects the **original** Pivoter executable in the neighboring `Pivoter/` directory. |
| `data/graphs/` | Synthetic and constructed graph inputs used in the reported experiments. |
| `results/` | Recorded histograms, runtimes, selector comparisons, and complete/limited-run statuses. |

The original authors' Pivoter implementation is **not** the `pivoter` mode in `balanced_mitm.cpp`. The paper's comparisons to **Original Pivoter** used the authors' separately supplied source/executable, available from <https://github.com/sjain12/Pivoter>. The complete Lab4 archive also places that source in a neighboring `source_code/Pivoter/` directory, not inside this repository.

## 3. Balanced / fixed-split MITM

### 3.1. Unique decomposition

Partition the vertices into disjoint sets $`L`$ and $`R`$, of sizes $`a`$ and $`b`$ respectively. For each clique $`A\subseteq L`$, compute the right-side compatibility set

```math
Q(A)=R\cap\bigcap_{v\in A}N(v), \qquad Q(\varnothing)=R.
```

Every global clique has **exactly one** representation $`A\cup T`$, where $`A`$ is a clique in $`L`$ and $`T`$ is a clique in the induced graph $`G[Q(A)]`$. Uniqueness follows because $`A=C\cap L`$ and $`T=C\cap R`$ are fixed by the original clique $`C`$.

Therefore,

```math
\boxed{P_G(x)=\sum_{\substack{A\subseteq L\\A\text{ clique}}}x^{|A|}P_{G[Q(A)]}(x).}
```

The benefit of preprocessing the right side is that many different left cliques share the same compatibility mask; the histogram for a mask is computed only once.

### 3.2. Subset DP: count cliques without enumerating every completion

For a right-side mask $`M\subseteq R`$, let $`F[M](x)=P_{G[M]}(x)`$. Choose any vertex $`v\in M`$ and write $`M'=M\setminus\{v\}`$. Every clique in $`M`$ either **excludes** $`v`$ or **includes** $`v`$. In the second case, all its remaining vertices must lie in $`M'\cap N(v)`$. Thus:

```math
 F[\varnothing](x)=1, \qquad
 \boxed{F[M](x)=F[M'](x)+xF[M'\cap N(v)](x).}
```

In coefficient form, if $`F[M](x)=\sum_j f[M,j]x^j`$, then

```math
 f[M,j]=f[M',j]+f[M'\cap N(v),j-1],
```

with invalid negative-degree terms treated as zero. The direct recurrence is implemented by `mitm-dp`. The table has $`2^b`$ masks and up to $`b+1`$ coefficient slots per mask.

The alternative `mitm` mode first assigns the monomial $`x^{|T|}`$ to each mask $`T\subseteq R`$ that is itself a clique (zero to other masks). A subset zeta transform then produces

```math
F[M](x)=\sum_{\substack{T\subseteq M\\T\text{ clique}}}x^{|T|}.
```

Both modes compute the same answers by different preprocessing routes. The choice is useful for comparing a recurrence-based DP with a dense SOS transform.

### 3.3. Enumerating left-side cliques exactly once

The implementation uses a fixed vertex order and a depth-first search. At a node representing clique $`A`$, its candidates are vertices that appear later than the last selected vertex and are adjacent to every vertex already chosen. Adding candidate $`v`$ replaces the candidate set with its intersection with $`N(v)`$, and the compatibility mask with $`Q(A)\cap N(v)`$.

At **every** visited left clique $`A`$, the solver adds the coefficients of $`x^{|A|}F[Q(A)](x)`$ to the answer. This includes intermediate, nonmaximal cliques. A conventional pivoted Bron–Kerbosch search designed only to list maximal cliques cannot simply be substituted here: skipping a nonmaximal left clique could miss global cliques whose right-side completion uses it.

An important difference from the block variant below: **the fixed-split search must not prune just because $`Q(A)`$ is empty**. Even then, $`F[\varnothing](x)=1`$, so $`A`$ itself is a valid clique and its descendants may also contribute their left-only cliques.

### 3.4. Complexity and maximum-clique relationship

Ignoring polynomial factors in $`n`$ and polynomial arithmetic, the preprocessing and search take

```math
T(n)=O^*(2^a+2^b), \qquad M(n)=O^*(2^b).
```

For $`a\approx b\approx n/2`$, this gives **$`O^*(2^{n/2})`$ time with exponential memory**. The fixed program limits $`b\leq21`$ for practicality. If $`n>42`$, this implementation cannot keep both halves near $`n/2`$; using a smaller DP half may make the left enumeration much more expensive. The theoretical balanced bound is not a promise that the current program efficiently handles all graphs up to its 63-vertex representation limit.

For optimization-only use, the related recurrence

```math
 F_{\max}[M]=\max\bigl(F_{\max}[M'],\ 1+F_{\max}[M'\cap N(v)]\bigr)
```

leads to $`\omega(G)=\max_{A\text{ clique in }L}(|A|+F_{\max}[Q(A)])`$. Tracking choices could reconstruct a witness, but **this repository's tested executables implement histogram counting**, not that witness-producing solver. Nor is $`O^*(2^{n/2})`$ the best known asymptotic guarantee for maximum-clique optimization alone.

### 3.5. A small worked example

Take $`L=\{a,b\}`$ and $`R=\{u,v\}`$, with edges $`ab`$, $`uv`$, $`au`$, $`av`$, and $`bu`$. The left cliques are $`\varnothing,\{a\},\{b\},\{a,b\}`$, with compatibility sets $`\{u,v\}`$, $`\{u,v\}`$, $`\{u\}`$, and $`\{u\}`$.

Since $`F[\{u,v\}]=1+2x+x^2`$ and $`F[\{u\}]=1+x`$,

```math
\begin{aligned}
P_G(x)&=(1+2x+x^2)+x(1+2x+x^2)\\
&\quad+x(1+x)+x^2(1+x)\\
&=1+4x+5x^2+2x^3.
\end{aligned}
```

Hence there are 4 singletons, 5 edges, 2 triangles, and a maximum clique size of 3.

## 4. Repeated support-block MITM

### 4.1. Motivation and counting invariant

A right-side DP over $`n/2`$ vertices quickly exceeds available memory. The block approach therefore uses a **small support budget** $`k`$ (for example, 16, 20, or 25), and repeats the computation over a shrinking residual graph.

At iteration $`t`$, let the current residual vertex set be $`V_t`$. Select a support block $`S_t\subseteq V_t`$ (up to $`k`$ vertices) and put $`U_t=V_t\setminus S_t`$. Build the clique-polynomial table for masks of $`S_t`$. For each clique $`A\subseteq U_t`$, define

```math
Q_t(A)=S_t\cap\bigcap_{v\in A}N(v).
```

At this iteration, count **only** cliques that intersect the support block: subtract the empty support completion from each table lookup.

```math
\boxed{P_{G[V_t]}(x)=P_{G[U_t]}(x)
+\sum_{\substack{A\subseteq U_t\\A\text{ clique}}}
 x^{|A|}\bigl(F_{S_t}[Q_t(A)](x)-1\bigr).}
```

After counting these cliques, remove $`S_t`$ and continue on $`U_t`$. The empty clique is initialized once, outside the block loop. Every **nonempty** clique belongs to exactly one iteration: the first removed block it intersects. Consequently, no clique is double-counted.

### 4.2. Why empty-mask pruning is correct here

If $`Q_t(A)=\varnothing`$, no support vertex can extend $`A`$. Any descendant $`A'\supseteq A`$ also satisfies $`Q_t(A')\subseteq Q_t(A)=\varnothing`$. Its contribution $`F_{S_t}[Q_t(A')]-1`$ is zero. The entire subtree may therefore be pruned **for this iteration**, because left-only cliques will be processed after removing $`S_t`$.

This exact pruning rule explains much of the variation with graph structure and support selection: the algorithm is fast when compatibility masks become empty early, but can visit an enormous number of states when they remain nonempty.

### 4.3. Support selection policies

The support block is **reselected on every residual graph**, not just once at the beginning.

| `--order` | Implemented rule | Main trade-off |
| --- | --- | --- |
| `high` | Take up to $`k`$ vertices of highest **residual** degree. | Cheap, but can leave many compatible left cliques alive. |
| `low` | Take up to $`k`$ vertices of lowest residual degree. | Cheap alternative; no general guarantee of improvement. |
| `peel` | Greedily construct the explicit side by repeatedly taking a vertex of low degree within the still-unselected pool; the remaining vertices become support. | Degree-based approximation to a favorable partition. |
| `clique3` | Greedily construct the explicit side $`U_t`$ to minimize how many *new small cliques* each candidate would create there; the remaining vertices become support. | More expensive selection, sometimes substantially fewer search states. |

For `clique3`, suppose $`U`$ is the explicit side selected so far. Adding $`v`$ creates a new $`j`$-clique for every $`(j-1)`$-clique in $`G[N(v)\cap U]`$. The score used by the block implementation is the **exact sum of the counts of cliques of sizes 0, 1, 2, and 3** in that neighborhood. In other words, it measures the increment in 1-, 2-, 3-, and 4-vertex cliques, with residual degree and vertex ID as deterministic tie-breakers. This is a local heuristic, not a globally optimal partitioning algorithm.

These policies are distinct. In particular, `clique3` does *not* simply select the highest-degree support vertices. A related `clique4` ordering exists in the fixed-split program, but not as an option in either block executable.

### 4.4. Time, space, and compact rows

For a block of size $`k`$, the standard implementation stores $`(k+1)2^k`$ 32-bit polynomial coefficients, even though a mask of size $`r`$ cannot contain a clique larger than $`r`$. The compact variant uses an offsets array and gives each mask only $`1+|M|`$ slots. Its total coefficient count is exactly

```math
\sum_{M\subseteq S}(1+|M|)=2^k+k2^{k-1}.
```

At $`k=25`$, this means 452,984,832 32-bit coefficient slots, or approximately **1.69 GiB for the coefficient array alone**, plus roughly 128 MiB of 32-bit offsets and other working data. This is consistent with the observed peak resident memory near 1.81 GiB.

This optimization reduces the memory constant but **not** the exponential dependence on $`k`$. Nor does repeated blocking retain the balanced method's $`O^*(2^{n/2})`$ worst-case time bound. On difficult dense graphs, an iteration can still explore nearly $`2^{|V_t|-k}`$ explicit-side states. `k=25` results additionally differ from `k=20` in storage layout, so a speedup between those executables cannot be attributed exclusively to the larger budget.

## 5. Building and running

Run from the `CliqueMITM/` repository root. The code uses GNU C++ extensions (`unsigned __int128` and built-in bit operations); **GCC / g++** is the tested compiler family. Commands below are for a POSIX shell; on Windows, compile under WSL/MinGW or adapt executable names and paths accordingly.

```sh
mkdir -p build

g++ -std=c++20 -O3 -march=native src/balanced_mitm.cpp -o build/balanced_mitm
g++ -std=c++20 -O3 -march=native src/block_mitm.cpp -o build/block_mitm
g++ -std=c++20 -O3 -march=native src/block_mitm_compact.cpp -o build/block_mitm_compact
```

Representative commands:

```sh
# Fixed split: 14 vertices on the DP side, direct polynomial recurrence.
./build/balanced_mitm data/graphs/MM39_shuffle29.edges mitm-dp --split 14 --order clique3

# Same fixed split, but preprocessing with clique flags + SOS transform.
./build/balanced_mitm data/graphs/MM39_shuffle29.edges mitm --split 14 --order clique3

# Count support-intersecting cliques block by block, with exact early pruning.
./build/block_mitm data/graphs/hamming6-4_constructed.edges \
  --budget 16 --order clique3 --seconds 5 --node-limit 100000000

# Memory-compact 25-support run (substantial RAM; see Section 4.4).
./build/block_mitm_compact data/graphs/rand_85_093_9301.edges \
  --budget 25 --order clique3 --seconds 300 --node-limit 0
```

**Options and output.** `balanced_mitm` accepts modes `enumerate`, `pivoter`, `mitm`, and `mitm-dp`, plus `--split` (number of right/DP vertices) and `--order` (`input`, `high`, `low`, `degen`, `clique3`, `clique4`). The block programs accept `--budget`, `--order` (`high`, `low`, `peel`, `clique3`), `--seconds` and `--node-limit`. A zero time or node limit disables that respective cutoff. The block runs print `status=COMPLETE` or `status=LIMIT`, `nodes`, `blocks`, `dp_states`, and `counts=`. The compact executable also provides `selection_ms`, `dp_ms`, and `dfs_ms`.

The `counts=` field contains coefficients in order from size 0 to size $`n`$. **Never treat the histogram from `status=LIMIT` as a certified result**: it reflects only partial progress. Search-node limits are checked periodically rather than at every recursion step.

### Input format

The provided `.edges` files are **zero-based** edge lists:

```text
4 5
0 1
0 2
0 3
1 2
2 3
```

The first line gives `n m`, followed by `m` undirected edges `(u, v)` with vertices numbered `0..n-1`. The parsers also recognize a DIMACS-like `p edge n m` header, but **do not convert 1-based DIMACS edge indices** to zero-based indices. Convert authentic DIMACS inputs before using them (including stripping any `e` prefixes from edge records), and record the conversion. The included Hamming6-4 file was reconstructed from its mathematical definition, not retrieved as the original DIMACS archive file.

## 6. Experiments and reproducibility

### 6.1. Small-scale correctness checks

Correctness is assessed by comparing **every coefficient**, not merely the maximum clique size or total number of cliques. The provided scripts exercise different graph sizes, densities, partitions, selectors, and support budgets:

```sh
python3 scripts/test_balanced.py
python3 scripts/test_blocks.py
python3 scripts/test_compact_blocks.py
```

These scripts expect executables in `build/` as compiled above. `test_blocks.py` writes generated inputs under `data/small_tests/`; the others use temporary directories. The balanced tests use exhaustive `enumerate` as their reference; the block tests independently brute-force clique subsets. The earlier selector-validation campaign described in the paper completed **448 configurations on 112 random graphs**. The scripts are also useful as *regression procedures*; their existence alone does not constitute proof that every possible graph is handled efficiently.

For structured cases, the mathematical answer gives another check. For a Moon–Moser complete multipartite graph with $`n=3r`$ and $`r`$ parts of size three,

```math
P_G(x)=(1+3x)^r.
```

For example, the 39-vertex instance should produce $`(1+3x)^{13}`$. The constructed Hamming6-4 graph has histogram $`1+64x+704x^2+960x^3+240x^4`$.

### 6.2. Short-run comparisons reported in the paper

These are **three-run median whole-process wall times in milliseconds** from the report, including startup and input handling. A dash means that configuration was not tested:

| Instance | $`n`$ | $`m`$ | Fixed MITM | Block MITM | Original Pivoter |
| --- | ---: | ---: | ---: | ---: | ---: |
| Shuffled Moon–Moser | 39 | 702 | 4.23 | 5.40 | 232.78 |
| Random $`G(40,p\approx0.93)`$ | 40 | 725 | 4.50 | 4.12 | 78.87 |
| Constructed Hamming6-4 | 64 | 704 | — | 4.75 | 78.97 |
| Random $`G(80,p\approx0.70)`$ | 80 | 2231 | — | 13.41 | 162.75 |

The corresponding complete histograms agreed coefficient-by-coefficient. For these easy inputs, process startup and parsing can dominate the original Pivoter executable, so the numbers should **not** be read as pure kernel-level algorithmic speedup factors. The bundled `results/short_runs.csv` contains measurements from a **different short-run campaign** (and hence not necessarily the exact wall times in this table); do not silently combine the two timing series.

### 6.3. Selector sensitivity

An earlier three-run selector campaign used $`k=16`$ and found substantial sensitivity to vertex order and graph structure:

| Input | `high` median | `clique3` median | Interpretation |
| --- | ---: | ---: | --- |
| Shuffled Moon–Moser 45 | 91.95 ms | 11.47 ms | Clique-aware partition reduced surviving search states. |
| Shuffled Moon–Moser 51 | 852.82 ms | 147.13 ms | Similar structural advantage. |
| Random 80-vertex, 70% dense | 10.60 ms | 12.06 ms | Clique-aware selector was slightly slower. |

The `results/support_benchmark.csv` file also contains `peel` runs and additional graph/order variations, including timeout statuses. The effect is not determined by edge density alone; it depends on which cross-block compatibility masks persist during recursion.

### 6.4. Larger stress tests and memory measurements

The long experiments used a **300-second cutoff**:

| Instance and solver | Budget | Reported result |
| --- | ---: | --- |
| Shuffled Moon–Moser 69, block `clique3` | 20 | Complete, 112.10 s |
| Shuffled Moon–Moser 69, block `high` | 20 | Timeout |
| Shuffled Moon–Moser 69, original Pivoter | — | Timeout |
| Shuffled Moon–Moser 75, block `clique3` | 20 | Timeout |
| Random 85, 93.28% dense, block `clique3` | 20 | Timeout |
| Random 85, 93.28% dense, **compact** block `clique3` | 25 | Complete, 75.65 s three-run median |
| Random 85, 93.28% dense, original Pivoter | — | Complete, 106.46 s (one run) |

For the 85-vertex graph with 3,330 edges, the successful compact runs took **77.46, 75.43, and 75.65 seconds**. The program and original Pivoter agreed on the entire histogram stored in `results/clique_histogram_85.csv`: **84,940,562,690,386 cliques including the empty clique**, with a maximum clique size of **33**. The compact block program used about **1.81 GiB peak RSS**, compared with roughly **4.3 MiB** for Pivoter. The cross-solver times differ in repetition count, and the memory gap is a central limitation rather than a minor implementation detail.

The 69-vertex multipartite result also agrees with $`(1+3x)^{23}`$. A timeout is **not** a negative correctness result, nor a complete histogram. Report the timeout as such; avoid extrapolating its partial counts.

### 6.5. Files and measurement discipline

- `results/short_runs.csv` — a separate short timing campaign and checked histograms.
- `results/support_benchmark.csv` — selector timing, node counts, block counts, and statuses.
- `results/long_runs_20_support.csv` — capped large experiments, including partially completed runs.
- `results/benchmark_25_support.csv` — the 85-vertex support-budget comparison, timing repeats, and peak RSS.
- `results/clique_histogram_85.csv` — side-by-side coefficient agreement for the completed 85-vertex instance.

For future comparisons, use identical graph bytes, compiler flags, machine conditions, time limits, and timing scopes. Distinguish wall-clock runtime from internally reported solver time, and record ordering/preprocessing overhead and memory. Benchmark files should retain their source/provenance and any conversion instructions.

## 7. Comparison with other exact clique methods

**Original Pivoter (Jain–Seshadhri, 2020)** is a particularly relevant comparator because it also counts clique-size distributions without enumerating every clique separately. It computes a degeneracy ordering and constructs a pivot-based Succinct Clique Tree. Its unique representation of a clique by mandatory (`hold`) vertices and a subset of optional (`pivot`) vertices allows a leaf with $`h`$ mandatory and $`p`$ optional vertices to contribute

```math
x^h(1+x)^p
```

at once. Its worst-case exponential factor is $`O^*(3^{n/3})`$, with polynomial memory. Our **balanced** MITM bound has a smaller exponential base for aggregate counting but pays for a large subset table. The **repeated-block** extension controls that table size but does not inherit the balanced worst-case bound. Neither method is uniformly faster across all instances.

**Bron–Kerbosch with pivoting** targets *maximal-clique enumeration*, so its outputs and costs cannot be equated with an all-clique histogram without additional reasoning. **MCQ, MCS, and BBMC** are primarily optimization approaches for finding a maximum clique, often using greedy coloring and bitsets, and are most relevant if the goal is only $`\omega(G)`$ rather than all $`c_j`$ values. These are task-matched comparisons, not interchangeable baselines.

## 8. Extensions and limitations

The graph decomposition can support other algebraic summaries when the contribution of the left and right parts can be combined efficiently: for instance, modular clique counts, suitable vertex-weighted counts, and maximum-weight objectives using a corresponding DP. A commutative-semiring formulation covers many such examples, but **not every arbitrary monoid statistic** or non-factorizable property. These extensions are conceptual here; they are **not all implemented or benchmarked** in the included executables.

Open work includes:

1. Implement and test a selector based on small-clique *participation* directly on the support side, rather than the present explicit-side construction.
2. Test local swaps or other **survival-aware** selectors that estimate how many nonempty compatibility states would remain.
3. Compare budgets and storage layouts separately so their performance effects are not confounded.
4. Validate original DIMACS instances such as `brock200_2` and `brock200_4`, preserving source hashes, correct 0-based conversion, and matched counting baselines. **No authentic brock run is claimed in the current paper.**
5. Extend experiments to vertex-weighted and other algebraically composable statistics, including memory and overflow analysis.

**Bottom line.** The fixed-split method offers a clean, exact $`O^*(2^{n/2})`$ counting framework when a balanced DP table fits memory. Repeated support blocks permit some larger instances and offer an exact zero-compatibility pruning rule, but their runtime depends strongly on partition selection and graph structure. The presented results include both competitive examples and severe time/memory limitations; the README and paper should be read together rather than as claims of uniform superiority.

## References

- C. Bron and J. Kerbosch, “Algorithm 457: Finding all cliques of an undirected graph,” *Communications of the ACM*, 1973.
- J. W. Moon and L. Moser, “On cliques in graphs,” *Israel Journal of Mathematics*, 1965.
- E. Tomita, A. Tanaka, and H. Takahashi, “The worst-case time complexity for generating all maximal cliques,” *Theoretical Computer Science*, 2006.
- P. Prosser, “Exact algorithms for maximum clique: A computational study,” *Algorithms*, 2012.
- P. San Segundo, D. Rodríguez-Losada, and A. Jiménez, “An exact bit-parallel algorithm for the maximum clique problem,” *Computers & Operations Research*, 2011.
- E. Tomita et al., “A simple and faster branch-and-bound algorithm for finding a maximum clique,” *IEICE Transactions*, 2013.
- J. M. Robson, “Finding a maximum independent set in time $`O(2^{n/4})`$,” technical report, 2001.
- S. Jain and C. Seshadhri, “The power of pivoting for exact clique counting,” *WSDM*, 2020. <https://doi.org/10.1145/3336191.3371839>
- Second DIMACS Implementation Challenge, “Cliques, Coloring and Satisfiability,” benchmark collection, 1992–1993. <https://archive.dimacs.rutgers.edu/Volumes/Vol26.html>
- The companion four-page report provides the precise bibliography and the distinction between published and informal MITM references.
