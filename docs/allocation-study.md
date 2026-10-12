# Fixed-candidate allocation study

## Frozen protocol (8 October 2026, before evaluation)

The versioned `fixtures/allocation_cases.csv` fixes three development shapes and
15 held-out cases. Development uses feature 0.04, zero translation and rotation.
Held-out cases isolate feature 0.027, translation, and rotation, then combine new
features 0.061/0.093 with new translations and rotations. Feature means tube
radius, box half-thickness, or sphere separation. Rotations are Rz Ry Rx in
degrees. No horizon is selected using held-out outcomes: report 1, 2, 4 and 8.

Each case gets its own field fitted on 4,096 samples (half uniform, half near
surface), using seeds 0/1/2. This tests an allocation algorithm on unseen shape
parameters, not transfer of a fitted field between shapes. Each seed has 4,096
independent uniform queries, 4,096 near-surface queries (sigma feature/2), 4,096
surface points, and deterministic transformed feature probes. Sample streams
use offsets 10,000/20,000/30,000. Export all coordinates and oracle labels, and
hash them; these exact exported points are the frozen run's sample definition.
C++ normal-distribution sampling can differ between standard libraries.

Use a full depth-four candidate octree (4,681 nodes), including empty regions.
Store float32 oracle corners, without sharing or caching corner calls. All
methods reuse the same candidate object, routed samples, and node losses. Empty
regions have zero fitting loss but remain eligible, unlike the original online
allocator. This deliberate difference keeps a fixed search space.

Budgets: 4,096 / 8,192 / 16,384 / 32,768 bytes. Count the field object plus every
36-byte query node, including internal nodes. Each split adds eight nodes.
Candidate tree, DP tables, queue, and oracle are temporary construction memory,
excluded from retained storage; report table states and transitions as work
diagnostics, not peak RSS. A bundle's loss reduction is divided by its added
bytes, with deterministic node-index and bundle-size tie breaking.

Report exact DP, greedy spending the cap, greedy's best prefix, and lookahead
2/4/8 with best-prefix retention. All greedy trajectories permit negative-gain
steps. Best-prefix retention can leave budget unused and prevents a harmful
last step from artificially inflating the main comparison. Bundles commit the
root split and up to h-1 descendant splits together; their internal allocation
is exact within the h-split limit. A shorter remaining budget can select smaller
bundles. No guarantee that a larger horizon improves every run.

## Objective and guarantee

For node v, L(v) is the sum of squared errors of its float32 trilinear corner
interpolant on fitting samples routed through v, divided by domain scale squared.
Do not use per-leaf means: sample counts differ. Labels/corners never change
during allocation, so loss is additive over the selected leaves.

Let D(v,k) be minimum loss for exactly k splits below and including v:

```
D(v,0) = L(v)
D(v,k) = min sum_i D(child_i,k_i), with sum_i k_i = k-1, k >= 1
```

Terminal nodes only admit k=0. Eight child convolutions solve the constrained
sum. The root selects the minimum over k <= the byte-derived cap, with fewer
splits on exact ties. Backtracking produces the actual packed query field;
the benchmark recomputes its fitting SSE independently. The combinatorial
guarantee is **optimal fitting loss among prunings of this fixed candidate
tree**, subject to floating-point arithmetic. It does not imply optimal corner
values, a globally optimal octree, test error, continuity or topology.

Exact tables require O(N K) states and at most O(N K^2) convolution work;
subtree-size bounds reduce actual work. Lookahead replaces K by h and uses a
priority queue of feasible bundles. Timings include table construction and
selection, then report packing separately. Both take zero new oracle calls.
Candidate generation costs 8N distance calls, plus separate fitting labels.
The shared full candidate is deliberately charged separately; this experiment
does not establish that precomputing it beats the existing lazy online builder.

The background is [Shukla et al., Rate-Distortion Optimized Tree Structured
Compression Algorithms for Piecewise Polynomial Images](https://www.commsp.ee.ic.ac.uk/~pld/publications/quadtree_coding.pdf),
particularly fixed-tree pruning with additive distortion. This implementation
uses a discrete budget DP; it does not implement that paper's prune-and-join
scheme or inherit its image compression guarantees.

## Measurements and reproduction

Run a Release build and CTest, then:

```powershell
.\build-core\Release\sdf_allocation_bench.exe --out runs\allocation-study-v1
python scripts\allocation_report.py runs\allocation-study-v1 --exe build-core\Release\sdf_allocation_bench.exe
python scripts\verify_allocation.py runs\allocation-study-v1
```

The executable refuses an existing output directory and flushes partial CSVs;
only successful completion writes `COMPLETE`. Raw points live in `samples/`.
The report pairs methods by case, seed and budget, checks optimality inequalities
and byte caps, and writes gaps, summaries, source/binary/data SHA-256 hashes and
a standalone HTML report. Timings retain all three serial repetitions; summaries
use their median. These are descriptive measurements on one machine, without
confidence intervals or claims about independent timing trials.

Primary comparisons: absolute fitting SSE gap, gap/optimum (undefined at zero
optimum), and (greedy gap - lookahead gap)/greedy gap (undefined at zero greedy
gap; negative means regression). Aggregate closure uses summed gaps to avoid
unstable ratios near zero. Also report per-pair closure and regressions. Evaluate
uniform/near RMSE, MAE, p99, sign error, and feature-probe inside/outside recall.
Surface MAE is |predicted SDF| on true surface samples, a one-sided residual;
it is not Chamfer distance, Hausdorff error, or a mesh-topology guarantee.

Held-out cases and queries are never used by the allocators. Any investigation
prompted by these results must be labelled post-evaluation analysis; changes to
horizons, distributions or case selection require a new protocol version.

## Results and interpretation

The frozen experiment found **zero fitting SSE gap** in all 216 paired
case/seed/budget combinations (180 held-out, 36 development). Greedy spending
the cap, greedy best-prefix, and all three lookahead horizons select the same
split sets as exact DP here. Gap closure is **undefined**, not 100%, because
there is no baseline gap. The independently enumerated synthetic negative-gain
fixture does distinguish the algorithms: greedy SSE 15, exact and two-split
lookahead SSE 9 at three splits. Synthetic losses demonstrate algorithmic
behaviour; they are not evidence of a natural-geometry improvement.

Post-evaluation diagnostics explain this result. Each split has a fixed marginal
gain `L(parent) - sum L(children)` and identical byte cost. Selected gains
telescope: tree loss equals root loss minus their sum. In every measured case,
greedy's selected gains are nonnegative and no smaller than any unselected
internal-node gain; it spends the split cap. Thus its feasible, ancestor-closed
selection already attains the unconstrained top-K gain upper bound. This is a
second optimality certificate, checked independently of the DP with a 1e-12
floating-point comparison tolerance. Of 31,590 internal candidates across 54
trees, only 53 have negative gain and only two parent/child edges have increasing
gain. Neither creates a hidden profitable branch at the measured budgets.

The implementation records every selected split and every candidate loss in
the evidence. Exact DP does 68,425 child-convolution transitions at 32 KiB;
one-split greedy does 4,680. This large difference in DP work produces only a
small timing difference on these small trees: allocation is roughly 0.6–0.8 ms,
while candidate generation averages roughly 2.6 ms and 37,448 oracle calls.
Fitting labels cost another 4,096 calls. All allocation policies make zero
additional calls. Queue, vector allocation, best-prefix copies, and traversal
overheads limit the practical savings; no large speedup is established.

At 32 KiB, held-out mean uniform RMSE is 0.01300, near-surface RMSE 0.02065,
surface residual MAE 0.01739 and mean feature-probe inside recall 65.99% for all
methods. These means weight shape/seed cases equally, rather than pooling probe
counts. The methods' held-out predictions are identical in this experiment.
Optimal allocation does not cure the remaining representation error or thin
feature failures. No sign-penalty or mesh results are inferred from this SSE
experiment. Larger trees, richer shapes, and different query distributions may
violate the gain-order certificate and need a separately frozen evaluation.

See [the full report](allocation-v1/report.html), [paired gaps](allocation-v1/paired_gaps.csv),
and [the frozen evidence archive](allocation-v1/evidence.zip). The ZIP preserves
all 270 fitting/evaluation point files, candidate losses, selected splits, raw
timings, metrics and source snapshots with a SHA-256 manifest. Verify without
unpacking using `python scripts/verify_allocation.py docs/allocation-v1/evidence.zip`.
The outer `SHA256SUMS` also covers the ZIP and readable report files. Reproduce
an archive by adding `--archive NEW_DIRECTORY` to the report command.

Validation: the existing core suite and 4,499 allocation checks pass in a local
MSVC Release build, including exhaustive enumeration of all 257 depth-two
prunings under 30 cost assignments, budget boundaries, a negative-gain unlocking
case, packed-field SSE reconstruction, and zero allocation oracle calls. CI
now runs the allocation benchmark and provenance verification; a remote CI run
has not been observed for this change.
