# Fixed-q8 tail allocation: staged protocol v1

Started 8 October 2026. Preserve precision-v1 as prior development evidence.
All stages use the same translated analytic wall, torus and sphere-gap oracles,
q8 decoder, byte accounting, sign probes and 32/64/96 surface checks as that study.
No topology claim follows from a tail-error improvement alone.

## Stage 1: resolution check (frozen before execution)

Fully subdivide at depths four and five, including training-empty cells. Test
features 0.02, 0.035, 0.06 and 0.10 on all three shapes (including failed shapes
and larger positive controls). Full trees are deterministic and do not depend
on the fitting seed. Use seed 20 and 100,000 independent uniform and 100,000
near-surface evaluation points, with offsets 10000 and 20000. Near noise sigma
is feature/2. This is a resolution diagnostic with unequal bytes, not a budget
comparison. Count object, branch descriptors and q8 leaf payloads. Save fields,
points and finest extracted meshes. Oracle controls determine where extraction
resolution still prevents an interpretable geometry comparison.

Proceed to deeper budget allocation if depth five improves at least one failed
feature under both final extraction resolutions and the existing probe criterion.
The criterion for a recovery threshold remains success at that feature and all
larger tested features, with oracle success at both final resolutions.

## Next stages

Stage 1 passes the progression rule: depth five recovers the 0.02 wall and 0.06
torus, while depth four fails. Both depths still fail smaller sphere gaps.
All subsequent allocations therefore use depth five (37,449 candidate nodes).

Development settings, frozen before execution: three shapes, features 0.03 and
0.08, fitting seed 101, 10,000 fitting points with 50/50 uniform/near mixture,
8/16/32 KiB caps. Compare grids of 9, 17 and 33 equally spaced thresholds on
[0, fitting CVaR99 of the feasible MSE field], always adding the MSE field's p99.
That extra threshold supplies a feasible baseline upper bound. The threshold
bracket is sufficient because F(t)>=t for nonnegative errors and the baseline
CVaR is feasible. Select the smallest grid whose fitting F objective is within
1% in mean and 5% in every case/cap of the 33-point grid. If neither passes,
use 33. Record all searches and total allocation times, including baseline DP.
No evaluation points or geometry checks select the search density.

MSE uses squared error divided by domain scale squared. Sign allocation adds
0.0025 per wrong sign (zero prediction is wrong; exclude reference |d|<=2e-9),
retaining the established baseline weight without retuning. Tail losses use
absolute distance in scene units. The fitting count is a multiple of 100 so
the variational CVaR99 and worst-1%-mean conventions agree exactly.

Choose the finite threshold-search procedure on separate development cases,
then record and freeze it before the equal-byte evaluation. For every threshold
t, the fitting objective is the additive sum of max(abs(error)-t,0). Exact tree
DP minimizes this loss at the byte cap. Compare the resulting value of
t + sum(excess)/(0.01*n) across thresholds. A finite search is exact only over
its enumerated thresholds and fixed candidate prunings, subject to floating point.
It is not a continuous-threshold optimum or a guarantee on evaluation CVaR99.

The final comparison will hold q8, candidate depth, fitting points and byte caps
fixed across MSE, MSE plus sign penalty and CVaR99. Use at least 100,000 independent
evaluation points per distribution. Report MAE, CVaR99, construction time and
minimum recoverable geometry together. Development selection and final settings
will be appended before running the final stage.

## Frozen final comparison

Development completed before final execution: the 9-point grid has mean F
excess 0.5907% and maximum 2.6784% versus 33 points across all 18 case/cap
pairs. It passes the predeclared 1%/5% rule, so **9 grid points plus baseline
p99** are frozen. The 17-point alternative has mean excess 0.1201%, maximum
1.2011%; these are development objective comparisons, not held-out gains.

Final settings: depth five, features 0.02/0.035/0.06/0.10, three shapes, fresh
fitting seeds 30 and 31, 10,000 fitting points, 8/16/32 KiB caps and all three
q8 policies (216 fields). Independent uniform/near evaluation streams contain
100,000 points each, offsets 10000/20000. Thresholds depend only on fitting
errors. Shapes were seen in resolution diagnostics: this is a controlled
comparison with fresh samples, not unseen-shape generalization. The same points
are paired across policies. Tail measurements contain exactly 1,000 errors.

Reuse all geometry criteria from precision-v1, including oracle pass at both
64/96 grids, probe recalls >=95%, both seeds, and success at every larger tested
feature. Keep failures and resolution disagreements. Report per-case paired
changes and equal-case means, without claiming population confidence from two
seeds. Report fitting-label, candidate, loss-cache, allocation and packing time.
Each policy is independently allocated at each cap; CVaR timing includes its
baseline DP and every threshold DP. Shared candidate time includes construction
of the existing three-encoding candidate container; retained fields use q8 only.

The excess-loss formulation follows [Rockafellar and Uryasev, Conditional
value-at-risk for general loss distributions](https://www.sciencedirect.com/science/article/pii/S0378426602002716).
Here exactness applies to each threshold's tree DP and the enumerated threshold
set. No continuous-threshold, held-out-error or topology optimality is claimed.
