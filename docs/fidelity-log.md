# SDF Budget Lab: bounded fidelity research round

Started 8 October 2026. At most three evidence-driven interventions, then a
frozen held-out evaluation and independent reproduction of the strongest choice.

## Evidence and protocol

The selective depth-six baseline uses K=32, q8 and exact MSE allocation.
Its 32,628 bytes recover torus radius .035 but widths remain 28--31% low;
sphere gap .060 fails because no proposed deeper splits survive allocation.
All earlier evaluation cases, including seeds 40/41, are now development.
No universal geometric fidelity or unseen-shape claim is intended.

Keep fitting 10,000 points, equal uniform/near mixture, feature/2 noise,
translation (.017,-.023,.011), 32,768-byte cap including query object, and
original probe recalls >=95%, closed manifold and component Euler values
[0] for torus / [2,2] for spheres at BOTH 64 and 96 extraction cells.
Require both seeds and every larger control: torus .035/.060/.100 and gap
.060/.100. Retain 32-cell diagnostics. Widths use the unchanged feature/2048
scan and 40 bisections; they supplement recovery without redefining it.

Development screening uses seeds 40/41 and 10,000 queries per distribution.
It saves fields, fitting/query/probe data and predictions; mesh checks are
deferred until a candidate is selected. Screen by worst target recall, then
dimensional error, subject to maintaining larger probe controls. Report all
error and construction tradeoffs. Full frozen evaluation uses fresh paired
seeds 140/141 and 100,000 queries per distribution. Those seeds will not guide
an intervention in this round. Independently rebuild the frozen source and
repeat the evaluation; compare fields and predictions exactly, timing separately.

## Intervention 1: make fitting sign failures costly

Hypothesis: the gap has candidate coverage but MSE undervalues its sign failures.
Keep residual-ranked K=32 and corners unchanged. Apply the existing exact
additive q8 loss SSE + lambda * wrong-sign fitting count, lambda in
{0.0001, 0.001, 0.01}, with the selective MSE baseline paired in each run.
No probes enter training or cell ranking. Inspect independent decoded metrics
before deciding whether increased candidate coverage or another intervention
is warranted. The development grid requires a reproducible harness and independent verification.

### Intervention 1 evidence and decision

All four C++ suites passed. The independent decoder verified the first two
penalties including every saved prediction, residual rank, fitting objective,
byte layout and cost counter. The .01 verification follows the same checks.
At .001/.01 the gap recall plateaus at 85.15625% on both development seeds;
the remaining 19 failing probes per seed terminate at depth five. Only 0/2
deeper gap parents are retained. Torus target retains all 32 candidates and
its widths barely change. Larger probe controls pass. Lambda .0001 leaves
one gap seed at 68.75%, so retain .001, the smaller plateau penalty.

## Intervention 2: relax the candidate bottleneck

Hypothesis: torus dimensional distortion persists because all 32 candidates
are exhausted, and broader candidate coverage may also help the gap allocation.
Keep lambda .001, q8, residual ranking and fitting samples fixed; compare
K=128 and K=512 on the same development cases. This tests candidate capacity
separately from changing the loss. Inspect recall, widths, error tails and
retained deeper cells before choosing the third and final intervention.

### Intervention 2 evidence and decision

Both K values pass independent prediction/objective/rank/byte/cost checks.
Torus target recall becomes 100% on both seeds. Median radial width becomes
.057097/.057185 and vertical .054911/.054959, against truth .070; improved
but still undersized. K=512 gives exactly the same target widths and gap
recalls as K=128. Gap recall remains 85.15625%, axis width .049409 vs .060.
Use K=128 for the next experiment. All larger probe controls remain passing.
A saved-output diagnostic finds 7 fitting points and no fitting sign errors
in the remaining failed gap leaf at seed 40 (seed 41 has 10 points, 2 errors).
Simply increasing sign cost cannot penalize an unobserved fitting failure.

## Intervention 3: deterministic allocation coverage

Keep K=128 and lambda .001. After ranking candidates using ONLY the original
10,000 samples, add a Cartesian grid of interior cell-centre oracle labels
inside each selected depth-five parent: q=2 or q=4 samples per axis. These
8 or 64 points per parent augment the allocation objective with unit weight;
they do not alter candidate ranking, corner values, evaluation, or probes.
Recompute the exact additive SSE/sign loss over original plus extra labels.
Count every extra oracle call and its construction time. Save all extra labels
and independently reconstruct their coordinates and truth. This is a changed
allocation-fitting distribution, not free oracle access or held-out data.
The candidate pool and queries remain analytic-shape agnostic. This is the
third and final intervention. After screening both q values, freeze the best
development choice, run full geometry checks, then fresh held-out evaluation.

### Intervention 3 evidence and selection

Both grids reach 100% target probe recall on both seeds. q=4 also reaches
100% on every larger probe control (q=2 leaves one larger gap at 95.3125%).
The independent decoder and coordinate reconstruction pass for both grids.
q=4 torus radial median widths are .057514/.057395; vertical .054955/.056016.
Mean absolute relative width errors over all angles are 18.91/19.21% and
19.26/19.00% respectively, smaller than q=2's roughly 19.8--21.7%.
Both q values give gap widths .066829/.063191 against .060: recovery improves,
but the dimensional error changes sign and remains 5.32--11.38%.
Choose K=128, lambda=.001, q=4. The selection script makes the dimensional
tie-break explicit as equal weight to six target case/axis mean absolute
relative errors; this was decided on development evidence, before fresh seeds.
The choice is also favored by all-controls recall, but incurs more oracle calls.

Full development confirmation saves every extraction mesh at 32/64/96 cells.
The C++ checks pass at required 64/96 for both targets and all larger controls.
At 32 cells the thin torus has five sphere-like components even for the oracle;
this diagnostic failure is retained and does not change the pre-existing rule.
Independent Python mesh connectivity, edge/vertex manifold and component Euler
checks must finish before freezing and starting seeds 140/141. No fourth
intervention will be made, regardless of the held-out result.

## Frozen evaluation

Independent development verification passed: 20 fields, 4,207,184 saved C++
queries, all 90 saved meshes rechecked independently, and both targets plus
larger controls recovered at the required resolutions. The selection and
computational source hashes are frozen in `runs/fidelity-frozen/selection.json`.
Its development score is 0.15514543 mean absolute relative width error across
six target case/axis groups. Fresh seeds 140/141 use 100,000 queries per
distribution. Three interventions are exhausted; evaluation will not alter
the method. Final evidence is verified and reproduced from saved source
in a new CMake build directory, including all C++ suites.

## Final results and completed reproduction

Fresh seeds 140/141 pass unchanged recovery for both targets and all larger
controls. Gap .060 improves from 68.75% outside recall to 100% on both seeds;
torus .035 improves 100%/99.21875% to 100%/100%. Actual bytes stay 32,628.
Mean absolute relative width errors fall from 30.13% to 19.12% radially,
30.96% to 19.17% vertically, and 27.55% to 5.32% for the gap. Gap width is
.06319079 vs .060; torus median thickness remains roughly 18--21% too small.

Near MAE improves 11.04% (torus) and 11.11% (gap), but near CVaR99 worsens
32.24% and 49.66%; uniform MAE worsens 21.59% and 14.16%. The retained
geometry benefit has an explicit error-tail cost. Oracle calls rise from
311,640 to 325,976. Initial target construction means: torus 73.45 to 94.81 ms,
gap 64.74 to 91.84 ms. Full metrics, controls and widths are in the final report.

The final verifier passed 20 fields, 4,207,184 saved C++ queries, all residual
ranks, augmented sample coordinates/objectives, byte/call accounting, and all
90 saved meshes including retained 32-cell failures. A fresh MSVC Release build
from `runs/fidelity-final-source`, in `runs/fidelity-rebuild`, passed all four
C++ suites. The saved source's own runner and independent verifier executed
`runs/fidelity-reproduction`; all 4,207,184 predictions and 90 meshes verified.
Reproduction skips redundant width scans; exact packed-field identity implies
identical width inputs and results. All 245 deterministic outputs match exactly:
20 fields, 50 samples, 80 predictions, 90 meshes, five CSV/JSON evidence files.
Timings are independently measured and excluded from exact matching.

Evidence is retained for 200 independently checked fields across seven screens,
development confirmation, final evaluation and reproduction: 16,871,840 saved
C++ query predictions in total. Earlier baseline buffers also match the previous
selective-depth-six evidence byte for byte on all 10 development cases.
The figure was visually inspected; report tables and source/export manifests
are checked. No remote publication or changes to earlier raw runs were made.

**Stopping decision:** the predeclared maximum of three interventions is reached,
and the strongest development choice is independently reproduced. This bounded
research round is complete. Future work
would need a new question/round, would treat these held-out cases as development
once used to guide an intervention, and would reserve fresh seeds. The remaining
torus dimensional error and worsened tails are explicit limits, not hidden
successes or a claim that further improvement is impossible.
