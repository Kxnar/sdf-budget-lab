# Joint spatial and precision allocation: protocol v1

Frozen before the first evaluation run, 8 October 2026. This is a new,
post-allocation-study experiment, not an extension of that study's held-out set.
The objective remains fitting SSE; tail-objective allocation is deferred.

## Representation and exact reference

Each node has a four-byte tagged offset. Branches point to eight consecutive
child descriptors and retain no values. Float32 leaves have 32 payload bytes;
local unsigned 16-bit leaves have two float32 parameters plus 16 value bytes;
local unsigned 8-bit leaves have two float32 parameters plus eight value bytes.
Every payload and descriptor is four-byte aligned. Total node costs are 36,
28, 20 bytes for leaves and four bytes for branches. Offsets include tags in
their low two bits. Count the actual C++ query object, including its domain,
pointer and allocation length, plus the exact byte-array allocation. Allocator
bookkeeping, oracle, fitting data, candidate tree and DP scratch are excluded.
Report DP states/transitions, not a claim of peak working-set memory.

Quantization uses minimum corner value as offset and (maximum-minimum)/Q as
scale, both stored as float32; Q=65535 or 255. Round to nearest integer with
clamping. A constant cell has zero scale. Decoding computes offset+scale*q in
double and rounds to float32 before the shared double trilinear interpolation.
Fit losses use precisely these decoded values. Numerical precision can change
the sign; there is no sign preservation or continuity constraint.

For each candidate v, exact-cost DP admits three leaf (cost, SSE) options or a
branch with cost four plus the sum of eight child costs. Child convolution and
Pareto pruning discard only states with at least as much cost and no less loss.
The root chooses minimum loss below the cap, preferring fewer bytes on ties.
Costs use exact four-byte units with floor((budget-object_bytes)/4); no lossy
cost binning, epsilon dominance or heuristic pruning. Guarantee: minimum SSE
among encoding choices and prunings of this fixed candidate, subject to floating
point arithmetic. It is not a guarantee of test accuracy or geometry.

Compare exact spatial allocation restricted to f32, q16, q8, and joint choice.
All use the same layout, candidate, samples, decoder and objective. Also time a
float32 copy of each selected topology separately as a decoding-cost diagnostic;
that diagnostic does not participate in equal-byte accuracy comparisons.

## Sweep and samples

Three shapes: thin box (half-thickness), torus (tube radius), double sphere
(surface gap). Features: 0.02, 0.035, 0.06, 0.10. Translate every shape by
(0.017,-0.023,0.011) to avoid alignment with octree/grid planes. No rotation in
this initial controlled sweep. Two fitting seeds: 0 and 1. Full depth-four
candidate (4,681 nodes). Fit 4,096 points, half uniform, half near-surface with
sigma feature/2. Budgets 8,192, 16,384, 32,768 bytes. Four methods give 288
fields. This finite sweep is descriptive, not a population confidence claim.

Evaluation uses independent streams: 32,768 uniform points (seed+10000), 32,768
near points (seed+20000), and transformed deterministic feature probes. Report
RMSE, mean absolute error, p99, CVaR99, sign error, inside/outside recall. No
evaluation point selects the representation. Export exact fitting/evaluation
points, all encoding losses, choices and packed buffers. Quantiles and CVaR
follow the existing metric convention (CVaR averages ceil(0.01*n) worst errors).
The roughly 328 tail observations per distribution are descriptive only.

Report smallest tested budget meeting near RMSE <=0.02 AND inside/outside probe
recall >=0.95 as a fixed-accuracy diagnostic; also report near-RMSE-only <=0.02.
These are observed crossings on the discrete sweep, not globally minimal bytes.
Do not assume monotonic held-out accuracy as the fitting budget increases.
Timing: report candidate construction, exact DP, backtracking/packing, and
three serial query batches after warmup. Same-topology f32 timings are paired;
use medians and retain raw repeats. One DP at the maximum budget supports all
three caps; its cost is shared, not silently divided between fields.

## Geometry success, defined before results

Extract a real triangle surface with six conforming marching tetrahedra per
cube at 32, 64 and 96 cells per axis, for the packed field and the oracle.
Use the domain [-1,1]^3; exact zero lattice values receive an outside perturbation
of 1e-14 times domain scale. Index shared lattice-edge intersections, check
edge incidence and single-cycle vertex links, count connected components and
Euler characteristic per component. Save OBJ meshes at 96 cells.

Wall succeeds only if the mesh is one closed manifold component with Euler
characteristic 2 and both probe recalls >=0.95. Torus succeeds only if one
closed manifold component has Euler characteristic 0 and both recalls >=0.95.
Gap succeeds only if exactly two closed manifold components each have Euler
characteristic 2 and both recalls >=0.95. Topology alone is insufficient to
pass. Report topology and probes separately. Success must hold at both 64 and
96, and oracle extraction must also pass at both. Report disagreements with 32
and between 64/96; these three grids do not prove convergence. The minimum
recoverable feature is the smallest tested feature with success at that and
every larger tested feature, for both seeds. If no such suffix exists, report
not reached. Box full thickness and torus full tube diameter equal 2*feature;
sphere gap equals feature.

Adaptive leaf boundaries may be discontinuous. Extraction samples the field
onto a conforming lattice, so mesh topology characterizes that sampled
reconstruction, not a proof about the continuous adaptive representation.
Oracle success with field failure diagnoses a representation/allocation failure
at the tested grids; oracle failure leaves extraction resolution unresolved.

Variable-precision storage is established work, including
[NanoVDB's Fp16/FpN construction](https://www.openvdb.org/documentation/doxygen/tools_2CreateNanoGrid_8h_source.html).
The contribution here is the explicit joint tree knapsack and controlled
equal-byte comparison, not invention of variable-precision volume storage.
