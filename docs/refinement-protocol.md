# Selective depth-six refinement v1

Frozen before development execution, 8 October 2026. Start with torus, then
apply exactly the same rule to sphere gaps. Fixed q8, decoded normalized SSE
objective, exact tree DP, 32,768 bytes including the query object. Preserve the
full depth-five candidate tree and all its options. Rank occupied depth-five
cells by descending sum of squared q8 fitting residuals (domain scale squared
normalization); break ties by original node index. Add all eight depth-six
children to the first K positive-loss cells. No evaluation/probe information,
sign penalty, shape-specific score, or changed fitting distribution is used.

Development: translated torus and double sphere, features 0.03 and 0.08,
fitting seeds 201/202, 10,000 samples, 50/50 uniform/near, noise feature/2.
Try K=32/64/128/256/512/1024, plus K=0 baseline. Choose the smallest positive K
whose final allocated fitting SSE is within 1% on average and 5% in every
development case of K=1024. If none smaller qualifies, use 1024. Equal weight
to each of the eight cases; no geometry tests or evaluation samples are
generated during selection. Record every search, cost, and selected-parent
count. A larger candidate set cannot worsen the exact fitting optimum.

Freeze K and development evidence hashes before final execution. Final: torus
radii 0.035/0.060/0.100, sphere gaps 0.060/0.100; fresh paired fitting seeds
40/41, same sampling parameters and shift (0.017,-0.023,0.011). Independently
build baseline and selective fields. Use 100,000 uniform and 100,000 near
evaluation points per case, offsets +10000/+20000. Same fit and evaluation
points for each paired comparison. Torus cases run first.

Unchanged recovery: both original sign-probe recalls >=95%, closed manifold,
expected component Euler values (torus [0], spheres [2,2]), at BOTH 64/96
extraction cells, oracle checks at both resolutions, both seeds, and success
at every larger tested feature. Keep 32-cell extraction diagnostic. Report
individual failures and controls. These are familiar shapes with fresh random
samples, not evidence of unseen-shape generalization. Two seeds do not support
population confidence intervals. Do not retune after evaluation.

Save fields, fit/evaluation/probe samples, all C++ predictions, selected-cell
rank/box/count/SSE and whether the allocation retains its split, meshes and
source/executable hashes. Independently decode bytes and verify saved queries,
losses, candidate ranking, counters, byte caps, and recovery aggregation.
Retain the existing width procedure: 256 torus angles, radial/vertical first
sign transitions around the analytic centreline; sphere axis gap interval;
step feature/2048 then 40 bisections. Wrong centre sign yields zero width;
missing brackets stay missing. Width does not replace recovery.

Report actual bytes, uniform/near MAE and CVaR99, width distributions, fitting
SSE, candidate and retained depth-six counts, oracle calls and construction
time. Total construction includes labels, full d5 generation/encoding/loss,
ranking and d6 generation/encoding/loss, DP and packing; excludes evaluation,
meshing and evidence writes. Count duplicate corner calls. Existing three-
encoding build container is retained, but only q8 may enter the packed field.
No temporary-memory or general runtime claim. Keep failures and interrupted
outputs; use new run directories. Success target: torus 0.035 and/or gap 0.060
under the cap with all larger controls passing.

## Frozen final selection

Development completed before final execution. K=32 passes: mean allocated SSE
excess 0.1713009911%, maximum 0.8653346701%, versus K=1024 across eight cases.
K>=64 matches the 1024-parent fitting optimum in all cases. Freeze **K=32**
under the predeclared smallest-qualifying rule. Development CSV SHA-256:
`ff188bf377364b374aa9a895611652be7de901ef8f38c5b62abf89ad5d4268a2`.
Evidence: `runs/refinement-v1-dev`, `runs/refinement-v1-selection/selection.json`.
No evaluation probes, geometry extraction or evaluation samples were used.
