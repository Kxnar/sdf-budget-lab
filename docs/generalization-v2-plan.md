# Representation diagnosis and bounded implementation plan

This records the development diagnosis established before v2 implementation. It is development evidence, not a final accuracy claim.

| Old qualified mesh | Wrong v1 feature witnesses | Rays with incorrect root count | Shortest selected interval |
|---|---:|---:|---:|
| Goathead | 55 | 42/96 | .000441829 |
| Spot | 46 | 24/96 | .000654012 |
| Brucewick | 73 | 47/96 | .000198213 |

Most diagnosed wrong witnesses lie in retained depth-five or depth-six cells even though the candidate tree permits depth seven. Depth-seven cell width is .015625. Query-local full depth-seven q8/unquantized checks still misclassify 45/16/56 diagnosed witnesses; full depth twelve classifies all previously wrong witnesses correctly. That is not a complete-ray, width, topology or joint-gate result. Increasing the old scan to 16,384 samples changes one/three/five ray counts, so a fixed scan is also an evaluator limitation. Penguin remains reference-unresolved.

Implementation order is the [v2 protocol](generalization-v2-protocol.md): correct the field-ray evaluator and independently reassess old evidence; quantify query-local depth-7/9/12 q8/f32/unquantized feasibility on every old selected probe and ray; run one resource pilot and the exactly five-policy development screen only after the feasibility checkpoint; confirm a single policy; then freeze its settings, fresh source split, reference inputs, evaluator and executable before fresh fitting. The bounded screen permits no implicit depth, sample-weight, budget or encoding expansion. Representation or resource infeasibility requires a quantified checkpoint and revised representation plan.

The focused method changes training coverage and allocation cost explicitly. Its fitting-only 32,768 anchors and paired multiscale offsets enter before the base tree is ranked, and normalized-max ranking applies at every stage, including shallow cells. Original, global and every-stage quadrature samples retain duplicates and declared weights. Exact tree allocation does not imply a globally optimal geometric representation. Construction starts with an empty bounded LRU, reports its extra labels and cold costs, and cannot load the large v1 cache.

Fresh planned inputs are fish, mushroom, springer and hammer at pinned upstream revision 835213277571d625722d25880c9c3010c6f45a5a. Expected individual licenses are public domain, CC BY-NC 4.0, GPL v2 and CC BY 4.0 respectively; these are to be verified from individual READMEs before acquisition/qualification. Stuffedtoy then boat are predeclared acquisition/license/input-validation fallbacks. Reference uncertainty and all representation failures remain in the report.
