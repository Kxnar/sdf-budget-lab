# Complete-ray feasibility checkpoint

No deeper-trilinear development screen or fresh held-out representation was fitted. The fixed depth-7/9/12 controls were evaluated on all old selected probes and all 96 rays for each qualified mesh, with independent field polynomials/root counts and corner-oracle checks. They are query-local virtual controls, not retained-memory contenders and not full-gate certification.

| Mesh, depth 12 | q8 incorrect-count rays | f32 incorrect-count rays | Unquantized incorrect-count rays | Probe/width interpretation |
|---|---:|---:|---:|---|
| Goathead | 1 | 0 | 0 | f32/unquantized pass probes and widths; q8 introduces extra transitions |
| Spot | 1 | 0 | 0 | f32/unquantized pass probes and widths; q8 introduces extra transitions |
| Brucewick | 3 | 1 | 1 | All pass the per-group 95% recall gate, but exact signed interval counts fail |

There is no recorded numerical ambiguity in these depth-12 failures. Brucewick ray 0 loses the two reference transitions at -0.037580215845532089 and -0.037382002810736802, an interval of length .000198213. The interval midpoint independently has signed distance -7.7642888269e-6, while full-depth-12 unquantized trilinear interpolation predicts +2.5642105855e-5. Q8 additionally introduces transitions on rays 6 and 11; ray 11 contains 36 transitions versus reference 6. This separates a resolution/interpolation failure from additional quantization sign instability. No gate or eligibility threshold was changed.

Independent brute-force triangle distances and ray-parity signs checked 144/144/152 diagnostic corners, including every corner of the missing-interval cell. Maximum label discrepancy was 1.74e-16, with no sign disagreement outside the 2e-9 boundary band. The raw supplementary JSON field name `sign_mismatches_away_2e9_band` is a labeling typo; its executable predicate uses **2e-9**, as the preserved source shows.

Native depth-12 diagnostic phases took 14.79/9.13/11.97 seconds, with 2,592,143/2,628,330/2,579,492 physical oracle evaluations and peak working sets below 50 MB. These are cache-reset phase times with an already constructed oracle and prior process warming; they include corner/probe labeling, streamed segment writes and all three encodings' root evaluation. Oracle setup and independent verification are excluded. No cold-process construction claim or speed ratio follows.

The 262,144-entry LRU stayed bounded and was emptied before each depth. No v1 persistent cache was loaded or appended. The checkpoint requires a quantitative representation plan; the five-policy depth-12-only screen cannot establish the declared interval gate on all three old meshes as currently specified.

See `feasibility-summary.csv`, `feasibility-CHECKPOINT.json`, `oracle-checks-COMPLETE.json`, the three cost CSVs and `CHECKPOINT-PROVENANCE.json`. Full source/executable and corner/segment evidence are under `runs/generalization-v2/`. The versioned [v1 reassessment](../generalization-v2-current-interpretation.md) preserves the rounded-box/capsule crossings but corrects fifteen higher-cap passes.
