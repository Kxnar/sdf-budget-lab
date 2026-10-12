# Current interpretation of the v1 result

The original [v1 report](generalization-v1/report.md), all original fields and native/evaluator flags remain preserved. A separately versioned complete-ray audit replaces the fixed 1,024-step scan with actual leaf/grid segments, cubic roots and one-sided interface signs. Native derivative-bracket roots are checked independently using power-basis polynomials, rational multiplicity/Sturm counts and companion roots.

All 243 previously joint-passing fields and all 116 reproduced fields were included, together with the three diagnosed mesh fields: 344 distinct case/seed/method/cap rows, 312 unique field buffers. Fifteen original passes fail the corrected root-count predicate: cylinder09 at 512 KiB, rounded box at 1 MiB, and spherical shell at 512 KiB, each across five fitting seeds. These are observed sign-transition failures, not unresolved numerical assessments. Other metrics and already extracted topology remain independently verified and unchanged.

The central rounded-box comparison remains valid at its tested 196,500-byte passing cap in all five seeds, versus 262,020 bytes for the strongest tested q8 comparator and 500,064 bytes for grid: 25.0% and 60.7% less retained field storage. Capsule06 remains at 130,980 bytes, with three q8 gains and two ties. The shell's smallest currently verified passing field increases from 524,244 to 1,048,548 bytes in all five seeds. This reinforces the requirement to check every cap independently: passing a lower budget does not imply passing a larger one.

This audit reassesses the previously passing and reproduction/diagnostic fields, not the entire originally rejected matrix. A newly passing ray assessment on a field whose topology was skipped does not become a joint pass. Numerical ambiguity is recorded separately from geometric failure; sampled distance bounds remain sampled rather than continuous certificates.

Compact reassessment data and evaluator provenance are in [generalization-v2](generalization-v2/README.md). Full ray evidence and the interrupted first evaluator version are preserved under `runs/generalization-v2/`. The ongoing bounded mesh intervention is governed by the [v2 protocol](generalization-v2-protocol.md) and [implementation plan](generalization-v2-plan.md); its subsequent compact branch is now evaluated on fresh held-out inputs.

<!-- v2-publication -->
## Completed compact representation comparison

On the three reference-qualified fresh meshes, Q24 used 9.64% less retained field storage at tested passing caps than the strongest passing compact control (f32); the per-mesh range was 9.22% to 10.18%. The declared held-out suite contained four meshes: Fish remained reference unresolved. All frozen sampled-distance, feature, signed-ray dimension and two-resolution extracted-topology gates were required.

The candidate stores explicit triangles from the source mesh with shared quantized vertices and its own BVH. This changes the representation family and input requirement; it is not a q8 allocation/interpolation breakthrough. The old-mesh diagnostic found a missed Brucewick interval even in depth-12 unquantized trilinear interpolation. Adaptive/grid comparators without a joint crossing through 1 MiB stay censored, with no invented percentage saving. Compact ordinary queries cost more than interpolation; batch sizes differ across families, so no cross-family speed ratio is claimed.

See the [final report](generalization-v2/final/report.md), [full matrix](generalization-v2/final/matrix.csv) and [current evidence commands](generalization-v2/commands.md). Historical sphere-gap gains, tail regressions and the separately versioned v1 ray reassessment remain intact.
