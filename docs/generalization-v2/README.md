# Generalization v2 evidence

See the [current interpretation](../generalization-v2-current-interpretation.md), [protocol](../generalization-v2-protocol.md) and [bounded implementation plan](../generalization-v2-plan.md). Reassessment is separate from the unchanged v1 study. The completed [fresh compact mesh study](final/report.md) is a separate representation-family comparison.

`reassessment-matrix.csv` contains every reassessed field's original and updated flags, numerical/geometry status and ray diagnostics. `reassessment-COMPLETE.json` binds evaluator and original freeze hashes. Large per-field ray evidence, source/executable snapshots and virtual full-depth diagnostic corner streams are kept in `runs/generalization-v2/` outside Git. Virtual full-depth controls are not storage contenders or full-gate certification.

The [feasibility checkpoint](feasibility-checkpoint.md) finds one remaining unquantized Brucewick ray failure and additional q8 transition failures. That checkpoint superseded the unrun deeper weighted-octree screen. The [compact mesh plan](../generalization-v2-compact-plan.md) selected Q24 using only old development meshes; the four fresh inputs and all final settings were frozen before held-out fields were built.

<!-- v2-publication -->
On the three reference-qualified fresh meshes, Q24 used 9.64% less retained field storage at tested passing caps than the strongest passing compact control (f32); the per-mesh range was 9.22% to 10.18%. The declared held-out suite contained four meshes: Fish remained reference unresolved. All frozen sampled-distance, feature, signed-ray dimension and two-resolution extracted-topology gates were required.

[Full final matrix](final/matrix.csv), [report and figures](final/report.md), [strict runner/verifier/rebuild/reproduction/report commands](commands.md), and [final provenance manifest](final/SHA256SUMS.json). Large binary evidence remains local under `runs/generalization-v2/fresh-v2`.
