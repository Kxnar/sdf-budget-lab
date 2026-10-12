# SDF Budget Lab

A C++ numerical research project: allocate a fixed storage budget to preserve
thin 3D geometry while controlling distance error and query latency.

This extends the question in [neural-geometry](https://github.com/Kxnar/neural-geometry):
what disappears when an implicit shape representation has limited capacity?
Here, capacity is a byte budget for a classical distance field. The analytic
shapes use the same dimensions and sign convention as that project. No neural
models are trained by this code.

[First development results and their limits](docs/development.md) include
multi-seed accuracy, feature diagnostics and visual comparisons.

[Fixed-candidate allocation study](docs/allocation-study.md) adds an exact
budget DP, greedy and bounded-lookahead comparison, and a frozen alignment
evaluation. In 216 tested case/seed/budget combinations, greedy matched the
exact optimum; exported split gains independently certify why. This guarantee
concerns fitting SSE within the candidate tree, not held-out geometry accuracy.

[Joint spatial and precision study](docs/precision-study.md) adds executable
f32/q16/q8 packed fields, exact byte-budget allocation, and 3D surface topology
sweeps. Across 288 fields, compression improves accuracy per byte, but joint
precision selection adds little over fixed q8 and does not recover the thinnest
features. The report includes decoding cost, failed surfaces and frozen evidence.

[Direct CVaR99 allocation with fixed q8](docs/tail-study.md) first checks full
depth-four/five fields, then compares MSE, sign penalty and tail allocation on
216 depth-five fields with 100,000 evaluation points per distribution. At
32 KiB, tail allocation reduces uniform CVaR99 by 10.6%, but near-surface MAE
rises 20.5% and minimum recovered geometry is unchanged between objectives.

[Geometry distortion diagnostic](docs/distortion-study.md) isolates the remaining
0.035 torus and 0.06 sphere-gap failures. Full depth-five f32 barely changes q8
retention; full depth six recovers both under the frozen criteria. Signed-error
profiles, measured widths and fitting coverage identify spatial interpolation
as the main limit, with an additional allocation effect in one gap seed. This
unequal-byte diagnostic motivates selective deeper q8 candidates.

[Selective depth-six refinement](docs/refinement-study.md) now recovers torus
radius **0.035 at 32,628 bytes**, versus the depth-five minimum tested radius
0.060. q8 and MSE stay fixed; 32 parents are chosen from fitting residuals using
a development-only rule. Both fresh seeds and every larger control pass the
unchanged recovery criteria. Near-surface MAE improves 19.2%, while uniform MAE
worsens 7.25%; widths remain too small. The 0.060 sphere gap still fails because
the MSE allocation retains none of its deeper candidates.

[Allocation coverage and geometric fidelity](docs/fidelity-study.md) completes
a three-intervention follow-up. A frozen sign-aware allocation with 128 deeper
candidates and deterministic interior fitting samples recovers the **0.060 sphere
gap at 32,628 bytes** on fresh seeds, preserving the 0.035 torus and every larger
control. Torus dimensional error falls to about 19%, but near-surface CVaR99
worsens 32% for the torus and 50% for the gap. The report includes unsuccessful
experiments, unchanged recovery criteria and independent saved-output checks.

[Frozen generalization study](docs/generalization-v1/report.md) compares retained
storage at stricter, joint geometric accuracy gates on 16 frozen instances,
including four novel analytic families and four imported meshes. The held-out
rounded box uses **25.0% fewer bytes than the strongest tested q8 comparator**
and **60.7% fewer than the grid** across five fitting seeds at tested passing
caps. The capsule has three q8 gains and two ties; all three reference-qualified
real meshes fail, and Penguin's reference topology remains unresolved. The
[protocol](docs/generalization-protocol.md) and
[compact provenance/reproduction guide](docs/generalization-v1/README.md)
retain censored cases, nonmonotonic failures, independent checks and costs.

The [current complete-ray interpretation](docs/generalization-v2-current-interpretation.md)
preserves those rounded-box and capsule comparisons but identifies fifteen
previously missed higher-cap transition failures and increases the shell's
verified crossing. Original v1 evidence remains intact. The separately declared
[completed v2 compact mesh study](docs/generalization-v2/final/report.md) compares Q24 with strong compact controls on four fresh meshes (three reference qualified). See its [full matrix](docs/generalization-v2/final/matrix.csv) and [current commands](docs/generalization-v2/commands.md).

## Allocation objectives

The central problem is constrained resource allocation, with numerical error
and rare geometry failures as competing objectives.

- **Allocation:** greedily choose the next octree split by estimated reduction
  in training loss per additional byte.
- **Risk:** report mean error, p99 and CVaR99 (mean of the worst 1% absolute errors),
  and track sign failures which can erase walls or fill gaps.
- **Experiments:** separate fitting and evaluation samples, pair evaluation
  seeds across methods, report run variability and stress a shifted query mix.
- **Trade-offs:** compare accuracy, storage and latency together; count oracle
  calls and construction time rather than treating preprocessing as free.

The two adaptive policies use the loss

`L(pred, truth) = ((pred - truth) / domain_scale)^2 + lambda * wrong_sign`.

`adaptive_mse` sets lambda to zero; `adaptive_risk` defaults to 0.0025. A split
adds eight nodes. Its priority is the training loss before minus after the split,
divided by those added bytes. This is a heuristic, not an optimal allocation
solver or a CVaR minimiser. Negative-gain candidates remain eligible to cross
refinement plateaus while spending the cap. Training-empty leaves are not split.

## Implemented scope

| Component | Behaviour |
|---|---|
| Uniform baseline | Cubic float32 grid with trilinear interpolation |
| Adaptive representation | Packed octree; eight float32 corners and a child index per node |
| Mesh reference | CGAL AABB tree with exact predicates/constructions; inside/outside sign |
| Input validation | Closed, triangulated, nonempty, nondegenerate meshes without self intersections |
| Analytic controls | Sphere, thin torus, thin box, two separated spheres |
| Outputs | Configuration, raw CSV, first-seed cross sections, leaf rectangles, completion marker |
| Analysis | Seed summaries, paired bootstrap intervals, plots and an offline HTML report |

Mesh vertices must be normalized into `[-1,1]^3`; queries outside that box fail.
The CGAL reference computes exact squared distance to the loaded triangles using
EPECK. Conversion to double and the final square root are rounded. Exact distance
to a tessellated torus is distance to that mesh, not to an ideal smooth torus.
Distance acceleration follows the [CGAL AABB-tree documentation](https://doc.cgal.org/latest/AABB_tree/).

All query nodes count toward storage, including internal nodes. Both methods use
the same maximum byte budget, counting the C++ object and its exact-size backing
array. Grid sizes and octree splits are discrete, so actual bytes can differ;
`representation_bytes` exposes utilization. This is storage of the representation,
not process RSS: allocator overhead, the oracle and temporary build structures
are excluded. Oracle rows have storage zero as a not-measured placeholder.
Build memory can substantially exceed retained storage.

Independent leaves can disagree at coarse/fine boundaries. This first version
does not balance the octree, share corners, guarantee continuity, certify
topology, or extract a 3D mesh. Probe retention and the cross section are practical
diagnostics. A mesh without custom probes has no feature-probe rows.

## Build and run

Requires CMake 3.21+, a C++17 compiler, Boost headers (including Multiprecision),
and CGAL for the mesh reference. The standalone compact field's exact predicates
require Boost even with `SDF_WITH_CGAL=OFF`; on Ubuntu install `libboost-dev`.
Use a release
build for timing. A pinned vcpkg manifest is provided.

```powershell
# vcpkg root may be anywhere; set it to your checkout.
$taskVcpkg = "C:\path\to\vcpkg"
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE="$taskVcpkg/scripts/buildsystems/vcpkg.cmake"
cmake --build build --config Release --parallel 4
ctest --test-dir build -C Release --output-on-failure
.\build\Release\sdf_bench.exe --shape torus --feature 0.04 --out runs\torus-first
```

For a classic-mode vcpkg checkout, disable manifest mode. GMP's dependency scripts
require a path without spaces; use a checkout or junction at such a path:

```powershell
$taskVcpkg = "C:/deps/vcpkg"
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE="$taskVcpkg/scripts/buildsystems/vcpkg.cmake" -DVCPKG_MANIFEST_MODE=OFF
cmake --build build --config Release --parallel 4
ctest --test-dir build -C Release --output-on-failure
```

For analytic-only development, use `-DSDF_WITH_CGAL=OFF`. On single-configuration
generators set `-DCMAKE_BUILD_TYPE=Release`; executables are directly in `build/`.

```powershell
# Four shapes, three fitting seeds, three byte caps; creates a new directory.
python scripts\sweep.py --exe build\Release\sdf_bench.exe --out runs\development

# Read-only input mesh plus optional known-point x,y,z CSV.
.\build\Release\sdf_bench.exe --mesh fixtures\thin_box.off --probes fixtures\thin_box_probes.csv --out runs\mesh-first

# Optional plotting environment, independent of the C++ executable.
python -m venv .venv
.\.venv\Scripts\python -m pip install -r requirements-analysis.txt
.\.venv\Scripts\python scripts\report.py runs\development
```

Open `runs/development/report.html` for figures and links to raw measurements.
Benchmark and sweep commands refuse existing output directories. Plotting may be
repeated to regenerate derived figures and summaries. `scripts/make_fixtures.py`
regenerates the tracked synthetic fixtures. `sdf_bench --help` lists all options.
The sweep driver also accepts `--mesh`, `--probes` and `--slice-plane` and records
input/executable hashes, source fingerprints and machine details for mesh runs.
Torus analytic slices default to the xy plane; other slices default to xz.

## Experiment protocol

For each shape and seed, fit on 8,192 reference samples: half uniform, half
area-uniform surface points plus Gaussian noise with sigma `feature/2`.
The adaptive methods receive the same fitting pool and exact corner queries.
Uniform grids query the reference at their vertices. Construction label counts
therefore differ and are reported; this comparison holds retained storage fixed,
not preprocessing compute or oracle access.

Each evaluation distribution uses 8,192 new samples. Uniform tests use seed
10,000 + seed; near-surface tests use 20,000 + seed. A shifted distribution uses
90% surface queries and four times the original noise, with seed 30,000 + seed.
This is distribution shift on the same geometry, not unseen-mesh generalization.
Deterministic interior/exterior probes diagnose specific thin features.

Query timing warms the field, then measures 11 serial CPU batches on the uniform
evaluation samples. `median_ns_per_query` and `p95_batch_ns_per_query` are batch
averages; the latter is not individual-query tail latency. A checksum prevents
discarding results. Oracle timing includes the C++ callable dispatch. Results
depend on compiler, machine load and cache state. Timing repeats do not provide
independent machine-level uncertainty estimates.

Distances are in normalized scene units. Signs exclude reference distances within
2e-9 of the surface; a zero prediction at a nonboundary point is a sign failure.
p99 uses nearest rank. CVaR99 averages the largest `ceil(0.01*n)` errors. Missing
interior/exterior classes give undefined recalls; counts accompany scores.
Bootstrap intervals resample paired seed-level differences. With three seeds,
interpret them as coarse development variability, not evidence across all meshes.
Seeded sampling is repeatable within the recorded toolchain. The C++ standard
library's normal distribution may generate different points on another compiler;
cross-platform bitwise reproduction is not promised.

## Research decisions and next milestone

Current runs are development evidence. Fix settings before a held-out comparison;
do not tune lambda on the final test matrix. The first hypothesis is that adaptive
allocation improves local feature retention at small storage caps, accepting a
possible query-speed or global-error cost. The second is that the sign penalty
improves retention beyond MSE allocation, while potentially worsening error tails.

Next, use development feature sizes 0.02, 0.04 and 0.08 to pick settings, then freeze
them and evaluate unseen thicknesses 0.015, 0.03 and 0.06 with fresh fitting seeds.
Include translated and rotated shapes to expose alignment dependence. Add stronger
surface/3D topology measurements and more real closed meshes before making broad
preservation claims. Subsequently add continuity constraints or quantized values
only if measured failures justify them. A neural baseline can join the same
accuracy-storage plots later, with tensor storage and CPU/GPU timing distinguished.

Evidence records what was actually run; no novelty or universal improvement is
claimed.
