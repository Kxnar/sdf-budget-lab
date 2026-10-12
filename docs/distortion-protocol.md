# Geometry distortion diagnostic v1

Frozen before execution, 8 October 2026. Diagnose the existing translated
0.035 torus and 0.06 sphere gap. Include every larger feature from the existing
tail-v1 sweep (torus 0.06/0.10, gap 0.10). No threshold is claimed below these
tested features. Keep the original shift (0.017,-0.023,0.011), feature probes,
>=95% recalls, manifold/component/Euler checks, oracle checks at BOTH 64/96,
both seeds and success at every larger tested feature unchanged. Retain 32
as a diagnostic extraction. Thickness/width are supplementary, not replacement
success criteria.

Compare fully subdivided q8/f32 at depths five and six, including empty cells.
A direct full builder uses the same codec, packed format and query arithmetic;
test its bytes against existing full candidate packing. No depth-six candidate
tree or exact allocation is constructed. Depth-five q8 MSE at 32,768 bytes is
the budgeted comparator. Fitting seeds 30/31, 10,000 points, 50/50 uniform/near;
100,000 uniform and near evaluation points each, seeds +10000/+20000, noise
feature/2. These are existing cases/sampling settings, not new held-out evidence.
Full fields do not depend on fitting; repeated builds measure two serial timings.

Save raw fields, samples, C++ query outputs, metrics, meshes and construction
costs. Full construction queries only terminal corners (duplicates counted),
so report this separately from candidate/DP construction. Independent Python
decoding must reproduce every saved C++ evaluation query and frozen metrics.

Localize error at 1,024 torus centreline angles and 1,025 points on the gap
axis over [-gap,+gap]. At the exact original probes record leaf depth, fitting
count and whether correct-sign fitting samples exist in the same selected
leaf. For full fields compare q8 minus f32, and f32 minus double-corner
trilinear interpolation, separately from interpolation minus oracle. These
are pointwise decompositions, not additive absolute-error metrics.

Measure torus radial and vertical contiguous negative intervals around the
analytic centreline at 256 angles, and the contiguous positive sphere gap
around its midpoint. Locate first sign transitions by dense scan (step <=
feature/2048) then 40 bisections. A centre with wrong sign has width zero;
retain missing/bracket failures rather than substituting an expected width.
Adaptive sign transitions can be discontinuities rather than true roots.

Choose the next intervention from this evidence. Local precision is justified
only if f32 resolves the failure; deeper spatial candidates if deeper full
fields resolve it; allocation/coverage if full fields pass and budgeted fail.
Any subsequent intervention must lower the minimum recovered feature at the
same 32 KiB cap with the frozen recovery rule, reporting uniform/near MAE,
CVaR99 and full construction cost. This diagnostic is not that equal-byte win.
