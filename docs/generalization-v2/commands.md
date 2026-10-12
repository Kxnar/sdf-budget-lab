# Current v2 commands

Run these from the repository root in PowerShell. The current study directory is `runs/generalization-v2/fresh-v2`. Its original `FROZEN.json`, source snapshot and executables remain authoritative; numbered procedural/reporting helpers preserve later driver and presentation corrections without changing fields or gates.

Strict completion checks and resumable serial execution use **resume2**, not the original frozen runner. This validates zero exit, COMPLETE, full job identity, manifests, the preparation declaration and development selection:

```powershell
.venv/Scripts/python.exe scripts/generalization_v2_resume2.py runs/generalization-v2/fresh-v2
```

Full independent verification uses the frozen evaluator and immutable per-job caches:

```powershell
.venv/Scripts/python.exe scripts/generalization_v2_verify.py runs/generalization-v2/fresh-v2
```

Rebuild the frozen source separately. Boost headers are required even when CGAL is OFF; this study's mesh validation/executable path uses CGAL ON and the existing classic-vcpkg installation:

```powershell
cmake -S runs/generalization-v2/fresh-v2/source -B runs/generalization-v2/fresh-v2/rebuilt -DSDF_WITH_CGAL=ON -DCMAKE_TOOLCHAIN_FILE=C:/deps/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_MANIFEST_MODE=OFF
cmake --build runs/generalization-v2/fresh-v2/rebuilt --config Release --parallel 1
ctest --test-dir runs/generalization-v2/fresh-v2/rebuilt -C Release --output-on-failure
# Run the unchanged, hashed v1 fixture supplement after the preserved CTest exception:
runs/generalization-v2/fresh-v2/rebuilt/Release/sdf_tests.exe runs/generalization-v2/fresh-v2/test-fixtures
```

The original frozen v2 source inventory omitted `fixtures/*.off`: its first CTest invocation passed 6/7 and the core test failed because its fixture path was absent. The failed `rebuilt/Testing/Temporary/LastTest.log` is retained. No frozen source changed. The separate unchanged v1 fixture supplement is bound by `REBUILD-TEST-FIXTURES.json`; the rebuilt core passed 10,081 checks against it, so all seven rebuilt executables passed. `REBUILD-PROVENANCE.json` binds this exception, supplemental success, source and binary hashes. Build output is explicitly partial, not a complete captured log.

With verification stopped, reproduce all declared method/mesh endpoints plus every contributing fitting-seed crossing and adjacent smaller tested cap. **reproduce_amended3** also saves independent complete baseline event/signed-witness replays; compact flags, exact events, witnesses, surfaces and decoded rejection outcomes are compared. Timing columns are retained separately:

```powershell
.venv/Scripts/python.exe scripts/generalization_v2_reproduce_amended3.py runs/generalization-v2/fresh-v2
```

Generate the current numerical report and scientific figures using **report_amended6** and **figures3**:

```powershell
.venv/Scripts/python.exe scripts/generalization_v2_report_amended6.py runs/generalization-v2/fresh-v2
.venv/Scripts/python.exe scripts/generalization_v2_figures3.py runs/generalization-v2/fresh-v2
.venv/Scripts/python.exe tests/test_generalization_v2_procedural.py
.venv/Scripts/python.exe tests/test_generalization_v2_replay.py
```

The generated report, full logical matrix and compact provenance live in `docs/generalization-v2/final`. Final presentation links/UTF-8 checks and bundle hashes are recorded separately after visual inspection. Large native buffers, prediction streams, reference meshes and rebuilt-source reproduction remain local under `runs/`, which is Git-ignored. Earlier numbered helpers and amendments remain available for provenance; they are not alternative evaluation settings.

After reproduction and figure QA, `scripts/generalization_v2_publish3.py` checks the exact prescribed completed reproduction selection, observed zero exits and fixed process guard, then publishes all twelve complete compact assessments, rebuild/fixture provenance, current interpretation/links, UTF-8 validation and final bundle hashes. Reporting amendment 10 binds these presentation helpers and this guide. Earlier figures2/publish2 and their amendments are preserved. The publisher preserves historical text and updates only its marked current-publication sections:

```powershell
.venv/Scripts/python.exe scripts/generalization_v2_publish3.py runs/generalization-v2/fresh-v2
```
