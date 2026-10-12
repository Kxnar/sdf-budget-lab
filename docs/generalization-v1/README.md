# Generalization v1 evidence

The protocol is `../generalization-protocol.md`. The machine-readable freeze, immutable source and inputs, paired reference samples, packed fields and verification records are in `../../runs/generalization-v1/`. Large binary evidence is intentionally excluded from Git by the existing `runs/` rule; the compact report does not replace it.

Run the following from the repository root in PowerShell. The frozen runner checks its own snapshot and all input hashes before resuming; it retains completed jobs and rejects incompatible data. Use the frozen runner rather than the edited live evaluator/report helpers to resume native fitting.

```powershell
.venv/Scripts/python.exe runs/generalization-v1/source/scripts/generalization_run.py runs/generalization-v1
.venv/Scripts/python.exe scripts/generalization_verify.py runs/generalization-v1 --merge-shards runs/generalization-v1/verification-shards
.venv/Scripts/python.exe scripts/generalization_reproduce.py runs/generalization-v1
.venv/Scripts/python.exe scripts/generalization_report.py runs/generalization-v1
.venv/Scripts/python.exe scripts/generalization_figures.py runs/generalization-v1
```

Stop this study's independent validation workers before fresh-process timing. The separately rebuilt native source remains frozen. Selected grid timings use a fresh process for each cap; adaptive timings construct the maximum allocation frontier without a shared on-disk cache. An otherwise idle host is not assumed.

`VERIFIER_CORRECTION.json` records a post-freeze correction enforcing the already-declared signed interval predicate. Original native width flags, extracted surfaces and source hashes remain available. `analysis-provenance.json` binds the corrected evaluator and report revision; `figure-provenance.json` binds plotted data and source meshes. The compact `SHA256SUMS.json` covers the final files in this directory.

Each crossing is the smallest tested actual retained-byte count passing every gate. No finite percentage is assigned to a censored comparator, and no exact continuous minimum is inferred between tested caps. Reference-resolution uncertainty is kept separate from a representation failure.
