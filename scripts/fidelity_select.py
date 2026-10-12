"""Freeze the completed development choice before any fresh-seed evaluation."""
import json
from pathlib import Path
import numpy as np
from tail_evidence import rows, hashes, sha, write_csv
from tail_run import ROOT

names=['i1-l0001','i1-l001','i1-l01','i2-k128','i2-k512','i3-q2','i3-q4']
summary=[]
for name in names:
    d=ROOT/'runs'/('fidelity-'+name)
    a=d.with_name(d.name+'-verified')
    hashes(d); hashes(a)
    assert (a/'VERIFIED.json').exists()
    m=[r for r in rows(d/'metrics.csv') if r['dataset']=='probes' and r['method']=='risk_q8']
    targets=[r for r in m if (r['shape']=='torus' and float(r['feature'])==.035) or (r['shape']=='double_sphere' and float(r['feature'])==.06)]
    controls=[r for r in m if r not in targets]
    w=[r for r in rows(a/'widths.csv') if r['method']=='risk_q8']
    keys={(r['case'],r['axis']) for r in w}
    assert len(keys)==6
    dimensional=np.mean([np.mean([abs(float(r['relative_error'])) for r in w if (r['case'],r['axis'])==key]) for key in keys])
    assert np.isfinite(dimensional)
    summary.append(dict(experiment=name,worst_target_recall=min(min(float(r['inside_recall']),float(r['outside_recall'])) for r in targets),controls_pass=all(min(float(r['inside_recall']),float(r['outside_recall']))>=.95 for r in controls),mean_axis_relative_error=dimensional,raw_sha256=sha(d/'SHA256SUMS.json'),verification_sha256=sha(a/'SHA256SUMS.json')))
chosen=min([r for r in summary if r['controls_pass']],key=lambda r:(-r['worst_target_recall'],r['mean_axis_relative_error']))
assert chosen['experiment']=='i3-q4'
confirm=ROOT/'runs/fidelity-dev-confirm-verified'
hashes(confirm)
assert all(r['recovered']=='True' for r in rows(confirm/'recovery.csv') if r['method']=='risk_q8')
out=ROOT/'runs/fidelity-frozen'
out.mkdir(exist_ok=False)
write_csv(out/'development_comparison.csv',summary)
files=[ROOT/'CMakeLists.txt']
for folder in ('src','include'):
    files += [f for f in (ROOT/folder).rglob('*') if f.is_file()]
frozen=dict(k=128,**{'lambda':.001},quadrature=4,development_choice=chosen,confirmation_sha256=sha(confirm/'SHA256SUMS.json'),held_out_seeds=[140,141],evaluation_queries_per_distribution=100000,source_sha256={f.relative_to(ROOT).as_posix():sha(f) for f in files},rule='Passing larger probe controls; maximize worst target recall, then minimize equal-weight mean absolute relative width error over six target case/axis groups. Full development recovery required before freeze.')
(out/'selection.json').write_text(json.dumps(frozen,indent=2))
(out/'SHA256SUMS.json').write_text(json.dumps({f.name:sha(f) for f in out.iterdir() if f.is_file()},indent=2))
print(json.dumps(frozen,indent=2))
