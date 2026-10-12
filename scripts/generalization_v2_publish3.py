"""Publish compact verified evidence and current interpretation; no field changes."""
import argparse,json,shutil

def base_document(text):
    return text.split('<!-- v2-publication -->',1)[0].rstrip()

from pathlib import Path
from generalization_common import *
from generalization_v2_reproduce_amended3 import verified_identity,selection

def completed_reproduction(d,f,r):
    expected=[(g['id'],s,m,caps) for g,s,m,caps in selection(d,f)]
    actual=[(x['geometry'],x['seed'],x['method'],x['caps']) for x in r['records']]
    assert r['jobs']==len(actual)==len(expected) and actual==expected
    assert r['helper_sha256']==sha(ROOT/'scripts/generalization_v2_reproduce_amended3.py')
    for x in r['records']:
        target=d/'reproduction-amended3'/f"{x['geometry']}_{x['seed']}_{x['method']}"
        proc=json.loads((target/'PROCESS.json').read_text())
        assert (target/'COMPLETE').exists() and proc['exit_code']==0 and not proc['resource_unresolved']
        assert proc['guard_bytes']==1610612736 and x['exact_field_and_predictions'] is True
    return dict(jobs=len(actual),observed_zero_exits=len(actual),selection_exact=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);a=ap.parse_args();d=a.directory.resolve();v=verified_identity(d);f=json.loads((d/'FROZEN.json').read_text());out=ROOT/'docs/generalization-v2/final';r=json.loads((d/'REPRODUCED-AMENDED.json').read_text());assert r['frozen_sha256']==sha(d/'FROZEN.json') and r['verification_manifest_sha256']==sha(d/'verification/SHA256SUMS.json');check_hashes(d/'reproduction-amended3',json.loads((d/'reproduction-amended3/SHA256SUMS.json').read_text()))
    completed=completed_reproduction(d,f,r)
    amendment=json.loads((d/'REPORTING-AMENDMENT-10.json').read_text())
    assert amendment['frozen_sha256']==sha(d/'FROZEN.json')
    for name,h in amendment['helpers'].items():assert sha(ROOT/'scripts'/name)==h
    assert amendment['commands_sha256']==sha(ROOT/'docs/generalization-v2/commands.md')
    assessments={}
    for g in f['geometries']:
        for m in ('compact_f64','compact_f32','compact_Q24'):
            p=d/'verification'/f"{g['id']}_0_{m}.json";assessments[p.stem]=dict(verified_cache_sha256=sha(p),assessment=json.loads(p.read_text())['assessment'])
    json_write(out/'compact-assessments.json',dict(frozen_sha256=sha(d/'FROZEN.json'),verification_manifest_sha256=sha(d/'verification/SHA256SUMS.json'),records=assessments))
    for name in ('REPORTING-AMENDMENT-6.json','REPORTING-AMENDMENT-7.json','REPORTING-AMENDMENT-8.json','REPORTING-AMENDMENT-9.json','REPORTING-AMENDMENT-10.json','REBUILD-PROVENANCE.json','REBUILD-TEST-FIXTURES.json','COMPACT-DECLARATION.json','REFERENCES_COMPLETE.json','BOUNDARY-RERUN-VERIFIED.json'):
        source=d/name if (d/name).exists() else ROOT/'runs/generalization-v2'/name
        assert source.exists(),f'Missing required publication provenance: {source}'
        shutil.copy2(source,out/name)
    shutil.copy2(ROOT/'runs/generalization-v2/compact-development-v2/SELECTION.json',out/'DEVELOPMENT-SELECTION.json')
    comparisons=rows(out/'paired_savings.csv');finite=[float(x['saving_percent']) for x in comparisons if x['comparator']=='compact_envelope' and x['saving_percent']!=''];assert len(finite)==3
    sizes={k:a['assessment']['actual_bytes'] for k,a in assessments.items()};median=float(np.median(finite));low=min(finite);high=max(finite)
    claim=f'On the three reference-qualified fresh meshes, Q24 used {median:.2f}% less retained field storage at tested passing caps than the strongest passing compact control (f32); the per-mesh range was {low:.2f}% to {high:.2f}%. The declared held-out suite contained four meshes: Fish remained reference unresolved. All frozen sampled-distance, feature, signed-ray dimension and two-resolution extracted-topology gates were required.'
    limits='The candidate stores explicit triangles from the source mesh with shared quantized vertices and its own BVH. This changes the representation family and input requirement; it is not a q8 allocation/interpolation breakthrough. The old-mesh diagnostic found a missed Brucewick interval even in depth-12 unquantized trilinear interpolation. Adaptive/grid comparators without a joint crossing through 1 MiB stay censored, with no invented percentage saving. Compact ordinary queries cost more than interpolation; batch sizes differ across families, so no cross-family speed ratio is claimed.'
    report=out/'report.md';text=base_document(report.read_text(encoding='utf-8'));text+='\n\n<!-- v2-publication -->\n'+claim+'\n\n'+limits+'\n\n[Complete compact metrics, feature recalls, dimensions, displacement and section accounting](compact-assessments.json) | [Current commands](../commands.md)\n\n![Tail and audit errors with five-fit ranges](tails.png)\n\n![Representative complete-ray failure](failed-ray.png)\n'
    for g in f['geometries']:
        if g['qualified_levels']:text+=f"\n![{g['id']} original and Q24 extraction](surface-{g['id']}.png)\n"
    batches=[];table=['','## Matched compact query variation','','Main sweep: microseconds/query across nine serial batches of the same 4,096 uniform sites per field. The range is repeated-batch variation on this host, not independent shape variation or an idle-machine confidence interval.','', '| Mesh | Field | Median | Minimum | Maximum |','|---|---|---:|---:|---:|']
    for g in f['geometries']:
        for m in ('compact_f32','compact_Q24'):
            for phase,folder in (('main_sweep','jobs'),('source_reproduction','reproduction-amended3')):
                source=d/folder/f"{g['id']}_0_{m}"/'timing.csv';rr=rows(source)
                assert len(rr)==9 and all(int(x['count'])==4096 for x in rr)
                us=[float(x['ns'])/1000 for x in rr]
                batches.extend(dict(phase=phase,geometry=g['id'],method=m,repeat=x['repeat'],count=x['count'],microseconds_per_query=float(x['ns'])/1000,source_sha256=sha(source)) for x in rr)
                if phase=='main_sweep':table.append(f"| {g['id']} | {m.removeprefix('compact_')} | {np.median(us):.2f} | {min(us):.2f} | {max(us):.2f} |")
    csv_write(out/'compact_query_batches.csv',batches)
    text+='\n'.join(table)+'\n\n[All 144 main-sweep and source-reproduction matched compact batches](compact_query_batches.csv) are labelled by phase and source-file hash. Source reproduction retains a separate fresh-process timing observation; deterministic equivalence does not imply equal host timings. Adaptive/grid batches use 100,000 sites, so these tables support no cross-family speed ratio.\n'
    text+='\n[Evidence completion and final QA](../COMPLETION-2026-10-12.md) | [Rebuild provenance and preserved test-input exception](REBUILD-PROVENANCE.json) | [Hashed unchanged v1 test-fixture supplement](REBUILD-TEST-FIXTURES.json)\n'
    report.write_text(text,encoding='utf-8')
    guide=ROOT/'docs/generalization-v2/README.md';old=base_document(guide.read_text(encoding='utf-8'));old=old.replace('No new held-out representation has been fitted.', 'The completed [fresh compact mesh study](final/report.md) is a separate representation-family comparison.').replace('The planned construction screen is held for a representation review; gates and fresh-data separation remain unchanged.', 'That checkpoint superseded the unrun deeper weighted-octree screen. The [compact mesh plan](../generalization-v2-compact-plan.md) selected Q24 using only old development meshes; the four fresh inputs and all final settings were frozen before held-out fields were built.')
    old+='\n\n<!-- v2-publication -->\n'+claim+'\n\n[Full final matrix](final/matrix.csv), [report and figures](final/report.md), [strict runner/verifier/rebuild/reproduction/report commands](commands.md), and [final provenance manifest](final/SHA256SUMS.json). Large binary evidence remains local under `runs/generalization-v2/fresh-v2`.\n';guide.write_text(old,encoding='utf-8')
    interpretation=ROOT/'docs/generalization-v2-current-interpretation.md';old=base_document(interpretation.read_text(encoding='utf-8'));old=old.replace('it has not yet fitted fresh held-out representations.', 'its subsequent compact branch is now evaluated on fresh held-out inputs.');old+='\n\n<!-- v2-publication -->\n## Completed compact representation comparison\n\n'+claim+'\n\n'+limits+'\n\nSee the [final report](generalization-v2/final/report.md), [full matrix](generalization-v2/final/matrix.csv) and [current evidence commands](generalization-v2/commands.md). Historical sphere-gap gains, tail regressions and the separately versioned v1 ray reassessment remain intact.\n';interpretation.write_text(old,encoding='utf-8')
    readme=ROOT/'README.md';old=readme.read_text(encoding='utf-8');old=old.replace('[v2 mesh intervention](docs/generalization-v2-protocol.md) is ongoing.', '[completed v2 compact mesh study](docs/generalization-v2/final/report.md) compares Q24 with strong compact controls on four fresh meshes (three reference qualified). See its [full matrix](docs/generalization-v2/final/matrix.csv) and [current commands](docs/generalization-v2/commands.md).');readme.write_text(old,encoding='utf-8')
    fidelity=ROOT/'docs/fidelity-study.md';old=base_document(fidelity.read_text(encoding='utf-8'));old+='\n\n<!-- v2-publication -->\n## Subsequent explicit mesh representation comparison\n\nThe [completed compact mesh study](generalization-v2/final/report.md) tests four fresh meshes, of which three have qualified references. '+claim+' '+limits+' The [current interpretation](generalization-v2-current-interpretation.md) preserves this earlier round and the corrected v1 ray assessment separately.\n';fidelity.write_text(old,encoding='utf-8')
    for p in out.rglob('*'):
        if p.is_file() and p.suffix in ('.md','.json','.csv','.svg','.txt'):p.read_bytes().decode('utf-8')
    json_write(out/'PUBLICATION-PROVENANCE.json',dict(helper_sha256=sha(Path(__file__)),frozen_sha256=sha(d/'FROZEN.json'),verified_manifest_sha256=sha(d/'verification/SHA256SUMS.json'),reproduction_manifest_sha256=sha(d/'reproduction-amended3/SHA256SUMS.json'),reproduction_seal_sha256=sha(d/'REPRODUCED-AMENDED.json'),completed_reproduction=completed,reporting_amendment_sha256=sha(d/'REPORTING-AMENDMENT-10.json'),report_sha256=sha(report),commands_sha256=sha(ROOT/'docs/generalization-v2/commands.md'),documentation_sha256={str(p.relative_to(ROOT)):sha(p) for p in (guide,interpretation,readme,fidelity)},utf8_strict_decode_pass=True))
    json_write(out/'SHA256SUMS.json',manifest(out))
if __name__=='__main__':main()
