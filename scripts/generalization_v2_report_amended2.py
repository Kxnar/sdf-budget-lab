"""Report only frozen verified rows; no interpolation across tested budgets."""
import argparse
from generalization_common import *
from generalization_v2_run import validate_freeze,FIT_SEEDS,COMPACT

def crossing(matrix,geometry,method,seed):
    passing=[r for r in matrix if r['geometry']==geometry and r['method']==method and int(r['seed'])==seed and str(r['pass'])=='True']
    return min(passing,key=lambda r:(int(r['bytes']),int(r['cap']))) if passing else None

def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);a=ap.parse_args();d=a.directory.resolve();f=validate_freeze(d,False);verification=json.loads((d/'verification/VERIFIED.json').read_text());assert verification['frozen_sha256']==sha(d/'FROZEN.json');check_hashes(d/'verification',json.loads((d/'verification/SHA256SUMS.json').read_text()));from generalization_v2_verify import evaluator_identity;assert verification['evaluator']==evaluator_identity();matrix=rows(d/'verification/matrix.csv');out=ROOT/'docs/generalization-v2/final';out.mkdir(parents=True,exist_ok=True)
    table=[];paired=[]
    for g in f['geometries']:
        for m in f['methods']:
            seeds=[0] if m=='grid' or m in COMPACT else FIT_SEEDS;cc=[crossing(matrix,g['id'],m,s) for s in seeds];sizes=[int(r['bytes']) for r in cc if r]
            table.append(dict(geometry=g['id'],method=m,reference_qualified=len(g['qualified_levels'])==2,passing_fits=len(sizes),fits=len(seeds),minimum_tested_actual_bytes=min(sizes) if sizes else '',median_tested_actual_bytes=float(np.median(sizes)) if sizes else '',maximum_tested_actual_bytes=max(sizes) if sizes else '',status='observed_crossing' if sizes else 'reference_unresolved' if not g['qualified_levels'] else 'no_observed_crossing_through_1MiB'))
        selected=crossing(matrix,g['id'],'compact_Q24',0);controls=[crossing(matrix,g['id'],m,0) for m in ('compact_f64','compact_f32')];controls=[r for r in controls if r]
        for seed in FIT_SEEDS:
            q8=[crossing(matrix,g['id'],m,seed) for m in METHODS if m!='grid'];q8=[r for r in q8 if r];grid=crossing(matrix,g['id'],'grid',0)
            for label,comparison in [('compact_envelope',min(controls,key=lambda r:int(r['bytes'])) if controls else None),('q8_envelope',min(q8,key=lambda r:int(r['bytes'])) if q8 else None),('grid',grid)]:
                if label!='q8_envelope' and seed!=340:continue
                saving=100*(1-int(selected['bytes'])/int(comparison['bytes'])) if selected and comparison else None
                paired.append(dict(geometry=g['id'],seed=seed if label=='q8_envelope' else 0,comparator=label,candidate_bytes=int(selected['bytes']) if selected else '',comparator_bytes=int(comparison['bytes']) if comparison else '',comparator_method=comparison['method'] if comparison else '',saving_percent=saving if saving is not None else '',status='finite_paired_comparison' if saving is not None else 'candidate_no_joint_crossing' if not selected else 'comparator_censored_no_finite_saving'))
    csv_write(out/'crossings.csv',table);csv_write(out/'paired_savings.csv',paired)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    labels=[g['id'].split('_',1)[1] for g in f['geometries']]
    fig,ax=plt.subplots(figsize=(9,4));x=np.arange(4);width=.25
    for j,m in enumerate(COMPACT):
        values=[]
        for g in f['geometries']:
            cache=json.loads((d/'verification'/f"{g['id']}_0_{m}.json").read_text());size=cache['assessment']['actual_bytes'];values.append(size/1024 if size else np.nan)
        ax.bar(x+(j-1)*width,values,width,label=m.replace('compact_',''))
    ax.set_xticks(x,labels);ax.set_ylabel('Actual retained field KiB');ax.legend();ax.set_title('Intrinsic compact storage; acceptance assessed separately');fig.tight_layout();fig.savefig(out/'storage.png',dpi=180);fig.savefig(out/'storage.svg');plt.close(fig)
    fig,ax=plt.subplots(figsize=(12,4));heat=[]
    for g in f['geometries']:
        heat.append([next(t['passing_fits']/t['fits'] for t in table if t['geometry']==g['id'] and t['method']==m) if g['qualified_levels'] else np.nan for m in f['methods']])
    ax.imshow(np.ma.masked_invalid(heat),vmin=0,vmax=1,cmap='RdYlGn',aspect='auto');ax.set_xticks(range(11),[m.replace('compact_','') for m in f['methods']],rotation=40,ha='right');ax.set_yticks(range(4),labels)
    for i,g in enumerate(f['geometries']):
        for j,m in enumerate(f['methods']):
            t=next(t for t in table if t['geometry']==g['id'] and t['method']==m);ax.text(j,i,f"{t['passing_fits']}/{t['fits']}" if g['qualified_levels'] else '?',ha='center',va='center',fontsize=9)
    ax.set_title('Fits with a joint crossing at some tested cap; ? reference unresolved');fig.tight_layout();fig.savefig(out/'joint_gates.png',dpi=180);fig.savefig(out/'joint_gates.svg');plt.close(fig)
    qualified=sum(len(g['qualified_levels'])==2 for g in f['geometries']);passing=sum(crossing(matrix,g['id'],'compact_Q24',0) is not None for g in f['geometries']);finite=[r for r in paired if r['comparator']=='compact_envelope' and r['seed']==0 and r['status']=='finite_paired_comparison']
    lines=[f'# Frozen held-out compact mesh study',f'Q24 has a joint crossing on {passing} of {qualified} reference-qualified held-out meshes; the declared suite contains {len(f["geometries"])} meshes. Reference-unresolved inputs remain unknown in the matrix.', 'The candidate stores shared quantized vertices, every triangle, codec metadata and its own conservative BVH. It queries only decoded geometry. This is a representation-family comparison, not improved q8 interpolation.',f'{len(matrix):,} logical cap rows and 156 physical jobs; 1,884 physical encoded-field records: twelve deterministic compact records, 52 grid records and 1,820 adaptive cap records. Different cap records may have identical binary content. Compact fields have an intrinsic size and are reused at every feasible cap; fitting variation applies only to the adaptive baselines.', 'All frozen scalar, feature, signed-ray dimension and two-resolution extraction gates remain mandatory. Sampled maxima are not continuous worst-case certification.', '', '| Held-out mesh | Reference levels | Q24 bytes at a tested passing cap | Strongest passing compact control | Saving |','|---|---|---:|---:|---:|']
    for g in f['geometries']:
        c=crossing(matrix,g['id'],'compact_Q24',0);p=next((r for r in finite if r['geometry']==g['id']),None)
        lines.append(f"| {g['id']} | {','.join(map(str,g['qualified_levels'])) or 'unresolved'} | {c['bytes'] if c else 'no joint crossing'} | {p['comparator_bytes'] if p else 'no observed crossing'} | {float(p['saving_percent']):.2f}% |" if p else f"| {g['id']} | {','.join(map(str,g['qualified_levels'])) or 'unresolved'} | {c['bytes'] if c else 'no joint crossing'} | no finite comparison |  -  |")
    lines+=['', 'The optimistic per-case q8 and compact envelopes are separate comparisons, not deployable selection policies. No finite or percentage saving is computed when either field has no observed joint crossing. The minimum is the smallest tested actual retained bytes; accuracy need not be monotonic across caps.', '', 'Retained storage includes the 88-byte owning query object and complete serialized buffer, including padding. Original meshes, CGAL oracles, temporary construction buffers, allocator overhead and query scratch are excluded; construction and total benchmark peak memory are recorded separately.', 'Nine serial warm compact batches use 4,096 uniform sites; adaptive/grid native batches use 100,000. No cross-family timing ratio is asserted from these different scopes. Source loading, original validation, encoding, decoded validation, extraction, oracle labels and ordinary queries are recorded separately. Host background load is uncontrolled.', '', 'Each accepted codec passes decoded orientation, cyclic vertex links, nondegeneracy, domain containment, component topology and EPECK self-intersection validation. Independent verification decodes every byte, checks exact f64/XOR geometry identity, lossy displacement, BVH coverage/bounds, all scalar metrics and every ray event; stratified brute-force decoded-triangle queries are a separate subset, not a second full 200,000-query oracle.', '', '[Full crossings](crossings.csv), [paired comparisons](paired_savings.csv), and [frozen provenance](FROZEN.json). The local binary evidence is retained under `runs/generalization-v2/fresh-v2`.', '', 'Input attribution: Fish - Keenan Crane, public domain; Mushroom - Josh Holinaty, CC BY-NC 4.0; Springer - Tim Edwards, GPLv2; Hammer - Andrea Spognetta, CC BY 4.0. Source revision 835213277571d625722d25880c9c3010c6f45a5a. Exact individual README text, URLs, hashes and transforms are frozen.']
    finite_values=[float(r['saving_percent']) for r in finite]
    lines+=['',f'Across {len(finite_values)} independent qualified meshes with finite compact comparisons, median retained-storage saving is {np.median(finite_values):.2f}% (range {min(finite_values):.2f}–{max(finite_values):.2f}%).' if finite_values else 'No finite compact comparison is observed.', 'These named held-out meshes are not a population estimate. Fish is reference-unresolved and excluded from joint-crossing denominators. Compact/grid comparisons are emitted once; q8 paired rows vary across five independent fitting seeds.', '', 'The preserved diagnosis explains the representation change: depth 12 unquantized trilinear interpolation loses a short Brucewick interval. Q24 explicitly retains triangles, requires a closed triangle-mesh input and carries a measured query-cost tradeoff. It is not a q8 allocation/interpolation breakthrough. Earlier sphere-gap and tail-regression results remain separate.']
    failures=[];costs=[];peaks=[];distinct=set()
    for g in f['geometries']:
        for m in f['methods']:
            rr=[r for r in matrix if r['geometry']==g['id'] and r['method']==m]
            failures.append(dict(geometry=g['id'],method=m,logical_rows=len(rr),passing_rows=sum(r['pass']=='True' for r in rr),scalar_failed=sum(r['scalar_pass']=='False' for r in rr),feature_failed=sum(r['probe_pass']=='False' for r in rr),width_failed=sum(r['width_pass']=='False' for r in rr),observed_topology_failed=sum(r.get('topology_observed_failure')=='True' for r in rr),storage_infeasible=sum(r['status']=='storage_infeasible' for r in rr),numerical_unresolved=sum(r['status']=='numerical_unresolved' for r in rr),reference_unresolved=sum(r['reference_qualified']=='False' for r in rr)))
            seeds=[0] if m=='grid' or m in COMPACT else FIT_SEEDS
            for seed in seeds:
                job=d/'jobs'/f"{g['id']}_{seed}_{m}";proc=json.loads((job/'PROCESS.json').read_text());timing=rows(job/'timing.csv') if (job/'timing.csv').exists() else []
                if m in COMPACT:
                    c=json.loads((job/'construction.json').read_text()) if (job/'construction.json').exists() else {};ns=[float(r['ns']) for r in timing];build=sum(c.get(k,0) for k in ('source_seconds','encoding_bvh_packing_seconds','decoded_validation_seconds'));peak=c.get('validated_construction_cumulative_peak_working_set_bytes');calls=c.get('distance_triangles','')
                    if (job/'field.bin').exists():distinct.add(sha(job/'field.bin'))
                    scope='4096 matched uniform sites, nine batches'
                else:
                    c={};ns=[float(r['ns_per_query']) for r in timing];ff=rows(job/'fields.csv');build=max(float(r['build_tree_seconds'])+float(r['loss_seconds'])+float(r['dp_seconds']) for r in ff)+sum(float(r['pack_seconds']) for r in ff);peak=None;calls=max(int(r['oracle_calls']) for r in ff);scope='100000 uniform sites, nine batches per cap; maximum frontier reused'
                    for cap in CAPS:distinct.add(sha(job/f'{cap}.bin'))
                costs.append(dict(geometry=g['id'],method=m,seed=seed,native_build_and_sweep_pack_seconds=build,query_median_ns=float(np.median(ns)) if ns else '',query_maximum_batch_ns=max(ns) if ns else '',query_scope=scope,logical_build_or_query_triangle_work=calls,owned_process_elapsed_seconds=proc['seconds']))
                peaks.append(dict(geometry=g['id'],method=m,seed=seed,construction_cumulative_peak_working_bytes=peak if peak else 'unavailable',total_native_peak_working_bytes=c.get('total_benchmark_peak_working_set_bytes','unavailable'),observed_total_process_peak_committed_or_working_bytes=proc['observed_peak_committed_or_working_bytes']))
    csv_write(out/'failure_gates.csv',failures);csv_write(out/'costs.csv',costs);csv_write(out/'peak_memory.csv',peaks)
    lines+=['',f'The 1,884 encoded-field records contain {len(distinct):,} distinct binary content hashes. Verification covers 156 physical jobs and 2,028 logical cap rows; independently decoded baseline predictions/objective queries total {verification["baseline_decoded_queries"]:,}.', '[Failure-gate counts](failure_gates.csv) distinguish skipped topology, observed failure, numerical uncertainty and unresolved reference. [Phase/query costs](costs.csv) and [peak memory](peak_memory.csv) keep construction and total process scopes separate; baseline construction-only peaks are unavailable.', '', 'Verification/reproduction identities are bound to the original frozen SHA and full verified manifests. Reproduction uses separately rebuilt frozen source with study verification stopped, retaining exact fields, predictions, flags, complete ray events/witnesses and surfaces or decoded-rejection outcomes.']
    lines+=['', '| Mesh | Field | Build/validate seconds | Warm query median microseconds | Construction cumulative peak MiB |','|---|---|---:|---:|---:|']
    for g in f['geometries']:
        for m in ('compact_f32','compact_Q24'):
            c=next(r for r in costs if r['geometry']==g['id'] and r['method']==m);p=next(r for r in peaks if r['geometry']==g['id'] and r['method']==m);peak=p['construction_cumulative_peak_working_bytes'];peak_text=f'{peak/1024**2:.2f}' if isinstance(peak,int) else str(peak)
            lines.append(f"| {g['id']} | {m.replace('compact_','')} | {c['native_build_and_sweep_pack_seconds']:.3f} | {float(c['query_median_ns'])/1000:.2f} | {peak_text} |")
    lines+=['', 'Both compact query columns use identical 4096-site, nine-batch scopes; construction includes source load, encoding/BVH/packing and mandatory decoded validation, excluding the evaluation oracle/extraction. Peaks are cumulative process high-water marks through decoded validation.', '', '![Retained storage](storage.png)', '![Joint crossings](joint_gates.png)']
    if (d/'REPRODUCED-AMENDED.json').exists():
        reproduction=json.loads((d/'REPRODUCED-AMENDED.json').read_text());assert reproduction['frozen_sha256']==sha(d/'FROZEN.json');lines.append(f'Frozen-source reproduction completed {reproduction["jobs"]} fresh-process jobs, including every contributing fitting-seed crossing and adjacent smaller tested cap or censored endpoint.')
    (out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    import shutil
    for name in ('FROZEN.json','FINAL_POLICY.txt','PREPARATION-DECLARATION.json'):shutil.copy2(d/name,out/name)
    shutil.copy2(d/'verification/VERIFIED.json',out/'VERIFIED.json')
    for name in ('PROCEDURAL-AMENDMENT-1.json','PROCEDURAL-AMENDMENT-2.json','FIELD-COUNT-CLARIFICATION.json','HARDWARE.json','REPRODUCED-AMENDED.json'):
        if (d/name).exists():shutil.copy2(d/name,out/name)
    json_write(out/'REPORT-PROVENANCE.json',dict(report_helper_sha256=sha(Path(__file__)),original_frozen_sha256=sha(d/'FROZEN.json'),verification_manifest_sha256=sha(d/'verification/SHA256SUMS.json'),frozen_source_hashes=f['source_sha256'],oracle_audit_summaries={g['id']:json.loads((d/'oracle-audit'/g['id']/'RESULT.json').read_text()) for g in f['geometries']}))
    json_write(out/'SHA256SUMS.json',manifest(out))
if __name__=='__main__':main()
