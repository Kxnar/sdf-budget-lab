"""Report only frozen verified rows; no interpolation across tested budgets."""
import argparse
from generalization_common import *
from generalization_v2_run import validate_freeze,FIT_SEEDS,COMPACT

def crossing(matrix,geometry,method,seed):
    passing=[r for r in matrix if r['geometry']==geometry and r['method']==method and int(r['seed'])==seed and str(r['pass'])=='True']
    return min(passing,key=lambda r:(int(r['bytes']),int(r['cap']))) if passing else None

def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);a=ap.parse_args();d=a.directory.resolve();f=validate_freeze(d,False);assert (d/'verification/VERIFIED.json').exists();matrix=rows(d/'verification/matrix.csv');out=ROOT/'docs/generalization-v2/final';out.mkdir(parents=True,exist_ok=True)
    table=[];paired=[]
    for g in f['geometries']:
        for m in f['methods']:
            seeds=[0] if m=='grid' or m in COMPACT else FIT_SEEDS;cc=[crossing(matrix,g['id'],m,s) for s in seeds];sizes=[int(r['bytes']) for r in cc if r]
            table.append(dict(geometry=g['id'],method=m,reference_qualified=len(g['qualified_levels'])==2,passing_fits=len(sizes),fits=len(seeds),minimum_tested_actual_bytes=min(sizes) if sizes else '',median_tested_actual_bytes=float(np.median(sizes)) if sizes else '',maximum_tested_actual_bytes=max(sizes) if sizes else '',status='observed_crossing' if sizes else 'reference_unresolved' if not g['qualified_levels'] else 'no_observed_crossing_through_1MiB'))
        selected=crossing(matrix,g['id'],'compact_Q24',0);controls=[crossing(matrix,g['id'],m,0) for m in ('compact_f64','compact_f32')];controls=[r for r in controls if r]
        for seed in FIT_SEEDS:
            q8=[crossing(matrix,g['id'],m,seed) for m in METHODS if m!='grid'];q8=[r for r in q8 if r];grid=crossing(matrix,g['id'],'grid',0)
            for label,comparison in [('compact_envelope',min(controls,key=lambda r:int(r['bytes'])) if controls else None),('q8_envelope',min(q8,key=lambda r:int(r['bytes'])) if q8 else None),('grid',grid)]:
                saving=100*(1-int(selected['bytes'])/int(comparison['bytes'])) if selected and comparison else None
                paired.append(dict(geometry=g['id'],seed=seed,comparator=label,candidate_bytes=int(selected['bytes']) if selected else '',comparator_bytes=int(comparison['bytes']) if comparison else '',comparator_method=comparison['method'] if comparison else '',saving_percent=saving if saving is not None else '',status='finite_paired_comparison' if saving is not None else 'candidate_no_joint_crossing' if not selected else 'comparator_censored_no_finite_saving'))
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
    qualified=sum(len(g['qualified_levels'])==2 for g in f['geometries']);passing=sum(crossing(matrix,g['id'],'compact_Q24',0) is not None for g in f['geometries']);finite=[r for r in paired if r['comparator']=='compact_envelope' and r['seed']==340 and r['status']=='finite_paired_comparison']
    lines=[f'# Frozen held-out compact mesh study',f'Q24 has a joint crossing on {passing} of {qualified} reference-qualified held-out meshes; the declared suite contains {len(f["geometries"])} meshes. Reference-unresolved inputs remain unknown in the matrix.', 'The candidate stores shared quantized vertices, every triangle, codec metadata and its own conservative BVH. It queries only decoded geometry. This is a representation-family comparison, not improved q8 interpolation.',f'{len(matrix):,} logical cap rows and156 physical jobs; twelve deterministic compact fields,52 grid fields and1,820 adaptive cap fields. Compact fields have an intrinsic size and are reused at every feasible cap; fitting variation applies only to the adaptive baselines.', 'All frozen scalar, feature, signed-ray dimension and two-resolution extraction gates remain mandatory. Sampled maxima are not continuous worst-case certification.', '', '| Held-out mesh | Reference levels | Q24 bytes at a tested passing cap | Strongest passing compact control | Saving |','|---|---|---:|---:|---:|']
    for g in f['geometries']:
        c=crossing(matrix,g['id'],'compact_Q24',0);p=next((r for r in finite if r['geometry']==g['id']),None)
        lines.append(f"| {g['id']} | {','.join(map(str,g['qualified_levels'])) or 'unresolved'} | {c['bytes'] if c else 'no joint crossing'} | {p['comparator_bytes'] if p else 'no observed crossing'} | {float(p['saving_percent']):.2f}% |" if p else f"| {g['id']} | {','.join(map(str,g['qualified_levels'])) or 'unresolved'} | {c['bytes'] if c else 'no joint crossing'} | no finite comparison | — |")
    lines+=['', 'The optimistic per-case q8 and compact envelopes are separate comparisons, not deployable selection policies. No finite or percentage saving is computed when either field has no observed joint crossing. The minimum is the smallest tested actual retained bytes; accuracy need not be monotonic across caps.', '', 'Retained storage includes the88-byte owning query object and complete serialized buffer, including padding. Original meshes, CGAL oracles, temporary construction buffers, allocator overhead and query scratch are excluded; construction and total benchmark peak memory are recorded separately.', 'Nine serial warm compact batches use4,096 uniform sites; adaptive/grid native batches use100,000. No cross-family timing ratio is asserted from these different scopes. Source loading, original validation, encoding, decoded validation, extraction, oracle labels and ordinary queries are recorded separately. Host background load is uncontrolled.', '', 'Each accepted codec passes decoded orientation, cyclic vertex links, nondegeneracy, domain containment, component topology and EPECK self-intersection validation. Independent verification decodes every byte, checks exact f64/XOR geometry identity, lossy displacement, BVH coverage/bounds, all scalar metrics and every ray event; stratified brute-force decoded-triangle queries are a separate subset, not a second full200,000-query oracle.', '', '[Full crossings](crossings.csv), [paired comparisons](paired_savings.csv), and [frozen provenance](FROZEN.json). The local binary evidence is retained under `runs/generalization-v2/fresh-v2`.', '', 'Input attribution: Fish—Keenan Crane, public domain; Mushroom—Josh Holinaty, CC BY-NC4.0; Springer—Tim Edwards, GPLv2; Hammer—Andrea Spognetta, CC BY4.0. Source revision835213277571d625722d25880c9c3010c6f45a5a. Exact individual README text, URLs, hashes and transforms are frozen.']
    (out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    import shutil
    for name in ('FROZEN.json','FINAL_POLICY.txt','PREPARATION-DECLARATION.json'):shutil.copy2(d/name,out/name)
    shutil.copy2(d/'verification/VERIFIED.json',out/'VERIFIED.json')
    json_write(out/'SHA256SUMS.json',manifest(out))
if __name__=='__main__':main()
