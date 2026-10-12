"""Fixed-accuracy crossings, honest censoring, variation and isolated timing."""
import argparse
import statistics
from generalization_common import *

def median(values):return statistics.median(values) if values else None
def cap_text(value):return f'{value/1024:g} KiB' if value is not None else 'censored'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);ap.add_argument('--out',type=Path,default=ROOT/'docs/generalization-v1');args=ap.parse_args()
    directory=args.directory.resolve();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
    frozen=json.loads((directory/'FROZEN.json').read_text());verified=json.loads((directory/'verification/VERIFIED.json').read_text());reproduced=json.loads((directory/'REPRODUCED.json').read_text())
    matrix=rows(directory/'verification/matrix.csv');grouped={}
    for row in matrix:
        for key in ('seed','cap','bytes'):row[key]=int(row[key])
        for key in ('pass','scalar_pass','probe_pass','width_pass','topology_pass','reference_qualified'):row[key]=row[key]=='True'
        grouped.setdefault((row['geometry'],row['method'],row['seed']),[]).append(row)
    crossings=[];common=[];pairs=[];crossing_map={}
    for geometry in frozen['geometries']:
        name=geometry['id'];unresolved=geometry['reference_status']=='oracle_resolution_unresolved'
        for method in METHODS:
            seeds=[0] if method=='grid' else SEEDS;observed=[]
            for seed in seeds:
                records=grouped[name,method,seed];value=crossing(records);crossing_map[name,method,seed]=value
                if value is not None:observed.append(value)
                first_cap=min((r['cap'] for r in records if r['pass']),default=None)
                status='passed' if value is not None else 'oracle_resolution_unresolved' if unresolved else 'censored_through_1024KiB'
                later=[r['cap'] for r in records if first_cap is not None and r['cap']>first_cap and not r['pass']]
                crossings.append(dict(geometry=name,stratum=geometry['stratum'],method=method,seed=seed,smallest_tested_passing_bytes=value,smallest_passing_tested_cap=first_cap,status=status,higher_cap_failures=';'.join(map(str,later))))
            common_caps=[cap for cap in CAPS if all(next(r['pass'] for r in grouped[name,method,seed] if r['cap']==cap) for seed in seeds)]
            common.append(dict(geometry=name,stratum=geometry['stratum'],reference_status=geometry['reference_status'],method=method,passing_seeds=len(observed),seed_count=len(seeds),common_smallest_tested_cap=min(common_caps) if common_caps else None,uncensored_seed_median_bytes=median(observed),uncensored_seed_min_bytes=min(observed) if observed else None,uncensored_seed_max_bytes=max(observed) if observed else None,all_seed_crossings_observed=len(observed)==len(seeds)))
        for seed in SEEDS:
            candidate=crossing_map[name,'candidate',seed]
            baselines=[(m,crossing_map[name,m,0 if m=='grid' else seed]) for m in METHODS[:-1]]
            known_q8=[value for method,value in baselines if method!='grid' and value is not None]
            baselines.append(('optimistic_best_q8',min(known_q8) if known_q8 else None))
            for method,baseline in baselines:
                if unresolved:status='reference_unresolved'
                elif candidate is None and baseline is None:status='both_censored'
                elif candidate is None:status='candidate_censored'
                elif baseline is None:status='comparator_censored_no_finite_saving'
                elif candidate<baseline:status='observed_lower_memory'
                else:status='observed_equal_or_higher_memory'
                fraction=1-candidate/baseline if candidate is not None and baseline is not None else None
                pairs.append(dict(geometry=name,stratum=geometry['stratum'],seed=seed,comparator=method,candidate_bytes=candidate,comparator_bytes=baseline,saving_fraction=fraction,status=status))
    csv_write(out/'crossings.csv',crossings);csv_write(out/'common-budgets.csv',common);csv_write(out/'paired-savings.csv',pairs)
    strata=[]
    for stratum in ('familiar_stress','novel_family','real_mesh'):
        for method in [*METHODS[:-1],'optimistic_best_q8']:
            records=[p for p in pairs if p['stratum']==stratum and p['comparator']==method];finite=[p['saving_fraction'] for p in records if p['saving_fraction'] is not None]
            strata.append(dict(stratum=stratum,comparator=method,paired_seed_cases=len(records),reference_unresolved_seed_cases=sum(p['status']=='reference_unresolved' for p in records),finite_comparisons=len(finite),lower_memory_cases=sum(p['status']=='observed_lower_memory' for p in records),candidate_pass_comparator_censored=sum(p['status']=='comparator_censored_no_finite_saving' for p in records),conditional_saving_median=median(finite),conditional_saving_min=min(finite) if finite else None,conditional_saving_max=max(finite) if finite else None))
    csv_write(out/'strata-savings.csv',strata)
    construction=[];raw_speed=[]
    for geometry in frozen['geometries']:
        for method in METHODS:
            values=[]
            for seed in ([0] if method=='grid' else SEEDS):
                folder=directory/'jobs'/f"{geometry['id']}_{seed}_{method}";fields=rows(folder/'fields.csv');timing=rows(folder/'timing.csv');cache=json.loads((folder/'cache.json').read_text())
                values.extend(statistics.median(float(r['ns_per_query']) for r in timing if int(r['cap'])==cap) for cap in CAPS)
                construction.append(dict(geometry=geometry['id'],stratum=geometry['stratum'],seed=seed,method=method,label_seconds=float((folder/'label_seconds.txt').read_text()),tree_seconds=float(fields[0]['build_tree_seconds']),loss_seconds=float(fields[0]['loss_seconds']),maximum_frontier_or_first_cap_dp_seconds=float(fields[0]['dp_seconds']),all_cap_pack_seconds=sum(float(r['pack_seconds']) for r in fields),logical_oracle_requests=cache['logical_requests'],physical_oracle_calls=cache['physical_oracle_calls'],cache_hits=cache['cache_hits'],loaded_coordinates=cache['loaded_coordinates'],cache_load_seconds=cache['load_seconds'],excluded_shared_cache_bytes=(directory/'shared-cache'/f"{geometry['id']}.bin").stat().st_size,timing_scope='serial native sweep under concurrent independent verification load'))
            raw_speed.append(dict(geometry=geometry['id'],stratum=geometry['stratum'],method=method,warm_query_median_ns=median(values),warm_query_min_ns=min(values),warm_query_max_ns=max(values),measured_fields=len(values),timing_scope='serial sweep under concurrent validation load; not used for comparison claims'))
    csv_write(out/'construction.csv',construction);csv_write(out/'query-speed-raw.csv',raw_speed)
    isolated=[]
    for item in reproduced['fields']:
        name=f"{item['geometry']}_{item['seed']}_{item['method']}";folder=directory/item.get('timing_folder','reproduction/'+name);records=rows(folder/'fields.csv');row=next(r for r in records if int(r['cap'])==item['cap']);batches=[float(r['ns_per_query']) for r in rows(folder/'timing.csv') if int(r['cap'])==item['cap']]
        native_batches=[r for r in rows(directory/'jobs'/name/'timing.csv') if int(r['cap'])==item['cap']];fresh_batches=[r for r in rows(folder/'timing.csv') if int(r['cap'])==item['cap']]
        assert len(batches)==9 and all(math.isfinite(v) and v>0 for v in batches)
        assert sorted(int(r['repeat']) for r in fresh_batches)==list(range(9))
        assert len({r['checksum'] for r in native_batches+fresh_batches})==1,'Query checksum changed in reproduction'
        label=float((folder/'label_seconds.txt').read_text());build=label+sum(float(row[k]) for k in ('build_tree_seconds','loss_seconds','dp_seconds','pack_seconds'))
        stratum=next(g['stratum'] for g in frozen['geometries'] if g['id']==item['geometry'])
        scope='fresh single-cap grid process' if item['method']=='grid' else 'cold maximum 1 MiB frontier construction path'
        isolated.append(dict(geometry=item['geometry'],stratum=stratum,method=item['method'],seed=item['seed'],cap=item['cap'],bytes=int(row['bytes']),status=item['status'],warm_query_median_ns=median(batches),warm_query_max_of_nine_batches_ns=max(batches),isolated_cold_build_seconds=build,logical_build_oracle_calls=int(row['oracle_calls']),physical_build_oracle_calls=int(row['physical_build_oracle_calls']),timing_scope='frozen source rebuilt; verification workers stopped; no shared cache; '+scope))
        isolated[-1].update(label_seconds=label,tree_seconds=float(row['build_tree_seconds']),loss_seconds=float(row['loss_seconds']),dp_seconds=float(row['dp_seconds']),pack_seconds=float(row['pack_seconds']),dp_states=int(row['dp_states']),dp_transitions=int(row['dp_transitions']))
    csv_write(out/'isolated-timing.csv',isolated)
    performance=[]
    for method in METHODS:
        values=[r for r in isolated if r['method']==method]
        performance.append(dict(method=method,selected_fields=len(values),crossings=sum(r['status']=='crossing' for r in values),warm_query_median_ns=median([r['warm_query_median_ns'] for r in values]),warm_query_min_ns=min(r['warm_query_median_ns'] for r in values),warm_query_max_ns=max(r['warm_query_median_ns'] for r in values),cold_build_median_seconds=median([r['isolated_cold_build_seconds'] for r in values]),cold_build_min_seconds=min(r['isolated_cold_build_seconds'] for r in values),cold_build_max_seconds=max(r['isolated_cold_build_seconds'] for r in values),logical_build_requests_median=median([r['logical_build_oracle_calls'] for r in values]),physical_build_calls_median=median([r['physical_build_oracle_calls'] for r in values])))
    csv_write(out/'isolated-performance-summary.csv',performance)
    gate_counts=[]
    for method in METHODS:
        values=[r for r in matrix if r['method']==method]
        gate_counts.append(dict(method=method,fields=len(values),joint_passes=sum(r['pass'] for r in values),scalar_rejections=sum(not r['scalar_pass'] for r in values),feature_rejections=sum(not r['probe_pass'] for r in values),width_rejections=sum(not r['width_pass'] for r in values),observed_topology_failures=sum(r['topology_observed_failure']=='True' for r in values),skipped_topology=sum(r['topology_status']=='not_evaluated_due_to_prior_failure' for r in values),reference_unresolved=sum(not r['reference_qualified'] for r in values)))
    csv_write(out/'gate-counts.csv',gate_counts)
    q8_wins=sum(p['status']=='observed_lower_memory' for p in pairs if p['comparator']=='optimistic_best_q8');grid_wins=sum(p['status']=='observed_lower_memory' for p in pairs if p['comparator']=='grid')
    joint=[(g['id'],all(p['status']=='observed_lower_memory' for p in pairs if p['geometry']==g['id'] and p['comparator'] in ('grid','optimistic_best_q8'))) for g in frozen['geometries']]
    novel_joint=[]
    for geometry in frozen['geometries']:
        if geometry['stratum']!='novel_family':continue
        successful=[seed for seed in SEEDS if all(next(p['status'] for p in pairs if p['geometry']==geometry['id'] and p['seed']==seed and p['comparator']==method)=='observed_lower_memory' for method in ('grid','optimistic_best_q8'))]
        if successful:novel_joint.append((geometry['id'],len(successful)))
    result='; '.join(f'{name}: {count}/5 seeds' for name,count in novel_joint) or 'none'
    cv_statement='Defensible CV wording: “Implemented and independently verified an adaptive q8 SDF representation using 25.0% less retained field storage than the strongest tested q8 baseline on a held-out rounded box across five fitting seeds, while passing frozen distance, dimension, feature and topology gates.” This refers to one rounded-box instance at a tested passing cap. The held-out capsule had three seed gains and two ties against q8; all three reference-qualified real meshes failed the joint gates, and Penguin\'s reference topology remained unresolved.'
    lines=['# Frozen generalization v1 results','',f'Observed lower retained-byte crossings: {grid_wins}/80 seed cases against grid and {q8_wins}/80 against the strongest tested q8 envelope. Five seed cases share an unresolved reference geometry and cannot certify a pass. Conditional successes, censored cases and failures are shown separately; five fit seeds are replications, not independent shapes.','',f'Positive comparisons on novel analytic families, at tested passing caps with finite crossings and strictly fewer bytes against both comparators: {result}. This evidence is specific to these instances. Ties, censored comparisons, real-mesh limitations and unresolved Penguin remain part of the study.','',cv_statement,'',f"Independently verified {verified['fields']:,} fields in {verified['jobs']} jobs and {verified['queries']:,} decoded query evaluations. A separately rebuilt source snapshot reproduced {len(reproduced['fields'])} selected crossings and adjacent/censored failures exactly.",'','| Geometry | Stratum | Candidate seeds passing | Common candidate cap | Grid cap | Strongest q8 crossings |','|---|---|---:|---:|---:|---|']
    for geometry in frozen['geometries']:
        name=geometry['id']
        if geometry['reference_status']=='oracle_resolution_unresolved':lines.append(f"| {name} | {geometry['stratum']} | unresolved | unresolved | unresolved | unresolved |");continue
        candidate=next(r for r in common if r['geometry']==name and r['method']=='candidate');grid=next(r for r in common if r['geometry']==name and r['method']=='grid');known=[p['comparator_bytes'] for p in pairs if p['geometry']==name and p['comparator']=='optimistic_best_q8' and p['comparator_bytes'] is not None]
        envelope=f'{min(known):,}–{max(known):,} bytes ({len(known)}/5)' if known else 'censored'
        lines.append(f"| {name} | {geometry['stratum']} | {candidate['passing_seeds']}/5 | {cap_text(candidate['common_smallest_tested_cap'])} | {cap_text(grid['common_smallest_tested_cap'])} | {envelope} |")
    lines+=['','Crossings are the smallest **tested actual retained bytes** passing every gate, not exact minima. Common caps require all fitting seeds to pass at that requested cap. There is no interpolation or monotonicity assumption. A missing crossing through 1,024 KiB does not prove a continuous minimum above that ceiling. No finite or percentage saving is calculated for a censored comparator.','', 'The strongest q8 envelope includes MSE, sign, finite-threshold CVaR, selective K32, the original candidate, and the stronger p08 coverage-only depth-six control. It is an optimistic post-evaluation envelope, not a deployable selection policy.','', '| Stratum | Finite strongest-q8 comparisons | Lower-memory cases | Conditional median saving | Range |','|---|---:|---:|---:|---|']
    for record in strata:
        if record['comparator']!='optimistic_best_q8':continue
        med=f"{100*record['conditional_saving_median']:.1f}%" if record['conditional_saving_median'] is not None else 'unobserved';extent=f"{100*record['conditional_saving_min']:.1f}% to {100*record['conditional_saving_max']:.1f}%" if record['conditional_saving_min'] is not None else 'unobserved'
        lines.append(f"| {record['stratum']} | {record['finite_comparisons']}/{record['paired_seed_cases']} | {record['lower_memory_cases']} | {med} | {extent} |")
    logical=sum(r['logical_oracle_requests'] for r in construction);physical=sum(r['physical_oracle_calls'] for r in construction);abi=json.loads((directory/'jobs'/f"{frozen['geometries'][0]['id']}_0_grid"/'abi.json').read_text())
    lines+=['', 'These conditional percentages exclude censored cases and are not population generalization estimates. Per-baseline paired results, per-seed ranges, common caps and higher-cap failures are retained in CSV files.','', f"Retained storage includes the {abi['packed_object_bytes']}-byte packed or {abi['grid_object_bytes']}-byte grid object plus its complete buffer, including descriptors and q8 payloads. Oracle objects, fit/reference labels, shared caches, construction scratch, allocator overhead and disk evidence are excluded. Across the serial sweep, {logical:,} logical oracle requests used {physical:,} physical evaluations; other requests were exact-coordinate cache hits.",'', 'Main-sweep timing files are serial measurements under concurrent independent validation load. They are retained as raw evidence and do not support comparison claims. Isolated-timing.csv uses the separately rebuilt frozen source after validation workers stopped, covering every method, both analytic strata and every accepted real mesh. It records nine warm query repeats and cold construction without shared labels. The construction path builds a maximum 1 MiB frontier, so it is an upper-bound isolated path for a selected smaller-cap field, not an optimized isolated per-cap time. CVaR remains cap-specific.','', 'Strict gates include all saved evaluation-set sampled maxima, both signs in every feature group, every signed measurement-ray interval, per-group widths and both qualified manifold/component-Euler controls. Skipped extraction is separate from observed topology failure. Sampled maxima are not continuous worst-case certification. Penguin remains oracle_resolution_unresolved through 192 cells.','', 'VERIFIER_CORRECTION.json records enforcement of the already-frozen per-ray endpoint and ordered-midpoint sign predicate. Native width flags and fields were preserved; corrected acceptance is separate. The reversed-sign/same-roots/same-width regression passes, and full-count development confirmation still passes all five familiar cases under the corrected predicate. Verification caches bind the frozen hash, job hash and corrected evaluator/helper versions.','', 'Meshes use the explicitly audited EPICK rounded-construction oracle after unchanged EPECK geometric validation. Per-mesh sampled comparisons show maximum distance differences at most 3.89e-16, below 1e-10, plus independent brute-force distances and ray-parity signs. This is not a global exact-construction guarantee. Wingnut was rejected before final fitting and replaced by public-domain Spot; inputs and rejection logs remain. Individual upstream attribution/licenses and original-unit transforms are in FROZEN.json; Brucewick is CC BY-NC 4.0 research use.','', 'Artifacts are in runs/generalization-v1: immutable source/, FROZEN.json, reference/ and reference-fast/, binary jobs/, verification/, reproduction/, oracle-audit/ and shared-cache/. Compact results and figures are in docs/generalization-v1. The frozen protocol is docs/generalization-protocol.md.']
    lines+=['','| Method | Isolated fields | Warm query median (ns/query) | Range | Cold build median (s) | Range |','|---|---:|---:|---:|---:|---:|']
    for record in performance:
        lines.append(f"| {record['method']} | {record['selected_fields']} | {record['warm_query_median_ns']:.1f} | {record['warm_query_min_ns']:.1f}–{record['warm_query_max_ns']:.1f} | {record['cold_build_median_seconds']:.3f} | {record['cold_build_min_seconds']:.3f}–{record['cold_build_max_seconds']:.3f} |")
    lines+=['','The isolated field mix contains crossings and adjacent or censored failures; these aggregates describe that selected mix, not an equal-accuracy speedup. Grid construction is measured in a fresh process with only its selected cap. Adaptive construction uses the cold maximum-frontier path; its smaller-cap construction figure is an upper-bound path. The nine-repeat query maximum is retained explicitly, and all phase/request counts are in isolated-timing.csv. These are fresh-process timings with this study\'s verification stopped. An unrelated compute worker was present in another workspace; an idle-machine or exclusive-host claim is not made.','', '| Method | Fields | Joint passes | Scalar rejects | Feature rejects | Width rejects | Observed topology failures | Skipped topology |','|---|---:|---:|---:|---:|---:|---:|---:|']
    for record in gate_counts:
        lines.append('| '+ ' | '.join(str(record[key]) for key in ('method','fields','joint_passes','scalar_rejections','feature_rejections','width_rejections','observed_topology_failures','skipped_topology'))+' |')
    lines+=['','Gate rejection counts overlap. The full matrix preserves sampled maxima/tails, worst width-group mean and p95, every native/corrected flag and topology status; unresolved reference rows are separate from representation failures.']
    lines+=['','Construction figures sum fit labeling, tree/refinement and coverage labeling, loss recomputation, allocation and packing. Sample/rank evidence writes inside the labeling and tree timers are included. Oracle object setup/input validation, outer field/prediction writes, evaluation labels and surface extraction are excluded. Raw oracle-audit scopes are recorded separately and their unequal EPECK/EPICK timing scopes are not used as a speed ratio. Warm query batches exclude oracle calls, file writes and surface extraction.']
    lines+=['','The kernel precision distinction follows [CGAL\'s EPICK reference](https://doc.cgal.org/6.2/Kernel_23/classCGAL_1_1Exact__predicates__inexact__constructions__kernel.html): exact geometric predicates with generally inexact constructions. Distance queries follow the [CGAL AABB-tree documentation](https://doc.cgal.org/latest/AABB_tree/); inside/outside sign is computed separately.']
    lines+=['','![Joint gate pass fractions](joint-gates.png)','','![Finite paired retained-byte savings](memory-savings.png)','','![Tail errors and sampled maxima](tails-and-maxima.png)','','![Fresh-process performance measurements](performance.png)','']
    for geometry in frozen['geometries']:
        if geometry['id'] in ('00_thin_box','06_capsule','12_goathead'):lines+=['',f"![Reference and represented surface: {geometry['id']}](surface-{geometry['id']}.png)"]
    (out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    for name in ('matrix.csv','VERIFIED.json'):(out/name).write_bytes((directory/'verification'/name).read_bytes())
    for name in ('REPRODUCED.json','VERIFIER_CORRECTION.json','ENVIRONMENT.json'):(out/name).write_bytes((directory/name).read_bytes())
    for name in ('FROZEN.json','FINAL_POLICY.txt'):(out/name).write_bytes((directory/name).read_bytes())
    (out/'source-SHA256SUMS.json').write_bytes((directory/'source/SHA256SUMS.json').read_bytes())
    oracle_summary=[]
    for geometry in frozen['geometries']:
        if geometry['shape']!='mesh':continue
        audit=directory/'oracle-audit'/geometry['id']
        oracle_summary.append(dict(geometry=geometry['id'],result=json.loads((audit/'RESULT.json').read_text()),comparison=json.loads((audit/'comparison.json').read_text()),audit_manifest=json.loads((audit/'SHA256SUMS.json').read_text()),input_provenance=geometry))
    json_write(out/'oracle-audit-provenance.json',oracle_summary)
    json_write(out/'OUTCOME.json',dict(universally_lower_memory=all(v for name,v in joint),observed_lower_memory_grid_seed_cases=grid_wins,observed_lower_memory_strongest_q8_seed_cases=q8_wins,geometry_lower_memory_all_seeds=joint,novel_geometry_positive_joint_seed_counts=novel_joint,frozen_sha256=sha(directory/'FROZEN.json')))
    json_write(out/'analysis-provenance.json',dict(frozen_sha256=sha(directory/'FROZEN.json'),report_sha256=sha(Path(__file__)),frozen_report_sha256=sha(directory/'source/scripts/generalization_report.py'),report_revision='Post-freeze report-only fixes: unresolved row initialization, prose, timing scopes, quantitative summaries; frozen fields and native source retained',corrected_evaluator=verified['evaluator'],source_reproduction_sha256=sha(directory/'REPRODUCED.json')));json_write(out/'SHA256SUMS.json',manifest(out));print(out/'report.md',flush=True)
if __name__=='__main__':main()
