"""Preserve original screen gates and independently audit all saved development fields."""
import argparse,shutil,time
from generalization_common import *
from generalization_verify import verify_job
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);a=ap.parse_args();directory=a.directory.resolve();declaration=json.loads((directory/'SCREEN_FROZEN.json').read_text());out=directory/'verification';out.mkdir(exist_ok=True);snapshot=directory/'source';snapshot.mkdir(exist_ok=True)
    # Recover the exact executed source's sole subsequently-added output block.
    path=snapshot/'generalization_benchmark.cpp';expected=next(h for p,h in declaration['source_at_declaration_sha256'].items() if p.replace('\\','/')=='src/generalization_benchmark.cpp');assert sha(path)==expected,'Executed source snapshot changed';check_hashes(snapshot,json.loads((snapshot/'SHA256SUMS.json').read_text()));all_rows=[];total=0
    for policy in range(12):
        base=directory/'policies'/f'p{policy:02d}'
        for ref in sorted((directory/'reference').iterdir()):
            g=json.loads((ref/'metadata.json').read_text())
            for seed in (40,41):
                name=f"p{policy:02d}_{g['id']}_{seed}";cache=out/(name+'.json')
                if cache.exists():result=json.loads(cache.read_text());rr,q=result['rows'],result['queries']
                else:
                    start=time.time();rr,q=verify_job(base,g,seed,'candidate');json_write(cache,dict(rows=rr,queries=q));print(name,'independently verified',round(time.time()-start,1),'s',flush=True)
                # New pre-final clarification is recorded separately, preserving
                # the original native screen results and classifications.
                metrics=rows(base/'jobs'/f"{g['id']}_{seed}_candidate"/'metrics.csv')
                for r in rr:
                    r['policy']=f'p{policy:02d}';m=next(m for m in metrics if int(m['cap'])==r['cap'] and m['dataset']=='probes');r['probe_max']=float(m['max_error']);r['probe_max_pass']=r['probe_max']<=.020;r['all_known_query_max_pass']=r['scalar_pass'] and r['probe_max_pass'];r['strict_probe_clarified_pass']=r['pass'] and r['probe_max_pass']
                all_rows.extend(rr);total+=q
    assert len(all_rows)==1320;csv_write(out/'matrix.csv',all_rows);json_write(out/'VERIFIED.json',dict(fields=len(all_rows),decoded_queries=total,original_screen_classifications_preserved=True,probe_max_clarification_separate=True,declaration_sha256=sha(directory/'SCREEN_FROZEN.json'),verifier_sha256=sha(Path(__file__))));json_write(out/'SHA256SUMS.json',manifest(out));print('All1320developmentfields independently verified',flush=True)
if __name__=='__main__':main()
