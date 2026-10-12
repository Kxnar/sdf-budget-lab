"""Independent diagnostic corner distances and Brucewick missing interval."""
from generalization_common import *
from generalization_v2_feasibility import virtual_query

old=ROOT/'runs/generalization-v1';base=ROOT/'runs/generalization-v2';out=base/'feasibility-oracle-checks';out.mkdir(exist_ok=True);frozen=json.loads((old/'FROZEN.json').read_text());summary=[]
for g in frozen['geometries']:
    if g['id'] not in ('12_goathead','14_spot','15_brucewick'):continue
    folder=base/'feasibility'/g['id'];samples=[]
    for depth in (7,9,12):
        cells=np.fromfile(folder/f'probe_cells_{depth}.bin',dtype='<f8').reshape(-1,16)
        for i in np.linspace(0,len(cells)-1,48,dtype=int):
            cell=cells[i];c=int(i%8);p=np.where(np.array([bool(c&(1<<a)) for a in range(3)]),cell[5:8],cell[2:5]);samples.append([*p,cell[8+c]])
    detail=None
    if g['id']=='15_brucewick':
        ref=reference_dir(old,g);ray=np.loadtxt(ref/'rays.txt')[0];truth=np.fromstring(rows(ref/'reference_roots.csv')[0]['roots'],sep=';');mid=(truth[2]+truth[3])/2;point=ray[1:4]+mid*ray[4:7];cells=np.fromfile(folder/'segments_12_0.bin',dtype='<f8').reshape(-1,16);index=np.searchsorted(cells[:,1],mid,side='right');cell=cells[index]
        for c in range(8):samples.append([*np.where(np.array([bool(c&(1<<a)) for a in range(3)]),cell[5:8],cell[2:5]),cell[8+c]])
        label=brute_mesh(g,np.array([point]))[0];predictions={e:float(virtual_query(cells,ray,e,np.array([mid]))[0]) for e in ('q8','f32','f64')};detail=dict(ray=0,reference_roots=truth.tolist(),missing_interval=truth[2:4].tolist(),midpoint_world=point.tolist(),independent_signed_distance=float(label),virtual_predictions=predictions,containing_cell=cell.tolist(),interpretation='Missing two transitions on full depth12, also in unquantized interpolation; no numerical ambiguity')
    samples=np.array(samples);pred=brute_mesh(g,samples[:,:3]);np.testing.assert_allclose(pred,samples[:,3],atol=2e-10,rtol=2e-9);mismatch=(pred*np.sign(samples[:,3])<=0)&(np.abs(samples[:,3])>2e-9);assert not mismatch.any();np.savez(out/f"{g['id']}.npz",points=samples[:,:3],oracle=samples[:,3],independent=pred)
    record=dict(geometry=g['id'],points=len(samples),maximum_difference=float(np.max(abs(pred-samples[:,3]))),sign_mismatches_away_2e9_band=int(mismatch.sum()),normalized_input_sha256=sha(Path(g['normalized_mesh'])),native_verified_sha256=sha(folder/'VERIFIED.json'),method='Independent brute-force triangle distances and ray parity; inherited reference hit deduplication convention',missing_interval_detail=detail)
    json_write(out/f"{g['id']}.json",record);summary.append(record);print(g['id'],len(samples),record['maximum_difference'],flush=True)
json_write(out/'COMPLETE.json',dict(results=summary,script_sha256=sha(Path(__file__)),v1_frozen_sha256=sha(old/'FROZEN.json')));json_write(out/'SHA256SUMS.json',manifest(out))
