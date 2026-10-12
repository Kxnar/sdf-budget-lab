"""Record the predeclared bounded intervention before any v2 fits."""
from generalization_common import *
import datetime
out=ROOT/'runs/generalization-v2';out.mkdir(exist_ok=True)
p=out/'DECLARATION.json'
assert not p.exists(),'Preserve the existing pre-fit declaration'
json_write(p,dict(protocol='generalization-v2',created_utc=datetime.datetime.now(datetime.UTC).isoformat(),protocol_sha256=sha(ROOT/'docs/generalization-v2-protocol.md'),old_frozen_sha256=sha(ROOT/'runs/generalization-v1/FROZEN.json'),policies=['old_d7','capacity4096','data4096','focused2048','focused4096'],depths=[7,9,12],virtual_encodings=['q8','f32','f64'],screen_meshes=['12_goathead','14_spot','15_brucewick'],screen_seeds=[240,241],caps=CAPS,fresh_fit_seeds=list(range(340,345)),fresh_query_namespace=1880001,fresh_revision='835213277571d625722d25880c9c3010c6f45a5a',fresh_files=['fish/fish_low_resolution.obj','mushroom/mushroom.obj','springer/springer.obj','hammer/hammer.obj'],fallback_order=['stuffedtoy','boat'],rank_scale_floor=1e-5,weight_scale=.01,multiscale_anchors=32768,multiscale_labels=557056,max_nodes=250441,max_dp_states=8000000,max_transitions=1000000000,cache_entries=262144,heavy_workers=1,reserved_free_disk_bytes=30*1024**3,target_new_evidence_bytes=40*1024**3))
print(p)
