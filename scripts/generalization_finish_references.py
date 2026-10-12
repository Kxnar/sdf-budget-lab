"""Complete reference-only qualification with explicitly audited fast meshes."""
import argparse,subprocess,time
from generalization_common import *
from generalization_prepare import create_reference
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);ap.add_argument('--exe',type=Path,default=ROOT/'build-generalization-cgal/Release/sdf_generalization_bench.exe');a=ap.parse_args();directory=a.directory.resolve();exe=a.exe.resolve();summary=json.loads((directory/'oracle-audit/SUMMARY.json').read_text());assert summary['all_passed'] and summary['chosen_oracle']=='EPICK-audited';assert not (directory/'FROZEN.json').exists();geometries=json.loads((directory/'qualification.json').read_text())
    for i,g in enumerate(geometries):
        if g['shape']!='mesh':g['oracle']='analytic-exact-formula';g['reference_subdir']='reference/'+g['id'];continue
        g['parameters'][1]=1;g['oracle']='CGAL:EPICK-audited';g['oracle_comparison_sha256']=sha(directory/'oracle-audit'/g['id']/'RESULT.json');g['epeck_validation_sha256']=g['oracle_comparison_sha256'];folder=directory/'reference-fast'/g['id'];g['reference_subdir']='reference-fast/'+g['id']
        if (folder/'metadata.json').exists():g=json.loads((folder/'metadata.json').read_text())
        else:
            assert not folder.exists(),'Preserve incomplete reference and choose a new qualification attempt directory';print('Qualifying audited fast reference',g['id'],flush=True);g=create_reference(g,folder,exe,seed=880001+i*1000);json_write(folder/'metadata.json',g)
        geometries[i]=g;json_write(directory/'qualification.json',geometries)
    json_write(directory/'qualification.json',geometries);json_write(directory/'REFERENCES_COMPLETE.json',dict(geometries=16,oracle='analytic exact formulas plus separately audited EPICK rounded construction',source_kernel_comparison_sha256=sha(directory/'oracle-audit/SUMMARY.json')))
if __name__=='__main__':main()
