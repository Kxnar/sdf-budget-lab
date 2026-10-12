"""Static scientific visuals from verified evidence; no fitting or gate changes."""
import argparse
from generalization_common import *
from generalization_v2_run import validate_freeze

def main():
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);a=ap.parse_args();d=a.directory.resolve();f=validate_freeze(d,False);verification=json.loads((d/'verification/VERIFIED.json').read_text());assert verification['frozen_sha256']==sha(d/'FROZEN.json');check_hashes(d/'verification',json.loads((d/'verification/SHA256SUMS.json').read_text()));out=ROOT/'docs/generalization-v2/final';out.mkdir(parents=True,exist_ok=True);matrix=rows(d/'verification/matrix.csv');labels=[g['id'].split('_',1)[1] for g in f['geometries']];sources={}
    fig,axes=plt.subplots(2,1,figsize=(9,7),sharex=True)
    for method in ('grid','cvar99','coverage_q8','candidate','compact_f32','compact_Q24'):
        for ax,ds,key in ((axes[0],'near','near_cvar99'),(axes[1],'boundary','boundary_max')):
            values=[];lower=[];upper=[]
            for g in f['geometries']:
                if method.startswith('compact_'):
                    cache=d/'verification'/f"{g['id']}_0_{method}.json";record=json.loads(cache.read_text())['assessment'];metric=next((m for m in record.get('metrics',[]) if m['set']==ds),{});values.append(metric.get('cvar99' if ds=='near' else 'maximum',np.nan));lower.append(values[-1]);upper.append(values[-1])
                else:
                    fits=[float(r[key]) for r in matrix if r['geometry']==g['id'] and r['method']==method and int(r['cap'])==CAPS[-1]];values.append(np.median(fits));lower.append(min(fits));upper.append(max(fits))
            line=ax.plot(range(4),values,marker='o',label=method.replace('compact_',''))[0]
            if method!='grid' and not method.startswith('compact_'):ax.fill_between(range(4),lower,upper,color=line.get_color(),alpha=.15)
            ax.set_yscale('log');ax.grid(axis='y',alpha=.2)
    axes[0].axhline(.010,color='black',ls='--',label='near gate');axes[1].axhline(.020,color='black',ls='--',label='sampled max gate');axes[0].set_ylabel('Near CVaR99');axes[1].set_ylabel('Boundary sampled maximum');axes[1].set_xticks(range(4),labels);axes[0].legend(ncol=3,fontsize=8);axes[1].legend(ncol=3,fontsize=8);fig.suptitle('At 1 MiB requested cap; adaptive median and five-fit min/max range\nIntrinsic compact outcomes shown even when reference topology is unknown',fontsize=11);fig.tight_layout();fig.savefig(out/'tails.png',dpi=180);fig.savefig(out/'tails.svg');plt.close(fig)
    for g in f['geometries']:
        if not g['qualified_levels']:continue
        level=g['qualified_levels'][0];ref=reference_dir(d,g);paths=[ref/f'oracle_0_{level}.obj',d/'jobs'/f"{g['id']}_0_compact_Q24"/f'surface_{level}.obj']
        if not all(p.exists() for p in paths):continue
        meshes=[obj(p) for p in paths];vv=np.concatenate([v for v,t in meshes]);center=(vv.min(axis=0)+vv.max(axis=0))/2;radius=.55*np.ptp(vv,axis=0).max();fig=plt.figure(figsize=(9,4))
        for i,((v,t),title) in enumerate(zip(meshes,('Original oracle extraction','Decoded Q24 extraction'))):
            ax=fig.add_subplot(1,2,i+1,projection='3d');ax.plot_trisurf(v[:,0],v[:,1],v[:,2],triangles=t,color='#628ba4',linewidth=0,antialiased=False);ax.set_box_aspect((1,1,1));ax.set(xlim=(center[0]-radius,center[0]+radius),ylim=(center[1]-radius,center[1]+radius),zlim=(center[2]-radius,center[2]+radius));ax.view_init(22,40);ax.set_title(title,fontsize=10);ax.set_axis_off()
        fig.suptitle(g['id']+f' at {level} cells; shared camera limits',fontsize=11);fig.tight_layout();fig.savefig(out/f"surface-{g['id']}.png",dpi=180);plt.close(fig)
        for p in paths:sources[str(p.relative_to(d))]=sha(p)
    g=next(g for g in f['geometries'] if g['qualified_levels']);path=d/'jobs'/f"{g['id']}_340_candidate"/'widths.csv';rr=[r for r in rows(path) if int(r['cap'])==CAPS[-1]];chosen=next((r for r in rr if len(np.fromstring(r['reference_roots'],sep=';'))!=len(np.fromstring(r['predicted_roots'],sep=';'))),next((r for r in rr if r['pass']=='0'),max(rr,key=lambda r:float(r['mean_relative_error']))));assessment='numerical unresolved' if chosen['numerical_unresolved']=='1' else 'observed geometric ray failure (native)' if chosen['pass']=='0' else 'dimension-error diagnostic; inspect group gate'
    truth=np.fromstring(chosen['reference_roots'],sep=';');pred=np.fromstring(chosen['predicted_roots'],sep=';');fig,ax=plt.subplots(figsize=(9,3))
    for y,values,label,color in ((1,truth,'Original reference transitions','#1f77b4'),(0,pred,'Depth 7 q8 candidate transitions','#ff7f0e')):
        ax.scatter(values,np.full(len(values),y),marker='|',s=200,label=label,color=color)
        for lo,hi in zip(values[::2],values[1::2]):ax.plot([lo,hi],[y,y],lw=7,alpha=.4,color=color)
    ax.set_yticks([0,1],['q8 candidate','reference']);ax.set_xlabel('Ray parameter (domain units)');ax.set_title(f"{g['id']}, ray {chosen['ray']}: complete signed transitions at 1 MiB\n{assessment}; numeric gates remain the acceptance evidence",fontsize=11);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'failed-ray.png',dpi=180);fig.savefig(out/'failed-ray.svg');plt.close(fig);sources[str(path.relative_to(d))]=sha(path)
    json_write(out/'FIGURE-PROVENANCE.json',dict(frozen_sha256=sha(d/'FROZEN.json'),verified_manifest_sha256=sha(d/'verification/SHA256SUMS.json'),helper_sha256=sha(Path(__file__)),sources=sources,no_geometry_repairs=True));json_write(out/'SHA256SUMS.json',manifest(out))
if __name__=='__main__':main()
