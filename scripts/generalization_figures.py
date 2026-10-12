"""Static scientific figures from verified data; no method selection or fitting."""
import argparse
from generalization_common import *
def main():
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);ap.add_argument('--out',type=Path,default=ROOT/'docs/generalization-v1');a=ap.parse_args();directory=a.directory.resolve();out=a.out.resolve();out.mkdir(exist_ok=True);frozen=json.loads((directory/'FROZEN.json').read_text());matrix=rows(directory/'verification/matrix.csv');geometry=frozen['geometries'];methods=frozen['methods'];values=np.zeros((len(geometry),len(methods)))
    for i,g in enumerate(geometry):
        for j,m in enumerate(methods):
            seeds=[0] if m=='grid' else SEEDS;values[i,j]=np.nan if g['reference_status']=='oracle_resolution_unresolved' else sum(any(r['geometry']==g['id'] and r['method']==m and int(r['seed'])==s and r['pass']=='True' for r in matrix) for s in seeds)/len(seeds)
    cmap=plt.get_cmap('Blues').copy();cmap.set_bad('#dddddd')
    fig,ax=plt.subplots(figsize=(12,7));image=ax.imshow(values,vmin=0,vmax=1,cmap=cmap,aspect='auto');ax.set_xticks(range(len(methods)),methods,rotation=30,ha='right');ax.set_yticks(range(len(geometry)),[g['id'] for g in geometry]);ax.set_title('Strict joint-gate pass fraction across fit seeds\n(any tested cap; annotations = passes/replications; grey = unresolved reference)')
    for i in range(len(geometry)):
        for j in range(len(methods)):
            count=1 if methods[j]=='grid' else len(SEEDS)
            ax.text(j,i,'?' if np.isnan(values[i,j]) else f'{round(values[i,j]*count)}/{count}',ha='center',va='center',color='white' if values[i,j]>.65 else'black',fontsize=9)
    fig.colorbar(image,ax=ax,label='Observed pass fraction');fig.tight_layout();fig.savefig(out/'joint-gates.png',dpi=180);fig.savefig(out/'joint-gates.svg');plt.close(fig)
    performance=rows(out/'isolated-performance-summary.csv')
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for ax,key,label in zip(axes,('warm_query_median_ns','cold_build_median_seconds'),('Warm query median (ns/query)','Fresh-process construction median (s)')):
        ax.bar(range(len(performance)),[float(r[key]) for r in performance],color='#5C83A0');ax.set_xticks(range(len(performance)),[r['method'] for r in performance],rotation=40,ha='right');ax.set_ylabel(label);ax.set_yscale('log');ax.grid(axis='y',alpha=.2)
    fig.suptitle('Selected reproduction fields: crossings and failures\nverification stopped; background workload possible; field mix differs by method',fontsize=11);fig.tight_layout();fig.savefig(out/'performance.png',dpi=180);fig.savefig(out/'performance.svg');plt.close(fig)
    pairs=rows(out/'paired-savings.csv');fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
    for ax,method in zip(axes,('grid','optimistic_best_q8')):
        for i,g in enumerate(geometry):
            values=[100*float(r['saving_fraction']) for r in pairs if r['geometry']==g['id'] and r['comparator']==method and r['saving_fraction']!='']
            if values:ax.errorbar(i,np.median(values),yerr=[[np.median(values)-min(values)],[max(values)-np.median(values)]],fmt='o',color='#33678A',capsize=4)
            else:ax.text(i,0,'?',ha='center',color='#888888')
        ax.axhline(0,color='black',lw=1);ax.set_ylabel('Retained-byte saving (%)');ax.set_title('Finite paired comparisons against '+method);ax.grid(axis='y',alpha=.2)
    axes[1].set_xticks(range(len(geometry)),[g['id'] for g in geometry],rotation=40,ha='right');fig.suptitle('Conditional median and range across fitting seeds\n? = no finite paired crossing; no percentage inferred',fontsize=11);fig.tight_layout();fig.savefig(out/'memory-savings.png',dpi=180);fig.savefig(out/'memory-savings.svg');plt.close(fig)
    # Error tails at the largest cap, explicitly showing the acceptance line.
    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
    for m in ('grid','mse','cvar99','candidate_original','coverage_q8','candidate'):
        for ax,key in zip(axes,('near_cvar99','boundary_max')):
            med=[]
            for g in geometry:med.append(np.median([float(r[key]) for r in matrix if r['geometry']==g['id'] and r['method']==m and int(r['cap'])==CAPS[-1]]))
            ax.plot(range(len(geometry)),med,marker='o',markersize=3,label=m)
    axes[0].axhline(.010,color='black',ls='--',lw=1,label='near CVaR gate');axes[1].axhline(.020,color='black',ls='--',lw=1,label='sampled max gate');axes[0].set_ylabel('Near CVaR99 (domain units)');axes[1].set_ylabel('Boundary sampled max');axes[0].set_title('Median across fit seeds at 1024 KiB requested cap');axes[1].set_xticks(range(len(geometry)),[g['id'] for g in geometry],rotation=40,ha='right');axes[0].legend(ncol=3,fontsize=8);axes[1].legend(ncol=3,fontsize=8);fig.tight_layout();fig.savefig(out/'tails-and-maxima.png',dpi=180);fig.savefig(out/'tails-and-maxima.svg');plt.close(fig)
    # Preserve exact source meshes; plotting changes neither reference nor field.
    selected=[];source_meshes={}
    for g in geometry:
        if g['stratum'] not in [s for _,s in selected]:selected.append((g,g['stratum']))
    for g,stratum in selected:
        folder=directory/'jobs'/f"{g['id']}_240_candidate";fields=[r for r in matrix if r['geometry']==g['id'] and r['method']=='candidate' and int(r['seed'])==240];passing=[r for r in fields if r['pass']=='True'];candidate_mesh=folder/'candidate_32768_diagnostic64.obj';label='Failed 32 KiB diagnostic'
        if passing:
            r=min(passing,key=lambda x:int(x['bytes']));n=g['qualified_levels'][0];candidate_mesh=folder/f"candidate_{r['cap']}_{n}.obj";label=f"Passing {int(r['bytes']):,} retained bytes"
        ref=reference_dir(directory,g);level=int(rows(ref/'reference_topology.csv')[0]['cells']);reference_mesh=ref/f'oracle_0_{level}.obj'
        if not candidate_mesh.exists():continue
        for p in (reference_mesh,candidate_mesh):source_meshes[str(p.relative_to(directory))]=sha(p)
        meshes=[obj(p) for p in (reference_mesh,candidate_mesh)];vertices=np.concatenate([v for v,f in meshes]);lo=vertices.min(axis=0);hi=vertices.max(axis=0);center=(lo+hi)/2;radius=.55*(hi-lo).max()
        fig=plt.figure(figsize=(9,4))
        for index,((vv,ff),title) in enumerate(zip(meshes,('Reference extraction',label))):
            ax=fig.add_subplot(1,2,index+1,projection='3d');ax.plot_trisurf(vv[:,0],vv[:,1],vv[:,2],triangles=ff,color='#5C83A0',linewidth=0,antialiased=False);ax.set_box_aspect((1,1,1));ax.set(xlim=(center[0]-radius,center[0]+radius),ylim=(center[1]-radius,center[1]+radius),zlim=(center[2]-radius,center[2]+radius));ax.view_init(elev=22,azim=40);ax.set_title(title,fontsize=10);ax.set_axis_off()
        fig.suptitle(g['id']+' — '+stratum,fontsize=11);fig.tight_layout();fig.savefig(out/f"surface-{g['id']}.png",dpi=180);plt.close(fig)
    json_write(out/'figure-provenance.json',dict(frozen_sha256=sha(directory/'FROZEN.json'),verified_matrix_sha256=sha(directory/'verification/matrix.csv'),script_sha256=sha(Path(__file__)),source_mesh_sha256=source_meshes,sources='saved oracle and represented extracted meshes, no repair or triangulation changes'))
if __name__=='__main__':main()
