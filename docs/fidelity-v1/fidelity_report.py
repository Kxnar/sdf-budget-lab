"""Publish compact tables and static figures from verified, frozen evidence."""
import argparse
import json
import shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tail_evidence import rows, hashes, sha, write_csv

p=argparse.ArgumentParser()
p.add_argument('directory',type=Path)
p.add_argument('analysis',type=Path)
p.add_argument('reproduction',type=Path)
p.add_argument('out',type=Path)
a=p.parse_args()
for directory in (a.directory,a.analysis,a.reproduction,a.reproduction.with_name(a.reproduction.name+'-verified')):
    hashes(directory)
assert not a.out.exists()
a.out.mkdir(parents=True)
metrics=rows(a.directory/'metrics.csv')
costs=rows(a.directory/'construction.csv')
widths=rows(a.analysis/'widths.csv')
recovery=rows(a.analysis/'recovery.csv')
summary=[]
for shape,features in [('torus',[.035,.06,.1]),('double_sphere',[.06,.1])]:
    for feature in features:
        for method in ('baseline_q8','risk_q8'):
            rs=[r for r in metrics if r['shape']==shape and float(r['feature'])==feature and r['method']==method]
            cases={r['case'] for r in rs}
            ds=dict(shape=shape,feature=feature,method=method,bytes_min=min(int(r['bytes']) for r in rs),bytes_max=max(int(r['bytes']) for r in rs))
            for distribution in ('uniform','near'):
                rr=[r for r in rs if r['dataset']==distribution]
                assert len(rr)==2
                for metric in ('mae','cvar99'):
                    ds[distribution+'_'+metric]=np.mean([float(r[metric]) for r in rr])
            probes=[r for r in rs if r['dataset']=='probes']
            ds['inside_recall_min']=min(float(r['inside_recall']) for r in probes)
            ds['outside_recall_min']=min(float(r['outside_recall']) for r in probes)
            cc=[r for r in costs if r['case'] in cases and r['method']==method]
            ds['oracle_calls']=np.mean([int(r['oracle_calls']) for r in cc])
            ds['construction_ms']=1000*np.mean([float(r['total_seconds']) for r in cc])
            ds['recovered']=next(r['recovered'] for r in recovery if r['shape']==shape and float(r['feature'])==feature and r['method']==method)
            summary.append(ds)
write_csv(a.out/'summary.csv',summary)
width_summary=[]
for shape,feature,method,axis in sorted({(r['shape'],float(r['feature']),r['method'],r['axis']) for r in widths}):
    rr=[r for r in widths if (r['shape'],float(r['feature']),r['method'],r['axis'])==(shape,feature,method,axis)]
    values=np.array([float(r['width']) for r in rr])
    relative=np.array([float(r['relative_error']) for r in rr])
    medians=[np.median([float(r['width']) for r in rr if r['case']==c]) for c in sorted({r['case'] for r in rr})]
    width_summary.append(dict(shape=shape,feature=feature,method=method,axis=axis,truth=rr[0]['expected'],median_seed_1=medians[0],median_seed_2=medians[1],mean_absolute_relative_error=np.mean(np.abs(relative)),p95_absolute_relative_error=np.quantile(np.abs(relative),.95),zero_count=int(np.sum(values==0)),missing_count=int(np.sum(~np.isfinite(values)))))
write_csv(a.out/'width_summary.csv',width_summary)

fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
for index,(method,label,color) in enumerate([('baseline_q8','Selective d6 baseline','#747b87'),('risk_q8','Frozen allocation + coverage','#167e93')]):
    groups=[('torus',.035,'radial'),('torus',.035,'vertical'),('double_sphere',.06,'gap')]
    vals=[]
    for shape,feature,axis in groups:
        r=next(r for r in width_summary if (r['shape'],r['feature'],r['axis'],r['method'])==(shape,feature,axis,method))
        vals.append(np.mean([r['median_seed_1'],r['median_seed_2']])/float(r['truth']))
    axs[0].bar(np.arange(3)+(index-.5)*.32,vals,.32,label=label,color=color)
axs[0].axhline(1,color='black',linestyle='--',linewidth=1)
axs[0].set(xticks=np.arange(3),xticklabels=['Torus radial','Torus vertical','Sphere gap'],ylabel='Median measured width / true width',title='Dimensional fidelity on fresh seeds')
axs[0].legend(fontsize=8,loc='lower right')
labels=[];ratios=[]
for shape,feature in [('torus',.035),('double_sphere',.06)]:
    base=next(r for r in summary if (r['shape'],r['feature'],r['method'])==(shape,feature,'baseline_q8'))
    new=next(r for r in summary if (r['shape'],r['feature'],r['method'])==(shape,feature,'risk_q8'))
    for metric in ('uniform_mae','uniform_cvar99','near_mae','near_cvar99'):
        labels.append(('Torus' if shape=='torus' else 'Gap')+' '+metric.replace('_',' '))
        ratios.append(100*(new[metric]/base[metric]-1))
axs[1].barh(labels,ratios,color=['#ad554e' if r>0 else '#167e93' for r in ratios])
axs[1].axvline(0,color='black',linewidth=.8)
axs[1].set(xlabel='Change from baseline (%)',title='Error tradeoffs (lower is better)')
fig.savefig(a.out/'fidelity.png',dpi=180)
plt.close(fig)

# Independent rebuild reproduction: compare every deterministic saved result.
comparisons=[]
for folder in ('fields','samples','predictions','meshes'):
    first={f.relative_to(a.directory).as_posix():sha(f) for f in (a.directory/folder).rglob('*') if f.is_file()}
    second={f.relative_to(a.reproduction).as_posix():sha(f) for f in (a.reproduction/folder).rglob('*') if f.is_file()}
    assert first==second,folder
    comparisons.append(dict(group=folder,files=len(first),identical=True))
for name in ('metrics.csv','allocation.csv','selection.csv','surfaces.csv','policy.json'):
    assert sha(a.directory/name)==sha(a.reproduction/name),name
    comparisons.append(dict(group=name,files=1,identical=True))
(a.out/'reproduction.json').write_text(json.dumps(dict(comparisons=comparisons,first_manifest=sha(a.directory/'SHA256SUMS.json'),reproduction_manifest=sha(a.reproduction/'SHA256SUMS.json'),timing='Construction times are independently measured and excluded from exact comparison.'),indent=2))
for src in [a.directory/'metrics.csv',a.directory/'construction.csv',a.directory/'allocation.csv',a.directory/'surfaces.csv',a.analysis/'recovery.csv',a.analysis/'VERIFIED.json',Path('runs/fidelity-frozen/selection.json'),Path('runs/fidelity-frozen/development_comparison.csv'),Path(__file__)]:
    shutil.copy2(src,a.out/src.name)
(a.out/'SHA256SUMS.json').write_text(json.dumps({f.name:sha(f) for f in a.out.iterdir() if f.is_file()},indent=2))
print(json.dumps(summary,indent=2))
print(json.dumps(width_summary,indent=2))
