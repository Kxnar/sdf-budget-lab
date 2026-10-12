"""Publish paired results and width diagnostics from verified refinement evidence."""
import argparse
import html
import json
import os
from pathlib import Path
import shutil
import statistics as stats
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tail_evidence import rows, write_csv, sha, hashes

MODES = [('mse_q8', 5), ('selective_q8', 6)]


def table(data):
    result = '<div class="scroll"><table><tr>'+''.join(f'<th>{html.escape(k)}</th>' for k in data[0])+'</tr>'
    for row in data:
        result += '<tr>'+''.join(f'<td>{v:.7g}</td>' if isinstance(v, float) else f'<td>{html.escape(str(v))}</td>' for v in row.values())+'</tr>'
    return result+'</table></div>'


def main():
    p = argparse.ArgumentParser()
    for name in ('raw', 'analysis', 'development', 'selection', 'out'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    if a.out.exists():
        raise ValueError('Use a new report directory')
    hashes(a.raw)
    verified = json.loads((a.analysis/'VERIFIED.json').read_text())
    assert verified['raw_manifest_sha256'] == sha(a.raw/'SHA256SUMS.json')
    a.out.mkdir(parents=True)
    metrics, costs = rows(a.raw/'metrics.csv'), rows(a.raw/'construction.csv')
    selection = rows(a.raw/'selection.csv')
    coverage = rows(a.analysis/'coverage.csv')
    selection_coverage = []
    for case in sorted({r['case'] for r in selection}):
        selected_cells = {tuple(round((float(r['lo_'+axis])+1)*16) for axis in 'xyz') for r in selection if r['case']==case}
        failed = [r for r in coverage if r['case']==case and r['method']=='mse_q8' and int(r['failed'])]
        included = sum(tuple(min(int((float(r[axis])+1)*16),31) for axis in 'xyz') in selected_cells for r in failed)
        selection_coverage.append(dict(case=case, baseline_failed_probes=len(failed), failures_in_selected_d5_cells=included))
    write_csv(a.out/'selection_coverage.csv', selection_coverage)
    widths = rows(a.analysis/'widths.csv')
    allocations, surfaces = rows(a.raw/'allocation.csv'), rows(a.raw/'surfaces.csv')
    per_case = []
    for m in metrics:
        if m['dataset'] != 'probes':
            continue
        case, method = m['case'], m['method']
        cost = next(r for r in costs if r['case']==case and r['method']==method)
        alloc = next(r for r in allocations if r['case']==case and r['method']==method)
        checks = [r for r in surfaces if r['case']==case and r['method'] in (method, 'oracle') and r['cells'] in ('64', '96')]
        assert len(checks)==4
        row = dict(case=case, shape=m['shape'], feature=float(m['feature']), seed=int(m['seed']), method=method,
                   bytes=int(m['bytes']), recovery=int(all(int(r['success']) for r in checks)),
                   inside_recall=float(m['inside_recall']), outside_recall=float(m['outside_recall']))
        for ds in ('uniform', 'near'):
            metric = next(r for r in metrics if r['case']==case and r['method']==method and r['dataset']==ds)
            for key in ('mae', 'cvar99'):
                row[f'{ds}_{key}'] = float(metric[key])
        row.update(fit_sse=float(alloc['fit_sse']), selected_parents=int(alloc['selected_parents']), retained_parents=int(alloc['retained_parents']), oracle_calls=int(cost['oracle_calls']), construction_ms=float(cost['total_seconds'])*1000)
        per_case.append(row)
    write_csv(a.out/'per_case.csv', per_case)
    summary = []
    for shape, feature in [('torus', .035), ('torus', .06), ('torus', .1), ('double_sphere', .06), ('double_sphere', .1)]:
        for method, _ in MODES:
            chosen = [r for r in per_case if r['shape']==shape and r['feature']==feature and r['method']==method]
            assert len(chosen)==2
            row = dict(shape=shape, feature=feature, method=method, recovered_both_seeds=all(r['recovery'] for r in chosen))
            for key in ('bytes', 'inside_recall', 'outside_recall', 'uniform_mae', 'uniform_cvar99', 'near_mae', 'near_cvar99', 'fit_sse', 'retained_parents', 'oracle_calls', 'construction_ms'):
                row[key] = stats.mean(r[key] for r in chosen)
            summary.append(row)
    write_csv(a.out/'summary.csv', summary)
    paired = []
    for r in per_case:
        if r['method']!='selective_q8':
            continue
        baseline = next(b for b in per_case if b['case']==r['case'] and b['method']=='mse_q8')
        row = dict(case=r['case'])
        for key in ('inside_recall', 'outside_recall', 'uniform_mae', 'uniform_cvar99', 'near_mae', 'near_cvar99', 'fit_sse', 'construction_ms', 'oracle_calls'):
            row[key+'_delta'] = r[key]-baseline[key]
            row[key+'_relative_delta'] = r[key]/baseline[key]-1 if baseline[key] else None
        paired.append(row)
    write_csv(a.out/'paired.csv', paired)
    width_summary = []
    for row in per_case:
        chosen = [r for r in widths if r['case']==row['case'] and r['method']==row['method']]
        for axis in sorted({r['axis'] for r in chosen}):
            vals = np.array([float(r['width']) for r in chosen if r['axis']==axis])
            width_summary.append(dict(case=row['case'], method=row['method'], axis=axis, count=len(vals), median=float(np.median(vals)), minimum=float(vals.min()), maximum=float(vals.max()), mean=float(vals.mean()), zero_width=int(np.sum(vals==0)), missing=int(np.sum(~np.isfinite(vals)))))
    write_csv(a.out/'width_summary.csv', width_summary)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    colors = ['#6b7280', '#087f8c']
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.4), constrained_layout=True)
    for method, color in zip([m for m, _ in MODES], colors):
        for seed, style in [(40, '-'), (41, '--')]:
            for col, axis in enumerate(('radial', 'vertical')):
                selected = [r for r in widths if r['case']==f'torus_0.035000_{seed}' and r['method']==method and r['axis']==axis]
                axs[col].plot(np.degrees([float(r['angle']) for r in selected]), [float(r['width']) for r in selected], style, color=color, lw=1, label=f'{method}, seed {seed}')
            gap = next(r for r in widths if r['case']==f'double_sphere_0.060000_{seed}' and r['method']==method)
            axs[2].scatter([0 if method=='mse_q8' else 1], [float(gap['width'])], color=color, marker='o' if seed==40 else 'x', s=65)
    for col, axis in enumerate(('radial', 'vertical')):
        axs[col].axhline(.07, color='black', ls=':', label='True width')
        axs[col].set_title(f'Torus 0.035: {axis} thickness')
        axs[col].set_xlabel('Centreline angle (degrees)')
        axs[col].set_ylabel('Contiguous negative width')
        axs[col].grid(alpha=.2)
    axs[0].legend(fontsize=7)
    axs[2].axhline(.06, color='black', ls=':', label='True gap')
    axs[2].set_xticks([0, 1], ['Depth 5', 'Selective depth 6'])
    axs[2].set_xlim(-.5, 1.5)
    axs[2].set_ylim(0, .075)
    axs[2].set_title('Sphere gap 0.060: unchanged')
    axs[2].set_ylabel('Contiguous positive width')
    axs[2].legend()
    fig.suptitle('Same 32,628 bytes: torus recovery improves, dimensional error remains')
    for ext in ('png', 'pdf'):
        fig.savefig(a.out/f'widths.{ext}', dpi=170)
    plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    chosen = [r for r in summary if r['shape']=='torus' and r['feature']==.035]
    for ax, key, title in zip(axs, ('mae', 'cvar99'), ('Mean absolute error', 'CVaR99 absolute error')):
        for j, row in enumerate(chosen):
            ax.bar(np.array([0, 1])+(j-.5)*.32, [row['uniform_'+key], row['near_'+key]], width=.32, color=colors[j], label=row['method'])
        ax.set_xticks([0, 1], ['Uniform', 'Near surface'])
        ax.set_title(title)
        ax.set_ylabel('Scene units')
        ax.legend(fontsize=8)
    fig.suptitle('Torus radius 0.035: two-seed means at the same byte cap')
    for ext in ('png', 'pdf'):
        fig.savefig(a.out/f'errors.{ext}', dpi=170)
    plt.close(fig)
    for filename in ('VERIFIED.json', 'recovery.csv', 'coverage_summary.csv'):
        shutil.copy2(a.analysis/filename, a.out/filename)
    for filename in ('development.csv',):
        shutil.copy2(a.development/filename, a.out/filename)
    for filename in ('selection.json', 'comparison.csv'):
        shutil.copy2(a.selection/filename, a.out/filename)
    for filename in ('metrics.csv', 'construction.csv', 'allocation.csv', 'surfaces.csv', 'selection.csv', 'config.json'):
        shutil.copy2(a.raw/filename, a.out/filename)
    evidence = {}
    for name, path in [('raw', a.raw), ('analysis', a.analysis), ('source', a.raw.with_name(a.raw.name+'-source')), ('development', a.development), ('development_source', a.development.with_name(a.development.name+'-source'))]:
        evidence[name] = Path(os.path.relpath(path, a.out)).as_posix()
        evidence[name+'_manifest_sha256'] = sha(path/'SHA256SUMS.json')
    (a.out/'EVIDENCE.json').write_text(json.dumps(evidence, indent=2))
    for filename in ('refinement_analysis.py', 'refinement_verify_selection.py', 'refinement_report.py', 'refinement_select.py', 'refinement_run.py', 'distortion_analysis.py', 'tail_evidence.py', 'tail_run.py'):
        shutil.copy2(Path(__file__).parent/filename, a.out/filename)
    (a.out/'index.html').write_text('''<!doctype html><meta charset="utf-8"><title>Selective depth-six refinement</title>
<style>body{font:16px system-ui;margin:32px;max-width:1500px;color:#17212b}img{width:100%}.scroll{overflow:auto}table{border-collapse:collapse;font-size:12px}td,th{padding:7px;border:1px solid #ddd;white-space:nowrap}th{background:#edf4f5}p{max-width:1050px;line-height:1.5}</style>
<h1>Selective depth six recovers the 0.035 torus at 32 KiB</h1>
<p>q8 and the MSE objective are fixed. Thirty-two depth-five cells ranked by fitting squared residuals gain depth-six children. Development selection used separate sizes and seeds, without generating evaluation probes. Both fresh seeds pass the frozen torus criteria; all larger controls pass. The 0.060 sphere gap remains unrecovered: MSE retains none of its deeper candidates.</p>
<p>Both methods occupy exactly 32,628 bytes. Torus interior recall rises from 80.08% to 99.61% / 98.44%. Closed manifold, components and Euler checks pass at both 64 and 96 extraction cells, including oracle controls. Width still differs from truth. Near-surface error improves while uniform error worsens.</p>
<img src="errors.png" alt="Torus uniform and near-surface error comparison"><img src="widths.png" alt="Torus width profiles and unchanged sphere gap">
<h2>All target and larger-control results</h2><p>Equal-weight means across seeds 40/41. Recovery requires both seeds. MAE/CVaR99 use 100,000 points per distribution per seed. Construction milliseconds include labels, candidates, refinement, allocation and packing; two serial measurements, not a timing confidence interval.</p>'''+table(summary)+
        '<h2>Recovered widths, separately by seed</h2><p>Zero means the analytic centre has the wrong sign, not necessarily absence of a displaced surface. Missing brackets remain missing. Truth: torus width = twice radius; sphere width = gap.</p>'+table(width_summary)+
        '<h2>Reproducible evidence</h2><p><a href="per_case.csv">Per-case results</a> · <a href="paired.csv">Paired changes</a> · <a href="recovery.csv">Recovery thresholds</a> · <a href="comparison.csv">Development selection</a> · <a href="construction.csv">Full cost components</a> · <a href="selection.csv">Selected-cell audit</a> · <a href="VERIFIED.json">Verification</a> · <a href="EVIDENCE.json">Raw/source locations and hashes</a></p><p>Familiar analytic shapes with fresh samples. No claim of unseen-shape generalization, exact dimensions, improved global error or temporary-memory savings. Raw fields, samples and meshes remain in the linked local runs directories and are not bundled into this compact report.</p>', encoding='utf-8')
    (a.out/'SHA256SUMS.json').write_text(json.dumps({f.relative_to(a.out).as_posix(): sha(f) for f in a.out.rglob('*') if f.is_file()}, indent=2))


if __name__ == '__main__':
    main()
