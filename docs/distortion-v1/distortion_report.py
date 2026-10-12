"""Create compact tables/figures with links to retained, hashed local raw evidence."""
import argparse
import html
import json
from pathlib import Path
import shutil
import statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from tail_evidence import rows, write_csv, sha

MODES = [('full_q8', 5), ('full_f32', 5), ('full_q8', 6), ('full_f32', 6), ('mse_q8', 5)]
COLORS = ['#dd7930', '#963d30', '#188579', '#224ea0', '#777777']


def main():
    p = argparse.ArgumentParser()
    p.add_argument('raw', type=Path)
    p.add_argument('analysis', type=Path)
    p.add_argument('out', type=Path)
    a = p.parse_args()
    if a.out.exists():
        raise ValueError('Use a new report directory')
    a.out.mkdir(parents=True)
    metrics, costs = rows(a.raw/'metrics.csv'), rows(a.raw/'construction.csv')
    widths, profiles = rows(a.analysis/'widths.csv'), rows(a.analysis/'profiles.csv')
    summary, width_summary = [], []
    for shape, feature in sorted({(r['shape'], float(r['feature'])) for r in metrics}):
        for method, depth in MODES:
            chosen = [r for r in metrics if r['shape']==shape and float(r['feature'])==feature and r['method']==method and int(r['depth'])==depth]
            cases = {r['case'] for r in chosen}
            row = dict(shape=shape, feature=feature, depth=depth, method=method, bytes=statistics.mean(int(r['bytes']) for r in chosen))
            for ds in ('uniform', 'near', 'probes'):
                for key in (('inside_recall', 'outside_recall') if ds=='probes' else ('mae', 'cvar99')):
                    row[f'{ds}_{key}'] = statistics.mean(float(r[key]) for r in chosen if r['dataset']==ds)
            for key in ('total_seconds', 'oracle_calls'):
                row[key] = statistics.mean(float(r[key]) for r in costs if r['case'] in cases and r['method']==method and int(r['depth'])==depth)
            summary.append(row)
            w = [r for r in widths if r['shape']==shape and float(r['feature'])==feature and r['method']==method and int(r['depth'])==depth]
            for axis in sorted({r['axis'] for r in w}):
                vals = np.array([float(r['width']) for r in w if r['axis']==axis])
                width_summary.append(dict(shape=shape, feature=feature, depth=depth, method=method, axis=axis, count=len(vals), minimum=float(vals.min()), median=float(np.median(vals)), mean=float(vals.mean()), maximum=float(vals.max()), zero_width=int(np.sum(vals==0)), missing=int(np.sum(~np.isfinite(vals)))))
    write_csv(a.out/'summary.csv', summary)
    write_csv(a.out/'width_summary.csv', width_summary)
    for name in ('VERIFIED.json', 'recovery.csv', 'coverage_summary.csv', 'decomposition.csv'):
        shutil.copy2(a.analysis/name, a.out/name)
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axs = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    for col, (shape, feature) in enumerate([('torus', .035), ('double_sphere', .06)]):
        for (method, depth), color in zip(MODES, COLORS):
            subset = [r for r in profiles if r['shape']==shape and float(r['feature'])==feature and r['method']==method and int(r['depth'])==depth and r['case'].endswith('_30')]
            x = np.array([float(r['coordinate']) for r in subset])
            if shape=='torus':
                x = np.degrees(x)
            label = f'{method}, d{depth}' + (' (32 KiB)' if method=='mse_q8' else '')
            style = '--' if method=='full_f32' else '-'
            axs[0, col].plot(x, [float(r['signed_error']) for r in subset], style, color=color, label=label, lw=1.2)
            axs[1, col].plot(x, [float(r['prediction']) for r in subset], style, color=color, lw=1.2)
        axs[0, col].axhline(0, color='black', lw=.7)
        axs[1, col].axhline(0, color='black', lw=.7)
        axs[1, col].plot(x, [float(r['truth']) for r in subset], color='black', ls=':', label='Oracle')
        for row in range(2):
            axs[row, col].set_xlabel('Centreline angle (degrees)' if shape=='torus' else 'Gap axis relative to midpoint')
            axs[row, col].grid(alpha=.2)
        axs[0, col].set_title('Torus radius 0.035' if shape=='torus' else 'Sphere gap 0.060')
    axs[0, 0].set_ylabel('Signed error: field minus oracle')
    axs[1, 0].set_ylabel('Signed distance (zero = boundary)')
    axs[0, 0].legend(fontsize=8, loc='upper right')
    axs[1, 1].legend(fontsize=8)
    fig.suptitle('Precision barely changes distortion; finer cells remove the sign failures\nTranslated analytic cases; seed 30 for the budgeted comparator', fontsize=13)
    fig.savefig(a.out/'signed_error.png', dpi=170)
    fig.savefig(a.out/'signed_error.pdf')
    plt.close(fig)
    fig, axs = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    for (method, depth), color in zip(MODES, COLORS):
        for col, axis in enumerate(('radial', 'vertical')):
            subset = [r for r in widths if r['shape']=='torus' and float(r['feature'])==.035 and r['axis']==axis and r['method']==method and int(r['depth'])==depth and r['case'].endswith('_30')]
            axs[col].plot(np.degrees([float(r['angle']) for r in subset]), [float(r['width']) for r in subset], '--' if method=='full_f32' else '-', color=color, lw=1.1, label=f'{method}, d{depth}')
            axs[col].set_title(f'Torus {axis} thickness (true 0.070)')
            axs[col].set_xlabel('Centreline angle (degrees)')
        gap = [r for r in widths if r['shape']=='double_sphere' and float(r['feature'])==.06 and r['method']==method and int(r['depth'])==depth]
        index = MODES.index((method, depth))
        axs[2].bar(index, statistics.mean(float(r['width']) for r in gap), color=color)
        axs[2].scatter([index]*len(gap), [float(r['width']) for r in gap], color='black', s=12)
    for ax in axs[:2]:
        ax.axhline(.07, color='black', ls=':')
        ax.set_ylabel('Contiguous negative width at centreline')
        ax.grid(alpha=.2)
    axs[0].legend(fontsize=7)
    axs[2].axhline(.06, color='black', ls=':')
    axs[2].set_title('Sphere gap width (true 0.060)')
    axs[2].set_xticks(range(5), ['q8 d5', 'f32 d5', 'q8 d6', 'f32 d6', '32 KiB\nq8 d5'], rotation=35)
    axs[2].set_ylabel('Contiguous positive width at midpoint')
    fig.savefig(a.out/'widths.png', dpi=170)
    fig.savefig(a.out/'widths.pdf')
    plt.close(fig)
    # Keep raw/analysis manifests here, with their exact local relative destinations.
    import os
    rawlink = Path(os.path.relpath(a.raw, a.out)).as_posix()
    analysislink = Path(os.path.relpath(a.analysis, a.out)).as_posix()
    source = a.raw.with_name(a.raw.name+'-source')
    sourcelink = Path(os.path.relpath(source, a.out)).as_posix()
    (a.out/'EVIDENCE.json').write_text(json.dumps(dict(raw=rawlink, analysis=analysislink, source=sourcelink, raw_manifest_sha256=sha(a.raw/'SHA256SUMS.json'), source_manifest_sha256=sha(source/'SHA256SUMS.json')), indent=2))
    table = '<table><tr>'+''.join(f'<th>{html.escape(k)}</th>' for k in summary[0])+'</tr>'
    for r in summary:
        table += '<tr>'+''.join(f'<td>{v:.6g}</td>' if isinstance(v, float) else f'<td>{html.escape(str(v))}</td>' for v in r.values())+'</tr>'
    table += '</table>'
    (a.out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Geometry distortion diagnostic</title><style>body{font:16px system-ui;margin:32px;max-width:1500px}img{width:100%}table{border-collapse:collapse;font-size:12px}td,th{padding:7px;border:1px solid #ddd}th{background:#eef3f5}</style><h1>Geometry distortion diagnostic</h1><p>Full depth six recovers the tested 0.035 torus and 0.060 gap; changing depth-five q8 to f32 does not. This is an unequal-byte resolution diagnostic, not a 32 KiB improvement.</p><p>Recovery criteria are unchanged. Widths supplement the original tests; they expose residual shape error even when sign and topology pass.</p><img src="signed_error.png"><img src="widths.png"><h2>Two-seed means, construction and retained bytes</h2>'+table+f'<p><a href="recovery.csv">Frozen recovery results</a> · <a href="coverage_summary.csv">Fitting coverage</a> · <a href="decomposition.csv">Error decomposition</a> · <a href="{rawlink}/SHA256SUMS.json">Raw evidence manifest</a> · <a href="{analysislink}/profiles.csv">Profiles</a> · <a href="{sourcelink}/SHA256SUMS.json">Source/executable manifest</a></p><p>Raw evidence is retained locally in the linked runs directories; it is not bundled into this compact report.</p>', encoding='utf-8')
    for filename in ('distortion_analysis.py', 'distortion_report.py', 'tail_evidence.py'):
        shutil.copy2(Path(__file__).parent/filename, a.out/filename)
    (a.out/'SHA256SUMS.json').write_text(json.dumps({f.relative_to(a.out).as_posix(): sha(f) for f in a.out.rglob('*') if f.is_file()}, indent=2))


if __name__ == '__main__':
    main()
