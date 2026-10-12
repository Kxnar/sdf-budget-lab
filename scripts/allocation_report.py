"""Validate paired allocation evidence and produce a portable report (stdlib only)."""
import argparse
import csv
import hashlib
import html
import json
import math
import platform
import statistics
import subprocess
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METHODS = ['exact', 'greedy_spend', 'greedy_best', 'lookahead_2', 'lookahead_4', 'lookahead_8']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def ratio(a, b):
    return a / b if b > 1e-12 else None


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--exe', required=True, type=Path)
    parser.add_argument('--archive', type=Path, help='New evidence directory with frozen-data ZIP and readable summaries')
    args = parser.parse_args()
    if args.archive and args.archive.exists():
        raise ValueError('Archive destination must not exist')
    run = args.run.resolve()
    if not (run / 'COMPLETE').exists():
        raise ValueError('Incomplete run; retain partial evidence but do not summarize it')
    config = json.loads((run / 'config.json').read_text())
    rows = read_csv(run / 'metrics.csv')
    groups = {}
    all_keys = set()
    for r in rows:
        key = (r['case'], int(r['seed']), int(r['budget']))
        full_key = (*key, r['method'], r['dataset'])
        if full_key in all_keys:
            raise ValueError('Duplicate metric row')
        all_keys.add(full_key)
        if int(r['bytes']) > int(r['budget']):
            raise ValueError('Byte cap exceeded')
        expected = config['field_object_bytes'] + (1 + 8 * int(r['splits'])) * config['node_bytes']
        if int(r['bytes']) != expected or int(r['allocation_oracle_calls']) != 0:
            raise ValueError('Storage or oracle accounting mismatch')
        if r['dataset'] == 'uniform':
            groups.setdefault(key, {})[r['method']] = r
    cases = read_csv(run / 'cases.csv')
    expected_keys = {(c['id'], s, b) for c in cases for s in config['seeds'] for b in config['budgets']}
    if set(groups) != expected_keys or len(rows) != len(expected_keys) * len(METHODS) * 4:
        raise ValueError('Missing cases or metrics')
    expected_full = {(*k, m, d) for k in expected_keys for m in METHODS for d in ['uniform', 'near', 'surface', 'probes']}
    if all_keys != expected_full:
        raise ValueError('Incomplete metric distributions')
    # Independently audit the selected sets using the telescoping gain identity.
    selections = read_csv(run / 'selections.csv')
    certificates = []
    diagnostics = []
    trees = {}
    for c in cases:
        for seed in config['seeds']:
            nodes = read_csv(run / 'candidates' / f'{c["id"]}_{seed}.csv')
            losses = [float(n['loss']) for n in nodes]
            children = [int(n['first_child']) for n in nodes]
            gains = {v: losses[v]-sum(losses[ch:ch+8]) for v, ch in enumerate(children) if ch}
            parents = {ch+i: v for v, ch in enumerate(children) if ch for i in range(8)}
            trees[c['id'], seed] = (losses, children, gains, parents)
            diagnostics.append(dict(case=c['id'], split=c['split'], seed=seed,
                internal_nodes=len(gains), negative_gain_nodes=sum(g < -1e-12 for g in gains.values()),
                increasing_gain_edges=sum(gains[v] > gains[p]+1e-12 for v, p in parents.items() if v in gains)))
    seen_selections = set()
    for row in selections:
        key = (row['case'], int(row['seed']), int(row['budget']))
        method = row['method']
        if (*key, method) in seen_selections:
            raise ValueError('Duplicate selection')
        seen_selections.add((*key, method))
        losses, children, gains, parents = trees[key[:2]]
        selected = set(int(v) for v in row['split_node_ids'].split(';') if v)
        if any(v not in gains or (v and parents[v] not in selected) for v in selected):
            raise ValueError('Invalid ancestor-closed split set')
        fit = losses[0] - sum(gains[v] for v in selected)
        result = groups[key][method]
        if len(selected) != int(result['splits']) or abs(fit-float(result['fit_sse'])) > 1e-9*(1+abs(fit)):
            raise ValueError('Independent telescoping SSE verification failed')
        minimum = min((gains[v] for v in selected), default=0.)
        maximum = max((g for v, g in gains.items() if v not in selected), default=0.)
        cap = min(len(gains), (key[2]-config['field_object_bytes']-config['node_bytes'])//(8*config['node_bytes']))
        certificate = (minimum >= -1e-12 and minimum >= maximum-1e-12 and
                       (len(selected) == cap or maximum <= 1e-12))
        certificates.append(dict(case=key[0], seed=key[1], budget=key[2], method=method,
                                 min_selected_gain=minimum, max_unselected_gain=maximum,
                                 top_k_certificate=int(certificate)))
    if seen_selections != {(*k, m) for k in expected_keys for m in METHODS}:
        raise ValueError('Missing selected split sets')
    write_csv(run / 'gain_diagnostics.csv', diagnostics)
    write_csv(run / 'certificates.csv', certificates)
    pairs = []
    for key, methods in groups.items():
        if set(methods) != set(METHODS):
            raise ValueError('Unpaired methods')
        opt = float(methods['exact']['fit_sse'])
        greedy = float(methods['greedy_best']['fit_sse'])
        for method, r in methods.items():
            fit = float(r['fit_sse'])
            if not math.isfinite(fit) or fit < opt - 1e-9 * (1 + opt):
                raise ValueError('Heuristic below exact optimum')
            pairs.append(dict(case=key[0], split=r['split'], seed=key[1], budget=key[2], method=method,
                              fit_sse=fit, optimum_sse=opt, gap_sse=max(0., fit-opt),
                              relative_gap=ratio(max(0., fit-opt), opt),
                              fraction_greedy_gap_closed=ratio(greedy-fit, greedy-opt),
                              allocation_seconds=float(r['allocation_seconds']),
                              total_construction_seconds=float(r['total_construction_seconds'])))
    write_csv(run / 'paired_gaps.csv', pairs)
    summaries = []
    for split in ['development', 'heldout']:
        for budget in config['budgets']:
            subset = [p for p in pairs if p['split'] == split and p['budget'] == budget]
            baseline = [p for p in subset if p['method'] == 'greedy_best']
            greedy_gap = sum(p['gap_sse'] for p in baseline)
            for method in METHODS:
                data = [p for p in subset if p['method'] == method]
                mean = lambda name: statistics.mean(p[name] for p in data)
                gap = sum(p['gap_sse'] for p in data)
                opt = sum(p['optimum_sse'] for p in data)
                selected = [r for r in rows if r['split'] == split and int(r['budget']) == budget and r['method'] == method]
                metric = lambda dataset, name: statistics.mean(float(r[name]) for r in selected if r['dataset'] == dataset)
                paired = [(p, groups[(p['case'], p['seed'], p['budget'])]['greedy_best']) for p in data]
                regressions = sum(p['fit_sse'] > float(g['fit_sse']) + 1e-10 for p, g in paired)
                summaries.append(dict(split=split, budget=budget, method=method, pairs=len(data),
                    mean_fit_sse=mean('fit_sse'), aggregate_relative_gap=ratio(gap, opt),
                    aggregate_fraction_gap_closed=ratio(greedy_gap-gap, greedy_gap),
                    worse_than_greedy_pairs=regressions,
                    mean_allocation_ms=1000*mean('allocation_seconds'),
                    mean_total_construction_ms=1000*mean('total_construction_seconds'),
                    uniform_rmse=metric('uniform', 'rmse'), near_rmse=metric('near', 'rmse'),
                    surface_residual_mae=metric('surface', 'mae'), near_sign_error=metric('near', 'sign_error'),
                    probe_inside_recall=metric('probes', 'inside_recall')))
    write_csv(run / 'summary.csv', summaries)
    (run / 'summary.json').write_text(json.dumps(summaries, indent=2, allow_nan=False) + '\n')
    costs = read_csv(run / 'candidate_costs.csv')
    if {(r['case'], int(r['seed'])) for r in costs} != {(c['id'], s) for c in cases for s in config['seeds']}:
        raise ValueError('Missing candidate costs')
    for r in costs:
        if int(r['candidate_oracle_calls']) != 8*int(r['nodes']):
            raise ValueError('Candidate oracle count mismatch')
    def display(value):
        if value is None:
            return 'undefined'
        return f'{value:.5g}' if isinstance(value, float) else str(value)
    columns = ['split', 'budget', 'method', 'aggregate_relative_gap', 'aggregate_fraction_gap_closed',
               'worse_than_greedy_pairs', 'mean_allocation_ms', 'mean_total_construction_ms',
               'uniform_rmse', 'near_rmse', 'surface_residual_mae', 'probe_inside_recall']
    table = '<table><thead><tr>' + ''.join(f'<th>{html.escape(c.replace("_", " "))}</th>' for c in columns) + '</tr></thead><tbody>'
    for row in summaries:
        table += '<tr>' + ''.join(f'<td>{html.escape(display(row[c]))}</td>' for c in columns) + '</tr>'
    table += '</tbody></table>'
    candidate_ms = statistics.mean(float(r['candidate_seconds']) for r in costs)*1000
    report = f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>SDF allocation study</title>
<style>body{{font:15px system-ui;margin:40px;color:#182b40;background:#fafbfc}}h1{{font-size:30px}}p{{max-width:1000px;line-height:1.6}}table{{border-collapse:collapse;font-size:12px}}th,td{{padding:8px;border-bottom:1px solid #cbd5e1;text-align:right}}th{{background:#e2e8f0;position:sticky;top:0}}tr:nth-child(even){{background:#f1f5f9}}td:nth-child(3){{text-align:left}}</style>
<h1>How much does greedy leave on the table?</h1>
<p>{len(cases)} frozen shape cases × 3 fitting seeds × 4 budgets; 15 cases are held out.
Exact means optimal fitting SSE within the fixed depth-{config['depth']} candidate tree.
Every method uses identical labels, float32 corners and packed-node accounting.
Candidate generation averages <b>{candidate_ms:.3f} ms</b> and is charged separately in candidate_costs.csv.
Allocation makes zero additional oracle calls. Total construction includes fitting labels, candidate generation, allocation and packing.</p>
<p>Relative gap and gap closure below are fractions, aggregated from summed paired SSE gaps.
Closure uses greedy_best; negative values mean regression. Undefined means a zero denominator.
Timing is the mean of three-repeat medians across shape/seed pairs, not an independent-machine study.
Surface residual is one-sided |predicted SDF| at true surface points, not a mesh distance guarantee.
Raw per-pair and distribution-specific data accompany this report.</p>{table}</html>'''
    (run / 'report.html').write_text(report, encoding='utf-8')
    source_paths = [ROOT/'CMakeLists.txt', ROOT/'fixtures/allocation_cases.csv', ROOT/'docs/allocation-study.md',
                    *sorted((ROOT/'src').glob('*.cpp')), *sorted((ROOT/'include/sdf').glob('*.hpp')),
                    ROOT/'tests/test_allocation.cpp', ROOT/'scripts/verify_allocation.py', Path(__file__).resolve()]
    # Preserve the exact source bytes behind hashes, including a dirty checkout.
    for p in source_paths:
        target = run/'source'/p.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and sha(target) != sha(p):
            raise ValueError('Source changed since snapshot; use a fresh run')
        if not target.exists():
            shutil.copyfile(p, target)
    provenance = dict(schema=1, platform=platform.platform(), machine=platform.machine(),
        python=platform.python_version(), git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        git_status=subprocess.check_output(['git', 'status', '--short'], cwd=ROOT, text=True),
        executable=dict(path=str(args.exe.resolve()), sha256=sha(args.exe)),
        source_sha256={str(p.relative_to(ROOT)).replace('\\', '/'): sha(p) for p in source_paths},
        artifact_sha256={str(p.relative_to(run)).replace('\\', '/'): sha(p) for p in sorted(run.rglob('*')) if p.is_file() and p.name != 'manifest.json'})
    (run/'manifest.json').write_text(json.dumps(provenance, indent=2) + '\n')
    if args.archive:
        args.archive.mkdir(parents=True)
        bundle = args.archive/'evidence.zip'
        with zipfile.ZipFile(bundle, 'x', zipfile.ZIP_DEFLATED, compresslevel=9) as z:
            for p in sorted(run.rglob('*')):
                if p.is_file():
                    z.write(p, p.relative_to(run).as_posix())
        for name in ['report.html', 'summary.csv', 'summary.json', 'paired_gaps.csv', 'candidate_costs.csv',
                     'gain_diagnostics.csv', 'certificates.csv', 'config.json', 'build.json']:
            shutil.copyfile(run/name, args.archive/name)
        (args.archive/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(args.archive.iterdir()) if p.is_file()))
    print(f'Validated {len(groups)} paired allocations and {len(rows)} metric rows; report: {run / "report.html"}')


if __name__ == '__main__':
    main()
