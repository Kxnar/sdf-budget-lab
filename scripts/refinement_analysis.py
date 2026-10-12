"""Independent verification of selective q8 refinement and frozen recovery."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from tail_evidence import rows, write_csv, hashes, sha, close
from distortion_analysis import SHIFT, decode, oracle, layout, intervals
from refinement_verify_selection import verify_selection


def analyze(directory, out):
    if out.exists():
        raise ValueError('Analysis output must be new')
    assert (directory/'COMPLETE').exists() and (directory/'EXIT_CODE').read_text() == '0'
    hashes(directory)
    hashes(directory.with_name(directory.name+'-source'))
    out.mkdir(parents=True)
    config = json.loads((directory/'config.json').read_text())
    metrics = rows(directory/'metrics.csv')
    fields = [r for r in metrics if r['dataset'] == 'probes']
    assert len(fields) == (4 if config['smoke'] else 20)
    summaries, coverage, profiles, widths = [], [], [], []
    sample_cache = {}
    verified = 0
    for row in fields:
        case, shape = row['case'], row['shape']
        feature, depth = float(row['feature']), int(row['depth'])
        method = row['method']
        name = f'{case}_{depth}_{method}'
        if case not in sample_cache:
            sample_cache.clear()
            sample_cache[case] = {s: np.loadtxt(directory/'samples'/f'{case}_{s}.csv', delimiter=',', skiprows=1, ndmin=2) for s in ('fit', 'uniform', 'near', 'probes')}
        samples = sample_cache[case]
        data = (directory/'fields'/f'{name}.bin').read_bytes()
        layout(data, method, depth)
        assert len(data)+config['field_object_bytes'] == int(row['bytes'])
        if method in ('mse_q8', 'selective_q8'):
            assert int(row['bytes']) <= 32768
        for ds in ('fit', 'uniform', 'near', 'probes'):
            sample = samples[ds]
            pred = decode(data, sample[:, :3])
            cpp = np.fromfile(directory/'predictions'/f'{name}_{ds}.bin', dtype='<f8')
            np.testing.assert_allclose(pred, cpp, rtol=2e-12, atol=2e-14)
            np.testing.assert_allclose(sample[:, 3], oracle(sample[:, :3], shape, feature), rtol=2e-12, atol=2e-14)
            errors = np.abs(pred-sample[:, 3])
            m = next(r for r in metrics if r['case'] == case and int(r['depth']) == depth and r['method'] == method and r['dataset'] == ds)
            assert len(sample) == int(m['count'])
            close(errors.mean(), m['mae'], 'MAE')
            close(np.sort(errors)[-math.ceil(.01*len(errors)):].mean(), m['cvar99'], 'CVaR99')
            for sign, key in ((-1, 'inside_recall'), (1, 'outside_recall')):
                mask = sample[:, 3]*sign > 2e-9
                close(np.mean(pred[mask]*sign > 0), m[key], key)
            valid = np.abs(sample[:, 3]) > 2e-9
            close(np.mean(pred[valid]*np.sign(sample[valid, 3]) <= 0), m['sign_error'], 'sign error')
            verified += len(pred)
        allocated = next(r for r in rows(directory/'allocation.csv') if r['case']==case and r['method']==method)
        fitting_prediction=decode(data, samples['fit'][:, :3])
        close(np.sum(((fitting_prediction-samples['fit'][:, 3])/2)**2), allocated['fit_sse'], 'allocated SSE')
        # Exact frozen probes, routed to the selected leaf, compared with fitting pool.
        probes, fit = samples['probes'], samples['fit']
        pred, nodes, depths = decode(data, probes[:, :3], True)
        _, fit_nodes, _ = decode(data, fit[:, :3], True)
        failed = pred*np.sign(probes[:, 3]) <= 0
        counts, same_counts = [], []
        for i in range(len(probes)):
            same_leaf = fit_nodes == nodes[i]
            count = int(same_leaf.sum())
            correct_sign = int(np.sum(same_leaf & (fit[:, 3]*np.sign(probes[i, 3]) > 2e-9)))
            counts.append(count); same_counts.append(correct_sign)
            coverage.append(dict(case=case, depth=depth, method=method, probe=i, x=probes[i, 0], y=probes[i, 1], z=probes[i, 2], truth=probes[i, 3], prediction=pred[i], failed=int(failed[i]), leaf_depth=depths[i], fit_count=count, same_sign_fit_count=correct_sign))
        summaries.append(dict(case=case, shape=shape, feature=feature, seed=row['seed'], depth=depth, method=method, failed_probes=int(failed.sum()), failing_with_fit=int(np.sum(failed & (np.array(counts)>0))), failing_with_same_sign_fit=int(np.sum(failed & (np.array(same_counts)>0))), failing_at_max_depth=int(np.sum(failed & (depths==depth)))))
        if shape == 'torus':
            u = np.arange(1024)*2*np.pi/1024
            points = SHIFT + np.column_stack([.55*np.cos(u), .55*np.sin(u), np.zeros_like(u)])
        else:
            u = np.linspace(-feature, feature, 1025)
            points = SHIFT + np.column_stack([u, np.zeros_like(u), np.zeros_like(u)])
        truth, pred = oracle(points, shape, feature), decode(data, points)
        quant = rounding = interpolation = np.full(len(points), np.nan)
        for i in range(len(points)):
            profiles.append(dict(case=case, shape=shape, feature=feature, depth=depth, method=method, coordinate=u[i], truth=truth[i], prediction=pred[i], signed_error=pred[i]-truth[i], quantization=quant[i], f32_rounding=rounding[i], interpolation=interpolation[i]))
        if shape == 'torus':
            a = np.arange(256)*2*np.pi/256
            radial = np.column_stack([np.cos(a), np.sin(a), np.zeros_like(a)])
            centers = SHIFT + .55*radial
            axes = [('radial', radial), ('vertical', np.tile([0., 0., 1.], (256, 1)))]
        else:
            a, centers = np.array([0.]), SHIFT[None, :]
            axes = [('gap', np.array([[1., 0., 0.]]))]
        for axis, directions in axes:
            if depth == 5 and method == 'mse_q8':
                ol, ore = intervals(lambda p: oracle(p, shape, feature), centers, directions, feature, shape != 'torus')
                np.testing.assert_allclose(ore-ol, feature if shape != 'torus' else 2*feature, atol=1e-12)
            left, right = intervals(data, centers, directions, feature, shape != 'torus')
            for i in range(len(a)):
                widths.append(dict(case=case, shape=shape, feature=feature, depth=depth, method=method, axis=axis, angle=a[i], left=left[i], right=right[i], width=right[i]-left[i], expected=feature if shape != 'torus' else 2*feature))
        print(name, 'verified and localized', flush=True)
    for filename, values in [('coverage.csv', coverage), ('coverage_summary.csv', summaries), ('profiles.csv', profiles), ('widths.csv', widths)]:
        write_csv(out/filename, values)
    verify_selection(directory)
    # Frozen recovery across both seeds, final resolutions, oracle and larger features.
    surface = rows(directory/'surfaces.csv')
    expected_cases = {(shape, feature, str(seed)) for shape, features in [('torus', [.035, .06, .1]), ('double_sphere', [.06, .1])] for feature in features for seed in [40, 41]}
    if not config['smoke']:
        assert {(r['shape'], float(r['feature']), r['seed']) for r in fields} == expected_cases
        assert len(surface) == 90
    for r in surface:
        if r['method'] == 'oracle':
            probe_pass = True
        else:
            m = next(m for m in fields if m['case']==r['case'] and m['method']==r['method'])
            probe_pass = float(m['inside_recall'])>=.95 and float(m['outside_recall'])>=.95
        assert bool(int(r['probe_pass'])) == probe_pass
        expected_topology = (1, 0) if r['shape']=='torus' else (2, 4)
        topology = (int(r['components']), int(r['euler'])) == expected_topology and all(int(r[k])==0 for k in ['boundary_edges', 'nonmanifold_edges', 'nonmanifold_vertices'])
        assert bool(int(r['topology_pass'])) == topology
        assert bool(int(r['success'])) == (topology and probe_pass)
    recovery = []
    for shape in ('torus', 'double_sphere'):
        features = sorted({float(r['feature']) for r in fields if r['shape'] == shape})
        for method, depth in [('mse_q8', 5), ('selective_q8', 6)]:
            passes = {}
            for f in features:
                selected = [r for r in surface if r['shape'] == shape and float(r['feature']) == f and int(r['cells']) in ((16, 24) if config['smoke'] else (64, 96)) and (r['method'] == 'oracle' or (r['method'] == method and int(r['depth']) == depth))]
                assert len(selected) == (4 if config['smoke'] else 8)
                passes[f] = all(int(r['success']) for r in selected)
            threshold = next((f for f in features if all(passes[g] for g in features if g>=f)), None)
            recovery.append(dict(shape=shape, method=method, depth=depth, minimum_tested_recovered_feature=threshold))
    write_csv(out/'recovery.csv', recovery)
    (out/'VERIFIED.json').write_text(json.dumps(dict(fields=len(fields), cpp_queries=verified, raw_manifest_sha256=sha(directory/'SHA256SUMS.json'), source_manifest_sha256=sha(directory.with_name(directory.name+'-source')/'SHA256SUMS.json'), verifier_sha256=sha(Path(__file__))), indent=2))
    print(f'Verified {len(fields)} fields and {verified} C++ queries', flush=True)
    (out/'SHA256SUMS.json').write_text(json.dumps({f.relative_to(out).as_posix(): sha(f) for f in out.rglob('*') if f.is_file()}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('directory', type=Path)
    p.add_argument('out', type=Path)
    a = p.parse_args()
    analyze(a.directory, a.out)
