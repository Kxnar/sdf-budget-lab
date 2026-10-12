"""Apply the predeclared fitting-only development selection rule."""
import argparse
import json
from pathlib import Path
from tail_evidence import rows, hashes, sha, write_csv


def main():
    p = argparse.ArgumentParser()
    p.add_argument('development', type=Path)
    p.add_argument('out', type=Path)
    a = p.parse_args()
    if a.out.exists():
        raise ValueError('Use new selection directory')
    hashes(a.development)
    hashes(a.development.with_name(a.development.name+'-source'))
    assert (a.development/'COMPLETE').exists()
    assert (a.development/'EXIT_CODE').read_text() == '0'
    data = rows(a.development/'development.csv')
    assert len(data) == 56
    keys = {(r['shape'], r['feature'], r['seed']) for r in data}
    assert len(keys) == 8
    comparisons = []
    for k in (32, 64, 128, 256, 512, 1024):
        excess = []
        for key in keys:
            case = [r for r in data if (r['shape'], r['feature'], r['seed']) == key]
            case.sort(key=lambda r: int(r['k']))
            assert [int(r['k']) for r in case] == [0, 32, 64, 128, 256, 512, 1024]
            assert all(float(b['fit_sse']) <= float(c['fit_sse'])+1e-12 for c, b in zip(case, case[1:]))
            chosen = next(r for r in case if int(r['k']) == k)
            reference = float(case[-1]['fit_sse'])
            excess.append(float(chosen['fit_sse'])/reference-1 if reference else 0.)
        comparisons.append(dict(k=k, mean_excess=sum(excess)/len(excess), max_excess=max(excess), qualifies=int(sum(excess)/len(excess)<=.01 and max(excess)<=.05)))
    k = next(r['k'] for r in comparisons if r['qualifies'])
    a.out.mkdir(parents=True)
    write_csv(a.out/'comparison.csv', comparisons)
    frozen = dict(k=k, development_csv=str((a.development/'development.csv').resolve()), development_sha256=sha(a.development/'development.csv'), rule='smallest K with mean SSE excess <=1% and every case <=5% versus K=1024')
    (a.out/'selection.json').write_text(json.dumps(frozen, indent=2))
    print(json.dumps(dict(selection=frozen, comparisons=comparisons), indent=2))


if __name__ == '__main__':
    main()
