"""Snapshot and execute a new, immutable fidelity experiment."""
import argparse
import json
import platform
import shutil
import subprocess
from pathlib import Path
from tail_run import ROOT, sha

p = argparse.ArgumentParser()
p.add_argument('out', type=Path)
p.add_argument('--exe', type=Path, default=Path('build-core/Release/sdf_fidelity_bench.exe'))
p.add_argument('--k', type=int, default=32)
p.add_argument('--penalty', type=float, default=0)
p.add_argument('--seed', type=int, default=40)
p.add_argument('--quick', action='store_true')
p.add_argument('--quadrature', type=int, default=0)
p.add_argument('--frozen', type=Path)
a = p.parse_args()
snapshot = a.out.with_name(a.out.name+'-source')
if a.out.exists() or snapshot.exists():
    raise ValueError('Use new directories')
if a.seed != 40:
    assert a.frozen, 'Fresh evaluation requires frozen selection'
if a.frozen:
    frozen = json.loads(a.frozen.read_text())
    assert (a.k, a.penalty) == (frozen['k'], frozen['lambda'])
    assert a.quadrature == frozen.get('quadrature',0)
    for file, expected in frozen.get('source_sha256',{}).items():
        assert sha(ROOT/file)==expected, f'Frozen source changed: {file}'
snapshot.mkdir(parents=True)
files = [ROOT/'CMakeLists.txt', ROOT/'docs/fidelity-log.md']
for folder in ('src', 'include', 'tests', 'scripts'):
    files += [f for f in (ROOT/folder).rglob('*') if f.suffix in ('.cpp','.hpp','.py')]
files += [f for f in (ROOT/'fixtures').rglob('*') if f.is_file()]
for src in files:
    dst = snapshot/src.relative_to(ROOT)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
shutil.copy2(a.exe, snapshot/a.exe.name)
if a.frozen:
    shutil.copy2(a.frozen, snapshot/'frozen.json')
command = [str(a.exe.resolve()), '--out', str(a.out.resolve()), '--k', str(a.k), '--lambda', str(a.penalty), '--seed-start', str(a.seed)]
command += ['--quadrature', str(a.quadrature)]
if a.quick:
    command += ['--quick']
(snapshot/'run.json').write_text(json.dumps(dict(command=command, platform=platform.platform()), indent=2))
def manifest(folder):
    (folder/'SHA256SUMS.json').write_text(json.dumps({f.relative_to(folder).as_posix():sha(f) for f in folder.rglob('*') if f.is_file()}, indent=2))
manifest(snapshot)
result = subprocess.run(command, check=False)
if a.out.exists():
    (a.out/'EXIT_CODE').write_text(str(result.returncode))
    manifest(a.out)
raise SystemExit(result.returncode)
