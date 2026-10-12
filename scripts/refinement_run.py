"""Freeze source/executable before a fresh refinement experiment; retain failures."""
import argparse
import json
from pathlib import Path
import platform
import shutil
import subprocess
from tail_run import sha, ROOT


def main():
    p = argparse.ArgumentParser()
    p.add_argument('out', type=Path)
    p.add_argument('--exe', type=Path, required=True)
    p.add_argument('--smoke', action='store_true')
    p.add_argument('--development', action='store_true')
    p.add_argument('--selection', type=Path)
    a = p.parse_args()
    snapshot = a.out.with_name(a.out.name + '-source')
    if a.out.exists() or snapshot.exists():
        raise ValueError('Use fresh output/source directories')
    if not a.development and not a.selection:
        raise ValueError('Final execution requires frozen development selection')
    snapshot.mkdir(parents=True)
    files = [ROOT / 'CMakeLists.txt', ROOT / 'docs/refinement-protocol.md']
    for folder in ('src', 'include', 'tests', 'scripts'):
        files += [f for f in (ROOT / folder).rglob('*') if f.suffix in ('.cpp', '.hpp', '.py')]
    for src in files:
        dst = snapshot / src.relative_to(ROOT)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    shutil.copy2(a.exe, snapshot / a.exe.name)
    command = [str(a.exe.resolve()), '--out', str(a.out.resolve())]
    if a.development:
        command += ['--development']
    else:
        frozen = json.loads(a.selection.read_text())
        assert sha(Path(frozen['development_csv'])) == frozen['development_sha256']
        command += ['--k', str(frozen['k'])]
        shutil.copy2(a.selection, snapshot / 'selection.json')
    if a.smoke:
        command += ['--smoke']
    (snapshot / 'run.json').write_text(json.dumps(dict(command=command, platform=platform.platform()), indent=2))
    (snapshot / 'SHA256SUMS.json').write_text(json.dumps({f.relative_to(snapshot).as_posix(): sha(f) for f in snapshot.rglob('*') if f.is_file()}, indent=2))
    result = subprocess.run(command, check=False)
    if a.out.exists():
        (a.out / 'EXIT_CODE').write_text(str(result.returncode))
        (a.out / 'SHA256SUMS.json').write_text(json.dumps({f.relative_to(a.out).as_posix(): sha(f) for f in a.out.rglob('*') if f.is_file()}, indent=2))
    raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
