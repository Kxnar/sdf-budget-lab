"""Run one frozen stage in a new directory, retaining partial failures and hashes."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("out", type=Path)
    p.add_argument("--exe", type=Path, required=True)
    p.add_argument("--stage", choices=["resolution", "development", "evaluation"], required=True)
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--grid", type=int)
    a = p.parse_args()
    snapshot = a.out.with_name(a.out.name + "-source")
    if a.out.exists() or snapshot.exists():
        raise ValueError("Use new output and source directories")
    snapshot.mkdir(parents=True)
    files = [ROOT / "CMakeLists.txt", ROOT / "docs/tail-protocol.md"]
    for folder in ("include", "src", "tests", "scripts"):
        files += [f for f in (ROOT / folder).rglob("*") if f.suffix in (".cpp", ".hpp", ".py")]
    for src in files:
        dst = snapshot / src.relative_to(ROOT)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    shutil.copy2(a.exe, snapshot / a.exe.name)
    command = [str(a.exe.resolve()), "--out", str(a.out.resolve()), "--stage", a.stage]
    if a.smoke:
        command.append("--smoke")
    if a.grid is not None:
        command += ["--grid", str(a.grid)]
    (snapshot / "run.json").write_text(json.dumps(dict(command=command, platform=platform.platform()), indent=2))
    (snapshot / "SHA256SUMS.json").write_text(json.dumps({f.relative_to(snapshot).as_posix(): sha(f) for f in snapshot.rglob("*") if f.is_file()}, indent=2))
    result = subprocess.run(command, check=False)
    if a.out.exists():
        (a.out / "EXIT_CODE").write_text(str(result.returncode))
        (a.out / "SHA256SUMS.json").write_text(json.dumps({f.relative_to(a.out).as_posix(): sha(f) for f in a.out.rglob("*") if f.is_file()}, indent=2))
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
