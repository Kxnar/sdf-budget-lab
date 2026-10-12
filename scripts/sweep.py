"""Run a fixed development matrix and fingerprint every input; never overwrite runs."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys


def fingerprint(root):
    digest = hashlib.sha256()
    for pattern in ("src/*.cpp", "include/sdf/*.hpp", "CMakeLists.txt", "scripts/*.py", "fixtures/*"):
        for path in sorted(root.glob(pattern)):
            digest.update(path.relative_to(root).as_posix().encode())
            # Git checkouts on Windows and Linux can use different line endings.
            digest.update((path.read_text(encoding="utf-8").rstrip() + "\n").encode("utf-8"))
    return digest.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--features", default="0.04")
    p.add_argument("--shapes", default="sphere,torus,thin_box,double_sphere")
    p.add_argument("--seeds", default="0,1,2")
    p.add_argument("--budgets", default="8192,32768,131072")
    p.add_argument("--train", type=int, default=8192)
    p.add_argument("--test", type=int, default=8192)
    p.add_argument("--slice-size", type=int, default=257)
    p.add_argument("--mesh", type=Path)
    p.add_argument("--probes", type=Path)
    p.add_argument("--slice-plane", choices=("xy", "xz"))
    args = p.parse_args()
    exe = args.exe.resolve()
    if not exe.is_file():
        p.error(f"Executable does not exist: {exe}")
    # Refuse before any writes or subprocess calls.
    if args.out.exists():
        p.error("Output already exists; choose a new directory")
    args.out.mkdir(parents=True)
    root = Path(__file__).resolve().parents[1]
    cpu = platform.processor()
    if sys.platform == "win32":
        cpu_query = subprocess.run(["powershell", "-NoProfile", "-Command",
                                    "(Get-CimInstance Win32_Processor).Name"],
                                   capture_output=True, text=True, check=True)
        cpu = cpu_query.stdout.strip()
    manifest = {"phase": "development", "source_sha256": fingerprint(root),
                "executable_sha256": hashlib.sha256(exe.read_bytes()).hexdigest(),
                "executable": str(exe), "python": sys.version, "platform": platform.platform(),
                "cpu": cpu, "arguments": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                "commands": []}
    manifest["inputs_sha256"] = {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
                                  for path in (args.mesh, args.probes) if path}
    for shape in (["mesh"] if args.mesh else args.shapes.split(",")):
        for feature in args.features.split(","):
            if shape == "sphere" and feature != args.features.split(",")[0]:
                continue  # Sphere radius does not depend on feature.
            target = args.out / f"{shape}-f{feature}"
            command = [str(exe), "--shape", shape, "--feature", feature, "--out", str(target),
                       "--seeds", args.seeds, "--budgets", args.budgets,
                       "--train", str(args.train), "--test", str(args.test),
                       "--slice-size", str(args.slice_size)]
            if args.mesh:
                command += ["--mesh", str(args.mesh.resolve())]
            if args.probes:
                command += ["--probes", str(args.probes.resolve())]
            if args.slice_plane:
                command += ["--slice-plane", args.slice_plane]
            manifest["commands"].append(command)
            (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
            subprocess.run(command, check=True)
    (args.out / "complete.json").write_text('{"status":"complete"}\n', encoding="utf-8")


if __name__ == "__main__":
    main()
