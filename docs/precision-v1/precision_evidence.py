"""Freeze source before running, verify packed evidence, and produce a compact report.

Standard library only. Use --run once; later invocations verify the archived run.
"""
import argparse
import csv
import hashlib
import html
import json
import math
from pathlib import Path
import shutil
import statistics
import struct
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def write_csv(path, values):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


def query(data, point):
    lo, hi, node = [-1.] * 3, [1.] * 3, 0
    while True:
        descriptor, = struct.unpack_from("<I", data, node)
        encoding, offset = descriptor & 3, descriptor & ~3
        if encoding:
            if encoding == 1:
                values = struct.unpack_from("<8f", data, offset)
            else:
                base, scale = struct.unpack_from("<2f", data, offset)
                integers = struct.unpack_from("<8H" if encoding == 2 else "<8B", data, offset + 8)
                values = [struct.unpack("<f", struct.pack("<f", base + scale * q))[0] for q in integers]
            t = [(point[i] - lo[i]) / (hi[i] - lo[i]) for i in range(3)]
            return sum(values[c] * math.prod(t[i] if c & (1 << i) else 1 - t[i] for i in range(3)) for c in range(8))
        mid = [(lo[i] + hi[i]) / 2 for i in range(3)]
        c = sum((1 << i) for i in range(3) if point[i] >= mid[i])
        for i in range(3):
            if c & (1 << i):
                lo[i] = mid[i]
            else:
                hi[i] = mid[i]
        node = offset + 4 * c


def verify_buffer(data, expected, candidate, choices):
    occupied = set()
    leaves = [0, 0, 0]
    branches = 0
    visited = set()

    def claim(offset, size):
        assert offset % 4 == 0 and 0 <= offset <= len(data) - size
        region = set(range(offset, offset + size, 4))
        assert not occupied.intersection(region), "Overlapping packed regions"
        occupied.update(region)

    def visit(node, v):
        nonlocal branches
        visited.add(v)
        claim(node, 4)
        descriptor, = struct.unpack_from("<I", data, node)
        encoding, offset = descriptor & 3, descriptor & ~3
        assert choices[v] == encoding
        if encoding:
            leaves[encoding - 1] += 1
            claim(offset, [0, 32, 24, 16][encoding])
        else:
            branches += 1
            first = int(candidate[v]["first_child"])
            assert first > v
            for c in range(8):
                visit(offset + 4 * c, first + c)

    visit(0, 0)
    assert visited == set(choices)
    assert len(occupied) * 4 == len(data), "Unaccounted padding/trailing storage"
    assert branches == int(expected["branches"])
    assert leaves == [int(expected[k]) for k in ("f32_leaves", "q16_leaves", "q8_leaves")]


def verify(directory):
    assert (directory / "COMPLETE").exists(), "Incomplete run"
    config = json.loads((directory / "config.json").read_text())
    metrics = rows(directory / "metrics.csv")
    fields = [r for r in metrics if r["dataset"] == "near"]
    expected_cases = 1 if config["smoke"] else 24
    assert len(fields) == expected_cases * len(config["budgets"]) * 4
    assert len(metrics) == len(fields) * 3
    assert len({(r["case"], r["budget"], r["method"], r["dataset"]) for r in metrics}) == len(metrics)
    cache = {}
    fits = {}
    grouped = {}
    for r in fields:
        case, method, budget = r["case"], r["method"], int(r["budget"])
        prefix = f"{case}_{budget}_{method}"
        if case not in cache:
            cache[case] = rows(directory / "candidates" / f"{case}.csv")
        candidate = cache[case]
        if case not in fits:
            fits[case] = [tuple(float(s[k]) for k in ("x", "y", "z", "distance")) for s in rows(directory / "samples" / f"{case}_fit.csv")]
        choices = {int(s["node"]): int(s["encoding"]) for s in rows(directory / "choices" / f"{prefix}.csv")}
        if method != "joint":
            assert all(e in (0, {"f32": 1, "q16": 2, "q8": 3}[method]) for e in choices.values())
        data = (directory / "fields" / f"{prefix}.bin").read_bytes()
        assert len(data) + config["field_object_bytes"] == int(r["bytes"]) <= budget
        verify_buffer(data, r, candidate, choices)
        selected_loss = sum(float(candidate[v][["", "f32_loss", "q16_loss", "q8_loss"][e]]) for v, e in choices.items() if e)
        packed_loss = sum(((query(data, s[:3]) - s[3]) / 2) ** 2 for s in fits[case])
        assert math.isclose(selected_loss, float(r["fit_sse"]), rel_tol=1e-9, abs_tol=1e-10)
        assert math.isclose(packed_loss, float(r["fit_sse"]), rel_tol=1e-9, abs_tol=1e-10), prefix
        grouped.setdefault((case, budget), {})[method] = float(r["fit_sse"])
    for values in grouped.values():
        assert set(values) == {"f32", "q16", "q8", "joint"}
        assert values["joint"] <= min(values.values()) + 1e-10
    surfaces = rows(directory / "surfaces.csv")
    assert len(surfaces) == (len(fields) + expected_cases) * len(config["resolutions"])
    assert len({(r["case"], r["budget"], r["method"], r["cells"]) for r in surfaces}) == len(surfaces)
    assert len(rows(directory / "timings.csv")) == 3 * len(fields)
    assert len(rows(directory / "construction.csv")) == 4 * expected_cases
    source_manifest = directory / "source/source_hashes.json"
    if source_manifest.exists():
        for name, digest in json.loads(source_manifest.read_text()).items():
            assert sha(directory / "source" / name) == digest, f"Frozen source changed: {name}"
    manifest = directory / "SHA256SUMS.json"
    if manifest.exists():
        for name, digest in json.loads(manifest.read_text()).items():
            assert sha(directory / name) == digest, f"Hash mismatch: {name}"
    print(f"Verified {len(fields)} fields: byte coverage, choices, independent packed SSE, caps, joint inequalities and evidence hashes.", flush=True)
    return config, metrics, surfaces


def report(directory, config, metrics, surfaces):
    timings = rows(directory / "timings.csv")
    costs = rows(directory / "construction.csv")
    near = [r for r in metrics if r["dataset"] == "near"]
    probes = {(r["case"], r["budget"], r["method"]): r for r in metrics if r["dataset"] == "probes"}
    uniform = {(r["case"], r["budget"], r["method"]): r for r in metrics if r["dataset"] == "uniform"}
    summary = []
    for budget in config["budgets"]:
        for method in ("f32", "q16", "q8", "joint"):
            selected = [r for r in near if int(r["budget"]) == budget and r["method"] == method]
            keys = [(r["case"], r["budget"], r["method"]) for r in selected]
            time = [r for r in timings if int(r["budget"]) == budget and r["method"] == method]
            cost = [r for r in costs if r["method"] == method]
            per_field = [[r for r in time if r["case"] == key[0]] for key in keys]
            summary.append(dict(budget=budget, method=method, cases=len(selected),
                mean_bytes=statistics.mean(float(r["bytes"]) for r in selected),
                uniform_rmse=statistics.mean(float(uniform[k]["rmse"]) for k in keys),
                near_rmse=statistics.mean(float(r["rmse"]) for r in selected),
                near_mae=statistics.mean(float(r["mae"]) for r in selected),
                near_cvar99=statistics.mean(float(r["cvar99"]) for r in selected),
                inside_recall=statistics.mean(float(probes[k]["inside_recall"]) for k in keys),
                outside_recall=statistics.mean(float(probes[k]["outside_recall"]) for k in keys),
                query_ns=statistics.mean(statistics.median(float(r["packed_ns"]) for r in group) for group in per_field),
                same_topology_f32_ns=statistics.mean(statistics.median(float(r["topology_f32_ns"]) for r in group) for group in per_field),
                dp_seconds=statistics.mean(float(r["dp_seconds"]) for r in cost),
                candidate_seconds=statistics.mean(float(r["candidate_seconds"]) for r in cost),
                pack_seconds=statistics.mean(float(group[0]["backtrack_pack_seconds"]) for group in per_field)))
    write_csv(directory / "summary.csv", summary)
    crossings = []
    for case in sorted({r["case"] for r in near}):
        for method in ("f32", "q16", "q8", "joint"):
            selected = sorted((r for r in near if r["case"] == case and r["method"] == method), key=lambda r: int(r["budget"]))
            rmse = [int(r["budget"]) for r in selected if float(r["rmse"]) <= .02]
            accuracy = [int(r["budget"]) for r in selected if float(r["rmse"]) <= .02 and all(float(probes[(case, r["budget"], method)][k]) >= .95 for k in ("inside_recall", "outside_recall"))]
            crossings.append(dict(case=case, method=method, rmse_budget=min(rmse) if rmse else "not reached", accuracy_budget=min(accuracy) if accuracy else "not reached"))
    write_csv(directory / "accuracy_crossings.csv", crossings)
    recovery = []
    convergence = []
    for shape in sorted({r["shape"] for r in near}):
        for budget in config["budgets"]:
            for method in ("f32", "q16", "q8", "joint"):
                success = []
                for feature in config["features"]:
                    selected = [r for r in surfaces if r["shape"] == shape and math.isclose(float(r["feature"]), feature) and ((int(r["budget"]) == budget and r["method"] == method) or r["method"] == "oracle")]
                    final = [r for r in selected if int(r["cells"]) in config["resolutions"][-2:]]
                    passed = len(final) == config["seeds"] * 4 and all(r["success"] == "1" for r in final)
                    success.append(passed)
                    convergence.append(dict(shape=shape, budget=budget, method=method, feature=feature, stable_success=int(passed), successes=";".join(f'{r["seed"]}:{r["method"]}:{r["cells"]}={r["success"]}' for r in selected)))
                minimum = next((f for i, f in enumerate(config["features"]) if all(success[i:])), None)
                recovery.append(dict(shape=shape, budget=budget, method=method, minimum_feature=minimum if minimum is not None else "not reached", physical_size=minimum * (1 if shape == "double_sphere" else 2) if minimum is not None else "not reached"))
    write_csv(directory / "recovery.csv", recovery)
    write_csv(directory / "convergence.csv", convergence)
    def table(data):
        headers = list(data[0])
        def fmt(x):
            return f"{x:.5g}" if isinstance(x, float) else str(x)
        return "<table><tr>" + "".join(f"<th>{html.escape(k)}</th>" for k in headers) + "</tr>" + "".join("<tr>" + "".join(f"<td>{html.escape(fmt(r[k]))}</td>" for k in headers) + "</tr>" for r in data) + "</table>"
    text = """<!doctype html><html><meta charset="utf-8"><title>Joint precision study</title>
<style>body{font:15px system-ui;margin:40px;color:#203044}table{border-collapse:collapse;font-size:12px;margin:24px 0}td,th{border:1px solid #ccd4df;padding:6px}th{background:#e8eff7}p{max-width:1000px}</style>
<h1>Joint spatial and numerical precision</h1>
<p>Exact fitting SSE tree knapsack; every baseline uses the same compact layout. Finite analytic sweep, two seeds, no universal topology guarantee. Query latency includes decoding. Same-topology float32 timing isolates the combined codec/payload-layout effect, but uses more memory. DP cost is shared across the three budgets. Means weight cases equally; timing uses within-field batch medians.</p>
<p><a href="summary.csv">Summary CSV</a> · <a href="recovery.csv">Minimum recovered features</a> · <a href="accuracy_crossings.csv">Accuracy budget crossings</a> · <a href="convergence.csv">Extraction convergence</a> · <a href="evidence.zip">Raw evidence and source archive</a> · Mesh archives: <a href="meshes-thin_box.zip">walls</a>, <a href="meshes-torus.zip">tori</a>, <a href="meshes-double_sphere.zip">sphere pairs</a></p>"""
    text += "<h2>Equal-byte results</h2>" + table(summary)
    text += "<h2>Minimum recoverable feature</h2><p>Success at both final extraction resolutions, both seeds, and every larger tested feature. Oracle controls must also succeed. Physical size is full wall thickness, tube diameter, or sphere gap. Not reached includes failures caused by the fitting objective.</p>" + table(recovery)
    text += "<h2>Fixed-accuracy crossings</h2><p>Smallest tested budget with near RMSE ≤0.02; combined accuracy additionally requires inside and outside probe recall ≥95%. These are discrete observed crossings, not exact storage minima.</p>" + table(crossings)
    text += "<p>Protocol and full source are frozen in the archive. Mesh extraction samples a potentially discontinuous field onto a conforming lattice; reported mesh topology is resolution-dependent.</p></html>"
    (directory / "report.html").write_text(text, encoding="utf-8")
    return summary


def verify_archives(directory):
    outer = json.loads((directory / "SHA256SUMS.json").read_text())
    for name, digest in outer.items():
        assert sha(directory / name) == digest, f"Publication hash mismatch: {name}"
    with zipfile.ZipFile(directory / "evidence.zip") as archive:
        manifest = json.loads(archive.read("SHA256SUMS.json"))
    found = set()
    for path in sorted(directory.glob("*.zip")):
        with zipfile.ZipFile(path) as archive:
            assert archive.testzip() is None
            for name in archive.namelist():
                if name == "SHA256SUMS.json":
                    continue
                assert name not in found, f"Duplicate archive member: {name}"
                found.add(name)
                assert hashlib.sha256(archive.read(name)).hexdigest() == manifest[name], name
    assert found == set(manifest), "Missing archived evidence"
    print(f"Verified {len(found)} archived files, CRCs, and SHA-256 hashes.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--exe", type=Path)
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--verify-archives", action="store_true")
    args = parser.parse_args()
    directory = args.directory
    if args.verify_archives:
        verify_archives(directory)
        return
    if args.run:
        assert args.exe and args.exe.is_file()
        assert not directory.exists()
        snapshot = directory.with_name(directory.name + "-source")
        snapshot.mkdir(parents=True, exist_ok=False)
        paths = [ROOT / "CMakeLists.txt", ROOT / "docs/precision-protocol.md"]
        for folder in ("src", "include", "tests", "scripts"):
            paths.extend(p for p in (ROOT / folder).rglob("*") if p.is_file() and "__pycache__" not in p.parts)
        for p in paths:
            target = snapshot / p.relative_to(ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(p, target)
        shutil.copyfile(args.exe, snapshot / args.exe.name)
        command = [str(args.exe.resolve()), "--out", str(directory)] + (["--smoke"] if args.smoke else [])
        (snapshot / "invocation.json").write_text(json.dumps({"command": command}, indent=2))
        (snapshot / "source_hashes.json").write_text(json.dumps({str(p.relative_to(snapshot)).replace("\\", "/"): sha(p) for p in snapshot.rglob("*") if p.is_file()}, indent=2))
        result = subprocess.run(command, cwd=ROOT)
        if directory.exists():
            shutil.copytree(snapshot, directory / "source")
        result.check_returncode()
    config, metrics, surfaces = verify(directory)
    summary = report(directory, config, metrics, surfaces)
    hashes = {p.relative_to(directory).as_posix(): sha(p) for p in directory.rglob("*") if p.is_file() and p.name != "SHA256SUMS.json"}
    (directory / "SHA256SUMS.json").write_text(json.dumps(hashes, indent=2))
    if args.archive:
        args.archive.mkdir(parents=True, exist_ok=False)
        for name in ("config.json", "summary.csv", "recovery.csv", "accuracy_crossings.csv", "convergence.csv", "report.html"):
            shutil.copyfile(directory / name, args.archive / name)
        groups = {"evidence": [], "meshes-thin_box": [], "meshes-torus": [], "meshes-double_sphere": []}
        for p in directory.rglob("*"):
            if p.is_file():
                if p.parent.name == "meshes":
                    group = next(k for k in groups if k.startswith("meshes-") and p.name.startswith(k[7:] + "_"))
                else:
                    group = "evidence"
                groups[group].append(p)
        for group, paths in groups.items():
            with zipfile.ZipFile(args.archive / f"{group}.zip", "w", zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
                for p in paths:
                    archive.write(p, p.relative_to(directory).as_posix())
        (args.archive / "SHA256SUMS.json").write_text(json.dumps({p.name: sha(p) for p in args.archive.iterdir() if p.is_file()}, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
