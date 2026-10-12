"""Independent NumPy decoder, full-sample verification, summaries and archives.

The allocation benchmark itself has no Python dependencies. This verifier uses
NumPy to recompute every saved fitting and evaluation query from raw q8 bytes.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics
import struct
import zipfile
import numpy as np


def rows(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, values):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(values[0]))
        w.writeheader()
        w.writerows(values)


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def hashes(directory):
    manifest = json.loads((directory / "SHA256SUMS.json").read_text())
    for name, digest in manifest.items():
        assert sha(directory / name.replace("\\", "/")) == digest, name


def decode(data, points):
    raw = np.frombuffer(data, dtype=np.uint8)
    words = np.frombuffer(data, dtype="<u4")
    floats = np.frombuffer(data, dtype="<f4")
    lo, hi = np.full_like(points, -1.), np.full_like(points, 1.)
    node = np.zeros(len(points), dtype=np.int64)
    active = np.ones(len(points), dtype=bool)
    for _ in range(7):
        descriptors = words[node // 4]
        branch = (descriptors & 3) == 0
        active &= branch
        if not active.any():
            break
        mid = (lo[active] + hi[active]) / 2
        side = points[active] >= mid
        octant = side @ np.array([1, 2, 4])
        node[active] = (descriptors[active] & ~np.uint32(3)).astype(np.int64) + 4 * octant
        lo[active] = np.where(side, mid, lo[active])
        hi[active] = np.where(side, hi[active], mid)
    assert not active.any(), "Unexpected depth/cycle"
    descriptors = words[node // 4]
    assert np.all((descriptors & 3) == 3), "Non-q8 leaf"
    offset = (descriptors & ~np.uint32(3)).astype(np.int64)
    base, scale = floats[offset // 4].astype(float), floats[offset // 4 + 1].astype(float)
    t = (points - lo) / (hi - lo)
    result = np.zeros(len(points))
    for c in range(8):
        value = (base + scale * raw[offset + 8 + c]).astype(np.float32).astype(float)
        weight = np.ones(len(points))
        for axis in range(3):
            weight *= t[:, axis] if c & (1 << axis) else 1 - t[:, axis]
        result += value * weight
    return result


def layout(data, row):
    used = set()
    branches = leaves = 0

    def claim(offset, size):
        assert offset % 4 == 0 and 0 <= offset <= len(data) - size
        region = set(range(offset, offset + size, 4))
        assert not used.intersection(region)
        used.update(region)

    def visit(node):
        nonlocal branches, leaves
        claim(node, 4)
        d, = struct.unpack_from("<I", data, node)
        e, offset = d & 3, d & ~3
        if e:
            assert e == 3
            claim(offset, 16)
            leaves += 1
        else:
            branches += 1
            for c in range(8):
                visit(offset + 4 * c)

    visit(0)
    assert len(used) * 4 == len(data)
    assert branches == int(row["branches"]) and leaves == int(row["q8_leaves"])


def close(actual, expected, name):
    assert math.isclose(actual, float(expected), rel_tol=2e-9, abs_tol=2e-10), (name, actual, expected)


def verify(directory):
    assert (directory / "COMPLETE").exists()
    assert (directory / "EXIT_CODE").read_text() == "0"
    hashes(directory)
    hashes(directory.with_name(directory.name + "-source"))
    config = json.loads((directory / "config.json").read_text())
    if config["stage"] == "development":
        data = rows(directory / "development.csv")
        assert len(data) == (18 if config["smoke"] else 54)
        reference = {(r["case"], r["budget"]): float(r["objective"]) for r in data if r["grid"] == "33"}
        selection = []
        for grid in (9, 17, 33):
            subset = [r for r in data if int(r["grid"]) == grid]
            excess = [float(r["objective"]) / reference[r["case"], r["budget"]] - 1 for r in subset]
            selection.append(dict(grid=grid, mean_excess=statistics.mean(excess), max_excess=max(excess), mean_seconds=statistics.mean(float(r["seconds"]) for r in subset)))
        write_csv(directory / "selection.csv", selection)
        selected = next(r["grid"] for r in selection if r["mean_excess"] <= .01 and r["max_excess"] <= .05)
        print(f"Development rule selects {selected} points", flush=True)
        return
    data = rows(directory / "metrics.csv")
    nfields = (2 if config["smoke"] else 24) if config["stage"] == "resolution" else (9 if config["smoke"] else 216)
    assert len(data) == nfields * 3
    assert len({(r["case"], r["depth"], r["budget"], r["method"], r["dataset"]) for r in data}) == len(data)
    samples = {}
    buffers = {}
    fit_errors = {}
    for r in data:
        case, dataset = r["case"], r["dataset"]
        if case not in samples:
            samples.clear()
            samples[case] = {s: np.loadtxt(directory / "samples" / f"{case}_{s}.csv", delimiter=",", skiprows=1, ndmin=2) for s in ("fit", "uniform", "near", "probes")}
        name = f"{case}_{r['depth']}_{r['budget']}_{r['method']}"
        if name not in buffers:
            buffers.clear()
            raw = (directory / "fields" / f"{name}.bin").read_bytes()
            buffers[name] = raw
            layout(raw, r)
            assert len(raw) + config["field_object_bytes"] == int(r["bytes"]) <= int(r["budget"])
            fit = samples[case]["fit"]
            pred = decode(raw, fit[:, :3])
            errors = np.abs(pred - fit[:, 3])
            if r["method"] == "cvar99":
                loss = np.maximum(errors - float(r["threshold"]), 0).sum()
                close(float(r["threshold"]) + loss / (.01 * len(fit)), r["fit_objective"], "variational objective")
            else:
                loss = np.sum((errors / 2) ** 2)
                if r["method"] == "sign":
                    valid = np.abs(fit[:, 3]) > 2e-9
                    wrong = np.where(fit[:, 3] < 0, pred >= 0, pred <= 0)
                    loss += .0025 * np.count_nonzero(valid & wrong)
            close(loss, r["fit_loss"], "packed fitting loss")
            if "fit_cvar99" in r:
                close(np.sort(errors)[-math.ceil(.01 * len(errors)):].mean(), r["fit_cvar99"], "fit CVaR")
                fit_errors[case, r["budget"], r["method"]] = (np.sort(errors), r)
            choices_path = directory / "choices" / f"{name}.csv"
            if choices_path.exists():
                choices = rows(choices_path)
                close(sum(float(c["loss"]) for c in choices if c["encoding"] == "3"), r["fit_loss"], "selected additive loss")
        sample = samples[case][dataset]
        assert len(sample) == int(r["count"])
        if dataset != "probes":
            assert len(sample) == config["test"] and (config["smoke"] or len(sample) >= 100000)
        prediction = decode(buffers[name], sample[:, :3])
        errors = np.sort(np.abs(prediction - sample[:, 3]))
        for key, value in dict(rmse=np.sqrt(np.mean(errors ** 2)), mae=errors.mean(), p99=errors[math.ceil(.99 * len(errors)) - 1], cvar99=errors[-math.ceil(.01 * len(errors)):].mean()).items():
            close(value, r[key], key)
        valid = np.abs(sample[:, 3]) > 2e-9
        inside, outside = valid & (sample[:, 3] < 0), valid & (sample[:, 3] > 0)
        for key, mask, correct in (("inside_recall", inside, prediction < 0), ("outside_recall", outside, prediction > 0)):
            if mask.any():
                close(np.mean(correct[mask]), r[key], key)
        if valid.any():
            close(np.mean(np.where(sample[valid, 3] < 0, prediction[valid] >= 0, prediction[valid] <= 0)), r["sign_error"], "sign_error")
    if config["stage"] == "evaluation":
        trace = rows(directory / "thresholds.csv")
        for (case, budget, method), (errors, row) in fit_errors.items():
            if method != "mse":
                continue
            upper = errors[-math.ceil(.01 * len(errors)):].mean()
            p99 = errors[math.ceil(.99 * len(errors)) - 1]
            expected = sorted(set([upper * i / (config["grid"] - 1) for i in range(config["grid"])] + [p99]))
            actual = [r for r in trace if r["case"] == case and r["budget"] == budget]
            assert len(actual) == len(expected)
            for a, e in zip(actual, expected):
                close(e, a["threshold"], "frozen grid")
                close(float(a["threshold"]) + float(a["excess"]) / (.01 * config["fit"]), a["objective"], "trace F")
            tail = fit_errors[case, budget, "cvar99"][1]
            close(min(float(a["objective"]) for a in actual), tail["fit_objective"], "finite search minimum")
            assert float(tail["fit_cvar99"]) <= float(tail["fit_objective"]) + 1e-10 <= upper + 2e-10
    surfaces = rows(directory / "surfaces.csv")
    cases = len({r["case"] for r in data})
    assert len(surfaces) == (nfields + cases) * (2 if config["smoke"] else 3)
    assert len(rows(directory / "construction.csv")) == nfields
    (directory / "VERIFIED.json").write_text(json.dumps(dict(fields=nfields, metric_rows=len(data), all_saved_points_recomputed=True, byte_coverage=True, fitting_losses=True, frozen_thresholds=config["stage"] == "evaluation"), indent=2))
    print(f"Verified {nfields} fields and every saved fitting/evaluation query", flush=True)
    summarize(directory, data, surfaces)


def summarize(directory, data, surfaces):
    construction = rows(directory / "construction.csv")
    summaries = []
    for depth, budget, method, dataset in sorted({(r["depth"], r["budget"], r["method"], r["dataset"]) for r in data}):
        selected = [r for r in data if (r["depth"], r["budget"], r["method"], r["dataset"]) == (depth, budget, method, dataset)]
        cost = [r for r in construction if (r["depth"], r["budget"], r["method"]) == (depth, budget, method)]
        summary = dict(depth=depth, budget=budget, method=method, dataset=dataset, cases=len(selected))
        summary.update({key: statistics.mean(float(r[key]) for r in selected) for key in ("bytes", "mae", "rmse", "cvar99", "inside_recall", "outside_recall")})
        for key in ("candidate_seconds", "allocation_seconds", "pack_seconds"):
            summary[key] = statistics.mean(float(r[key]) for r in cost)
        summary["construction_seconds"] = statistics.mean(sum(float(r.get(key, 0)) for key in ("label_seconds", "candidate_seconds", "loss_cache_seconds", "allocation_seconds", "pack_seconds")) for r in cost)
        summaries.append(summary)
    write_csv(directory / "summary.csv", summaries)
    paired = []
    for r in data:
        if r["method"] not in ("sign", "cvar99"):
            continue
        base = next(x for x in data if (x["case"], x["budget"], x["dataset"], x["method"]) == (r["case"], r["budget"], r["dataset"], "mse"))
        paired.append(dict(case=r["case"], shape=r["shape"], feature=r["feature"], budget=r["budget"], method=r["method"], dataset=r["dataset"], mae_change=float(r["mae"]) / float(base["mae"]) - 1, cvar99_change=float(r["cvar99"]) / float(base["cvar99"]) - 1))
    if paired:
        write_csv(directory / "paired.csv", paired)
    recovery = []
    for shape, depth, budget, method in sorted({(r["shape"], r["depth"], r["budget"], r["method"]) for r in data}):
        features = sorted({float(r["feature"]) for r in data if r["shape"] == shape})
        passes = []
        for feature in features:
            field = [r for r in surfaces if (r["shape"], r["depth"], r["budget"], r["method"]) == (shape, depth, budget, method) and float(r["feature"]) == feature and int(r["cells"]) in (64, 96)]
            oracle = [r for r in surfaces if r["shape"] == shape and float(r["feature"]) == feature and r["method"] == "oracle" and int(r["cells"]) in (64, 96)]
            passes.append(bool(field and oracle) and all(r["success"] == "1" for r in field + oracle))
        minimum = next((f for i, f in enumerate(features) if all(passes[i:])), "not reached")
        recovery.append(dict(shape=shape, depth=depth, budget=budget, method=method, minimum_feature=minimum, full_size=minimum if minimum == "not reached" or shape == "double_sphere" else 2 * minimum))
    write_csv(directory / "recovery.csv", recovery)


def publish(directories, destination):
    if destination.exists():
        raise ValueError("Publication must be new")
    destination.mkdir(parents=True)
    archives = []
    for directory in directories:
        for name in ("summary.csv", "recovery.csv", "paired.csv", "selection.csv", "VERIFIED.json", "config.json"):
            if (directory / name).exists():
                shutil.copy2(directory / name, destination / f"{directory.name}-{name}")
        files = sorted(f for folder in (directory, directory.with_name(directory.name + "-source")) for f in folder.rglob("*") if f.is_file())
        # Split by uncompressed size, keeping each ZIP comfortably below 100 MB.
        groups, current, size = [], [], 0
        for f in files:
            if current and size + f.stat().st_size > 80_000_000:
                groups.append(current)
                current, size = [], 0
            current.append(f)
            size += f.stat().st_size
        if current:
            groups.append(current)
        for index, group in enumerate(groups, 1):
            path = destination / f"{directory.name}-{index:02}.zip"
            with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
                for f in group:
                    z.write(f, f.relative_to(directory.parent))
            archives.append(dict(file=path.name, sha256=sha(path), files={f.relative_to(directory.parent).as_posix(): sha(f) for f in group}))
    for name in ("tail_evidence.py", "tail_run.py"):
        shutil.copy2(Path(__file__).with_name(name), destination / name)
    (destination / "ARCHIVES.json").write_text(json.dumps(archives, indent=2))
    verify_archives(destination)


def verify_archives(directory):
    if (directory / "SHA256SUMS.json").exists():
        hashes(directory)
    for archive in json.loads((directory / "ARCHIVES.json").read_text()):
        path = directory / archive["file"]
        assert sha(path) == archive["sha256"]
        with zipfile.ZipFile(path) as z:
            assert set(z.namelist()) == set(archive["files"])
            for name, digest in archive["files"].items():
                with z.open(name) as f:
                    assert hashlib.file_digest(f, "sha256").hexdigest() == digest, name
    print("All archive/member hashes verified", flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("directories", nargs="+", type=Path)
    p.add_argument("--publish", type=Path)
    p.add_argument("--verify-archives", action="store_true")
    a = p.parse_args()
    if a.verify_archives:
        for directory in a.directories:
            verify_archives(directory)
    else:
        for directory in a.directories:
            verify(directory)
        if a.publish:
            publish(a.directories, a.publish)
