"""Aggregate seed-paired results and render scientific comparison figures."""
import argparse
from collections import defaultdict
import csv
import html
import hashlib
import json
import math
from pathlib import Path
import random
import statistics
import sys

# Optional task-local dependencies, never required by the C++ core.
LOCAL = Path(__file__).resolve().parents[1] / ".deps" / "plotting"
if LOCAL.is_dir():
    sys.path.insert(0, str(LOCAL))

METHODS = ("uniform", "adaptive_mse", "adaptive_risk")
COLORS = ("#4361ee", "#f48c06", "#007f5f")
METRICS = ("mae", "p99", "cvar99", "sign_error", "inside_recall", "outside_recall", "median_ns_per_query")


def interval(values):
    """Descriptive paired bootstrap of seed means, NOT independent query CIs."""
    if not values:
        return [None, None, None]
    mean = statistics.mean(values)
    if len(values) < 2:
        return [mean, None, None]
    rng = random.Random(101)
    boot = sorted(sum(rng.choices(values, k=len(values)))/len(values) for _ in range(5000))
    return [mean, boot[124], boot[4874]]


def analyse(directory):
    config = json.loads((directory / "config.json").read_text())
    with (directory / "metrics.csv").open(newline="") as file:
        rows = list(csv.DictReader(file))
    groups, matched = defaultdict(list), {}
    for row in rows:
        if row["method"] not in METHODS:
            continue
        key = (int(row["budget_bytes"]), row["dataset"], row["method"])
        groups[key].append(row)
        matched[(int(row["budget_bytes"]), row["dataset"], int(row["seed"]), row["method"])] = row
    result = []
    for (budget, dataset, method), group in sorted(groups.items()):
        summary = {"budget_bytes": budget, "dataset": dataset, "method": method,
                   "seeds": len(group), "representation_bytes_mean": statistics.mean(float(r["representation_bytes"]) for r in group),
                   "metrics": {}, "delta_vs_uniform": {}}
        for metric in METRICS:
            values = [float(r[metric]) for r in group if math.isfinite(float(r[metric]))]
            summary["metrics"][metric] = {"mean": statistics.mean(values) if values else None,
                                          "seed_sd": statistics.stdev(values) if len(values) > 1 else None}
            diffs = []
            for row in group:
                baseline = matched[(budget, dataset, int(row["seed"]), "uniform")]
                a, b = float(row[metric]), float(baseline[metric])
                if math.isfinite(a) and math.isfinite(b):
                    diffs.append(a-b)
            mean, low, high = interval(diffs)
            summary["delta_vs_uniform"][metric] = {"mean": mean, "bootstrap95_low": low, "bootstrap95_high": high}
        result.append(summary)
    payload = {"config": config, "summaries": result,
               "analysis_source_sha256": hashlib.sha256((Path(__file__).read_text(encoding="utf-8").rstrip()+"\n").encode("utf-8")).hexdigest(),
               "limits": "Development results. Bootstrap resamples paired seeds, not queries. Three seeds give a coarse variability estimate. Feature probes and 2D contours do not certify 3D topology."}
    (directory / "summary.json").write_text(json.dumps(payload, indent=2, allow_nan=False), encoding="utf-8")
    return payload


def plots(directory, payload):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    cfg = payload["config"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7), constrained_layout=True)
    for axis, dataset, metric, label in zip(axes, ("near", "near", "probes"),
                                           ("mae", "cvar99", "inside_recall"),
                                           ("Near-surface MAE", "Mean worst 1% error", "Interior feature retention")):
        for method, color in zip(METHODS, COLORS):
            selected = sorted((r for r in payload["summaries"] if r["method"] == method and r["dataset"] == dataset), key=lambda r:r["budget_bytes"])
            valid = [r for r in selected if r["metrics"][metric]["mean"] is not None]
            axis.plot([r["budget_bytes"]/1024 for r in valid], [r["metrics"][metric]["mean"] for r in valid],
                      "o-", color=color, label=method.replace("_", " "))
        axis.set(xlabel="Storage cap (KiB)", ylabel=label, xscale="log", title=label)
        axis.grid(alpha=0.2)
        if metric == "inside_recall":
            axis.set_ylim(-0.03, 1.03)
    axes[0].legend(fontsize=8)
    fig.suptitle(f"{cfg['shape']} · feature {cfg['feature']} · {len(cfg['seeds'])} development seeds")
    fig.savefig(directory / "frontier.png", dpi=180)
    plt.close(fig)

    # Middle budget gives a legible first visual; export all budgets for inspection.
    for budget in cfg["budgets"]:
        data = np.genfromtxt(directory / f"slice_{budget}.csv", delimiter=",", names=True)
        n = cfg["slice_size"]
        plane = cfg.get("slice_plane", "xz")
        second_axis = plane[1]
        x, z = data["x"].reshape(n, n), data[second_axis].reshape(n, n)
        truth = data["reference"].reshape(n, n)
        fig, axes = plt.subplots(1, 4, figsize=(13, 3.5), constrained_layout=True)
        for axis, name, title in zip(axes, ("reference", *METHODS), ("Reference", "Uniform grid", "Adaptive MSE", "Adaptive sign penalty")):
            values = data[name].reshape(n, n)
            axis.contourf(x, z, values < 0, levels=[-0.5, 0.5, 1.5], colors=["#faf9f6", "#94d2bd"])
            if truth.min() < 0 < truth.max():
                axis.contour(x, z, truth, levels=[0], colors="#222", linewidths=1, linestyles="dashed")
            if values.min() < 0 < values.max():
                axis.contour(x, z, values, levels=[0], colors="#006d77", linewidths=1.2)
            if name == "adaptive_risk":
                leaf_data = np.atleast_1d(np.genfromtxt(directory / f"leaves_{budget}.csv", delimiter=",", names=True))
                lines = []
                for b in leaf_data:
                    lines += [[(b["x0"],b["z0"]),(b["x1"],b["z0"])], [(b["x0"],b["z1"]),(b["x1"],b["z1"])],
                              [(b["x0"],b["z0"]),(b["x0"],b["z1"])], [(b["x1"],b["z0"]),(b["x1"],b["z1"])]]
                axis.add_collection(LineCollection(lines, colors="#555", linewidths=0.3, alpha=0.3))
            axis.set(title=title, xlabel="x", ylabel=second_axis, aspect="equal", xlim=(-0.9,0.9), ylim=(-0.9,0.9))
        fixed_axis = "z" if plane == "xy" else "y"
        fig.suptitle(f"{cfg['shape']} · {budget/1024:g} KiB cap · {fixed_axis}=0 · seed {cfg['seeds'][0]}\nDashed black: reference boundary; green: predicted interior")
        fig.savefig(directory / f"comparison_{budget}.png", dpi=180)
        plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path)
    args = p.parse_args()
    directories = [args.run] if (args.run / "config.json").is_file() else sorted(args.run.glob("*/"))
    cards = []
    for directory in directories:
        if not (directory / "complete.json").is_file():
            continue
        payload = analyse(directory)
        plots(directory, payload)
        cfg = payload["config"]
        budget = sorted(cfg["budgets"])[len(cfg["budgets"])//2]
        rel = directory.relative_to(args.run).as_posix()
        prefix = "" if rel == "." else rel + "/"
        cards.append(f'<section><h2>{html.escape(cfg["shape"])} · feature {cfg["feature"]}</h2>'
                     f'<img src="{prefix}comparison_{budget}.png" alt="Distance field cross section">'
                     f'<img src="{prefix}frontier.png" alt="Accuracy and feature retention versus storage">'
                     f'<p><a href="{prefix}summary.json">Seed summaries and paired intervals</a> · '
                     f'<a href="{prefix}metrics.csv">Raw measurements</a></p></section>')
    if not cards:
        p.error("No completed runs found")
    page = '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">'
    page += '<title>SDF Budget Lab · Development report</title><style>body{font:16px/1.6 system-ui;background:#faf9f6;color:#172a3a;max-width:1200px;margin:40px auto;padding:0 24px}h1{font-size:36px}section{border-top:1px solid #bbb;margin-top:36px;padding-top:12px}img{width:100%;height:auto}a{color:#006d77}</style>'
    page += '<h1>SDF Budget Lab</h1><p>Development evidence for adaptive distance-field compression under a shared byte cap.</p>'
    page += '<p>CVaR99 is the mean of the worst 1% absolute errors. Feature scores are known-point diagnostics. Timings are warm serial CPU batch averages; storage excludes build scratch space and allocator overhead. These plots do not certify topology.</p>'
    page += ''.join(cards) + '</html>'
    (args.run / "report.html").write_text(page, encoding="utf-8")
    print(f"Report: {(args.run / 'report.html').resolve()}")


if __name__ == "__main__":
    main()
