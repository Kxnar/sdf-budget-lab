"""Static scientific figures from verified tail-study results."""
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def rows(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("resolution", type=Path)
    p.add_argument("evaluation", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    data = rows(a.evaluation / "summary.csv")
    methods = [("mse", "MSE", "#236dba"), ("sign", "MSE + sign", "#cd8638"), ("cvar99", "CVaR99", "#168779")]
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), layout="constrained")
    for row, dataset in enumerate(("uniform", "near")):
        for col, metric in enumerate(("mae", "cvar99")):
            ax = axes[row, col]
            for method, label, color in methods:
                selected = sorted((r for r in data if r["method"] == method and r["dataset"] == dataset), key=lambda r: int(r["budget"]))
                ax.plot([int(r["budget"]) / 1024 for r in selected], [float(r[metric]) for r in selected], "o-", color=color, label=label)
            ax.set(title=f"{'Near surface' if row else 'Uniform'} · {'mean absolute error' if col == 0 else 'CVaR99'}", xlabel="Budget cap (KiB)", ylabel="Distance error", xticks=[8, 16, 32])
            ax.grid(alpha=.2)
    ax = axes[0, 2]
    for method, label, color in methods:
        selected = sorted((r for r in data if r["method"] == method and r["dataset"] == "near"), key=lambda r: int(r["budget"]))
        ax.plot([int(r["budget"]) / 1024 for r in selected], [1000 * float(r["construction_seconds"]) for r in selected], "o-", color=color, label=label)
    ax.set(title="Total construction", xlabel="Budget cap (KiB)", ylabel="Milliseconds", xticks=[8, 16, 32])
    ax.grid(alpha=.2)
    ax = axes[1, 2]
    ax.axis("off")
    ax.text(0, .98, "Geometry thresholds are identical", weight="bold", fontsize=12, va="top")
    ax.text(0, .83, "Minimum tested feature passing all checks:\n\nWall half-thickness     0.035 at all caps\nTorus tube radius       0.10 at 8 KiB\n                                   0.06 at 16/32 KiB\nSphere surface gap     0.10 at all caps\n\nRequires both seeds, 64/96 extraction grids,\noracle success and ≥95% sign-probe recall.\n\nLower tail error did not lower these thresholds.", va="top", linespacing=1.5)
    axes[0, 0].legend(frameon=False)
    fig.suptitle("Fixed q8, depth five · mean of 24 cases · 100,000 queries per distribution", fontsize=14)
    fig.savefig(a.out / "tradeoffs.png", dpi=170)
    fig.savefig(a.out / "tradeoffs.pdf")
    plt.close(fig)
    # Resolve the full-tree byte counts from actual evidence rather than ABI assumptions.
    views = [(next((a.resolution / "meshes").glob(f"torus_0.060000_20_{depth}_*_full_q8.obj")), f"Full depth {depth}") for depth in (4, 5)]
    views += [(a.evaluation / f"meshes/torus_0.060000_30_5_16384_{method}.obj", f"16 KiB · {label}") for method, label, _ in methods]
    fig = plt.figure(figsize=(15, 3.3), layout="constrained")
    for index, (path, title) in enumerate(views, 1):
        vertices, triangles = [], []
        for line in path.read_text().splitlines():
            if line.startswith("v "):
                vertices.append([float(x) for x in line.split()[1:]])
            elif line.startswith("f "):
                triangles.append([int(x) - 1 for x in line.split()[1:]])
        vertices, triangles = np.asarray(vertices), np.asarray(triangles)
        ax = fig.add_subplot(1, 5, index, projection="3d")
        ax.add_collection3d(Poly3DCollection(vertices[triangles], facecolor="#368b93", edgecolor="none", alpha=1))
        ax.set(xlim=(-.8, .8), ylim=(-.8, .8), zlim=(-.25, .25), title=title)
        ax.set_box_aspect((1, 1, .35))
        ax.view_init(58, -65)
        ax.set_axis_off()
    fig.suptitle("Torus radius 0.06 · actual meshes extracted at 96 cells per axis", fontsize=13)
    fig.savefig(a.out / "surfaces.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
