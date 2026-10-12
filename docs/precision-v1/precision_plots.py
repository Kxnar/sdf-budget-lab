"""Static scientific figures from the frozen precision run (requires matplotlib)."""
import argparse
import csv
from pathlib import Path
import sys

LOCAL = Path(__file__).resolve().parents[1] / ".deps" / "plotting"
if LOCAL.is_dir():
    sys.path.insert(0, str(LOCAL))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

METHODS = ["f32", "q16", "q8", "joint"]
COLORS = ["#4267ad", "#bc772b", "#588754", "#9b3d80"]


def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    summary = read(args.run / "summary.csv")
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), constrained_layout=True)
    for method, color in zip(METHODS, COLORS):
        data = [r for r in summary if r["method"] == method]
        x = [int(r["budget"]) / 1024 for r in data]
        for ax, metric, label in zip(axes, ["near_rmse", "inside_recall", "query_ns"], ["Near-surface RMSE", "Feature inside recall", "Query latency (ns)"]):
            ax.plot(x, [float(r[metric]) for r in data], "o-", color=color, label=method)
            ax.set(xlabel="Storage cap (KiB)", ylabel=label, xticks=[8, 16, 32])
            ax.grid(alpha=.2)
    axes[0].legend()
    axes[1].set_ylim(0, 1.04)
    fig.suptitle("Equal storage caps · means over 24 analytic cases · actual packed queries")
    fig.savefig(args.out / "tradeoffs.png", dpi=180)
    plt.close(fig)

    metrics = read(args.run / "metrics.csv")
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), constrained_layout=True)
    for ax, shape in zip(axes, ["thin_box", "torus", "double_sphere"]):
        for method, color in zip(METHODS, COLORS):
            data = [r for r in metrics if r["dataset"] == "probes" and r["method"] == method and r["shape"] == shape and r["budget"] == "32768"]
            features = sorted({float(r["feature"]) for r in data})
            values = [np.mean([float(r["inside_recall"]) for r in data if float(r["feature"]) == f]) for f in features]
            ax.plot(features, values, "o-", color=color, label=method)
        ax.set(title=shape, xlabel="Half-thickness / tube radius / gap", ylabel="Inside recall", ylim=(0, 1.04))
        ax.grid(alpha=.2)
    axes[0].legend()
    fig.suptitle("32 KiB feature probes · two-seed means · not a topology measurement")
    fig.savefig(args.out / "feature_sweep.png", dpi=180)
    plt.close(fig)

    surfaces = read(args.run / "surfaces.csv")
    fig = plt.figure(figsize=(15, 8), layout="constrained")
    for row, feature in enumerate([.035, .06, .10]):
        case = f"torus_{feature:.6f}_0"
        for col, method in enumerate(["oracle"] + METHODS):
            budget = 0 if method == "oracle" else 16384
            path = args.run / "meshes" / f"{case}_{budget}_{method}.obj"
            vertices, faces = [], []
            with path.open() as f:
                for line in f:
                    values = line.split()
                    if values[0] == "v":
                        vertices.append([float(v) for v in values[1:]])
                    elif values[0] == "f":
                        faces.append([int(v) - 1 for v in values[1:]])
            ax = fig.add_subplot(3, 5, 5 * row + col + 1, projection="3d")
            if faces:
                triangles = np.asarray(vertices)[np.asarray(faces)]
                normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
                light = np.array([.25, -.4, 1.])
                light /= np.linalg.norm(light)
                intensity = .3 + .7 * np.maximum(0, normals @ light / np.maximum(np.linalg.norm(normals, axis=1), 1e-30))
                colors = np.array(matplotlib.colors.to_rgb("#779fbe" if method == "oracle" else COLORS[col-1]))[None, :] * intensity[:, None]
                ax.add_collection3d(Poly3DCollection(triangles, facecolors=colors, linewidths=0, rasterized=True))
            data = next(r for r in surfaces if r["case"] == case and int(r["budget"]) == budget and r["method"] == method and r["cells"] == "96")
            ax.set_title(f'{method} · {data["components"]} components\nχ={data["euler"]} · radius {feature:g}', fontsize=10)
            ax.set(xlim=(-.75, .75), ylim=(-.75, .75), zlim=(-.22, .22))
            ax.set_box_aspect((1.5, 1.5, .44))
            ax.view_init(elev=55, azim=-65)
            ax.set_axis_off()
    fig.suptitle("Extracted 3D torus surfaces · 16 KiB · seed 0 · 96 cells per axis\nClosed ring requires one component with χ=0; fragmentation remains visible", fontsize=14)
    fig.savefig(args.out / "surfaces.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
