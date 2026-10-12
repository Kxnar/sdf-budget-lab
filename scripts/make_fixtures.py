"""Generate small closed, triangulated fixtures without external dependencies."""
from pathlib import Path
import math

ROOT = Path(__file__).resolve().parents[1] / "fixtures"
ROOT.mkdir(exist_ok=True)


def off(name, vertices, faces):
    path = ROOT / name
    path.write_text("OFF\n" + f"{len(vertices)} {len(faces)} 0\n" +
                    "".join(" ".join(f"{x:.17g}" for x in v) + "\n" for v in vertices) +
                    "".join("3 " + " ".join(map(str, f)) + "\n" for f in faces), encoding="utf-8")


def box(h):
    vertices = [(h[0]*(1 if i & 1 else -1), h[1]*(1 if i & 2 else -1),
                 h[2]*(1 if i & 4 else -1)) for i in range(8)]
    faces = [(0, 2, 3), (0, 3, 1), (4, 5, 7), (4, 7, 6), (0, 1, 5), (0, 5, 4),
             (2, 6, 7), (2, 7, 3), (0, 4, 6), (0, 6, 2), (1, 3, 7), (1, 7, 5)]
    return vertices, faces


vertices, faces = box((0.5, 0.5, 0.5))
off("cube.off", vertices, faces)
off("open.off", vertices, faces[:-1])
off("outside.off", [(x*3, y*3, z*3) for x, y, z in vertices], faces)
off("degenerate.off", [(0, 0, 0), (0.5, 0, 0), (0.8, 0, 0), (0, 0.5, 0.5)],
    [(0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)])
off("intersecting.off", vertices + [(x+0.2, y+0.2, z+0.2) for x, y, z in vertices],
    faces + [tuple(i+8 for i in f) for f in faces])
off("thin_box.off", *box((0.6, 0.6, 0.04)))
points, triangles = [], []
major, minor = 64, 16
for i in range(major):
    a = 2*math.pi*i/major
    for j in range(minor):
        b = 2*math.pi*j/minor
        r = 0.55+0.04*math.cos(b)
        points.append((r*math.cos(a), r*math.sin(a), 0.04*math.sin(b)))
for i in range(major):
    for j in range(minor):
        a, b = i*minor+j, ((i+1) % major)*minor+j
        c, d = ((i+1) % major)*minor+(j+1) % minor, i*minor+(j+1) % minor
        triangles += [(a, b, c), (a, c, d)]
off("torus.off", points, triangles)
(ROOT / "thin_box_probes.csv").write_text(
    "x,y,z\n" + "".join(f"{x},{y},{z}\n" for x in (-0.4, 0, 0.4)
                        for y in (-0.4, 0, 0.4) for z in (0, -0.06, 0.06)), encoding="utf-8")
(ROOT / "torus_probes.csv").write_text(
    "x,y,z\n" + "".join(f"{r*math.cos(2*math.pi*i/256):.17g},{r*math.sin(2*math.pi*i/256):.17g},0\n"
                        for i in range(256) for r in (0.55, 0.61)), encoding="utf-8")
print(f"Fixtures written to {ROOT}")
