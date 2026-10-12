"""Add readable figures/index to an already verified tail evidence package."""
import argparse
import hashlib
import html
import json
from pathlib import Path
import shutil


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("publication", type=Path)
    p.add_argument("figures", type=Path)
    a = p.parse_args()
    for name in ("tradeoffs.png", "tradeoffs.pdf", "surfaces.png"):
        shutil.copy2(a.figures / name, a.publication / name)
    for name in ("tail_run.py", "tail_evidence.py", "tail_plots.py", "tail_publication.py"):
        shutil.copy2(Path(__file__).with_name(name), a.publication / name)
    archives = json.loads((a.publication / "ARCHIVES.json").read_text())
    links = []
    for f in sorted(a.publication.glob("*.csv")):
        links.append(f'<li><a href="{f.name}">{html.escape(f.name)}</a></li>')
    archive_links = [f'<li><a href="{r["file"]}">{r["file"]}</a> ({(a.publication / r["file"]).stat().st_size / 1e6:.1f} MB)</li>' for r in archives]
    content = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Fixed-q8 CVaR99 allocation evidence</title>
<style>body{max-width:1100px;margin:40px auto;padding:0 24px;font:17px/1.6 system-ui;color:#20333c;background:#fcfdfc}h1,h2{line-height:1.2}img{max-width:100%;height:auto}a{color:#096b80}figure{margin:30px 0}.lead{font-size:21px}li{margin:5px 0}code{font-size:14px}table{border-collapse:collapse}td,th{padding:8px 20px;border-bottom:1px solid #ccd8dc;text-align:left}</style>
<h1>Direct CVaR99 allocation with fixed q8</h1>
<p>8 October 2026 · 216 equal-budget fields · depth five · 100,000 independent queries per distribution</p>
<p class="lead">At 32 KiB, uniform tail error falls 10.6%, but near-surface tail error falls only 3.0%, with a 20.5% increase in near-surface mean error. No allocation objective lowers the minimum recovered feature sizes.</p>
<p><a href="../tail-study.md">Full study and limitations</a> · <a href="../tail-protocol.md">Staged protocol</a></p>
<figure><img src="tradeoffs.png" alt="Mean error, tail error, construction time and geometric thresholds across equal budget caps"><figcaption>Equal-case means across 24 shape/feature/seed cases. CVaR construction includes the baseline solve and every threshold solve.</figcaption></figure>
<h2>Resolution comes first</h2><p>Full depth five recovers a 0.02 wall half-thickness and 0.06 torus radius where full depth four fails. Full-tree storage rises from 84,324 to 674,148 bytes. Smaller sphere gaps remain unresolved.</p>
<figure><img src="surfaces.png" alt="Actual torus meshes: fragmented at full depth four, connected at depth five and under all three 16 KiB objectives"><figcaption>Actual torus meshes at radius 0.06, extracted at 96 cells per axis. Geometry success also requires the 64-cell extraction and independent sign probes.</figcaption></figure>
<h2>What was verified</h2><p>All four C++ test suites pass, including 4,747 tail checks. An independent NumPy decoder checks every saved fitting/evaluation point, q8 byte coverage, storage caps, fitting losses and finite-search selection. Raw failures and finest meshes are retained. Remote CI has not been observed.</p>
<p>The DP is exact for each threshold and the search selects the minimum over its frozen finite set. This is not continuous-threshold optimality or a feature-preservation guarantee.</p>
<h2>Readable measurements</h2><ul>''' + "".join(links) + '''</ul><h2>Frozen evidence archives</h2><p>Each archive is below 100 MB. Together they preserve samples, buffers, choices, threshold traces, meshes, source/executable snapshots and hashes. Source snapshots are the versions taken before each stage ran.</p><ul>''' + "".join(archive_links) + '''</ul><p><a href="ARCHIVES.json">Archive and member SHA-256 manifest</a> · <a href="SHA256SUMS.json">Publication SHA-256 manifest</a></p><p>Verify with <code>python scripts/tail_evidence.py docs/tail-v1 --verify-archives</code> in a NumPy environment.</p></html>'''
    (a.publication / "index.html").write_text(content, encoding="utf-8")
    (a.publication / "SHA256SUMS.json").write_text(json.dumps({f.relative_to(a.publication).as_posix(): sha(f) for f in a.publication.rglob("*") if f.is_file() and f.name != "SHA256SUMS.json"}, indent=2))


if __name__ == "__main__":
    main()
