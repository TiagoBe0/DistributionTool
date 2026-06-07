#!/usr/bin/env python3
"""
PKA energy sweep figures for the Clementina FeCrNi cascades (1-8 keV).

Reads the per-energy reconciliation CSVs produced on the cluster by
distool_runs/run_energy.sh (columns:
  step, ws_vac, ws_int, hybrid_vac, hybrid_filtered, hybrid_only,
  grid_clusters, grid_vol_A3)
and makes two figures:

  figures/clementina_evolution.png  — WS vs hybrid-robust damage evolution,
                                       one panel per PKA energy.
  figures/clementina_survivors.png  — WS peak vs hybrid-robust-at-peak vs final
                                       survivors, as a function of PKA energy.

Usage:
    python3 scripts/plot_clementina_sweep.py [csv_dir]
Default csv_dir: clementina/sweep_out
"""
import os
import sys
import glob
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))


def energy_kev(label):
    # "5kev", "8kev_real" -> 5, 8
    import re
    m = re.match(r"(\d+)", label)
    return int(m.group(1)) if m else 0


def load(csv_dir):
    data = {}
    for path in sorted(glob.glob(os.path.join(csv_dir, "*_hybrid.csv"))):
        label = os.path.basename(path).replace("_hybrid.csv", "")
        steps, ws, hy, gc = [], [], [], []
        with open(path) as f:
            for row in csv.DictReader(f):
                def num(k):
                    v = row.get(k, "NA")
                    try: return float(v)
                    except (ValueError, TypeError): return np.nan
                steps.append(num("step")); ws.append(num("ws_vac"))
                hy.append(num("hybrid_vac")); gc.append(num("grid_clusters"))
        if steps:
            order = np.argsort(steps)
            data[label] = {
                "step": np.array(steps)[order], "ws": np.array(ws)[order],
                "hy": np.array(hy)[order], "gc": np.array(gc)[order],
                "E": energy_kev(label),
            }
    return dict(sorted(data.items(), key=lambda kv: kv[1]["E"]))


def fig_evolution(data, figdir):
    n = len(data)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.2 * nrow),
                             squeeze=False)
    for ax in axes.flat:
        ax.set_visible(False)
    for i, (label, d) in enumerate(data.items()):
        ax = axes.flat[i]; ax.set_visible(True)
        ax.plot(d["step"], d["ws"], "o-", color="#c0392b", ms=4, lw=1.6,
                label="Wigner-Seitz")
        ax.plot(d["step"], d["hy"], "s-", color="#1f6fb2", ms=4, lw=1.6,
                label="Hybrid `robust`")
        ax.set_title(f"{d['E']} keV", fontsize=10)
        ax.grid(alpha=0.25)
        ax.set_xlabel("MD step"); ax.set_ylabel("Frenkel pairs")
        if i == 0:
            ax.legend(frameon=False, fontsize=8)
    fig.suptitle("FeCrNi PKA cascades — WS over-counts the ballistic peak; "
                 "hybrid `robust` estimates surviving defects", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(figdir, "clementina_evolution.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
    print("wrote", os.path.relpath(out, ROOT))


def fig_survivors(data, figdir):
    E, ws_peak, hy_peak, survivors = [], [], [], []
    for label, d in data.items():
        ws = d["ws"]; hy = d["hy"]
        if np.all(np.isnan(ws)):
            continue
        ipk = int(np.nanargmax(ws))               # ballistic peak (max WS)
        E.append(d["E"])
        ws_peak.append(ws[ipk])
        hy_peak.append(hy[ipk])
        survivors.append(ws[np.where(~np.isnan(ws))[0][-1]])  # last valid = relaxed
    order = np.argsort(E)
    E = np.array(E)[order]; ws_peak = np.array(ws_peak)[order]
    hy_peak = np.array(hy_peak)[order]; survivors = np.array(survivors)[order]

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(E, ws_peak, "o-", color="#c0392b", lw=2, ms=8, label="WS peak (ballistic)")
    ax.plot(E, hy_peak, "s-", color="#1f6fb2", lw=2, ms=8, label="Hybrid `robust` at peak")
    ax.plot(E, survivors, "^-", color="#27ae60", lw=2, ms=8, label="Final survivors (relaxed)")
    ax.set_xlabel("PKA energy [keV]")
    ax.set_ylabel("Frenkel-pair count")
    ax.set_title("WS peak over-count vs hybrid estimate vs survivors")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    out = os.path.join(figdir, "clementina_survivors.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
    print("wrote", os.path.relpath(out, ROOT))


def main():
    csv_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "clementina", "sweep_out")
    data = load(csv_dir)
    if not data:
        sys.exit(f"no *_hybrid.csv found in {csv_dir}")
    print("energies:", ", ".join(f"{d['E']}keV(n={len(d['step'])})" for d in data.values()))
    figdir = os.path.join(ROOT, "figures")
    os.makedirs(figdir, exist_ok=True)
    fig_evolution(data, figdir)
    fig_survivors(data, figdir)


if __name__ == "__main__":
    main()
