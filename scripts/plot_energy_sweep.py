#!/usr/bin/env python3
"""
Damage-evolution figures for the FeCrNi PKA cascades, read straight from the
per-energy reconciliation CSVs produced by run_energy_sweep_local.sh:

    results/<label>/<label>_hybrid.csv
    columns: step,ws_vac,ws_int,hybrid_vac,hybrid_filtered,hybrid_only,
             grid_clusters,grid_vol_A3

Produces, for each energy, a damage-evolution curve (WS vs hybrid vs grid), and
a cross-energy comparison panel. The last simulated frame is treated as the
"cold" state ONLY if the count has plateaued (WS == hybrid, i.e. no transit
atoms); otherwise it is annotated as not-yet-relaxed.

Usage:
    python3 scripts/plot_energy_sweep.py 9kev 10kev
Outputs:
    figures/damage_evolution_<label>.png
    figures/energy_comparison.png
"""
import os
import sys
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _figstyle import setup, save, C, ROOT, FIGDIR


def load(label):
    path = os.path.join(ROOT, "results", label, f"{label}_hybrid.csv")
    rows = {"step": [], "ws": [], "hyb": [], "filt": [], "grid": [], "vol": []}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            try:
                rows["step"].append(int(r["step"]))
                rows["ws"].append(int(r["ws_vac"]))
                rows["hyb"].append(int(r["hybrid_vac"]))
                rows["filt"].append(int(r["hybrid_filtered"]))
                rows["grid"].append(int(r["grid_clusters"]))
                rows["vol"].append(float(r["grid_vol_A3"]))
            except (ValueError, KeyError):
                continue
    return rows


def plot_one(label, d):
    ts = d["step"]
    peak_i = d["ws"].index(max(d["ws"]))
    cold = d["ws"][-1] == d["hyb"][-1]

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(ts, d["ws"], "o-", color=C["ws"], lw=2, ms=6,
            label="Wigner-Seitz (topological Frenkel pairs)")
    ax.plot(ts, d["hyb"], "s-", color=C["hybrid"], lw=2, ms=6,
            label="Hybrid `robust` (stable defects)")
    ax.plot(ts, d["grid"], "^--", color=C["grid"], lw=1.5, ms=5,
            label="Grid open-void clusters")

    ax.annotate(f"peak WS = {d['ws'][peak_i]}\n@ t={ts[peak_i]}",
                xy=(ts[peak_i], d["ws"][peak_i]),
                xytext=(ts[peak_i] + 0.12 * (ts[-1] - ts[0]), d["ws"][peak_i] * 0.85),
                arrowprops=dict(arrowstyle="->", color=C["muted"]), fontsize=9)
    tag = (f"cold: {d['ws'][-1]} survivors" if cold
           else f"NOT relaxed\nWS={d['ws'][-1]}, hyb={d['hyb'][-1]}")
    ax.annotate(tag, xy=(ts[-1], d["ws"][-1]),
                xytext=(ts[-1] - 0.18 * (ts[-1] - ts[0]), max(d["ws"]) * 0.25),
                arrowprops=dict(arrowstyle="->", color=C["muted"]), fontsize=9,
                color=C["grid"] if cold else C["ws"])

    ax.set_xlabel("MD step")
    ax.set_ylabel("defect count")
    ax.set_title(f"FeCrNi {label} PKA cascade — damage evolution")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out = save(fig, f"damage_evolution_{label}")
    plt.close(fig)
    return out, peak_i, cold


def plot_compare(data):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
    colors = {"5kev": "#e08a1e", "9kev": "#1f6fb2", "10kev": "#c0392b"}
    for label, d in data.items():
        c = colors.get(label, None)
        a1.plot(d["step"], d["ws"], "o-", ms=4, lw=1.8, color=c, label=f"{label} WS")
        a2.plot(d["step"], d["hyb"], "s-", ms=4, lw=1.8, color=c, label=f"{label} hybrid")
    a1.set_title("Wigner-Seitz Frenkel pairs")
    a2.set_title("Hybrid `robust` stable defects")
    for ax in (a1, a2):
        ax.set_xlabel("MD step")
        ax.set_ylabel("defect count")
        ax.legend()
        ax.grid(alpha=0.25)
    fig.suptitle("FeCrNi cascades — energy comparison")
    fig.tight_layout()
    out = save(fig, "energy_comparison")
    plt.close(fig)
    return out


def main():
    setup()
    os.makedirs(FIGDIR, exist_ok=True)
    labels = sys.argv[1:] or ["9kev", "10kev"]
    data = {}
    for label in labels:
        d = load(label)
        data[label] = d
        out, peak_i, cold = plot_one(label, d)
        print(f"{label}: peak WS={d['ws'][peak_i]} @t={d['step'][peak_i]}  "
              f"final(t={d['step'][-1]})=WS:{d['ws'][-1]}/hyb:{d['hyb'][-1]}  "
              f"{'COLD' if cold else 'NOT-RELAXED'} -> {out}")
    print("comparison ->", plot_compare(data))


if __name__ == "__main__":
    main()
