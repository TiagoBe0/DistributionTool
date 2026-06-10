#!/usr/bin/env python3
"""
Compare the DistributionTool (SOAP + Wigner-Seitz + hybrid) sweep against the
pre-computed Wigner-Seitz ground truth, across all PKA energies.

Inputs
  results/distool_sweep/{label}_hybrid.csv   (this tool: step,ws_vac,ws_int,
        hybrid_vac,hybrid_filtered,hybrid_only,grid_clusters,grid_vol_A3)
  results/{label}.csv                        (ground-truth WS time series:
        step,n_atoms_def,n_vacancies,n_interstitials,elapsed_s)

Outputs (results/distool_sweep/)
  damage_evolution_all.png   per-energy WS vs hybrid vs ground-truth WS
  survivors_vs_energy.png     plateau survivors vs energy (WS vs hybrid vs GT)
  ws_crosscheck.png           distool-WS vs ground-truth-WS on shared steps
  summary.csv                 plateau numbers per energy

Usage: python3 scripts/analyze_distool_vs_ws.py
"""
import os, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Works both on the cluster (/home/santi/pka_acero) and from the local sshfs mount.
_CANDS = ["/home/santi/pka_acero",
          "/run/user/1000/gvfs/sftp:host=192.168.220.25/home/santi/pka_acero"]
BASE = next((p for p in _CANDS if os.path.isdir(p)), _CANDS[0])
SWEEP = os.path.join(BASE, "results", "distool_sweep")
GT = os.path.join(BASE, "results")

# label -> (ground-truth csv, pretty energy keV)
ENERGIES = [
    ("1kev",      "1kev.csv",         1),
    ("2kev",      "2kev.csv",         2),
    ("3kev",      "3kev_results.csv", 3),
    ("4kev",      "4kev.csv",         4),
    ("5kev",      "5kev.csv",         5),
    ("6kev",      "6kev.csv",         6),
    ("7kev",      "7kev.csv",         7),
    ("8kev_real", "8kev_real.csv",    8),
]


def load_sweep(label):
    p = os.path.join(SWEEP, f"{label}_hybrid.csv")
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p)
    for c in ("ws_vac", "ws_int", "hybrid_vac", "hybrid_filtered",
              "hybrid_only", "grid_clusters", "grid_vol_A3"):
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.sort_values("step").reset_index(drop=True)


def load_gt(csv):
    p = os.path.join(GT, csv)
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p)
    return df.sort_values("step").reset_index(drop=True)


def plateau(df, col, ntail=3):
    """Median of the last `ntail` finite values (stable-defect estimate)."""
    if df is None or col not in df:
        return np.nan
    v = df[col].dropna()
    if v.empty:
        return np.nan
    return float(np.median(v.tail(ntail)))


def main():
    data = {}
    for label, gtcsv, keV in ENERGIES:
        data[label] = (keV, load_sweep(label), load_gt(gtcsv))

    # ---- Figure 1: per-energy damage evolution -----------------------------
    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    for ax, (label, gtcsv, keV) in zip(axes.flat, ENERGIES):
        _, sw, gt = data[label]
        if sw is not None:
            ax.plot(sw.step, sw.ws_vac, "-o", ms=3, color="tab:red",
                    label="WS (distool)")
            ax.plot(sw.step, sw.hybrid_vac, "-s", ms=3, color="tab:blue",
                    label="hybrid robust")
        if gt is not None:
            ax.plot(gt.step, gt.n_vacancies, "--", color="0.4", lw=1,
                    label="WS ground-truth")
        ax.set_title(f"{keV} keV")
        ax.set_xlabel("MD step")
        ax.set_ylabel("vacancies / Frenkel pairs")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("FeCrNi PKA cascades — DistributionTool WS vs hybrid vs "
                 "ground-truth Wigner-Seitz", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    f1 = os.path.join(SWEEP, "damage_evolution_all.png")
    fig.savefig(f1, dpi=130); fig.savefig(f1[:-4] + ".pdf")
    plt.close(fig)

    # ---- summary table + Figure 2: survivors vs energy ---------------------
    rows = []
    for label, gtcsv, keV in ENERGIES:
        _, sw, gt = data[label]
        rows.append(dict(
            label=label, keV=keV,
            ws_distool=plateau(sw, "ws_vac"),
            hybrid=plateau(sw, "hybrid_vac"),
            ws_groundtruth=plateau(gt, "n_vacancies"),
            ws_peak=(float(np.nanmax(sw.ws_vac)) if sw is not None else np.nan),
        ))
    summ = pd.DataFrame(rows).sort_values("keV")
    summ.to_csv(os.path.join(SWEEP, "summary.csv"), index=False)
    print(summ.to_string(index=False))

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(summ.keV, summ.ws_distool, "-o", color="tab:red",
            label="WS (distool, plateau)")
    ax.plot(summ.keV, summ.hybrid, "-s", color="tab:blue",
            label="hybrid robust (plateau)")
    ax.plot(summ.keV, summ.ws_groundtruth, "--^", color="0.4",
            label="WS ground-truth (plateau)")
    ax.set_xlabel("PKA energy (keV)")
    ax.set_ylabel("surviving Frenkel pairs")
    ax.set_title("Surviving defects vs PKA energy — FeCrNi")
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout()
    f2 = os.path.join(SWEEP, "survivors_vs_energy.png")
    fig.savefig(f2, dpi=130); fig.savefig(f2[:-4] + ".pdf")
    plt.close(fig)

    # ---- Figure 3: WS cross-check (distool vs ground truth on shared steps) -
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    allx, ally = [], []
    for label, gtcsv, keV in ENERGIES:
        _, sw, gt = data[label]
        if sw is None or gt is None:
            continue
        m = pd.merge(sw[["step", "ws_vac"]], gt[["step", "n_vacancies"]],
                     on="step", how="inner").dropna()
        if m.empty:
            continue
        ax.scatter(m.n_vacancies, m.ws_vac, s=18, alpha=0.6, label=f"{keV} keV")
        allx += list(m.n_vacancies); ally += list(m.ws_vac)
    if allx:
        hi = max(max(allx), max(ally)) * 1.05
        ax.plot([0, hi], [0, hi], "k--", lw=1, label="y = x")
        ax.set_xlim(0, hi); ax.set_ylim(0, hi)
    ax.set_xlabel("WS vacancies — ground truth (OVITO/LAMMPS)")
    ax.set_ylabel("WS vacancies — distool")
    ax.set_title("Wigner-Seitz cross-check (shared steps)")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout()
    f3 = os.path.join(SWEEP, "ws_crosscheck.png")
    fig.savefig(f3, dpi=130); fig.savefig(f3[:-4] + ".pdf")
    plt.close(fig)

    print(f"\nWrote:\n  {f1}\n  {f2}\n  {f3}\n  {os.path.join(SWEEP,'summary.csv')}")


if __name__ == "__main__":
    main()
