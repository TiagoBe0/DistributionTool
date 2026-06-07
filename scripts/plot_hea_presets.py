#!/usr/bin/env python3
"""
HEA preset-comparison figure for the FeCrNi `dump-finalCool.160000` frame
(README §3.9 — the complement of the 5 keV cascade story).

This frame is a COOLED / relaxed high-entropy alloy: the MD has already resolved
recombination, so every surviving Frenkel pair is a stable real defect. Here the
`robust` preset is the WRONG choice — its ballistic transit/Frenkel penalties fire
on stable defects and make it UNDER-count, whereas WS and the `relaxed` preset agree
exactly.

Panel A: the four estimators on this single relaxed frame.
Panel B: why `robust` rejects — mean per-signal value, accepted vs rejected
         candidates. The positive base signals are high for BOTH groups (all are
         genuine WS-empty sites); only the transit and Frenkel penalties separate
         them, proving the filtering is penalty-driven.

This script is fully self-reproducible: it parses the numbers directly from the run
logs and the per-candidate breakdown CSV — nothing is hard-coded.

Required inputs (produced by distool):
    results/hea_finalCool/robust.log
    results/hea_finalCool/relaxed.log
    results/hea_finalCool/hybrid_vacancies_robust.csv
Regenerate those with:
    ./build/distool --hybrid --hybrid-preset robust  --output results/hea_finalCool/robust.csv  dump-finalCool.0 dump-finalCool.160000
    ./build/distool --hybrid --hybrid-preset relaxed --output results/hea_finalCool/relaxed.csv dump-finalCool.0 dump-finalCool.160000

Usage:
    python3 scripts/plot_hea_presets.py
Outputs:
    figures/hea_presets.png  (+ .pdf)
"""
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE   = os.path.dirname(os.path.abspath(__file__))
ROOT   = os.path.normpath(os.path.join(HERE, ".."))
RESDIR = os.path.join(ROOT, "results", "hea_finalCool")

# Per-signal columns to plot in panel B and their CSV header names.
#   header: # x y z consensus_score ws soft_ws vor dens soap topo
#           transit_pen frenkel_pen transit_filtered in_recomb ref_site_idx accepted
PANEL_B = [("ws", "ws"), ("soft_ws", "soft_ws"), ("dens", "dens"),
           ("soap", "soap"), ("topo", "topo"),
           ("transit\npen.", "transit_pen"), ("frenkel\npen.", "frenkel_pen")]
PENALTY_COLS = {"transit_pen", "frenkel_pen"}


def parse_reconciliation(log_path):
    """Return dict with ws, hybrid, grid_clusters from a distool log."""
    out = {}
    with open(log_path) as f:
        txt = f.read()
    m = re.search(r"topological \(WS\)\s*:\s*(\d+)\s+vacancies", txt)
    if m: out["ws"] = int(m.group(1))
    m = re.search(r"topological \(hybrid\)\s*:\s*(\d+)", txt)
    if m: out["hybrid"] = int(m.group(1))
    m = re.search(r"open void \(grid\)\s*:\s*(\d+)\s+clusters", txt)
    if m: out["grid"] = int(m.group(1))
    return out


def parse_breakdown(csv_path):
    """Return (header_names, mean_accepted, mean_rejected, n_acc, n_rej)."""
    with open(csv_path) as f:
        header = f.readline().lstrip("#").split()
        idx = {name: i for i, name in enumerate(header)}
        acc_idx = idx["accepted"]
        rows_acc, rows_rej = [], []
        for line in f:
            parts = line.split()
            if len(parts) < len(header):
                continue
            vals = [float(p) for p in parts]
            (rows_acc if vals[acc_idx] >= 0.5 else rows_rej).append(vals)
    acc = np.array(rows_acc)
    rej = np.array(rows_rej)
    means_acc = {name: acc[:, idx[name]].mean() for _, name in PANEL_B}
    means_rej = {name: rej[:, idx[name]].mean() for _, name in PANEL_B}
    return means_acc, means_rej, len(rows_acc), len(rows_rej)


def main():
    rob_log = os.path.join(RESDIR, "robust.log")
    rel_log = os.path.join(RESDIR, "relaxed.log")
    csv     = os.path.join(RESDIR, "hybrid_vacancies_robust.csv")
    for p in (rob_log, rel_log, csv):
        if not os.path.exists(p):
            sys.exit(f"missing input: {p}\nRun the two distool commands in the "
                     f"module docstring first.")

    rob = parse_reconciliation(rob_log)
    rel = parse_reconciliation(rel_log)
    means_acc, means_rej, n_acc, n_rej = parse_breakdown(csv)

    ws_truth   = rob["ws"]
    ESTIMATORS = ["Wigner-Seitz", "Hybrid\n`relaxed`", "Hybrid\n`robust`", "Grid\n(open void)"]
    COUNTS     = [rob["ws"], rel["hybrid"], rob["hybrid"], rob["grid"]]
    COLORS     = ["#c0392b", "#1f6fb2", "#7fb0d6", "#27ae60"]
    n_filtered = rob["ws"] - rob["hybrid"]

    figdir = os.path.join(ROOT, "figures")
    os.makedirs(figdir, exist_ok=True)

    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(11.5, 4.7), gridspec_kw={"width_ratios": [1.0, 1.35]})

    # ===== Panel A =====
    x = np.arange(len(ESTIMATORS))
    axA.bar(x, COUNTS, color=COLORS, width=0.62, edgecolor="0.3")
    ytop = max(COUNTS) * 1.15
    for xi, c in zip(x, COUNTS):
        axA.text(xi, c + ytop * 0.015, str(c), ha="center", va="bottom",
                 fontsize=10, fontweight="bold")
    axA.axhline(ws_truth, color="0.4", ls=":", lw=1.3)
    axA.set_xticks(x)
    axA.set_xticklabels(ESTIMATORS, fontsize=9)
    axA.set_ylabel("Vacancy count")
    axA.set_ylim(0, ytop)
    axA.set_title("Relaxed FeCrNi HEA (finalCool, t=160000)")
    axA.annotate("WS = relaxed\n(stable defects)", xy=(0.5, ws_truth),
                 xytext=(0.5, ytop * 0.95), ha="center", fontsize=8.5, color="#1a5276")
    axA.annotate("wrong preset:\npenalties fire on\nstable defects →\nunder-counts",
                 xy=(2, rob["hybrid"]), xytext=(2.55, ytop * 0.56), fontsize=8.2,
                 color="#922b21", ha="center",
                 arrowprops=dict(arrowstyle="->", color="#922b21"))
    axA.text(3, ytop * 0.20, "(different\nquantity)", ha="center",
             fontsize=7.5, color="#1e7d4f")

    # ===== Panel B =====
    labels = [lbl for lbl, _ in PANEL_B]
    rej_vals = [means_rej[name] for _, name in PANEL_B]
    acc_vals = [means_acc[name] for _, name in PANEL_B]
    xb = np.arange(len(PANEL_B))
    w = 0.4
    axB.bar(xb - w/2, rej_vals, w, label=f"rejected (N={n_rej})",
            color="#e08283", edgecolor="0.3")
    axB.bar(xb + w/2, acc_vals, w, label=f"accepted (N={n_acc})",
            color="#5499c7", edgecolor="0.3")
    # shade the penalty region (first penalty column onward)
    first_pen = next(i for i, (_, name) in enumerate(PANEL_B) if name in PENALTY_COLS)
    axB.axvspan(first_pen - 0.5, len(PANEL_B) - 0.5, color="#f9e79f", alpha=0.35)
    axB.text((first_pen - 0.5 + len(PANEL_B) - 0.5) / 2, 1.06,
             "penalties\n(the discriminator)", ha="center", va="bottom",
             fontsize=8.5, color="#7d6608")
    axB.text((first_pen - 0.5) / 2, 1.06, "positive base signals\n(high for BOTH groups)",
             ha="center", va="bottom", fontsize=8.5, color="0.35")
    axB.set_xticks(xb)
    axB.set_xticklabels(labels, fontsize=8.5)
    axB.set_ylabel("Mean signal value")
    axB.set_ylim(0, 1.25)
    axB.set_title(f"Why `robust` rejects {n_filtered} stable WS vacancies")
    axB.legend(frameon=False, fontsize=8.5, loc="center right")

    fig.tight_layout()
    png = os.path.join(figdir, "hea_presets.png")
    pdf = os.path.join(figdir, "hea_presets.pdf")
    fig.savefig(png, dpi=150, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    print(f"WS={rob['ws']} relaxed={rel['hybrid']} robust={rob['hybrid']} "
          f"grid={rob['grid']} | filtered={n_filtered} (acc={n_acc}, rej={n_rej})")
    print("wrote", os.path.relpath(png, ROOT))
    print("wrote", os.path.relpath(pdf, ROOT))


if __name__ == "__main__":
    main()
