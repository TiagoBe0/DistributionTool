#!/usr/bin/env python3
"""
Single-cascade scatter at 6 keV: three independent FeNiCr cascades.

The 1-cascade-per-energy survivors-vs-energy curve shows a 6>7 keV "inversion"
(10 vs 9) that is statistical, not physical.  We have THREE independent 6 keV
cascades, so we can show the real spread and that the energy trend is monotonic
once the 6 keV point carries an error bar.

  run 1  -> results/6kev.csv        (WS ground truth)
  run 2  -> results/6kev_it1.csv    (WS ground truth)
  run 3  -> results/6kev_it2.csv    (WS from the validated distool sweep)

run 3 additionally has a hybrid count of 14 vs WS 18 at the plateau: the hybrid
discounts 3-4 close (unrecombined) Frenkel pairs (high frenkel_pen, in_recomb=1).

Output: results/physics/sixkev_statistics.{png,pdf,csv}
Usage:  python3 scripts/analyze_6kev_statistics.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_CANDS = ["/home/santi/pka_acero",
          "/run/user/1000/gvfs/sftp:host=192.168.220.25/home/santi/pka_acero"]
BASE = next((p for p in _CANDS if os.path.isdir(p)), _CANDS[0])
GT = os.path.join(BASE, "results")
OUT = os.path.join(GT, "physics")
os.makedirs(OUT, exist_ok=True)

# full single-cascade series (survivors) for the energy trend
SERIES = [("1kev", "1kev.csv", 1), ("2kev", "2kev.csv", 2),
          ("3kev", "3kev_results.csv", 3), ("4kev", "4kev.csv", 4),
          ("5kev", "5kev.csv", 5), ("7kev", "7kev.csv", 7),
          ("8kev_real", "8kev_real.csv", 8)]
SIX = [("run 1", "6kev.csv"), ("run 2", "6kev_it1.csv"), ("run 3", "6kev_it2.csv")]
HYB6 = {"run 3": 14}  # hybrid (close-pair-corrected) where it differs from WS


def plateau(csv, ntail=5):
    df = pd.read_csv(os.path.join(GT, csv))
    v = pd.to_numeric(df["n_vacancies"], errors="coerce").dropna()
    return float(np.median(v.tail(ntail)))


def main():
    six_vals = [(name, plateau(csv)) for name, csv in SIX]
    sv = np.array([v for _, v in six_vals])
    mean6, sd6 = float(sv.mean()), float(sv.std(ddof=1))
    pd.DataFrame(dict(run=[n for n, _ in six_vals], ws_survivors=sv)).to_csv(
        os.path.join(OUT, "sixkev_statistics.csv"), index=False)
    print("6 keV WS survivors:", {n: v for n, v in six_vals},
          f"\n  mean = {mean6:.1f} ± {sd6:.1f}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.4))

    # ---- panel 1: the three 6 keV runs ----
    x = np.arange(len(six_vals))
    ax1.bar(x, sv, 0.55, color="tab:blue", label="WS survivors")
    ax1.axhspan(mean6 - sd6, mean6 + sd6, color="tab:blue", alpha=0.15)
    ax1.axhline(mean6, color="tab:blue", ls="--", lw=1.4,
                label=f"mean = {mean6:.1f} ± {sd6:.1f}")
    for xi, (name, _) in zip(x, six_vals):
        if name in HYB6:
            ax1.plot(xi, HYB6[name], "D", ms=9, color="tab:orange",
                     label="hybrid (close-pair corrected)")
    ax1.set_xticks(x); ax1.set_xticklabels([n for n, _ in six_vals])
    ax1.set_ylabel("surviving Frenkel pairs")
    ax1.set_title("Three independent 6 keV cascades — single-cascade scatter")
    ax1.grid(alpha=0.3, axis="y"); ax1.legend(fontsize=9)

    # ---- panel 2: energy trend with the 6 keV error bar ----
    es = [(keV, plateau(csv)) for _, csv, keV in SERIES]
    ek = [e for e, _ in es]; ev = [v for _, v in es]
    ax2.plot(ek, ev, "o", ms=8, color="0.4", label="single cascade (WS)")
    ax2.errorbar([6], [mean6], yerr=[sd6], fmt="s", ms=10, color="tab:blue",
                 capsize=5, lw=1.8, label="6 keV: 3 cascades (mean ± sd)")
    allk = sorted(ek + [6]); allv = dict(zip(ek, ev)); allv[6] = mean6
    ax2.plot(allk, [allv[k] for k in allk], "-", color="0.7", lw=1, zorder=0)
    ax2.set_xlabel("PKA energy (keV)")
    ax2.set_ylabel("surviving Frenkel pairs")
    ax2.set_title("Energy trend — the 6>7 keV 'inversion' is statistical")
    ax2.set_xticks(range(1, 9))
    ax2.grid(alpha=0.3); ax2.legend(fontsize=9)

    fig.suptitle("FeNiCr cascades — 6 keV multi-cascade statistics", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"sixkev_statistics.{ext}"), dpi=140)
    print(f"\nWrote: {os.path.join(OUT, 'sixkev_statistics.png')} (+pdf, +csv)")


if __name__ == "__main__":
    main()
