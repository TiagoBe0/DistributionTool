#!/usr/bin/env python3
"""
Why the hybrid estimator is useful: during the ballistic peak and the thermal-
spike cooling, raw Wigner-Seitz counts every transiently displaced atom as a
Frenkel pair and massively over-counts the damage.  The hybrid consensus filters
those transients, so it tracks the *surviving* defect level far earlier in the
cascade than raw WS does.

For every energy whose sweep captured the ballistic phase (1,4,5,6,7,8 keV) we
read results/distool_sweep/<label>_hybrid.csv and report:
  - ballistic peak WS count and the hybrid count at the same frame,
  - the surviving level (median of the last 3 hybrid frames),
  - the MD step at which WS vs hybrid first settle within tolerance of survivors.

Output: results/physics/hybrid_vs_ws_transient.{png,pdf,csv}
Usage:  python3 scripts/analyze_hybrid_at_peak.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/run/user/1000/gvfs/sftp:host=192.168.220.25/home/santi/pka_acero"
SWEEP = os.path.join(BASE, "results", "distool_sweep")
OUT = os.path.join(BASE, "results", "physics")
os.makedirs(OUT, exist_ok=True)

# only energies whose sweep includes the early ballistic frames
ENERGIES = [("1kev", 1), ("4kev", 4), ("5kev", 5),
            ("6kev", 6), ("7kev", 7), ("8kev_real", 8)]


def load(label):
    df = pd.read_csv(os.path.join(SWEEP, f"{label}_hybrid.csv"))
    for c in ("ws_vac", "hybrid_vac", "hybrid_filtered"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.dropna(subset=["ws_vac", "hybrid_vac"]).sort_values("step").reset_index(drop=True)


def settle_step(df, col, surv, peak_step):
    """First post-peak step from which the count stays within tol of survivors."""
    tol = max(1.0, round(0.25 * surv))
    d = df[df.step >= peak_step].reset_index(drop=True)
    within = (np.abs(d[col] - surv) <= tol).values
    for i in range(len(within)):
        if within[i:].all():
            return int(d.step[i])
    return np.nan


def main():
    rows, frames = [], {}
    for label, keV in ENERGIES:
        df = load(label)
        frames[keV] = df
        surv = float(np.median(df.hybrid_vac.tail(3)))
        ipk = int(df.ws_vac.idxmax())
        pk_step = int(df.step[ipk])
        ws_oc = float(df.ws_vac[ipk]) / surv
        hyb_oc = float(df.hybrid_vac[ipk]) / surv
        rows.append(dict(
            keV=keV, peak_step=pk_step,
            ws_peak=float(df.ws_vac[ipk]),
            hybrid_at_peak=float(df.hybrid_vac[ipk]),
            survivors=surv,
            ws_overcount=ws_oc,
            hyb_overcount=hyb_oc,
            hyb_closer_x=ws_oc / hyb_oc,            # how many times closer hybrid is
            settle_step=settle_step(df, "hybrid_vac", surv, pk_step),
        ))
    tab = pd.DataFrame(rows).sort_values("keV")
    tab.to_csv(os.path.join(OUT, "hybrid_vs_ws_transient.csv"), index=False)
    print(tab.to_string(index=False))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.6))

    # ---- panel 1: representative trajectory (5 keV) ----
    keV_demo = 5
    d = frames[keV_demo]
    surv = float(np.median(d.hybrid_vac.tail(3)))
    ax1.plot(d.step / 1000, d.ws_vac, "-o", ms=4, color="tab:red",
             label="raw Wigner-Seitz")
    ax1.plot(d.step / 1000, d.hybrid_vac, "-s", ms=4, color="tab:blue",
             label="hybrid (transients filtered)")
    ax1.axhline(surv, color="0.4", ls="--", lw=1.2,
                label=f"surviving level = {surv:.0f}")
    ax1.fill_between(d.step / 1000, d.hybrid_vac, d.ws_vac,
                     color="tab:orange", alpha=0.18, label="filtered transients")
    ax1.set_yscale("log")
    ax1.set_xlabel("MD step (×1000)")
    ax1.set_ylabel("Frenkel pairs (log)")
    ax1.set_title(f"{keV_demo} keV cascade — WS over-counts the transient")
    ax1.grid(alpha=0.3, which="both"); ax1.legend(fontsize=9)

    # ---- panel 2: peak over-count factor per energy ----
    x = np.arange(len(tab)); w = 0.38
    ax2.bar(x - w/2, tab.ws_overcount, w, color="tab:red",
            label="raw WS / survivors")
    ax2.bar(x + w/2, tab.hyb_overcount, w, color="tab:blue",
            label="hybrid / survivors")
    ax2.axhline(1, color="0.3", ls=":", lw=1)
    ax2.set_yscale("log")
    ax2.set_xticks(x); ax2.set_xticklabels([f"{k}" for k in tab.keV])
    ax2.set_xlabel("PKA energy (keV)")
    ax2.set_ylabel("over-count factor at ballistic peak (log)")
    ax2.set_title("Over-counting at the peak frame")
    ax2.grid(alpha=0.3, axis="y", which="both"); ax2.legend(fontsize=9)

    fig.suptitle("FeCrNi cascades — hybrid filtering of ballistic transients "
                 "vs raw Wigner-Seitz", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"hybrid_vs_ws_transient.{ext}"), dpi=140)
    print(f"\nWrote: {os.path.join(OUT, 'hybrid_vs_ws_transient.png')} (+pdf, +csv)")


if __name__ == "__main__":
    main()
