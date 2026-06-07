#!/usr/bin/env python3
"""
Damage-evolution figure for the FeCrNi 5 keV cascade (README §3.9).

Data extracted from the `Vacancy reconciliation` blocks of the run logs in
results/5kev_20_05/ (ball_robust.log, evol/t*.log, final_robust.log).
Reproduces README §3.9: WS over-counts transient Frenkel pairs in the ballistic
regime while the hybrid `robust` preset filters them; both converge to the 6
surviving pairs once the cascade relaxes (t >= 15000).

Usage:
    python3 scripts/plot_damage_evolution.py
Outputs:
    figures/damage_evolution.png  (+ .pdf)
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# timestep : (WS_vac, hybrid_robust_accepted, WS_filtered, grid_voids, grid_volume_A3)
DATA = {
    5000:   (424, 79, 345, 57, 297.0),
    10000:  (152, 20, 132,  1,   3.4),
    15000:  (  6,  6,   0,  0,   0.0),
    20000:  (  6,  6,   0,  0,   0.0),
    25000:  (  6,  6,   0,  0,   0.0),
    30000:  (  6,  6,   0,  0,   0.0),
    35000:  (  6,  6,   0,  0,   0.0),
    40000:  (  6,  6,   0,  0,   0.0),
    45000:  (  6,  6,   0,  0,   0.0),
    100000: (  6,  6,   0,  0,   0.0),
}

SURVIVORS = 6  # stable Frenkel pairs after relaxation

def main():
    ts       = sorted(DATA)
    ws       = [DATA[t][0] for t in ts]
    hybrid   = [DATA[t][1] for t in ts]
    filtered = [DATA[t][2] for t in ts]
    grid     = [DATA[t][3] for t in ts]

    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(here, "..", "figures")
    os.makedirs(figdir, exist_ok=True)

    fig, (ax, axz) = plt.subplots(
        1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [2.0, 1.0]})

    # ---- main panel: full damage evolution -------------------------------
    ax.plot(ts, ws,     "o-",  color="#c0392b", lw=2, ms=7,
            label="Wigner-Seitz (topological count)")
    ax.plot(ts, hybrid, "s-",  color="#1f6fb2", lw=2, ms=7,
            label="Hybrid  `robust`  (stable defects)")
    ax.plot(ts, grid,   "^--", color="#27ae60", lw=1.6, ms=6,
            label="Grid open-void clusters")
    ax.axhline(SURVIVORS, color="0.4", ls=":", lw=1.4)
    ax.annotate(f"{SURVIVORS} surviving Frenkel pairs",
                xy=(45000, SURVIVORS), xytext=(46000, 70),
                fontsize=9, color="0.3",
                arrowprops=dict(arrowstyle="->", color="0.5"))

    # shade the ballistic regime
    ax.axvspan(0, 12500, color="#f4d03f", alpha=0.18)
    ax.text(13500, 360,
            "ballistic / peak-damage\n(WS over-counts transients)",
            fontsize=8.5, color="#7d6608", ha="left", va="top")

    ax.set_xlabel("MD timestep after PKA")
    ax.set_ylabel("Defect count")
    ax.set_title("FeCrNi 5 keV cascade — vacancy estimators vs. regime")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.grid(alpha=0.25)
    ax.margins(x=0.02)

    # ---- inset/zoom panel: relaxed plateau -------------------------------
    mask = [t >= 15000 for t in ts]
    tz   = [t for t, m in zip(ts, mask) if m]
    axz.plot(tz, [v for v, m in zip(ws, mask) if m],     "o-", color="#c0392b", lw=2, ms=7)
    axz.plot(tz, [v for v, m in zip(hybrid, mask) if m], "s-", color="#1f6fb2", lw=2, ms=6)
    axz.axhline(SURVIVORS, color="0.4", ls=":", lw=1.4)
    axz.set_ylim(0, 12)
    axz.set_xlabel("MD timestep")
    axz.set_title("Relaxed plateau:\nWS = hybrid = 6 (exact)", fontsize=10)
    axz.grid(alpha=0.25)
    axz.text(0.5, 0.10, "penalties don't fire\n(no atoms in transit)",
             transform=axz.transAxes, ha="center", fontsize=8.5, color="0.35")

    fig.tight_layout()
    png = os.path.join(figdir, "damage_evolution.png")
    pdf = os.path.join(figdir, "damage_evolution.pdf")
    fig.savefig(png, dpi=150, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    print("wrote", os.path.normpath(png))
    print("wrote", os.path.normpath(pdf))

if __name__ == "__main__":
    main()
