#!/usr/bin/env python3
"""
Defect-production efficiency of the FeCrNi PKA cascades vs the NRT and arc-dpa
displacement models.

  N_NRT(T)  = 0.8 * T / (2 E_d)                       (Norgett-Robinson-Torrens)
  xi_MD(T)  = N_surv(T) / N_NRT(T)                     (our MD efficiency)
  xi_arc(T) = (1 - c) * (0.8 T/(2 E_d))^b + c          (Nordlund et al. 2018)

For classical cascade MD with no electronic stopping the damage energy equals the
PKA energy, so T = E_PKA.  E_d = 40 eV (ASTM standard for Fe); arc-dpa Fe params
b = -0.568, c = 0.286 (Nordlund 2018, Nat. Commun. 9:1084).

Surviving Frenkel-pair counts are read from the ground-truth WS time series
(results/<label>.csv, plateau = median of the last 5 frames).  6 keV has two
independent cascades -> plotted as two points plus their mean.

Output: results/physics/dpa_efficiency.{png,pdf,csv}
Usage:  python3 scripts/analyze_dpa_efficiency.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/run/user/1000/gvfs/sftp:host=192.168.220.25/home/santi/pka_acero"
GT = os.path.join(BASE, "results")
OUT = os.path.join(BASE, "results", "physics")
os.makedirs(OUT, exist_ok=True)

E_D = 40.0            # eV, threshold displacement energy (Fe, ASTM)
B_ARC, C_ARC = -0.568, 0.286   # Nordlund 2018, Fe

# label -> (ground-truth csv, energy keV).  6 keV: two independent runs.
RUNS = [
    ("1 keV",        "1kev.csv",         1),
    ("2 keV",        "2kev.csv",         2),
    ("3 keV",        "3kev_results.csv", 3),
    ("4 keV",        "4kev.csv",         4),
    ("5 keV",        "5kev.csv",         5),
    ("6 keV (run 1)", "6kev.csv",        6),
    ("6 keV (run 2)", "6kev_it1.csv",    6),
    ("6 keV (run 3)", "6kev_it2.csv",    6),
    ("7 keV",        "7kev.csv",         7),
    ("8 keV",        "8kev_real.csv",    8),
]


def plateau(csv, ntail=5):
    df = pd.read_csv(os.path.join(GT, csv))
    v = pd.to_numeric(df["n_vacancies"], errors="coerce").dropna()
    return float(np.median(v.tail(ntail)))


def n_nrt(T_eV):
    return 0.8 * T_eV / (2.0 * E_D)


def xi_arc(T_eV):
    return (1.0 - C_ARC) * (n_nrt(T_eV)) ** B_ARC + C_ARC


def main():
    rows = []
    for name, csv, keV in RUNS:
        nsurv = plateau(csv)
        T = keV * 1000.0
        nnrt = n_nrt(T)
        rows.append(dict(run=name, keV=keV, T_eV=T, n_surv=nsurv,
                         n_nrt=nnrt, xi_md=nsurv / nnrt, xi_arc=xi_arc(T)))
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(OUT, "dpa_efficiency.csv"), index=False)
    print(tab.to_string(index=False,
          formatters={"T_eV": "{:.0f}".format, "n_nrt": "{:.1f}".format,
                      "xi_md": "{:.3f}".format, "xi_arc": "{:.3f}".format}))

    # mean over the two 6 keV runs for a guide-to-the-eye trend
    agg = tab.groupby("keV").agg(n_surv=("n_surv", "mean"),
                                 n_surv_sd=("n_surv", "std"),
                                 xi_md=("xi_md", "mean")).reset_index()
    agg["n_surv_sd"] = agg["n_surv_sd"].fillna(0.0)

    Egrid = np.linspace(0.5, 8.5, 400) * 1000.0

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.4))

    # ---- panel 1: defect counts vs energy ----
    ax1.plot(Egrid / 1000, n_nrt(Egrid), "-", color="0.3", lw=1.6,
             label="NRT  (0.8 T / 2$E_d$)")
    ax1.plot(Egrid / 1000, n_nrt(Egrid) * xi_arc(Egrid), "--", color="tab:green",
             lw=1.8, label="arc-dpa (Nordlund 2018, Fe)")
    ax1.errorbar(agg.keV, agg.n_surv, yerr=agg.n_surv_sd, fmt="o", ms=8,
                 color="tab:blue", capsize=4, lw=1.5, label="MD survivors (WS)")
    ax1.set_xlabel("PKA energy (keV)")
    ax1.set_ylabel("surviving Frenkel pairs $N_d$")
    ax1.set_title("Defect production vs displacement models")
    ax1.grid(alpha=0.3); ax1.legend()

    # ---- panel 2: efficiency vs energy ----
    ax2.plot(Egrid / 1000, xi_arc(Egrid), "--", color="tab:green", lw=1.8,
             label="arc-dpa efficiency (Fe)")
    ax2.axhline(C_ARC, color="0.6", ls=":", lw=1,
                label=f"arc-dpa high-E limit $c$={C_ARC}")
    for _, r in tab.iterrows():
        ax2.plot(r.keV, r.xi_md, "o", ms=8, color="tab:blue")
    ax2.plot([], [], "o", color="tab:blue", label="MD efficiency $\\xi=N_d/N_{NRT}$")
    ax2.plot(agg.keV, agg.xi_md, "-", color="tab:blue", lw=1, alpha=0.5)
    ax2.set_xlabel("PKA energy (keV)")
    ax2.set_ylabel("defect-production efficiency $\\xi$")
    ax2.set_ylim(0, 1.05)
    ax2.set_title(f"Efficiency vs NRT  ($E_d$={E_D:.0f} eV)")
    ax2.grid(alpha=0.3); ax2.legend()

    fig.suptitle("FeCrNi PKA cascades — defect-production efficiency "
                 "(NRT / arc-dpa)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"dpa_efficiency.{ext}"), dpi=140)
    print(f"\nWrote: {os.path.join(OUT, 'dpa_efficiency.png')} (+pdf, +csv)")


if __name__ == "__main__":
    main()
