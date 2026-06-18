#!/usr/bin/env python3
"""
plot_survival_overview.py — Figura resumen de la dinámica de supervivencia de
vacancias en las cascadas pka_acero trackeadas.

Paneles:
  (a) población de vacancias WS vivas vs tiempo, por cascada (semilog)
  (b) conteo en el pico vs sobrevivientes finales en función de la energía
      del PKA, con la eficiencia de supervivencia anotada
  (c) curva de supervivencia de la población del pico (fracción aún viva
      vs Δt desde el pico), por cascada

Requiere los outputs de run_pka_tracking.py.
Salida: figures/survival_overview.png
"""
import glob
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

from _figstyle import setup, save, C

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))

ENERGY_KEV = {"1kev": 1, "2kev": 2, "3kev": 3, "4kev": 4, "5kev": 5,
              "6kev": 6, "6kev_it1": 6, "6kev_it2": 6, "7kev": 7,
              "8kev_real": 8}
RESOLVED = ["4kev", "5kev", "6kev", "6kev_it1", "6kev_it2", "7kev", "8kev_real"]


def main():
    tracking = os.path.join(ROOT, "results/pka_tracking")
    data = {}
    for en in RESOLVED:
        f_tr = os.path.join(tracking, en, "tracks.csv")
        f_pt = os.path.join(tracking, en, "tracks_points.csv")
        if not (os.path.exists(f_tr) and os.path.exists(f_pt)):
            continue
        data[en] = (pd.read_csv(f_tr), pd.read_csv(f_pt))

    colors = {en: cm.viridis(i / max(1, len(data) - 1))
              for i, en in enumerate(data)}

    setup()
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))

    # (a) población viva vs tiempo
    ax = axes[0]
    for en, (tr, pts) in data.items():
        pop = pts.groupby("step").size()
        ax.plot(pop.index / 1000, pop.values, "o-", ms=3.5, lw=1.4,
                color=colors[en], label=en)
    ax.set_yscale("log")
    ax.set_xlabel("timestep [×1000]")
    ax.set_ylabel("vacancias WS vivas")
    ax.set_title("(a) Población de vacancias por cascada")
    ax.legend(fontsize=7, frameon=False, ncol=2)
    ax.grid(alpha=0.25)

    # (b) pico vs final vs energía
    ax = axes[1]
    rows = []
    for en, (tr, pts) in data.items():
        peak_step = int(pts.groupby("step").size().idxmax())
        n_peak = (pts.step == peak_step).sum()
        n_fin = int(tr.survived.sum())
        rows.append((ENERGY_KEV[en], n_peak, n_fin, en))
    rows.sort()
    E = [r[0] for r in rows]
    ax.plot(E, [r[1] for r in rows], "o", ms=8, color=C["transient"],
            label="pico de daño")
    ax.plot(E, [r[2] for r in rows], "s", ms=8, color=C["survive"],
            label="sobreviven al final")
    for e, npk, nfn, _ in rows:
        ax.annotate(f"{100*nfn/npk:.0f}%", (e, nfn), textcoords="offset points",
                    xytext=(0, -14), ha="center", fontsize=7, color=C["survive"])
    ax.set_yscale("log")
    ax.set_xlabel("energía del PKA [keV]")
    ax.set_ylabel("vacancias WS")
    ax.set_title("(b) Pico vs supervivientes (eficiencia anotada)")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25)

    # (c) curva de supervivencia de la población del pico
    ax = axes[2]
    for en, (tr, pts) in data.items():
        peak_step = int(pts.groupby("step").size().idxmax())
        peak_tracks = set(pts[pts.step == peak_step].track_id)
        steps = sorted(s for s in pts.step.unique() if s >= peak_step)
        frac = [len(peak_tracks & set(pts[pts.step == s].track_id))
                / len(peak_tracks) for s in steps]
        ax.plot((np.array(steps) - peak_step) / 1000, frac, "o-", ms=3.5,
                lw=1.4, color=colors[en], label=en)
    ax.set_yscale("log")
    ax.set_xlabel("Δt desde el pico [×1000 pasos]")
    ax.set_ylabel("fracción de la población del pico aún viva")
    ax.set_title("(c) Recombinación de la población del pico")
    ax.grid(alpha=0.25)

    fig.suptitle("Supervivencia de vacancias en cascadas FeCrNi (tracking WS temporal)",
                 y=1.03, fontsize=13)
    fig.tight_layout()
    out = save(fig, "survival_overview")
    print("Figura:", os.path.relpath(out, ROOT))


if __name__ == "__main__":
    main()
