#!/usr/bin/env python3
"""
plot_cascade_anatomy.py — Anatomía espacial de una cascada con el destino de
cada vacancia del pico, más las distribuciones pooled de las dos features
espaciales que separan sobrevivientes de transitorios:

  (a) proyección 2D del pico de daño: vacancias transitorias (mueren),
      sobrevivientes con su trayectoria hasta la posición final, e
      intersticiales WS del pico — muestra la estructura core-shell
  (b) distancia al centroide de la nube de vacancias (surv vs doomed, pooled)
  (c) densidad local de vacancias a <10 Å (ídem)

Requiere los outputs de run_pka_tracking.py.
Uso:
    python3 scripts/plot_cascade_anatomy.py [--energy 8kev_real]
Salida: figures/cascade_anatomy_<energy>.png
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _figstyle import setup, save, C

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))

BOX_L = 243.74
RESOLVED = ["4kev", "5kev", "6kev", "6kev_it1", "6kev_it2", "7kev", "8kev_real"]
WS_SITES_COLS = ["ref_id", "x", "y", "z", "occupancy", "min_dist", "site_type"]


def unwrap(pos, ref):
    d = pos - ref
    d -= BOX_L * np.round(d / BOX_L)
    return ref + d


def spatial_features(pts, tracks, peak):
    """dist-al-centroide y densidad local (<10 Å) por vacancia del pico."""
    P = pts[pts.step == peak]
    y = np.array([tracks.loc[t].survived for t in P.track_id], float)
    pos = unwrap(P[["x", "y", "z"]].values, P[["x", "y", "z"]].values[0])
    cent = pos.mean(0)
    dcent = np.sqrt(((pos - cent) ** 2).sum(1))
    dm = np.sqrt(((pos[:, None, :] - pos[None, :, :]) ** 2).sum(2))
    locdens = (dm < 10).sum(1) - 1
    return y, dcent, locdens, pos, cent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--energy", default="8kev_real")
    args = ap.parse_args()
    tracking = os.path.join(ROOT, "results/pka_tracking")

    # ── datos de la cascada protagonista ─────────────────────────────────────
    en = args.energy
    tracks = pd.read_csv(f"{tracking}/{en}/tracks.csv").set_index("track_id")
    pts = pd.read_csv(f"{tracking}/{en}/tracks_points.csv")
    peak = int(pts.groupby("step").size().idxmax())
    y, dcent, locdens, pos, cent = spatial_features(pts, tracks, peak)
    P = pts[pts.step == peak]

    # intersticiales WS del pico
    sites = pd.read_csv(f"{tracking}/{en}/ws_sites_t{peak:06d}.csv",
                        sep=r"\s+", skiprows=1, names=WS_SITES_COLS)
    ints = sites[sites.occupancy >= 2][["x", "y", "z"]].values
    ints = unwrap(ints, pos[0])

    # ejes de proyección: los dos de mayor varianza de la nube de vacancias
    var = pos.var(0)
    ax1, ax2 = np.argsort(var)[::-1][:2]
    lbl = "xyz"

    setup()
    fig = plt.figure(figsize=(15.5, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1, 1])

    # (a) proyección
    ax = fig.add_subplot(gs[0])
    ax.scatter(ints[:, ax1], ints[:, ax2], marker="+", s=14, c="#7fa8d9",
               alpha=0.5, label=f"intersticiales WS ({len(ints)})")
    doom = pos[y == 0]
    ax.scatter(doom[:, ax1], doom[:, ax2], s=12, c=C["transient"], alpha=0.6,
               label=f"vacancias transitorias ({int((y==0).sum())})")
    surv_ids = [t for t, s in zip(P.track_id, y) if s == 1]
    first = True
    for tid in surv_ids:
        tp = pts[pts.track_id == tid].sort_values("step")
        tr_pos = unwrap(tp[["x", "y", "z"]].values, pos[0])
        ax.plot(tr_pos[:, ax1], tr_pos[:, ax2], "-", lw=1.4, color=C["traj"],
                alpha=0.8)
        ax.scatter(tr_pos[0, ax1], tr_pos[0, ax2], marker="*", s=160,
                   c=C["survive"], edgecolors="k", linewidths=0.6, zorder=5,
                   label=f"sobreviven ({len(surv_ids)}) + trayectoria"
                         if first else None)
        ax.scatter(tr_pos[-1, ax1], tr_pos[-1, ax2], marker="o", s=36,
                   facecolors="none", edgecolors=C["traj"], zorder=5)
        first = False
    ax.scatter(*[[cent[ax1]], [cent[ax2]]], marker="x", s=90, c="k",
               label="centroide de la nube")
    ax.set_xlabel(f"{lbl[ax1]} [Å]")
    ax.set_ylabel(f"{lbl[ax2]} [Å]")
    ax.set_title(f"(a) Pico de daño de {en} (t={peak}) — destino de cada vacancia")
    ax.legend(fontsize=8, frameon=False, loc="best")
    ax.set_aspect("equal")

    # ── features pooled sobre todas las cascadas resueltas ───────────────────
    DC, LD, YY = [], [], []
    for e2 in RESOLVED:
        f_tr = f"{tracking}/{e2}/tracks.csv"
        f_pt = f"{tracking}/{e2}/tracks_points.csv"
        if not (os.path.exists(f_tr) and os.path.exists(f_pt)):
            continue
        tr2 = pd.read_csv(f_tr).set_index("track_id")
        pt2 = pd.read_csv(f_pt)
        pk2 = int(pt2.groupby("step").size().idxmax())
        y2, dc2, ld2, _, _ = spatial_features(pt2, tr2, pk2)
        YY.append(y2); DC.append(dc2); LD.append(ld2)
    YY = np.concatenate(YY); DC = np.concatenate(DC); LD = np.concatenate(LD)

    def panel(ax, vals, xlabel, title, bins):
        ax.hist(vals[YY == 0], bins=bins, density=True, alpha=0.6,
                color=C["transient"], label="transitorias")
        ax.hist(vals[YY == 1], bins=bins, density=True, alpha=0.75,
                color=C["survive"], label="sobreviven")
        ax.axvline(np.median(vals[YY == 0]), color="#b23a22", ls="--", lw=1.2)
        ax.axvline(np.median(vals[YY == 1]), color=C["traj"], ls="--", lw=1.2)
        ax.set_xlabel(xlabel)
        ax.set_title(title)
        ax.legend(fontsize=8, frameon=False)

    panel(fig.add_subplot(gs[1]), DC,
          "distancia al centroide de la nube [Å]",
          "(b) Sobrevivientes: en el núcleo del daño",
          np.linspace(0, np.percentile(DC, 99), 30))
    panel(fig.add_subplot(gs[2]), LD.astype(float),
          "vacancias vecinas a < 10 Å",
          "(c) Sobrevivientes: en zonas densas",
          np.arange(0, np.percentile(LD, 99.5), 8))

    fig.suptitle("Anatomía de cascada: el daño que sobrevive vive en el núcleo "
                 "denso de la nube de vacancias (FeCrNi 4–8 keV)",
                 y=1.02, fontsize=12.5)
    fig.tight_layout()
    out = save(fig, f"cascade_anatomy_{en}")
    print("Figura:", os.path.relpath(out, ROOT))


if __name__ == "__main__":
    main()
