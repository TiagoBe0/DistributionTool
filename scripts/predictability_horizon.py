#!/usr/bin/env python3
"""
predictability_horizon.py — ¿Cuándo "se decide" el daño de una cascada?

Para cada cascada trackeada y cada frame t ≥ pico, calcula el AUC de predecir
la supervivencia final de las vacancias VIVAS en t usando (a) la señal de
topología y (b) el score del consenso robust, ambos del breakdown híbrido de
ese frame. El resultado es la curva de "horizonte de predictibilidad":
qué tan temprano puede saberse qué defectos van a sobrevivir.

Requiere los outputs de run_pka_tracking.py.

Uso:
    python3 scripts/predictability_horizon.py [--tracking results/pka_tracking]
Salidas:
    results/pka_tracking/predictability_horizon.csv
    figures/predictability_horizon.png
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))

BOX_L = 243.74   # lado del box cúbico compartido por todas las cascadas pka_acero
SIGNALS = ["ws", "soft_ws", "vor", "dens", "soap", "topo",
           "transit_pen", "frenkel_pen"]
EXCLUDE = {"1kev", "2kev", "3kev"}   # sin fase de recombinación resuelta


def read_hv(path):
    """Lee un hybrid_vacancies CSV usando su header comentado (tolera tanto el
    formato viejo de 16 columnas como el nuevo con centrality/cloud_dens)."""
    with open(path) as f:
        names = f.readline().lstrip("#").split()
    df = pd.read_csv(path, sep=r"\s+", skiprows=1, names=names)
    return df.rename(columns={"consensus_score": "score"})


def auc(y, s):
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    return float((pos[:, None] > neg[None, :]).mean()
                 + 0.5 * (pos[:, None] == neg[None, :]).mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tracking", default=os.path.join(ROOT, "results/pka_tracking"))
    args = ap.parse_args()

    rows = []
    for endir in sorted(glob.glob(os.path.join(args.tracking, "*"))):
        en = os.path.basename(endir)
        if en in EXCLUDE or not os.path.isdir(endir):
            continue
        f_tracks = os.path.join(endir, "tracks.csv")
        f_points = os.path.join(endir, "tracks_points.csv")
        if not (os.path.exists(f_tracks) and os.path.exists(f_points)):
            continue
        tracks = pd.read_csv(f_tracks).set_index("track_id")
        pts = pd.read_csv(f_points)
        peak = int(pts.groupby("step").size().idxmax())

        for step in sorted(pts.step.unique()):
            if step < peak:
                continue
            f_hv = os.path.join(endir, f"hybrid_vacancies_t{step:06d}.csv")
            if not os.path.exists(f_hv):
                continue
            hv = read_hv(f_hv)
            hv = hv[hv.ref_site_idx >= 0].reset_index(drop=True)
            P = pts[pts.step == step]
            if len(P) < 4:
                continue
            hp = hv[["x", "y", "z"]].values
            ys, topo, score = [], [], []
            for _, p in P.iterrows():
                d = hp - [p.x, p.y, p.z]
                d -= BOX_L * np.round(d / BOX_L)
                r2 = (d ** 2).sum(1)
                i = int(r2.argmin())
                if r2[i] > 0.01:
                    continue
                ys.append(int(tracks.loc[p.track_id].survived))
                topo.append(hv.iloc[i].topo)
                score.append(hv.iloc[i].score)
            ys = np.array(ys, float)
            if not (0 < ys.sum() < len(ys)):
                continue   # una sola clase viva: AUC indefinido
            rows.append(dict(energy=en, step=step, dt=step - peak, n=len(ys),
                             pos=int(ys.sum()),
                             auc_topo=auc(ys, np.array(topo)),
                             auc_score=auc(ys, np.array(score))))

    df = pd.DataFrame(rows)
    out_csv = os.path.join(args.tracking, "predictability_horizon.csv")
    df.to_csv(out_csv, index=False)
    agg = df.groupby("dt").agg(n_casc=("energy", "size"),
                               topo=("auc_topo", "mean"),
                               score=("auc_score", "mean"))
    print(agg.round(3).to_string())

    fig, ax = plt.subplots(figsize=(8, 5))
    for en, g in df.groupby("energy"):
        ax.plot(g.dt / 1000, g.auc_score, "o", ms=4, alpha=0.35, color="#c0392b")
        ax.plot(g.dt / 1000, g.auc_topo, "s", ms=4, alpha=0.35, color="#2a9d8f")
    ax.plot(agg.index / 1000, agg.score, "-", lw=2.5, color="#c0392b",
            label="score consenso (robust)")
    ax.plot(agg.index / 1000, agg.topo, "-", lw=2.5, color="#2a9d8f",
            label="topología sola")
    ax.axhline(0.5, ls="--", color="0.5", lw=1, label="azar")
    ax.set_xlabel("Δt desde el pico de daño [×1000 pasos]")
    ax.set_ylabel("AUC — predicción de supervivencia de vacancias vivas")
    ax.set_ylim(0, 1.05)
    ax.set_title("Horizonte de predictibilidad del daño superviviente\n"
                 "(cascadas FeCrNi 4–8 keV; puntos = cascadas, línea = media)")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out_fig = os.path.join(ROOT, "figures", "predictability_horizon.png")
    fig.savefig(out_fig, dpi=150, bbox_inches="tight")
    print(f"\nFigura: {os.path.relpath(out_fig, ROOT)}")
    print(f"CSV:    {os.path.relpath(out_csv, ROOT)}")


if __name__ == "__main__":
    main()
