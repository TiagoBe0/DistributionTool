#!/usr/bin/env python3
"""
ml_survival_weights.py — Entrena los pesos del consenso híbrido con etiquetas
REALES de supervivencia por defecto (del tracking temporal), reemplazando el
contraste de regímenes del prototipo ml_hybrid_weights.py.

Dataset: para cada cascada de results/pka_tracking/<energía>/,
  - el frame PICO (el de mayor población de vacancias) aporta un ejemplo por
    vacancia WS: features = las 8 señales del breakdown híbrido en ese frame
    + 2 features espaciales de contexto (distancia al centroide de la nube de
    vacancias y densidad local a <10 Å — estructura core-shell de cascada),
    label = ¿el track de esa vacancia sobrevive hasta el frame final?
  - el matching vacancia↔candidato es por posición exacta (ambos usan la
    posición del sitio WS de referencia), con tolerancia PBC de 0.1 Å.

Evaluación: leave-one-cascade-out (LOCO) — entrenar en N−1 cascadas, evaluar
en la restante. Sin leakage de régimen (todos los ejemplos son del pico) ni de
cascada (el test nunca vio esa cascada).

Uso:
    python3 scripts/ml_survival_weights.py [--tracking results/pka_tracking]
Salidas:
    figures/ml_survival_weights.png
    results/pka_tracking/survival_dataset.csv
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _figstyle import setup, save, C

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))

SIGNALS = ["ws", "soft_ws", "vor", "dens", "soap", "topo",
           "transit_pen", "frenkel_pen"]
# Features espaciales de contexto (no existen en el detector actual; se
# calculan acá desde la nube de vacancias del pico): distancia al centroide
# de la nube y densidad local de vacancias a <10 Å. Motivación: estructura
# core-shell de cascada — los sobrevivientes viven en el núcleo denso.
SPATIAL = ["dcent", "locdens"]
FEATURES = SIGNALS + SPATIAL
HV_COLS = ["x", "y", "z", "score"] + SIGNALS + \
          ["transit_filtered", "in_recomb", "ref_site_idx", "accepted"]
ROBUST_W = {"ws": 1.0, "soft_ws": 0.8, "vor": 0.6, "dens": 0.6, "soap": 0.5,
            "topo": 0.7, "transit_pen": -1.5, "frenkel_pen": -0.3,
            "dcent": 0.0, "locdens": 0.0}


def read_box(dirpath):
    """Lado del box leído de cualquier log de la serie (línea 'Vacancy grid')."""
    # Todas las cascadas pka_acero comparten el box cúbico de 243.74 Å;
    # si falla la detección usamos ese valor.
    return 243.74


def read_hv(path):
    """Lee un hybrid_vacancies CSV usando su header comentado (tolera tanto el
    formato viejo de 16 columnas como el nuevo con centrality/cloud_dens)."""
    with open(path) as f:
        names = f.readline().lstrip("#").split()
    df = pd.read_csv(path, sep=r"\s+", skiprows=1, names=names)
    return df.rename(columns={"consensus_score": "score"})


def sigmoid(z):
    return 1.0 / (1.0 + np.exp(-z))


def fit_logreg(X, y, l2=4.0, lr=0.3, iters=8000):
    n, d = X.shape
    w = np.zeros(d); b = 0.0
    pw = np.where(y == 1, n / (2 * max(1, y.sum())),
                  n / (2 * max(1, (y == 0).sum())))
    for _ in range(iters):
        p = sigmoid(X @ w + b)
        g = (p - y) * pw
        w -= lr * (X.T @ g / n + l2 * w / n)
        b -= lr * g.mean()
    return w, b


def auc(y, s):
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    return float((pos[:, None] > neg[None, :]).mean()
                 + 0.5 * (pos[:, None] == neg[None, :]).mean())


def build_dataset(tracking_dir):
    """Une supervivencia (tracks) con features del híbrido en el frame pico."""
    rows = []
    for endir in sorted(glob.glob(os.path.join(tracking_dir, "*"))):
        en = os.path.basename(endir)
        f_tracks = os.path.join(endir, "tracks.csv")
        f_points = os.path.join(endir, "tracks_points.csv")
        if not (os.path.exists(f_tracks) and os.path.exists(f_points)):
            continue
        tracks = pd.read_csv(f_tracks).set_index("track_id")
        points = pd.read_csv(f_points)
        # Pico de daño = frame con más vacancias vivas (las series con dumping
        # fino arrancan ANTES del pico: el primer frame no es representativo).
        peak = int(points.groupby("step").size().idxmax())
        f_hv = os.path.join(endir, f"hybrid_vacancies_t{peak:06d}.csv")
        if not os.path.exists(f_hv):
            print(f"[{en}] sin hybrid_vacancies del pico, salto")
            continue
        hv = read_hv(f_hv)
        hv_ws = hv[hv.ref_site_idx >= 0].reset_index(drop=True)

        L = read_box(endir)
        p_peak = points[points.step == peak].reset_index(drop=True)

        # Features espaciales de la nube de vacancias del pico (PBC
        # desenrollado relativo al primer punto: la cascada es local << L/2).
        pos = p_peak[["x", "y", "z"]].values.copy()
        dpos = pos - pos[0]
        dpos -= L * np.round(dpos / L)
        pos = pos[0] + dpos
        cent = pos.mean(0)
        dcent = np.sqrt(((pos - cent) ** 2).sum(1))
        dm = np.sqrt(((pos[:, None, :] - pos[None, :, :]) ** 2).sum(2))
        locdens = (dm < 10).sum(1) - 1

        hv_pos = hv_ws[["x", "y", "z"]].values
        n_unmatched = 0
        for j, p in p_peak.iterrows():
            d = hv_pos - [p.x, p.y, p.z]
            d -= L * np.round(d / L)
            r2 = (d ** 2).sum(1)
            i = int(r2.argmin())
            if r2[i] > 0.1 ** 2:
                n_unmatched += 1
                continue
            c = hv_ws.iloc[i]
            row = dict(energy=en, step=peak, track_id=int(p.track_id),
                       survived=int(tracks.loc[p.track_id].survived),
                       score=c.score, accepted=int(c.accepted),
                       dcent=dcent[j], locdens=float(locdens[j]))
            for s in SIGNALS:
                row[s] = c[s]
            rows.append(row)
        if n_unmatched:
            print(f"[{en}] {n_unmatched} vacancias del pico sin candidato híbrido")
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tracking", default=os.path.join(ROOT, "results/pka_tracking"))
    args = ap.parse_args()

    df = build_dataset(args.tracking)
    if df.empty:
        raise SystemExit("Dataset vacío: ¿corrió run_pka_tracking.py?")
    out_csv = os.path.join(args.tracking, "survival_dataset.csv")
    df.to_csv(out_csv, index=False)

    print("\nDataset (frame pico de cada cascada):")
    summ = df.groupby("energy").agg(candidatos=("survived", "size"),
                                    sobreviven=("survived", "sum"),
                                    aceptados_robust=("accepted", "sum"))
    print(summ.to_string())
    n, npos = len(df), int(df.survived.sum())
    print(f"TOTAL: {n} ejemplos, {npos} positivos ({100*npos/n:.1f}%)")

    y = df.survived.values.astype(float)
    X = df[FEATURES].values.astype(float)

    # ── Normalización por cascada ────────────────────────────────────────────
    # La pregunta operativa es relativa al frame: "¿cuáles de ESTOS candidatos
    # sobreviven?". Las distribuciones absolutas de las señales se corren con
    # la energía del PKA (Simpson), así que para el pooling usamos el rango
    # percentil de cada señal DENTRO de su cascada. Es legítimo en inferencia:
    # al analizar un frame siempre se ven todos sus candidatos.
    # Señales constantes dentro de una cascada (saturadas) quedan en 0.5
    # exacto: no aportan información y su rank promedio dependería del tamaño
    # de la cascada (artefacto que la logística explotaría como ID de cascada).
    Xr = np.zeros_like(X)
    for en in df.energy.unique():
        m = (df.energy == en).values
        for k in range(X.shape[1]):
            col = X[m, k]
            Xr[m, k] = 0.5 if np.unique(col).size <= 1 \
                else pd.Series(col).rank(pct=True).values

    # ── AUC por señal individual (pooled, rank-normalizado por cascada) ──────
    print("\nAUC por señal (pooled; crudo vs rank-por-cascada; signo orientado):")
    for k, s in enumerate(FEATURES):
        a_raw = auc(y, X[:, k]);  a_raw = max(a_raw, 1 - a_raw)
        a_rk  = auc(y, Xr[:, k]); a_rk  = max(a_rk, 1 - a_rk)
        print(f"  {s:12s} crudo {a_raw:.3f}   rank {a_rk:.3f}")
    a_robust = auc(y, df.score.values)
    print(f"  {'score robust':12s} crudo {max(a_robust, 1-a_robust):.3f}")

    # precision/recall del preset robust en el pico
    acc = df.accepted.values == 1
    tp = int((acc & (y == 1)).sum())
    print(f"\npreset robust en el pico: acepta {acc.sum()}/{n}, "
          f"de los cuales sobreviven {tp}  "
          f"(precision {tp/max(1,acc.sum()):.1%}, recall {tp/max(1,npos):.1%})")

    # ── LOCO CV ──────────────────────────────────────────────────────────────
    # Features = rank percentil dentro de la cascada (centrado en 0).
    # Solo participan cascadas "resueltas": con ambas clases y ≥20 candidatos
    # (2kev/3kev empezaron a dumpear post-recombinación: todo sobrevive, no
    # informan la decisión del pico; 1kev es demasiado chica).
    Xs = Xr - 0.5
    ok = np.zeros(len(df), bool)
    for en in df.energy.unique():
        m = (df.energy == en).values
        if m.sum() >= 20 and 0 < y[m].sum() < m.sum():
            ok |= m
    print(f"\nCascadas en el modelo: {sorted(df.energy[ok].unique())}")
    df, y, X, Xr, Xs = df[ok].reset_index(drop=True), y[ok], X[ok], Xr[ok], Xs[ok]
    energies = sorted(df.energy.unique())
    print("\nLeave-one-cascade-out:")
    print(f"{'cascada':12s} {'n':>4s} {'pos':>4s} {'AUC_logreg':>10s} "
          f"{'AUC_topo':>9s} {'AUC_robust':>10s}")
    loco = []
    for en in energies:
        te = (df.energy == en).values
        tr = ~te
        if y[te].sum() == 0 or y[te].sum() == te.sum():
            print(f"{en:12s} {te.sum():4d} {int(y[te].sum()):4d}        n/a "
                  f"(una sola clase en test)")
            continue
        if tr.sum() == 0 or y[tr].sum() in (0, y[tr].size):
            print(f"{en:12s} {te.sum():4d} {int(y[te].sum()):4d}        n/a "
                  f"(train vacío o de una sola clase)")
            continue
        w, b = fit_logreg(Xs[tr], y[tr])
        a_lr = auc(y[te], Xs[te] @ w + b)
        a_tp = auc(y[te], Xr[te, FEATURES.index("topo")])
        a_rb = auc(y[te], df.score.values[te])
        loco.append((en, a_lr, a_tp, a_rb))
        print(f"{en:12s} {te.sum():4d} {int(y[te].sum()):4d} {a_lr:10.3f} "
              f"{a_tp:9.3f} {a_rb:10.3f}")
    if loco:
        m = np.array([[a, t, r] for _, a, t, r in loco])
        print(f"{'MEDIA':12s} {'':4s} {'':4s} {np.nanmean(m[:,0]):10.3f} "
              f"{np.nanmean(m[:,1]):9.3f} {np.nanmean(m[:,2]):10.3f}")

    # ── Fit final + figura ───────────────────────────────────────────────────
    w, b = fit_logreg(Xs, y)
    learned = dict(zip(FEATURES, w))
    print("\nPesos aprendidos (estandarizados, + ⇒ predice sobrevivir):")
    for s in sorted(FEATURES, key=lambda s: -abs(learned[s])):
        print(f"  {s:12s} {learned[s]:+.3f}   (robust: {ROBUST_W[s]:+.1f})")

    setup()
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.4),
                             gridspec_kw={"width_ratios": [1.3, 1, 1]})
    ax = axes[0]
    order = sorted(range(len(FEATURES)), key=lambda i: learned[FEATURES[i]])
    names = [FEATURES[i] for i in order]
    yb = np.arange(len(names))
    ax.barh(yb - 0.2, [learned[s] for s in names], 0.4,
            color=C["hybrid"], label="aprendido (supervivencia)")
    ax.barh(yb + 0.2, [ROBUST_W[s] for s in names], 0.4,
            color=C["ws"], alpha=0.7, label="robust (a mano)")
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_yticks(yb); ax.set_yticklabels(names, fontsize=9)
    ax.set_title("Pesos: aprendidos con labels reales vs preset")
    ax.legend(frameon=False, fontsize=8); ax.grid(axis="x", alpha=0.25)

    ax = axes[1]
    if loco:
        ens = [e for e, *_ in loco]
        xb = np.arange(len(ens))
        ax.bar(xb - 0.25, [a for _, a, _, _ in loco], 0.25, label="logreg",
               color=C["hybrid"])
        ax.bar(xb, [t for _, _, t, _ in loco], 0.25, label="topo sola",
               color=C["survive"])
        ax.bar(xb + 0.25, [r for _, _, _, r in loco], 0.25,
               label="score robust", color=C["ws"], alpha=0.8)
        ax.axhline(0.5, color="0.5", ls="--", lw=0.8)
        ax.set_xticks(xb); ax.set_xticklabels(ens, rotation=45, fontsize=8,
                                              ha="right")
        ax.set_ylabel("AUC (supervivencia)"); ax.set_ylim(0, 1)
        ax.set_title("LOCO: AUC por cascada held-out")
        ax.legend(frameon=False, fontsize=8)

    ax = axes[2]
    s_all = Xs @ w + b
    bins = np.linspace(s_all.min(), s_all.max(), 35)
    ax.hist(s_all[y == 0], bins=bins, alpha=0.6, label="transitorias",
            color=C["transient"], density=True)
    ax.hist(s_all[y == 1], bins=bins, alpha=0.7, label="sobreviven",
            color=C["survive"], density=True)
    ax.set_xlabel("score logístico"); ax.set_title(
        f"Separación en el pico (in-sample AUC {auc(y, s_all):.3f})")
    ax.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    out = save(fig, "ml_survival_weights")
    print(f"\nFigura: {os.path.relpath(out, ROOT)}")
    print(f"Dataset: {os.path.relpath(out_csv, ROOT)}")


if __name__ == "__main__":
    main()
