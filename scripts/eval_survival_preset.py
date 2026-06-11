#!/usr/bin/env python3
"""
eval_survival_preset.py — Evalúa el preset nativo 'survival' de distool contra
las etiquetas de supervivencia del tracking temporal, cascada por cascada.

Requiere, por cada cascada en results/pka_tracking/<energía>/:
  - tracks.csv / tracks_points.csv         (de run_pka_tracking.py)
  - hybrid_vacancies_survival_peak.csv     (distool --hybrid-preset survival
                                            sobre el frame pico, p.ej.:
    build/distool --no-dump --ws --hybrid --hybrid-preset survival \
        --grid-spacing 2.0 --output results/pka_tracking/<en>/survival_peak.csv \
        <ref.0> <frame_pico>)

Reporta por cascada: AUC del score nativo, AUC de la señal SOAP (si la corrida
la computó), AUC de centralidad, candidatos aceptados y recall de
sobrevivientes. Nota: los pesos del preset se derivaron de ESTAS cascadas
(redondeados a mano del ajuste pooled), así que esto no es una evaluación
totalmente held-out — para eso hacen falta cascadas nuevas.

Uso:
    python3 scripts/eval_survival_preset.py [--tracking results/pka_tracking]
"""
import argparse
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))

BOX_L = 243.74
PEAKS = {"4kev": 5000, "5kev": 5000, "6kev": 5000, "6kev_it1": 5000,
         "6kev_it2": 5000, "7kev": 10000, "8kev_real": 10000}


def read_hv(path):
    with open(path) as f:
        names = f.readline().lstrip("#").split()
    return pd.read_csv(path, sep=r"\s+", skiprows=1, names=names) \
             .rename(columns={"consensus_score": "score"})


def auc(y, s):
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    return float((pos[:, None] > neg[None, :]).mean()
                 + 0.5 * (pos[:, None] == neg[None, :]).mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tracking", default=os.path.join(ROOT, "results/pka_tracking"))
    args = ap.parse_args()

    print(f"{'cascada':10s} {'n':>4s} {'pos':>3s} {'AUC_score':>9s} {'AUC_soap':>8s} "
          f"{'AUC_centr':>9s} {'acepta':>6s} {'recall':>7s}")
    aucs, tot_acc, tot_tp, tot_pos, tot_n = [], 0, 0, 0, 0
    for en, peak in PEAKS.items():
        endir = os.path.join(args.tracking, en)
        f_hv = os.path.join(endir, "hybrid_vacancies_survival_peak.csv")
        if not os.path.exists(f_hv):
            print(f"{en:10s} — falta {os.path.basename(f_hv)}")
            continue
        tracks = pd.read_csv(os.path.join(endir, "tracks.csv")).set_index("track_id")
        pts = pd.read_csv(os.path.join(endir, "tracks_points.csv"))
        P = pts[pts.step == peak]
        hv = read_hv(f_hv)
        hv = hv[hv.ref_site_idx >= 0].reset_index(drop=True)
        hp = hv[["x", "y", "z"]].values

        y, idx, unmatched = [], [], 0
        for _, p in P.iterrows():
            d = hp - [p.x, p.y, p.z]
            d -= BOX_L * np.round(d / BOX_L)
            r2 = (d ** 2).sum(1)
            i = int(r2.argmin())
            if r2[i] > 0.01:
                unmatched += 1
                continue
            y.append(int(tracks.loc[p.track_id].survived))
            idx.append(i)
        if unmatched > len(P) * 0.05:
            print(f"{en:10s} [!] {unmatched}/{len(P)} sin match — ¿referencia "
                  f"distinta entre el tracking y la corrida survival_peak?")
            continue
        y = np.array(y, float)
        h = hv.iloc[idx].reset_index(drop=True)
        a_sc = auc(y, h.score.values)
        a_so = auc(y, h.soap.values)
        a_ce = auc(y, h.centrality.values)
        acc = h.accepted.values == 1
        tp = int((acc & (y == 1)).sum())
        aucs.append((a_sc, a_so, a_ce))
        tot_acc += int(acc.sum()); tot_tp += tp
        tot_pos += int(y.sum());  tot_n += len(y)
        print(f"{en:10s} {len(y):4d} {int(y.sum()):3d} {a_sc:9.3f} {a_so:8.3f} "
              f"{a_ce:9.3f} {int(acc.sum()):6d} {tp}/{int(y.sum()):>2d}")

    if aucs:
        m = np.array(aucs)
        print(f"{'MEDIA':10s} {'':4s} {'':3s} {np.nanmean(m[:,0]):9.3f} "
              f"{np.nanmean(m[:,1]):8.3f} {np.nanmean(m[:,2]):9.3f}")
        print(f"\nTotal: acepta {tot_acc}/{tot_n} candidatos de pico; "
              f"recall global {tot_tp}/{tot_pos} sobrevivientes "
              f"({100*tot_tp/max(1,tot_pos):.0f}%)")


if __name__ == "__main__":
    main()
