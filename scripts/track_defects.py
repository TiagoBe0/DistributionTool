#!/usr/bin/env python3
"""
track_defects.py — Tracking temporal de vacancias Wigner-Seitz a lo largo de
una cascada.

Toma la serie de ws_sites_t*.csv que produce scripts/run_ws_series.sh (todos
los frames clasificados contra la MISMA referencia prístina, así que los
ref_id de los sitios son comparables entre frames) y enlaza las vacancias
(occupancy == 0) entre frames consecutivos en trayectorias ("tracks"):

  1. PERSIST : el mismo sitio (ref_id) sigue vacante en el frame siguiente.
  2. HOP     : el sitio se llenó pero hay una vacancia nueva a < max_hop Å
               (la vacancia migró a un sitio vecino). Matching greedy por
               distancia mínima con imagen mínima (PBC).
  3. DEATH   : sin match → la vacancia se recombinó (aniquilada).
  4. BIRTH   : vacancia nueva sin antecesor → defecto recién creado.

Cada track recibe la etiqueta `survived` si está vivo en el ÚLTIMO frame de
la serie (que debe ser el frame relajado final). Esas etiquetas por-defecto
son la verdad de terreno para entrenar/evaluar el score del detector híbrido
en el pico balístico.

Limitaciones (V1, documentadas a propósito):
  * Tracking a nivel de sitio: una divacancia son 2 tracks adyacentes.
  * En el núcleo balístico denso (t≈5000) el matching greedy puede enlazar
    mal identidades individuales; las estadísticas agregadas (births/deaths/
    supervivencia) son robustas, las identidades individuales tempranas no.
  * Gaps: un track que muere en t y "renace" en t+2Δt cuenta como dos tracks.

Uso:
    python3 scripts/track_defects.py \
        --sites-dir results/5kev_20_05/tracking \
        --box-from data/5kev_20_05/dump.ballistic.FeCrNi_5kev_s12345.0 \
        --max-hop 4.0 \
        --out-prefix results/5kev_20_05/tracking/tracks \
        --fig figures/defect_tracking_5kev.png
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WS_SITES_COLS = ["ref_id", "x", "y", "z", "occupancy", "min_dist", "site_type"]


# ─────────────────────────────────────────────────────────────────────────────
# IO
# ─────────────────────────────────────────────────────────────────────────────

def read_box(dump_path):
    """Lee Lx, Ly, Lz y flags periódicos del header de un dump LAMMPS."""
    bounds, periodic = [], [True, True, True]
    with open(dump_path) as f:
        lines = [next(f) for _ in range(9)]
    for i, line in enumerate(lines):
        if line.startswith("ITEM: BOX BOUNDS"):
            toks = line.split()
            if len(toks) >= 3:
                flags = toks[-3:]
                periodic = [t in ("pp", "p") for t in flags]
            for j in range(i + 1, i + 4):
                lo, hi = lines[j].split()[:2]
                bounds.append(float(hi) - float(lo))
            break
    if len(bounds) != 3:
        raise ValueError(f"No pude leer BOX BOUNDS de {dump_path}")
    return np.array(bounds), periodic


def load_sites_series(sites_dir):
    """Devuelve [(step, DataFrame de sitios vacantes)], ordenado por step."""
    files = sorted(glob.glob(os.path.join(sites_dir, "ws_sites_t*.csv")))
    if not files:
        sys.exit(f"No hay ws_sites_t*.csv en {sites_dir}")
    series = []
    for path in files:
        m = re.search(r"t(\d+)", os.path.basename(path))
        if not m:
            continue
        step = int(m.group(1))
        df = pd.read_csv(path, sep=r"\s+", skiprows=1, names=WS_SITES_COLS)
        vac = df[df.occupancy == 0][["ref_id", "x", "y", "z"]].reset_index(drop=True)
        series.append((step, vac))
    series.sort(key=lambda t: t[0])
    return series


# ─────────────────────────────────────────────────────────────────────────────
# Matching
# ─────────────────────────────────────────────────────────────────────────────

def pbc_dist_matrix(pos_a, pos_b, box, periodic):
    """Matriz de distancias |a_i - b_j| con imagen mínima."""
    d = pos_a[:, None, :] - pos_b[None, :, :]          # (na, nb, 3)
    for k in range(3):
        if periodic[k] and box[k] > 0:
            d[:, :, k] -= box[k] * np.round(d[:, :, k] / box[k])
    return np.sqrt((d ** 2).sum(axis=2))


def greedy_match(pos_old, pos_new, box, periodic, max_hop):
    """Matching bipartito greedy por distancia mínima creciente.

    Devuelve lista de (i_old, j_new, dist) con dist <= max_hop.
    Greedy es suficiente aquí: los pares reales (hops de ~1 NN) dominan
    claramente sobre los matches espurios a > max_hop.
    """
    if len(pos_old) == 0 or len(pos_new) == 0:
        return []
    dm = pbc_dist_matrix(pos_old, pos_new, box, periodic)
    pairs = []
    cand = np.argwhere(dm <= max_hop)
    order = np.argsort(dm[cand[:, 0], cand[:, 1]])
    used_old, used_new = set(), set()
    for k in order:
        i, j = cand[k]
        if i in used_old or j in used_new:
            continue
        used_old.add(i)
        used_new.add(j)
        pairs.append((int(i), int(j), float(dm[i, j])))
    return pairs


# ─────────────────────────────────────────────────────────────────────────────
# Tracking
# ─────────────────────────────────────────────────────────────────────────────

def track(series, box, periodic, max_hop, max_hop_final=None):
    """Enlaza vacancias entre frames. Devuelve (tracks, points, per_frame).

    tracks: dict track_id -> dict(birth_step, death_step, hops, hop_dist)
    points: lista de dicts (track_id, step, ref_id, x, y, z, how)
    per_frame: lista de dicts con estadísticas por transición

    max_hop_final: radio de matching para la ÚLTIMA transición (típicamente el
    salto grande balístico→relajado, donde la población es rala y los defectos
    pueden haber migrado más que entre frames balísticos consecutivos).
    """
    if max_hop_final is None:
        max_hop_final = max_hop
    last_step = series[-1][0]
    tracks, points, per_frame = {}, [], []
    next_id = 0

    step0, vac0 = series[0]
    active = {}  # ref_id -> track_id
    for _, row in vac0.iterrows():
        tracks[next_id] = dict(birth_step=step0, death_step=None,
                               hops=0, hop_dist=0.0)
        points.append(dict(track_id=next_id, step=step0, ref_id=int(row.ref_id),
                           x=row.x, y=row.y, z=row.z, how="birth"))
        active[int(row.ref_id)] = next_id
        next_id += 1
    per_frame.append(dict(step=step0, n_vac=len(vac0), persisted=0,
                          hopped=0, births=len(vac0), deaths=0))

    prev_vac = vac0
    for step, vac in series[1:]:
        old_ids = set(active.keys())
        new_ids = set(vac.ref_id.astype(int))

        persisted = old_ids & new_ids
        gone      = old_ids - new_ids   # candidatos a hop o muerte
        fresh     = new_ids - old_ids   # candidatos a hop o nacimiento

        new_active = {}
        vac_by_id = {int(r.ref_id): r for _, r in vac.iterrows()}
        prev_by_id = {int(r.ref_id): r for _, r in prev_vac.iterrows()}

        for rid in persisted:
            tid = active[rid]
            r = vac_by_id[rid]
            points.append(dict(track_id=tid, step=step, ref_id=rid,
                               x=r.x, y=r.y, z=r.z, how="persist"))
            new_active[rid] = tid

        # Matching de hops entre los que desaparecieron y los nuevos
        gone_l, fresh_l = sorted(gone), sorted(fresh)
        pos_old = np.array([[prev_by_id[r].x, prev_by_id[r].y, prev_by_id[r].z]
                            for r in gone_l]).reshape(-1, 3)
        pos_new = np.array([[vac_by_id[r].x, vac_by_id[r].y, vac_by_id[r].z]
                            for r in fresh_l]).reshape(-1, 3)
        hop_r = max_hop_final if step == last_step else max_hop
        pairs = greedy_match(pos_old, pos_new, box, periodic, hop_r)

        matched_old, matched_new = set(), set()
        for i, j, dist in pairs:
            rid_old, rid_new = gone_l[i], fresh_l[j]
            tid = active[rid_old]
            tracks[tid]["hops"] += 1
            tracks[tid]["hop_dist"] += dist
            r = vac_by_id[rid_new]
            points.append(dict(track_id=tid, step=step, ref_id=rid_new,
                               x=r.x, y=r.y, z=r.z, how="hop"))
            new_active[rid_new] = tid
            matched_old.add(rid_old)
            matched_new.add(rid_new)

        deaths = 0
        for rid in gone:
            if rid in matched_old:
                continue
            tracks[active[rid]]["death_step"] = step
            deaths += 1

        births = 0
        for rid in fresh:
            if rid in matched_new:
                continue
            r = vac_by_id[rid]
            tracks[next_id] = dict(birth_step=step, death_step=None,
                                   hops=0, hop_dist=0.0)
            points.append(dict(track_id=next_id, step=step, ref_id=rid,
                               x=r.x, y=r.y, z=r.z, how="birth"))
            new_active[rid] = next_id
            next_id += 1
            births += 1

        per_frame.append(dict(step=step, n_vac=len(vac),
                              persisted=len(persisted), hopped=len(pairs),
                              births=births, deaths=deaths))
        active = new_active
        prev_vac = vac

    survivors = set(active.values())
    return tracks, points, per_frame, survivors


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites-dir", required=True,
                    help="directorio con los ws_sites_t*.csv de la serie")
    ap.add_argument("--box-from", required=True,
                    help="dump LAMMPS del que leer el box (la referencia sirve)")
    ap.add_argument("--max-hop", type=float, default=4.0,
                    help="distancia máxima de salto entre frames [Å] (def: 4.0)")
    ap.add_argument("--max-hop-final", type=float, default=None,
                    help="radio de matching para la última transición "
                         "(salto balístico→relajado) [Å] (def: = max-hop)")
    ap.add_argument("--out-prefix", required=True,
                    help="prefijo de salida: <prefix>.csv y <prefix>_points.csv")
    ap.add_argument("--fig", default=None, help="ruta del PNG de la figura resumen")
    args = ap.parse_args()

    box, periodic = read_box(args.box_from)
    series = load_sites_series(args.sites_dir)
    steps = [s for s, _ in series]
    print(f"Frames: {steps}")
    print(f"Box: {box.round(2)} Å  periodic={periodic}  max_hop={args.max_hop} Å")

    tracks, points, per_frame, survivors = track(series, box, periodic,
                                                 args.max_hop, args.max_hop_final)

    # ── Tablas de salida ─────────────────────────────────────────────────────
    pts = pd.DataFrame(points)
    rows = []
    for tid, t in tracks.items():
        tp = pts[pts.track_id == tid]
        rows.append(dict(
            track_id=tid,
            birth_step=t["birth_step"],
            death_step=t["death_step"] if t["death_step"] is not None else -1,
            n_frames=len(tp),
            hops=t["hops"],
            hop_dist_total=round(t["hop_dist"], 3),
            last_step=int(tp.step.max()),
            last_ref_id=int(tp.sort_values("step").iloc[-1].ref_id),
            survived=int(tid in survivors),
        ))
    tr = pd.DataFrame(rows).sort_values(["survived", "birth_step"],
                                        ascending=[False, True])
    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)
    tr.to_csv(args.out_prefix + ".csv", index=False)
    pts.to_csv(args.out_prefix + "_points.csv", index=False)

    # ── Resumen en consola ───────────────────────────────────────────────────
    pf = pd.DataFrame(per_frame)
    print("\nPor frame (transición desde el frame anterior):")
    print(pf.to_string(index=False))

    n_total = len(tracks)
    n_surv = len(survivors)
    print(f"\nTracks totales: {n_total}   sobreviven al frame final: {n_surv}")
    print(f"Vida media de los que mueren: "
          f"{tr[tr.survived == 0].n_frames.mean():.2f} frames")

    surv = tr[tr.survived == 1]
    if not surv.empty:
        print("\nSobrevivientes (lineage):")
        cols = ["track_id", "birth_step", "n_frames", "hops", "hop_dist_total"]
        print(surv[cols].to_string(index=False))

    # ── Figura ───────────────────────────────────────────────────────────────
    if args.fig:
        # Vacancias vivas en cada frame, divididas en "van a sobrevivir" vs no
        alive_surv, alive_doom = [], []
        for s in steps:
            at_s = pts[pts.step == s].track_id.unique()
            ns = sum(1 for t in at_s if t in survivors)
            alive_surv.append(ns)
            alive_doom.append(len(at_s) - ns)

        fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
        ax = axes[0]
        ax.stackplot(steps, alive_surv, alive_doom,
                     labels=["sobreviven al final", "transitorias (mueren)"],
                     colors=["#2a9d8f", "#e76f51"], alpha=0.85)
        ax.set_yscale("symlog")
        ax.set_xlabel("timestep")
        ax.set_ylabel("vacancias WS vivas")
        ax.set_title("Vacancias por destino final")
        ax.legend(loc="upper right", fontsize=9)

        ax = axes[1]
        ax.plot(pf.step, pf.births, "o-", label="nacimientos", color="#264653")
        ax.plot(pf.step, pf.deaths, "s-", label="muertes (recombinación)",
                color="#e76f51")
        ax.plot(pf.step, pf.hopped, "^-", label="hops", color="#2a9d8f")
        ax.set_yscale("symlog")
        ax.set_xlabel("timestep")
        ax.set_title("Eventos por transición")
        ax.legend(fontsize=9)

        ax = axes[2]
        for tid in sorted(survivors):
            tp = pts[pts.track_id == tid].sort_values("step")
            ax.plot(tp.step, [tid] * len(tp), "-", lw=1.0, color="#bbbbbb")
            ax.scatter(tp.step, [tid] * len(tp), s=18,
                       c=["#264653" if h == "persist" else
                          "#e9c46a" if h == "hop" else "#2a9d8f"
                          for h in tp.how], zorder=3)
        ax.set_xlabel("timestep")
        ax.set_ylabel("track_id")
        ax.set_title("World-lines de los sobrevivientes\n"
                     "(verde=nace, amarillo=hop, azul=persiste)")

        fig.suptitle(f"Tracking temporal de vacancias WS — max_hop={args.max_hop} Å",
                     y=1.02)
        fig.tight_layout()
        fig.savefig(args.fig, dpi=160, bbox_inches="tight")
        print(f"\nFigura: {args.fig}")

    print(f"Tracks:  {args.out_prefix}.csv")
    print(f"Puntos:  {args.out_prefix}_points.csv")


if __name__ == "__main__":
    main()
