#!/usr/bin/env python3
"""
run_pka_tracking.py — Orquesta el tracking temporal sobre todas las cascadas
del dataset pka_acero (1–8 keV).

Por cada energía:
  1. Descubre la serie de dumps (prefijo balístico + frames relax) y elige
     frames a Δt≈5000: todos los pasos múltiplos de 5000, más el primer paso
     disponible (2kev/3kev arrancan en 12000/16000), el último balístico y
     los frames relax.
  2. Corre scripts/run_ws_series.sh (distool --no-soap --ws --hybrid) contra
     la referencia t=0 — la propia, o la de 1kev para las cascadas que no
     guardaron el frame 0 (mismo cristal termalizado: box idéntico bit a bit).
  3. Corre scripts/track_defects.py → tracks.csv con etiquetas de
     supervivencia + figura por energía.

Es resumible: salta los frames cuyo ws_sites ya existe y las energías con
tracks.csv ya generado (usar --force para rehacer).

Uso:
    python3 scripts/run_pka_tracking.py [--base /home/santi/pka_acero]
        [--out results/pka_tracking] [--energies 1kev,2kev,...] [--force]
"""

import argparse
import glob
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (dir, prefijo balístico, prefijo relax)  — None = autodetectar
ENERGIES = ["1kev", "2kev", "3kev", "4kev", "5kev",
            "6kev", "6kev_it1", "6kev_it2", "7kev", "8kev_real"]

FALLBACK_REF = "1kev"   # energías sin frame .0 usan la referencia de 1kev


def discover(dirpath):
    """Devuelve (ref_path|None, [(step, path)] balísticos, [(step, path)] relax)."""
    ball, relax = {}, {}
    for path in glob.glob(os.path.join(dirpath, "dump.*")):
        m = re.match(r".*\.(\d+)$", path)
        if not m:
            continue
        step = int(m.group(1))
        base = os.path.basename(path)
        if base.startswith("dump.relax"):
            relax[step] = path
        elif base.startswith("dump.ballistic") or base.startswith("dump.final"):
            ball[step] = path
    ref = ball.pop(0, None)
    return ref, sorted(ball.items()), sorted(relax.items())


def select_steps(ball, relax):
    """Frames a Δt≈5000 + primero + último balístico + relax posteriores."""
    if not ball:
        return []
    steps = [s for s, _ in ball]
    chosen = {s for s in steps if s % 5000 == 0}
    chosen.add(steps[0])
    chosen.add(steps[-1])
    frames = [(s, p) for s, p in ball if s in chosen]
    last_ball = steps[-1]
    # Relax: solo pasos posteriores al último balístico usado (evita el
    # duplicado relax.100000 == ballistic.100000 de las corridas 6–8 keV).
    frames += [(s, p) for s, p in relax if s > last_ball]
    return frames


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="/home/santi/pka_acero")
    ap.add_argument("--out", default="results/pka_tracking")
    ap.add_argument("--energies", default=",".join(ENERGIES))
    ap.add_argument("--grid-spacing", default="2.0")
    ap.add_argument("--max-hop", default="4.0")
    ap.add_argument("--max-hop-final", default="10.0")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    energies = args.energies.split(",")
    fallback_ref, _, _ = discover(os.path.join(args.base, FALLBACK_REF))
    runner = os.path.join(REPO, "scripts", "run_ws_series.sh")
    tracker = os.path.join(REPO, "scripts", "track_defects.py")

    for en in energies:
        d = os.path.join(args.base, en)
        if not os.path.isdir(d):
            print(f"[{en}] no existe, salto", flush=True)
            continue
        outdir = os.path.join(args.out, en)
        tracks_csv = os.path.join(outdir, "tracks.csv")
        if os.path.exists(tracks_csv) and not args.force:
            print(f"[{en}] tracks.csv ya existe, salto", flush=True)
            continue

        ref, ball, relax = discover(d)
        if ref is None:
            ref = fallback_ref
            ref_note = f"(referencia prestada de {FALLBACK_REF})"
        else:
            ref_note = ""
        frames = select_steps(ball, relax)
        if not frames:
            print(f"[{en}] sin frames, salto", flush=True)
            continue

        # Resumibilidad a nivel frame
        pending = [p for s, p in frames
                   if not os.path.exists(os.path.join(outdir, f"ws_sites_t{s:06d}.csv"))]
        print(f"[{en}] {len(frames)} frames ({len(pending)} pendientes), "
              f"pasos {frames[0][0]}–{frames[-1][0]} {ref_note}", flush=True)

        if pending:
            env = dict(os.environ,
                       DISTOOL=os.path.join(REPO, "build", "distool"),
                       GRID_SPACING=args.grid_spacing)
            r = subprocess.run([runner, ref, outdir] + pending, env=env,
                               cwd=REPO, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True)
            tail = "\n".join(r.stdout.strip().splitlines()[-3:])
            if r.returncode != 0:
                print(f"[{en}] FALLÓ la serie:\n{tail}", flush=True)
                continue

        r = subprocess.run(
            [sys.executable, tracker,
             "--sites-dir", outdir,
             "--box-from", ref,
             "--max-hop", args.max_hop,
             "--max-hop-final", args.max_hop_final,
             "--out-prefix", os.path.join(outdir, "tracks"),
             "--fig", os.path.join(REPO, "figures", f"tracking_{en}.png")],
            cwd=REPO, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if r.returncode != 0:
            print(f"[{en}] FALLÓ el tracker:\n{r.stdout[-800:]}", flush=True)
            continue
        # Resumen compacto: tracks totales y sobrevivientes
        for line in r.stdout.splitlines():
            if line.startswith("Tracks totales") or line.startswith("Por frame"):
                print(f"[{en}] {line}", flush=True)
        print(f"[{en}] LISTO → {tracks_csv}", flush=True)

    print("Barrido completo.", flush=True)


if __name__ == "__main__":
    main()
