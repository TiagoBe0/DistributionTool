#!/usr/bin/env python3
"""
Wigner-Seitz defect tracking across a cascade simulation.

Toma el frame con el step más bajo del run como referencia (estado pre-cascada)
y corre WS estándar (simple_ws.wigner_seitz) en cada frame, escribiendo un CSV
con vacancias e intersticiales vs step.

Uso típico:
    python analyze_cascade.py --run-dir ../2kev --output-csv results/2kev.csv \
        --nproc 16
"""

import argparse
import csv
import glob
import os
import re
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wigner_seitz import LAMMPSDumpReader  # reutilizamos el reader
from simple_ws import Box, wigner_seitz


STEP_RE = re.compile(r"\.(\d+)$")


def find_ballistic_frames(run_dir):
    # Soporta dos convenciones de naming:
    #   dump.ballistic_<label>.<step>  (corridas viejas: 1kev, 2kev, 3kev, 5kev)
    #   dump.ballistic.<label>.<step>  (corridas nuevas: 4kev, 6kev, 7kev, 8kev_real)
    candidates = (glob.glob(os.path.join(run_dir, "dump.ballistic_*")) +
                  glob.glob(os.path.join(run_dir, "dump.ballistic.*")))
    frames = []
    for f in candidates:
        m = STEP_RE.search(os.path.basename(f))
        if m:
            frames.append((int(m.group(1)), f))
    frames.sort()
    return frames


_ref_pos = None
_box = None


def _init_worker(ref_file):
    global _ref_pos, _box
    ref = LAMMPSDumpReader.read(ref_file)
    _ref_pos = ref.positions
    _box = Box(ref.box.xlo, ref.box.xhi,
               ref.box.ylo, ref.box.yhi,
               ref.box.zlo, ref.box.zhi)


def _analyze_one(args):
    step, frame_file = args
    t0 = time.time()
    df = LAMMPSDumpReader.read(frame_file)
    res = wigner_seitz(_ref_pos, df.positions, _box)
    return {
        "step": step,
        "n_atoms_def": int(len(df.positions)),
        "n_vacancies": res["n_vacancies"],
        "n_interstitials": res["n_interstitials"],
        "elapsed_s": round(time.time() - t0, 2),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", required=True)
    p.add_argument("--reference", default=None,
                   help="Dump de referencia (default: frame con step más bajo)")
    p.add_argument("--output-csv", required=True)
    p.add_argument("--every", type=int, default=1, help="Submuestrear cada N frames")
    p.add_argument("--max-frames", type=int, default=None)
    p.add_argument("--nproc", type=int, default=1)
    p.add_argument("--skip-existing", action="store_true",
                   help="No re-procesar steps ya presentes en el CSV")
    args = p.parse_args()

    frames = find_ballistic_frames(args.run_dir)
    if not frames:
        sys.exit(f"No se encontraron dumps en {args.run_dir}")

    ref_file = args.reference or frames[0][1]

    frames = frames[::args.every]
    if args.max_frames:
        frames = frames[:args.max_frames]

    existing = set()
    if args.skip_existing and os.path.exists(args.output_csv):
        with open(args.output_csv) as f:
            for row in csv.DictReader(f):
                existing.add(int(row["step"]))
        frames = [(s, f) for (s, f) in frames if s not in existing]

    print(f"Reference: {ref_file}")
    print(f"Frames a analizar: {len(frames)} (every={args.every}, max={args.max_frames})")
    print(f"Workers: {args.nproc}")

    if not frames:
        print("Nada para hacer."); return

    os.makedirs(os.path.dirname(os.path.abspath(args.output_csv)) or ".", exist_ok=True)

    t_total = time.time()
    if args.nproc > 1:
        with Pool(args.nproc, initializer=_init_worker, initargs=(ref_file,)) as pool:
            results = []
            for i, r in enumerate(pool.imap_unordered(_analyze_one, frames), 1):
                results.append(r)
                if i % 10 == 0 or i == len(frames):
                    print(f"  [{i}/{len(frames)}] step={r['step']} "
                          f"vac={r['n_vacancies']} int={r['n_interstitials']} "
                          f"({r['elapsed_s']}s)")
    else:
        _init_worker(ref_file)
        results = []
        for i, fr in enumerate(frames, 1):
            r = _analyze_one(fr)
            results.append(r)
            print(f"  [{i}/{len(frames)}] step={r['step']} "
                  f"vac={r['n_vacancies']} int={r['n_interstitials']} "
                  f"({r['elapsed_s']}s)")

    results.sort(key=lambda r: r["step"])

    fields = ["step", "n_atoms_def", "n_vacancies", "n_interstitials", "elapsed_s"]
    mode = "a" if (args.skip_existing and existing) else "w"
    write_header = mode == "w" or not os.path.exists(args.output_csv)
    with open(args.output_csv, mode, newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        if write_header:
            w.writeheader()
        w.writerows(results)

    elapsed = time.time() - t_total
    print(f"\nWrote {args.output_csv} ({len(results)} new rows) in {elapsed:.1f}s "
          f"({elapsed/max(1,len(results)):.2f}s/frame avg)")


if __name__ == "__main__":
    main()
