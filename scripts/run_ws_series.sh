#!/bin/bash
# run_ws_series.sh — Corre distool en modo topológico (--no-soap) sobre una
# serie temporal de frames contra una misma referencia, generando los
# ws_sites_*.csv por frame que consume scripts/track_defects.py.
#
# Uso:
#   scripts/run_ws_series.sh <ref.dump> <outdir> <frame1> [frame2 ...]
#
# Ejemplo (cascada 5 keV):
#   scripts/run_ws_series.sh \
#       data/5kev_20_05/dump.ballistic.FeCrNi_5kev_s12345.0 \
#       results/5kev_20_05/tracking \
#       data/5kev_20_05/dump.ballistic.FeCrNi_5kev_s12345.{5000,10000,...} \
#       data/5kev_20_05/dump.relax.FeCrNi_5kev_s12345.100000
#
# El timestep de cada frame se toma del sufijo numérico final del nombre de
# archivo (convención LAMMPS dump.<algo>.<timestep>).
set -euo pipefail

DISTOOL="${DISTOOL:-./build/distool}"
GRID_SPACING="${GRID_SPACING:-1.0}"
PRESET="${PRESET:-robust}"

[[ $# -ge 3 ]] || { echo "Uso: $0 <ref.dump> <outdir> <frame...>" >&2; exit 1; }
REF="$1"; OUTDIR="$2"; shift 2
[[ -f "$REF" ]] || { echo "Referencia no encontrada: $REF" >&2; exit 1; }
[[ -x "$DISTOOL" ]] || { echo "distool no encontrado: $DISTOOL (set DISTOOL=...)" >&2; exit 1; }

mkdir -p "$OUTDIR"

for FRAME in "$@"; do
    [[ -f "$FRAME" ]] || { echo "[!] Frame no encontrado, salto: $FRAME" >&2; continue; }
    STEP="${FRAME##*.}"
    [[ "$STEP" =~ ^[0-9]+$ ]] || { echo "[!] Sin timestep numérico, salto: $FRAME" >&2; continue; }
    TAG=$(printf "t%06d" "$STEP")
    echo "── $TAG ← $FRAME"
    "$DISTOOL" --no-soap --no-dump --ws --hybrid --hybrid-preset "$PRESET" \
        --grid-spacing "$GRID_SPACING" \
        --output "$OUTDIR/$TAG.csv" \
        "$REF" "$FRAME" > "$OUTDIR/$TAG.log" 2>&1 \
        || { echo "[✗] distool falló en $TAG (ver $OUTDIR/$TAG.log)" >&2; exit 1; }
    grep -E "topological \(WS\)|topological \(hybrid\)" "$OUTDIR/$TAG.log" | sed 's/^/    /'
done

echo "Listo: $(ls "$OUTDIR"/ws_sites_t*.csv 2>/dev/null | wc -l) frames en $OUTDIR"
