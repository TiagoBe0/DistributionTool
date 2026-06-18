#!/bin/bash
# Local sweep: run distool (hybrid robust) over every ballistic frame of one PKA
# energy in data/<DIR>/ and write a WS-vs-hybrid reconciliation CSV + per-step
# hybrid breakdowns. Mirrors run_energy_local.sh but uses the repo-local data/
# layout (ref = <PREFIX>.0) and keeps EVERY frame (fine 1000-step sampling).
#   Args: LABEL  DATA_DIR  DUMP_PREFIX
set -u
LABEL=$1; DIR=$2; PREFIX=$3

ROOT=/home/santi/DistributionTool
DISTOOL="$ROOT/build/distool"
OUTDIR="$ROOT/results/$LABEL"
BDDIR="$OUTDIR/breakdowns"
mkdir -p "$OUTDIR" "$BDDIR"

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$(nproc)}

REF="$ROOT/$DIR/${PREFIX}.0"
CSV="$OUTDIR/${LABEL}_hybrid.csv"
if [ ! -f "$REF" ]; then echo "MISSING REF: $REF"; exit 1; fi

SCRATCH=$(mktemp -d "/tmp/distool_${LABEL}_XXXXXX")
trap 'rm -rf "$SCRATCH"' EXIT
ln -sf "$REF" "$SCRATCH/ref.dump"

echo "step,ws_vac,ws_int,hybrid_vac,hybrid_filtered,hybrid_only,grid_clusters,grid_vol_A3" > "$CSV"

ALL=$(ls "$ROOT/$DIR/${PREFIX}."* 2>/dev/null | sed -E 's/.*\.([0-9]+)$/\1/' | grep -E '^[0-9]+$' | sort -n | uniq)

for s in $ALL; do
  [ "$s" = "0" ] && continue
  DMG="$ROOT/$DIR/${PREFIX}.$s"
  [ -f "$DMG" ] || continue
  ln -sf "$DMG" "$SCRATCH/dmg.dump"
  LOG="$SCRATCH/log.txt"
  "$DISTOOL" --n-max 4 --l-max 4 --hybrid --hybrid-preset robust --grid-spacing 1.0 \
       --output "$SCRATCH/out.csv" "$SCRATCH/ref.dump" "$SCRATCH/dmg.dump" > "$LOG" 2>&1 || true

  wsv=$(grep -oP 'topological \(WS\)\s*:\s*\K[0-9]+' "$LOG" | head -1)
  wsi=$(grep -oP 'vacancies / \K[0-9]+' "$LOG" | head -1)
  hyl=$(grep -oP 'topological \(hybrid\)\s*:\s*\K[0-9]+' "$LOG" | head -1)
  filt=$(grep -oP 'filters \K[0-9]+' "$LOG" | head -1)
  hyo=$(grep -oP '\+\K[0-9]+(?= non-WS)' "$LOG" | head -1)
  gc=$(grep -oP 'open void \(grid\)\s*:\s*\K[0-9]+' "$LOG" | head -1)
  gv=$(grep -oP 'clusters,\s+\K[0-9.]+' "$LOG" | head -1)
  echo "$s,${wsv:-NA},${wsi:-NA},${hyl:-NA},${filt:-0},${hyo:-0},${gc:-NA},${gv:-NA}" >> "$CSV"
  echo "  [$LABEL] step=$s  WS=${wsv:-NA}  hybrid=${hyl:-NA}  (filtered ${filt:-0})  grid=${gc:-NA}/${gv:-NA}A3"
  if [ -f "$SCRATCH/hybrid_vacancies_out.csv" ]; then
    cp "$SCRATCH/hybrid_vacancies_out.csv" "$BDDIR/bd_${s}.csv"
  fi
  rm -f "$SCRATCH"/out* "$SCRATCH"/analyzed_* "$SCRATCH"/vacancies_* \
        "$SCRATCH"/ws_* "$SCRATCH"/hybrid_* "$SCRATCH"/interstitials_* "$SCRATCH"/hist_*
done
echo "DONE $LABEL -> $CSV"
