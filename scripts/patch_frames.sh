#!/bin/bash
# Re-run specific (label,step) frames that came out NA and patch the CSV in place.
#   Args: LABEL  ENERGY_DIR  DUMP_PREFIX  STEP [STEP ...]
set -u
LABEL=$1; DIR=$2; PREFIX=$3; shift 3
BASE=/home/santi/pka_acero
DISTOOL=/home/santi/DistributionTool/build/distool
OUTDIR=/home/santi/pka_acero/results/distool_sweep
CSV="$OUTDIR/${LABEL}_hybrid.csv"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$(nproc)}

SCRATCH=$(mktemp -d "/tmp/distpatch_${LABEL}_XXXXXX")
trap 'rm -rf "$SCRATCH"' EXIT
ln -sf "$BASE/$DIR/${PREFIX}.0" "$SCRATCH/ref.dump"

for s in "$@"; do
  DMG="$BASE/$DIR/${PREFIX}.$s"
  [ -f "$DMG" ] || { echo "  skip $s (no dump)"; continue; }
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
  row="$s,${wsv:-NA},${wsi:-NA},${hyl:-NA},${filt:-0},${hyo:-0},${gc:-NA},${gv:-NA}"
  # replace existing line starting with "$s," else append
  if grep -q "^$s," "$CSV"; then
    sed -i "s|^$s,.*|$row|" "$CSV"
  else
    echo "$row" >> "$CSV"
  fi
  echo "  [$LABEL] step=$s  WS=${wsv:-NA}  hybrid=${hyl:-NA}  (patched)"
  rm -f "$SCRATCH"/out* "$SCRATCH"/analyzed_* "$SCRATCH"/vacancies_* \
        "$SCRATCH"/ws_* "$SCRATCH"/hybrid_* "$SCRATCH"/interstitials_* "$SCRATCH"/hist_*
done
# keep CSV sorted by step (header first)
{ head -1 "$CSV"; tail -n +2 "$CSV" | sort -t, -k1 -n; } > "$CSV.tmp" && mv "$CSV.tmp" "$CSV"
echo "PATCHED $LABEL"
