#!/bin/bash
# Variant of run_energy_local.sh for PKA energies that LACK a step-0 reference
# dump (2/3 keV). They share the same equilibrated box as every other cascade
# run, so we reuse the 1 keV step-0 dump as the pristine reference (validated:
# WS=5/hybrid=5 on 2keV@12000, matching the live WS ground truth).
#
# Runs on the LOCAL workstation, reading the sampled ballistic frames over the
# sshfs mount; distool itself does not link libtorch so it executes locally.
# All heavy side-outputs go to a LOCAL /tmp scratch; only the per-energy CSV and
# the tiny per-frame breakdowns are written back to the mounted results dir.
#   Args: LABEL  ENERGY_DIR  DUMP_PREFIX  REF_PATH
set -u
LABEL=$1; DIR=$2; PREFIX=$3; REF=$4

MNT="/run/user/1000/gvfs/sftp:host=192.168.220.25/home/santi/pka_acero"
DISTOOL="/run/user/1000/gvfs/sftp:host=192.168.220.25/home/santi/DistributionTool/build/distool"
OUTDIR="$MNT/results/distool_sweep"
BDDIR="$OUTDIR/breakdowns_${LABEL}"
mkdir -p "$OUTDIR" "$BDDIR"

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$(nproc)}

CSV="$OUTDIR/${LABEL}_hybrid.csv"
if [ ! -f "$REF" ]; then echo "MISSING REF: $REF"; exit 1; fi

SCRATCH=$(mktemp -d "/tmp/distool_${LABEL}_XXXXXX")
trap 'rm -rf "$SCRATCH"' EXIT
cp "$REF" "$SCRATCH/ref.dump"   # copied once; read every frame

echo "step,ws_vac,ws_int,hybrid_vac,hybrid_filtered,hybrid_only,grid_clusters,grid_vol_A3" > "$CSV"

ALL=$(ls "$MNT/$DIR/${PREFIX}."* 2>/dev/null | sed -E 's/.*\.([0-9]+)$/\1/' | grep -E '^[0-9]+$' | sort -n | uniq)

for s in $ALL; do
  [ "$s" = "0" ] && continue
  # sampling: keep every frame up to 20000, then every 5000 (matches sweep)
  if [ "$s" -gt 20000 ] && [ $((s % 5000)) -ne 0 ]; then continue; fi
  SRC="$MNT/$DIR/${PREFIX}.$s"
  [ -f "$SRC" ] || continue
  cp "$SRC" "$SCRATCH/dmg.dump"
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
  echo "  [$LABEL] step=$s  WS=${wsv:-NA}  hybrid=${hyl:-NA}  (filtered ${filt:-0})"
  if [ -f "$SCRATCH/hybrid_vacancies_out.csv" ]; then
    cp "$SCRATCH/hybrid_vacancies_out.csv" "$BDDIR/bd_${s}.csv"
  fi
  rm -f "$SCRATCH"/out* "$SCRATCH"/analyzed_* "$SCRATCH"/vacancies_* \
        "$SCRATCH"/ws_* "$SCRATCH"/hybrid_* "$SCRATCH"/interstitials_* "$SCRATCH"/hist_*
done
echo "DONE $LABEL -> $CSV"
