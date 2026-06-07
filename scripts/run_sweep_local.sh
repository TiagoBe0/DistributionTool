#!/bin/bash
# Drive run_energy_local.sh over every PKA energy that has a pristine step-0
# reference (1,4,5,6,7,8 keV; 2/3 keV lack a step-0 dump). Sequential — each
# distool run already uses all cores via OpenMP.
set -u
HERE=$(dirname "$(readlink -f "$0")")
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-$(nproc)}

LABELS=(1kev    4kev    5kev    6kev    7kev    8kev_real)
DIRS=(  1kev    4kev    5kev    6kev    7kev    8kev_real)
PREFIXES=("dump.ballistic_1kev.FeCrNi_1kev" \
          "dump.ballistic.FeCrNi_4keV" \
          "dump.ballistic_5kev.FeCrNi_5kev_v2" \
          "dump.ballistic.FeCrNi_6keV" \
          "dump.ballistic.FeCrNi_7keV" \
          "dump.ballistic.FeCrNi_8keV")

t0=$(date +%s)
for i in "${!LABELS[@]}"; do
  echo "==== [$((i+1))/${#LABELS[@]}] ${LABELS[$i]} ===="
  bash "$HERE/run_energy_local.sh" "${LABELS[$i]}" "${DIRS[$i]}" "${PREFIXES[$i]}"
done
echo "ALL DONE in $(( ($(date +%s)-t0)/60 )) min -> /home/santi/pka_acero/results/distool_sweep/"
