#!/bin/bash
#SBATCH --job-name=distool_sweep
#SBATCH --account=pad_70
#SBATCH --partition=cpunode
#SBATCH --array=0-5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

LABELS=(1kev    4kev                       5kev                             6kev                     7kev                     8kev_real)
DIRS=(1kev      4kev                       5kev                             6kev                     7kev                     8kev_real)
PREFIXES=("dump.ballistic_1kev.FeCrNi_1kev" "dump.ballistic.FeCrNi_4keV" "dump.ballistic_5kev.FeCrNi_5kev_v2" "dump.ballistic.FeCrNi_6keV" "dump.ballistic.FeCrNi_7keV" "dump.ballistic.FeCrNi_8keV")

i=$SLURM_ARRAY_TASK_ID
echo "Task $i -> ${LABELS[$i]} (dir ${DIRS[$i]}, prefix ${PREFIXES[$i]}) on $(hostname)"
bash "$HOME/distool_runs/run_energy.sh" "${LABELS[$i]}" "${DIRS[$i]}" "${PREFIXES[$i]}"
