#!/bin/bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

start_time=$(date +%s)

# Configurable parameters
BOND_THRESHOLD=${1:-1.0}    # Minimum bond length threshold (Angstrom)

{{ENV_SETUP}}


echo "Starting bond length check with threshold: $BOND_THRESHOLD A..."
python -m gewum.src.common.postprocess.bond_check --threshold $BOND_THRESHOLD
echo "Bond length check completed."

echo "Starting energy check..."
python -m gewum.src.PTworkflows.energy_post_pt
echo "Energy check completed."


end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Total runtime: $hours hours, $minutes minutes, $seconds seconds" > Time_post_processing.log
