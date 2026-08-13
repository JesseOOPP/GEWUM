#!/bin/bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

echo "=============================================="
echo "GEWUM MD Post-processing"
echo "=============================================="
echo "Start time: $(date)"

{{ENV_SETUP}}

python -m gewum.src.MDworkflows.md_post --output-dir ./md_plots

echo ""
echo "=============================================="
echo "MD Post-processing Complete"
echo "=============================================="
echo "End time: $(date)"
echo "Output files:"
echo "  - md_plots/*.png (Energy-Temperature plots)"
echo "  - md_plots/md_summary.csv (Summary statistics)"
