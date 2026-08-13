#!/bin/bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

# ============================================================
# GEWUM Structure Refinement Script
# Second-stage optimization: relax cell parameters
# Processes all CIF files in current directory
# ============================================================

start_time=$(date +%s)

# Configurable parameters
MODE=2             # Always use mode 2 for cell relaxation
FMAX=${1:-0.01}    # Force convergence threshold (eV/AA)
MAX_STEPS=${2:-300} # Maximum optimization steps

echo "============================================================"
echo "GEWUM Structure Refinement (Cell Relaxation)"
echo "============================================================"
echo "  MODE: $MODE (atoms + cell)"
echo "  FMAX: $FMAX eV/AA"
echo "  MAX_STEPS: $MAX_STEPS"
echo "============================================================"

{{ENV_SETUP}}

top_dir=$(pwd)
echo "Working directory: $top_dir"

run_refinement() {
    CIF_FILE=$1
    REFINED_DIR=$2
    OPT_MODE=$3
    OPT_FMAX=$4
    OPT_MAX_STEPS=$5

    echo "Refining: $CIF_FILE"
    
    if [ -z "$CIF_FILE" ] || [ -z "$REFINED_DIR" ]; then
        echo "Error: CIF_FILE or REFINED_DIR is empty. Skipping."
        return 1
    fi
    
    mkdir -p "$REFINED_DIR"
    
    python -m gewum.src.common.relaxation.umlip_relax "$CIF_FILE" "$REFINED_DIR" --mode "$OPT_MODE" --fmax "$OPT_FMAX" --max-steps "$OPT_MAX_STEPS" &>> refine.out
}

export -f run_refinement

tasks_file="refine_tasks.txt"
> "$tasks_file"

for cif_file in $(ls "$top_dir" | grep '\.cif$'); do
    cif_file_path="$top_dir/$cif_file"
    refined_dir="$top_dir/refined"
    
    if [ -f "$cif_file_path" ]; then
        echo "$cif_file_path $refined_dir" >> "$tasks_file"
    fi
done

task_count=$(wc -l < "$tasks_file")
echo "Found $task_count CIF files to refine"

if [ ! -s "$tasks_file" ]; then
    echo "No CIF files found in current directory."
    exit 1
fi

TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-64}
CORES_PER_TASK=1
NUM_TASKS=$((TOTAL_CPUS / CORES_PER_TASK))

export OMP_NUM_THREADS=$CORES_PER_TASK

parallel -j $NUM_TASKS  -a "$tasks_file" --colsep ' ' run_refinement {1} {2} $MODE $FMAX $MAX_STEPS

grep -c "Optimized" refine.out >> TOT_refined 2>/dev/null

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo "============================================================"
echo "Refinement completed!"
echo "Refined structures saved to: refined/"
echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"
echo "============================================================" | tee Time_refine.log
