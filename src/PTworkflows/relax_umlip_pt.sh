#!/bin/bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

start_time=$(date +%s)

# Configurable parameters (can be overridden via command line)
MODE=${1:-1}           # 1=atoms only, 2=atoms+cell (default: 1 for PT)
FMAX=${2:-0.05}        # Force convergence threshold (eV/AA)
MAX_STEPS=${3:-200}    # Maximum optimization steps

echo "Optimization settings:"
echo "  MODE: $MODE"
echo "  FMAX: $FMAX eV/AA"
echo "  MAX_STEPS: $MAX_STEPS"

{{ENV_SETUP}}

top_dir=$(pwd)
echo "Top directory: $top_dir"

if [ ! -d "$top_dir" ]; then
    echo "Top directory does not exist: $top_dir"
    exit 1
fi

run_optimization() {
    CIF_FILE=$1
    RELAXED_DIR=$2
    OPT_MODE=$3
    OPT_FMAX=$4
    OPT_MAX_STEPS=$5
    
    echo "CIF_FILE: $CIF_FILE"
    echo "RELAXED_DIR: $RELAXED_DIR"
    
    if [ -z "$CIF_FILE" ] || [ -z "$RELAXED_DIR" ]; then
        echo "Error: CIF_FILE or RELAXED_DIR is empty. Skipping this task."
        return 1
    fi
    
    mkdir -p "$RELAXED_DIR"
    
    python -m gewum.src.common.relaxation.umlip_relax "$CIF_FILE" "$RELAXED_DIR" --mode "$OPT_MODE" --fmax "$OPT_FMAX" --max-steps "$OPT_MAX_STEPS" &>> cif_relax.out
    
    last_bfgs_line=$(grep '^BFGS:' cif_relax.out | tail -n 1)
    
    if [ -n "$last_bfgs_line" ]; then
        echo "$last_bfgs_line" >> cif_relax.out
    fi
}

export -f run_optimization

tasks_file="tasks.txt"
> "$tasks_file" 

for cif_file in $(ls "$top_dir" | grep '\.cif$'); do
    cif_file_path="$top_dir/$cif_file"
    relaxed_dir="$top_dir/relaxed"
    
    echo "Adding task: CIF_FILE=$cif_file_path, RELAXED_DIR=$relaxed_dir"
    
    if [ -n "$cif_file_path" ] && [ -n "$relaxed_dir" ]; then
        echo "$cif_file_path $relaxed_dir" >> "$tasks_file"
    else
        echo "Error: CIF_FILE or RELAXED_DIR is empty. Skipping this task."
    fi
done

if [ ! -s "$tasks_file" ]; then
    echo "No tasks found."
    exit 1
fi

TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-64}
CORES_PER_TASK=1
NUM_TASKS=$((TOTAL_CPUS / CORES_PER_TASK))

export OMP_NUM_THREADS=$CORES_PER_TASK  

parallel -j $NUM_TASKS  -a "$tasks_file" --colsep ' ' run_optimization {1} {2} $MODE $FMAX $MAX_STEPS

grep -c "Optimized" cif_relax.out  >> TOT_relaxed

rm  cif_relax.out 
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Total runtime: $hours hours, $minutes minutes, $seconds seconds" > Time_cif_relax.log