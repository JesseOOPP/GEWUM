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
MODE=${1:-1}           # 1=atoms only, 2=atoms+cell (default: 2 for RD)
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

run_optimization_db() {
    DB_PATH=$1
    TMP_DIR=$2
    SG=$3
    NAME=$4
    OPT_MODE=$5
    OPT_FMAX=$6
    OPT_MAX_STEPS=$7

    if [ -z "$DB_PATH" ] || [ -z "$SG" ] || [ -z "$NAME" ]; then
        echo "Error: missing DB/SG/NAME. Skipping."
        return 1
    fi

    python -m gewum.src.common.relaxation.relax_db_io worker \
        --db "$DB_PATH" --sg "$SG" --name "$NAME" \
        --tmp-dir "$TMP_DIR" \
        --mode "$OPT_MODE" --fmax "$OPT_FMAX" --max-steps "$OPT_MAX_STEPS" \
        &>> cif_relax.out
}
export -f run_optimization_db

run_optimization() {
    CIF_FILE=$1
    RELAXED_DIR=$2
    OPT_MODE=$3
    OPT_FMAX=$4
    OPT_MAX_STEPS=$5

    if [ -z "$CIF_FILE" ] || [ -z "$RELAXED_DIR" ]; then
        echo "Error: CIF_FILE or RELAXED_DIR is empty. Skipping this task."
        return 1
    fi

    mkdir -p "$RELAXED_DIR"
    python -m gewum.src.common.relaxation.umlip_relax "$CIF_FILE" "$RELAXED_DIR" \
        --mode "$OPT_MODE" --fmax "$OPT_FMAX" --max-steps "$OPT_MAX_STEPS" \
        &>> cif_relax.out
}
export -f run_optimization

TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-64}
CORES_PER_TASK=1
NUM_TASKS=$((TOTAL_CPUS / CORES_PER_TASK))
export OMP_NUM_THREADS=$CORES_PER_TASK

formula_dirs=()
for entry in "$top_dir"/*/; do
    [ ! -d "$entry" ] && continue
    dirname=$(basename "$entry")
    [[ "$dirname" == 0_* ]] && continue
    [[ "$dirname" == "final_cifs" ]] && continue
    [[ "$dirname" == ".energy_temp" ]] && continue
    formula_dirs+=("$entry")
done

echo "Found ${#formula_dirs[@]} formula directories to process"

total_relaxed=0
total_skipped=0
total_new=0
total_relaxed_db=0  

db_count_relaxed() {
    local n
    n=$(sqlite3 "$1" "SELECT COUNT(*) FROM structures WHERE stage='relaxed';" 2>/dev/null)
    if [ -z "$n" ]; then
        n=$(python -c "
import sqlite3, sys
try:
    c = sqlite3.connect(sys.argv[1])
    print(c.execute(\"SELECT COUNT(*) FROM structures WHERE stage='relaxed'\").fetchone()[0])
except Exception:
    print(0)
" "$1" 2>/dev/null)
    fi
    echo "${n:-0}"
}

for formula_dir in "${formula_dirs[@]}"; do
    formula_name=$(basename "$formula_dir")
    formula_dir_clean="${formula_dir%/}"
    echo ""
    echo "============================================================"
    echo "Processing formula: $formula_name"
    echo "============================================================"

    DB_PATH="${formula_dir_clean}/structures.db"
    TMP_DIR="${formula_dir_clean}/.relax_tmp"
    CSV_PATH="${formula_dir_clean}/energy_results.csv"
    tasks_file="${formula_dir_clean}/tasks.txt"
    > "$tasks_file"

    if [ -f "$DB_PATH" ]; then
        # --- DB MODE ---
        echo "  Mode: DB ($DB_PATH)"
        rm -rf "$TMP_DIR"
        mkdir -p "$TMP_DIR"
        python -c "
from gewum.src.common.cif_db import CifDatabase
with CifDatabase('${formula_dir_clean}') as db:
    for sg, name, _ in db.query_initial_tasks():
        print(sg, name)
" > "$tasks_file"

        new_tasks=$(wc -l < "$tasks_file")
        echo "  Initial tasks in DB: $new_tasks"
        total_new=$((total_new + new_tasks))

        if [ -s "$tasks_file" ]; then
            parallel -j $NUM_TASKS -a "$tasks_file" --colsep ' ' \
                run_optimization_db "$DB_PATH" "$TMP_DIR" \
                {1} {2} $MODE $FMAX $MAX_STEPS
        fi

        python -m gewum.src.common.relaxation.relax_db_io commit \
            --db "$DB_PATH" --tmp-dir "$TMP_DIR" --csv "$CSV_PATH" --cleanup

        formula_relaxed=$(db_count_relaxed "$DB_PATH")
        echo "  Cumulative relaxed (DB, this formula): $formula_relaxed"
        total_relaxed_db=$((total_relaxed_db + formula_relaxed))
    else
        echo "  Mode: disk (no structures.db)"

        all_cifs=$(mktemp /tmp/all_cifs.XXXXXX)
        find "$formula_dir_clean" -mindepth 2 -maxdepth 2 -name '*.cif' -not -path '*/relaxed/*' > "$all_cifs"
        formula_cif_count=$(wc -l < "$all_cifs")
        echo "  CIF files found: $formula_cif_count"

        relaxed_set=$(mktemp /tmp/relaxed_set.XXXXXX)
        find "$formula_dir_clean" -mindepth 3 -maxdepth 3 -name '*_relaxed.cif' -path '*/relaxed/*' > "$relaxed_set"

        skip_count=0
        while IFS= read -r cif_file_path; do
            dir_path=$(dirname "$cif_file_path")
            cif_file=$(basename "$cif_file_path")
            base_name="${cif_file%.cif}"
            relaxed_file="$dir_path/relaxed/${base_name}_relaxed.cif"

            if grep -qF "$relaxed_file" "$relaxed_set"; then
                skip_count=$((skip_count + 1))
                continue
            fi

            echo "$cif_file_path $dir_path/relaxed" >> "$tasks_file"
        done < "$all_cifs"

        rm -f "$all_cifs" "$relaxed_set"

        new_tasks=$(wc -l < "$tasks_file")
        echo "  Already relaxed (skipped): $skip_count"
        echo "  New tasks: $new_tasks"
        total_skipped=$((total_skipped + skip_count))
        total_new=$((total_new + new_tasks))

        if [ -s "$tasks_file" ]; then
            parallel -j $NUM_TASKS -a "$tasks_file" --colsep ' ' \
                run_optimization {1} {2} $MODE $FMAX $MAX_STEPS
        fi
    fi

    rm -f "$tasks_file"
    echo "  Processing complete for $formula_name."
done

[ -f cif_relax.out ] && rm cif_relax.out
echo "$total_relaxed_db" >> TOT_relaxed

echo ""
echo "============================================================"
echo "Summary: ${#formula_dirs[@]} formulas, $total_new new tasks, $total_skipped skipped"
echo "============================================================"

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Total runtime: $hours hours, $minutes minutes, $seconds seconds" > Time_cif_relax.log
