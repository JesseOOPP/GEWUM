#!/bin/bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

# ============ Post-processing Configuration ============
# Bond length threshold in Angstrom
BOND_THRESHOLD=${1:-1.0}

# Minimum energy difference between selected structures (eV/atom)
ENERGY_GAP=${2:-0.001}

# Maximum number of structures to select per chemistry
MAX_N=${3:-100}
# ========================================================

start_time=$(date +%s)

{{ENV_SETUP}}

top_dir=$(pwd)
echo "Top directory: $top_dir"
echo "Bond threshold: $BOND_THRESHOLD Angstrom"
echo "Energy gap: $ENERGY_GAP eV/atom"
echo "Max structures: $MAX_N"

if [ ! -d "$top_dir" ]; then
    echo "Top directory does not exist: $top_dir"
    exit 1
fi

run_post_processing() {
    local chem_dir="$1"
    local gap="$2"
    local max_n="$3"
    
    if [ ! -d "$chem_dir" ]; then
        echo "Error: Chemistry directory does not exist: $chem_dir"
        return 1
    fi
    
    (cd "$chem_dir" && python -m gewum.src.RDworkflows.energy_post --gap "$gap" --max-n "$max_n")
}

export -f run_post_processing

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

for formula_dir in "${formula_dirs[@]}"; do
    formula_name=$(basename "$formula_dir")
    echo ""
    echo "============================================================"
    echo "Post-processing formula: $formula_name"
    echo "============================================================"
    echo "  Running bond length check..."
    python -m gewum.src.common.postprocess.bond_check --threshold $BOND_THRESHOLD --base-dir "$formula_dir"
    run_post_processing "$formula_dir" $ENERGY_GAP $MAX_N

done

echo "Chemical_Formula,CIF_Base_Name,Total_Energy_eV,Energy_per_Atom_eV,Relaxed_CIF_Path,SG_ori" > 0_final_result_tot.txt

for formula_dir in "${formula_dirs[@]}"; do
    result_file="$formula_dir/0_final_results.txt"
    if [ -f "$result_file" ]; then
        tail -n +2 "$result_file" >> 0_final_result_tot.txt
    fi
done

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Total runtime: $hours hours, $minutes minutes, $seconds seconds" > Time_post_processing.log
