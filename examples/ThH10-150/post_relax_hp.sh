#!/bin/bash

#SBATCH --job-name=post_relax_hp
#SBATCH --output=post_relax_hp.out
#SBATCH --error=post_relax_hp.err
#SBATCH --time=2400:00:00
#SBATCH --cpus-per-task=64
#SBATCH -p <partition>
#SBATCH -N 1

start_time=$(date +%s)

module purge
module load cmake/<version>
module load gcc/<version>
module load intel/<version>
module load mpi/<version>
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate <env_name>
export PATH=/path/to/your/tools/bin:$PATH   # adapt to your cluster

# ============ High Pressure Configuration ============
# Target pressure in GPa
TARGET_PRESSURE=${1:-150.0}

# Pressure tolerance in GPa (structures outside this range will be filtered)
PRESSURE_TOLERANCE=${2:-10.0}

# Bond length threshold in Angstrom
BOND_THRESHOLD=${3:-0.8}

# Minimum enthalpy difference between selected structures (eV/atom)
ENTHALPY_GAP=${4:-0.001}

# Maximum number of structures to select per chemistry
MAX_N=${5:-100}
# ======================================================

top_dir=$(pwd)
echo "=== High Pressure Post-processing ==="
echo "Top directory: $top_dir"
echo "Target pressure: $TARGET_PRESSURE GPa (tolerance: +/- $PRESSURE_TOLERANCE GPa)"
echo "Bond threshold: $BOND_THRESHOLD Angstrom"
echo "Enthalpy gap: $ENTHALPY_GAP eV/atom"
echo "Max structures: $MAX_N"

if [ ! -d "$top_dir" ]; then
    echo "Top directory does not exist: $top_dir"
    exit 1
fi

run_post_processing() {
    local chem_dir="$1"
    local pressure="$2"
    local tolerance="$3"
    local gap="$4"
    local max_n="$5"
    
    if [ ! -d "$chem_dir" ]; then
        echo "Error: Chemistry directory does not exist: $chem_dir"
        return 1
    fi
    
    (cd "$chem_dir" && python -m gewum.src.RDworkflows.energy_post_hp --pressure "$pressure" --tolerance "$tolerance" --gap "$gap" --max-n "$max_n")
}

export -f run_post_processing

# ============================================================
# Per-formula processing loop: unpack -> post-process -> pack
# Keeps peak inode usage to ONE formula directory at a time
# ============================================================

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

    # Post-processing reads files directly (no unpack needed with DB)
    echo "  Running bond length check..."
    python -m gewum.src.common.postprocess.bond_check --threshold $BOND_THRESHOLD --base-dir "$formula_dir"

    # Step 3: Energy post-processing
    run_post_processing "$formula_dir" "$TARGET_PRESSURE" "$PRESSURE_TOLERANCE" "$ENTHALPY_GAP" "$MAX_N"
done

# Collect all results
echo "Chemical_Formula,CIF_Base_Name,Total_Energy_eV,Energy_per_Atom_eV,Final_Pressure_GPa,Enthalpy_per_Atom_eV,Corrected_Enthalpy_per_Atom_eV,Relaxed_CIF_Path,SG_ori" > 0_final_result_tot.txt

for formula_dir in "${formula_dirs[@]}"; do
    result_file="$formula_dir/0_final_results.txt"
    if [ -f "$result_file" ]; then
        tail -n +2 "$result_file" >> 0_final_result_tot.txt
    fi
done

mkdir -p final_cifs
cp ./*/0_cif_final/*.cif ./final_cifs 2>/dev/null


end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Total runtime: $hours hours, $minutes minutes, $seconds seconds" > Time_post_processing_hp.log

echo "=== High Pressure Post-processing Complete ==="
echo "Output files:"
echo "  - 0_final_result_tot.txt (all selected structures)"
echo "  - */0_final_results.txt (per-chemistry results)"
echo "  - */0_cif_final/ (selected CIF files)"
