#!/bin/bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

start_time=$(date +%s)

# ============ Configurable Parameters ============
TEMPERATURE=${1:-300}      # Temperature in Kelvin (default: 300)
STEPS=${2:-10000}          # Number of MD steps (default: 10000)
TIMESTEP=${3:-1.0}         # Timestep in fs (default: 1.0)
SUPERCELL_X=${4:-2}        # Supercell expansion in X (default: 2)
SUPERCELL_Y=${5:-2}        # Supercell expansion in Y (default: 2)
SUPERCELL_Z=${6:-1}        # Supercell expansion in Z (default: 1)
DUMP_INTERVAL=${7:-100}    # Trajectory dump interval (default: 100)
FRICTION=${8:-0.01}        # Langevin friction (default: 0.01)

echo "=============================================="
echo "GEWUM NVT Molecular Dynamics"
echo "=============================================="
echo "Temperature: $TEMPERATURE K"
echo "Steps: $STEPS"
echo "Timestep: $TIMESTEP fs"
echo "Supercell: ${SUPERCELL_X}x${SUPERCELL_Y}x${SUPERCELL_Z}"
echo "Dump interval: $DUMP_INTERVAL"
echo "Friction: $FRICTION"
echo "=============================================="

{{ENV_SETUP}}

cif_files=(*.cif)

if [ ${#cif_files[@]} -eq 0 ] || [ "${cif_files[0]}" == "*.cif" ]; then
    echo "Error: No .cif files found in current directory!"
    exit 1
fi

echo "Found ${#cif_files[@]} CIF file(s)"

for cif_file in "${cif_files[@]}"; do
    base_name="${cif_file%.cif}"
    output_dir="md_output_${base_name}_T${TEMPERATURE}K"
    
    echo ""
    echo "Processing: $cif_file -> $output_dir"
    
    python -m gewum.src.MDworkflows.md_nvt "$cif_file" "$output_dir" \
        --temperature $TEMPERATURE \
        --steps $STEPS \
        --timestep $TIMESTEP \
        --supercell $SUPERCELL_X $SUPERCELL_Y $SUPERCELL_Z \
        --dump-interval $DUMP_INTERVAL \
        --friction $FRICTION
    
    if [ $? -eq 0 ]; then
        echo "[OK] Completed: $cif_file"
    else
        echo "[FAIL] Failed: $cif_file"
    fi
done

echo ""
echo "=============================================="
echo "All MD simulations completed"
echo "=============================================="

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Total runtime: $hours hours, $minutes minutes, $seconds seconds" | tee Time_md_nvt.log
