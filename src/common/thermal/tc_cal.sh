#!/bin/bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

# ============================================================
# GEWUM Thermal Conductivity - FC3 Calculation
# Computes third-order force constants for all CIF files
# ============================================================

unset DISPLAY
start_time=$(date +%s)

# Configurable parameters
SUPERCELL=${1:-"2 2 2"}
DEVICE=${2:-cpu}

echo "============================================================"
echo "GEWUM Thermal Conductivity - FC3 Calculation"
echo "============================================================"
echo "  Supercell: $SUPERCELL"
echo "  Device: $DEVICE"
echo "============================================================"

{{ENV_SETUP}}

for cif_file in *.cif; do
    if [ ! -f "$cif_file" ]; then
        echo "No CIF files found"
        exit 1
    fi
    
    prefix="${cif_file%.cif}"
    echo "Processing: $cif_file -> $prefix/"
    
    mkdir -p "$prefix"
    cp "$cif_file" "$prefix/"
    
    pushd "$prefix" > /dev/null
    python -m gewum.src.common.thermal.tc_fc3 --supercell $SUPERCELL --device $DEVICE
    popd > /dev/null
    
    echo "Completed: $prefix"
done

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo "============================================================"
echo "FC3 calculation completed!"
echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"
echo "============================================================" | tee Time_tc_fc3.log
