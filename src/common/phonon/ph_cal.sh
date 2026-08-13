#!/bin/bash
# GEWUM Phonon Calculation Script
# Calculates phonon dispersion using MLIP

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

start_time=$(date +%s)

{{ENV_SETUP}}

# ============ Phonon Configuration ============
# Supercell dimensions (default: 2 2 2)
SUPERCELL="2 2 2"

# Band path mode: auto / manual / default
# - auto: Automatically detect high-symmetry path from crystal structure
# - manual: Use custom band path specified in BAND_PATH
# - default: Use default cubic path (G-X-M-R-G)
PATH_MODE="auto"

# Custom band path (only used when PATH_MODE=manual)
# Format: "G:0,0,0|X:0.5,0,0|M:0.5,0.5,0|G:0,0,0"
BAND_PATH=""

# Band dispersion line color (e.g., purple, blue, red, green, black, #FF5733)
BAND_COLOR="purple"
# ===============================================

PH_ARGS=(--supercell $SUPERCELL)

if [ "$PATH_MODE" = "manual" ] && [ -n "$BAND_PATH" ]; then
    PH_ARGS+=(--band-path "$BAND_PATH")
elif [ "$PATH_MODE" = "default" ]; then
    PH_ARGS+=(--no-auto-path)
fi

if [ -n "$BAND_COLOR" ]; then
    PH_ARGS+=(--band-color "$BAND_COLOR")
fi

python -m gewum.src.common.phonon.ph_uMLIP "${PH_ARGS[@]}" >> ph.log

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Total runtime: $hours hours, $minutes minutes, $seconds seconds" > Time_phonon.log
