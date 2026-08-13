#!/bin/bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

# ============================================================
# GEWUM PT viz - Structure Analysis Visualization
# Space group Sankey, energy violin, RDF comparison
# ============================================================

# Configuration - modify these parameters as needed

# CIF directories (PT work directories containing *.cif and relaxed/)
CIF_DIRS="."                    # Space-separated for multiple: "./dir1 ./dir2"

# Plot selection
PLOT="all"                      # Which plot(s): all, sankey, violin, rdf

# Symmetry analysis
SYMPREC=0.1                     # Symmetry precision for spglib (default: 0.1)

# RDF parameters
RDF_RMAX=10.0                   # RDF maximum radius in Angstrom (default: 10.0)
RDF_BIN=0.1                     # RDF bin size in Angstrom (default: 0.1)

# Sankey / Violin parameters
IGNORE_P1=false                 # Exclude P1 flows: true or false (default: false)
VIOLIN_TOP_N=10                 # Number of top space groups shown in violin plot
SANKEY_CMAP_LEFT="autumn"             # Colormap for left column (initial SG) in sankey_spacegroup.png (e.g. viridis)
SANKEY_CMAP_RIGHT="autumn_r"            # Colormap for right column (relaxed SG) in sankey_spacegroup.png (e.g. viridis_r)

# Font settings
FONT_FAMILY="Liberation Sans"  # Font family for figures (default: Liberation Sans)
FONT_SIZE=""                    # Global font size for figures (empty = matplotlib default)

# Structure sampling
MAX_STRUCTURES=200              # Max structures per stage for RDF (default: 200)

# Parallel processing
WORKERS="64"                      # Number of parallel workers (empty = default: min(cpu_count, 4))

# Output settings
OUTPUT_DIR="./viz_pt_output"    # Output directory
DPI=600                         # Figure DPI

# ============================================================

start_time=$(date +%s)

{{ENV_SETUP}}

echo "============================================================"
echo "GEWUM PT viz - Structure Analysis Visualization"
echo "============================================================"
echo "Plot type: $PLOT"
echo "CIF dirs: $CIF_DIRS"
echo "Symprec: $SYMPREC"
echo "RDF rmax: $RDF_RMAX, bin: $RDF_BIN"
echo "Ignore P1: $IGNORE_P1"
echo "Violin Top-N: $VIOLIN_TOP_N"
if [ -n "$SANKEY_CMAP_LEFT" ]; then
    echo "Sankey cmap left: $SANKEY_CMAP_LEFT"
fi
if [ -n "$SANKEY_CMAP_RIGHT" ]; then
    echo "Sankey cmap right: $SANKEY_CMAP_RIGHT"
fi
echo "Max structures (RDF): $MAX_STRUCTURES"
echo "Output: $OUTPUT_DIR"
echo "DPI: $DPI"
echo "============================================================"

CMD="python -m gewum.src.PTworkflows.viz_pt"

for DIR in $CIF_DIRS; do
    CMD="$CMD --cif-dir $DIR"
done

CMD="$CMD --plot $PLOT"
CMD="$CMD --symprec $SYMPREC"
CMD="$CMD --rdf-rmax $RDF_RMAX"
CMD="$CMD --rdf-bin $RDF_BIN"

if [ "$IGNORE_P1" = true ]; then
    CMD="$CMD --ignore-p1"
fi

CMD="$CMD --max-structures $MAX_STRUCTURES"
CMD="$CMD -o $OUTPUT_DIR"
CMD="$CMD --dpi $DPI"

if [ -n "$FONT_FAMILY" ]; then
    CMD="$CMD --font-family \"$FONT_FAMILY\""
fi
if [ -n "$FONT_SIZE" ]; then
    CMD="$CMD --font-size $FONT_SIZE"
fi

if [ -n "$VIOLIN_TOP_N" ]; then
    CMD="$CMD --violin-top-n $VIOLIN_TOP_N"
fi

if [ -n "$SANKEY_CMAP_LEFT" ]; then
    CMD="$CMD --sankey-cmap-left $SANKEY_CMAP_LEFT"
fi
if [ -n "$SANKEY_CMAP_RIGHT" ]; then
    CMD="$CMD --sankey-cmap-right $SANKEY_CMAP_RIGHT"
fi

if [ -n "$WORKERS" ]; then
    CMD="$CMD --workers $WORKERS"
fi

echo "Running: $CMD"
eval $CMD

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo "============================================================"
echo "PT viz analysis completed!"
echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"
echo "============================================================"
