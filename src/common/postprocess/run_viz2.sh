#!/bin/bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

# ============================================================
# GEWUM viz2 - Structure Analysis Visualization
# Space group Sankey, RDF comparison, structure funnel charts
# ============================================================

# Configuration - modify these parameters as needed
DIM="3d"                        # Structure dimension: 0d, 1d, 2d, 3d

# CIF directory
CIF_DIRS="./Na1Cl2"            # Space-separated for multiple: "./Na2Cl2 ./Na2Cl3 ./Na2Cl5"

# Plot selection
PLOT="all"                      # Which plot(s): all, sankey, chord, rdf, funnel, ehull, violin

# Symmetry analysis
SYMPREC=0.1                     # Symmetry precision for spglib (default: 0.1)

# RDF parameters
RDF_RMAX=10.0                   # RDF maximum radius in Angstrom (default: 10.0)
RDF_BIN=0.1                     # RDF bin size in Angstrom (default: 0.1)

# Sankey parameters
TOP_N=229                        # Top-N flows to show in Sankey diagram (default: 15)
IGNORE_P1=true                 # Exclude P1 flows from Sankey: true or false (default: false)
SANKEY_CMAP_LEFT="viridis"             # Colormap for left column in sankey_spacegroup (e.g. viridis)
SANKEY_CMAP_RIGHT="viridis_r"            # Colormap for right column in sankey_spacegroup (e.g. viridis_r)

# Violin parameters
VIOLIN_TOP_N=20                 # Number of top space groups shown in violin plot

# Ehull parameters
MP_DATA=""                      # Path to offline MP JSON file (e.g., /path/to/MPtrj_2022.9_full.json)
API_KEY=""                      # Materials Project API key (alternative to MP_DATA)
EHULL_COMPAT=false              # Apply MP2020 compatibility corrections: true or false (default: false)
EHULL_CMAP="viridis"                    # Colormap for sankey_ehull.png (e.g. viridis).

# Font settings
FONT_FAMILY="Liberation Sans"  # Font family for figures (default: Liberation Sans)
FONT_SIZE=""                    # Global font size for figures (empty = matplotlib default)

# Structure sampling
MAX_STRUCTURES=500              # Max structures per stage for RDF (default: 200)

# Parallel processing
WORKERS="64"                      # Number of parallel workers (empty = default: min(cpu_count, 4))

# Output settings
OUTPUT_DIR="./viz2_output"      # Output directory
DPI=600                         # Figure DPI

# ============================================================

start_time=$(date +%s)

echo "============================================================"
echo "GEWUM viz2 - Structure Analysis Visualization"
echo "============================================================"
echo "Dimension: $DIM"
echo "Plot type: $PLOT"
echo "CIF dirs: $CIF_DIRS"
echo "Symprec: $SYMPREC"
echo "RDF rmax: $RDF_RMAX, bin: $RDF_BIN"
echo "Top-N (Sankey): $TOP_N"
echo "Violin Top-N: $VIOLIN_TOP_N"
echo "Max structures (RDF): $MAX_STRUCTURES"
echo "Output: $OUTPUT_DIR"
echo "DPI: $DPI"
echo "Ehull compat: $EHULL_COMPAT"
if [ -n "$EHULL_CMAP" ]; then
    echo "Ehull cmap: $EHULL_CMAP"
fi
if [ -n "$MP_DATA" ]; then
    echo "MP data: $MP_DATA"
elif [ -n "$API_KEY" ]; then
    echo "MP API: (key provided)"
else
    echo "MP data: (not configured - Ehull plots will be skipped)"
fi
echo "============================================================"

CMD="python -m gewum.src.common.postprocess.viz2_analysis"

for DIR in $CIF_DIRS; do
    CMD="$CMD --cif-dir $DIR"
done

CMD="$CMD --dim $DIM"
CMD="$CMD --plot $PLOT"
CMD="$CMD --symprec $SYMPREC"
CMD="$CMD --rdf-rmax $RDF_RMAX"
CMD="$CMD --rdf-bin $RDF_BIN"
CMD="$CMD --top-n $TOP_N"

if [ "$IGNORE_P1" = true ]; then
    CMD="$CMD --ignore-p1"
fi

CMD="$CMD --max-structures $MAX_STRUCTURES"
CMD="$CMD -o $OUTPUT_DIR"
CMD="$CMD --dpi $DPI"

if [ -n "$MP_DATA" ]; then
    CMD="$CMD --mp-data $MP_DATA"
fi

if [ -n "$API_KEY" ]; then
    CMD="$CMD --api-key $API_KEY"
fi

if [ "$EHULL_COMPAT" = true ]; then
    CMD="$CMD --ehull-compat"
fi

if [ -n "$EHULL_CMAP" ]; then
    CMD="$CMD --ehull-cmap $EHULL_CMAP"
fi

if [ -n "$FONT_FAMILY" ]; then
    CMD="$CMD --font-family \"$FONT_FAMILY\""
fi
if [ -n "$FONT_SIZE" ]; then
    CMD="$CMD --font-size $FONT_SIZE"
fi

if [ -n "$VIOLIN_TOP_N" ]; then
    CMD="$CMD --violin-top-n $VIOLIN_TOP_N"
fi

if [ -n "$WORKERS" ]; then
    CMD="$CMD --workers $WORKERS"
fi

if [ -n "$SANKEY_CMAP_LEFT" ]; then
    CMD="$CMD --sankey-cmap-left $SANKEY_CMAP_LEFT"
fi
if [ -n "$SANKEY_CMAP_RIGHT" ]; then
    CMD="$CMD --sankey-cmap-right $SANKEY_CMAP_RIGHT"
fi

echo "Running: $CMD"
eval $CMD

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo "============================================================"
echo "viz2 analysis completed!"
echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"
echo "============================================================"
