#!/bin/bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

# ============================================================
# GEWUM Structure Visualization Script
# Feature extraction, dimensionality reduction and comparison
# ============================================================

# Configuration - modify these parameters as needed
DIM="3d"                        # Options: 0d, 1d, 2d, 3d
DESCRIPTOR="simple"             # Options: simple, soap, coulomb

# Option 1: Auto-detect mode (recommended)
CIF_DIRS="./Na2Cl2"            # Space-separated for multiple: "./Na2Cl2 ./Na2Cl3 ./Na2Cl5"

# Option 2: Manual mode (for fine-grained control)
# Leave CIF_DIRS empty and set individual directories below
TOTAL_DIR=""                    # Total/reference CIF directory
SELECT_DIR=""                   # Selected CIF directory (optional)
RELAX_DIR=""                    # Relaxed CIF directory (optional)


REDUCTION="umap"                # Reduction method: umap, tsne
OUTPUT_DIR="./visualization"    # Output directory
LOAD_FEATURES=""                # Load saved .npz feature file (optional)

# SOAP parameters (used when DESCRIPTOR="soap")
SOAP_R_CUT=""                   # SOAP cutoff radius (auto per dimension if empty)
SOAP_N_MAX=4                    # SOAP radial basis n_max
SOAP_L_MAX=4                    # SOAP angular basis l_max
PCA_VARIANCE=0.95               # PCA variance ratio to preserve

# Plot settings
DPI=600                         # Figure DPI
TITLE=""                        # Custom plot title (optional)
NO_DENSITY=""                   # Set to "true" to skip density heatmap
FONT_FAMILY="Liberation Sans"  # Font family for plots
FONT_SIZE=""                    # Font size (optional, uses default if empty)
WORKERS="64"                      # Parallel workers (optional)

# Word cloud settings
CIFGEN_INP="./cifgen.inp"                   # Path to cifgen.inp for composition word cloud
WC_CMAP=""                      # Word cloud colormap (optional, default: pink-blue)

# ============================================================
# Dimension-specific defaults
# ============================================================
# 3D: descriptor=simple/soap,            r_cut=6.0
# 2D: descriptor=simple/soap,            r_cut=5.0
# 1D: descriptor=soap only,              r_cut=4.0
# 0D: descriptor=simple/soap/coulomb,    r_cut=5.0
# ============================================================

start_time=$(date +%s)

echo "============================================================"
echo "GEWUM Structure Visualization"
echo "============================================================"
echo "Dimension: $DIM"
echo "Descriptor: $DESCRIPTOR"
echo "Reduction: $REDUCTION"
echo "Output: $OUTPUT_DIR"
if [ -n "$TOTAL_DIR" ] || [ -n "$SELECT_DIR" ] || [ -n "$RELAX_DIR" ]; then
    echo "Mode: Manual (--total-dir / --select-dir / --relax-dir)"
fi
if [ -n "$CIF_DIRS" ]; then
    echo "Mode: Auto-detect"
    echo "CIF dirs: $CIF_DIRS"
fi
if [ -n "$TOTAL_DIR" ]; then
    echo "Total dir: $TOTAL_DIR"
fi
if [ -n "$SELECT_DIR" ]; then
    echo "Select dir: $SELECT_DIR"
fi
if [ -n "$RELAX_DIR" ]; then
    echo "Relax dir: $RELAX_DIR"
fi
if [ -n "$LOAD_FEATURES" ]; then
    echo "Load features: $LOAD_FEATURES"
fi
if [ "$DESCRIPTOR" = "soap" ]; then
    echo "SOAP params: r_cut=${SOAP_R_CUT:-auto}, n_max=$SOAP_N_MAX, l_max=$SOAP_L_MAX"
    echo "PCA variance: $PCA_VARIANCE"
fi
echo "DPI: $DPI"
echo "============================================================"

CMD="python -m gewum.src.common.postprocess.visualization"
CMD="$CMD --dim $DIM --descriptor $DESCRIPTOR"
CMD="$CMD --reduction $REDUCTION"
CMD="$CMD -o $OUTPUT_DIR"

if [ -n "$LOAD_FEATURES" ]; then
    CMD="$CMD --load-features $LOAD_FEATURES"
elif [ -n "$CIF_DIRS" ]; then
    for DIR in $CIF_DIRS; do
        CMD="$CMD --cif-dir $DIR"
    done
else
    if [ -n "$TOTAL_DIR" ]; then
        CMD="$CMD --total-dir $TOTAL_DIR"
    fi
    if [ -n "$SELECT_DIR" ]; then
        CMD="$CMD --select-dir $SELECT_DIR"
    fi
    if [ -n "$RELAX_DIR" ]; then
        CMD="$CMD --relax-dir $RELAX_DIR"
    fi
fi

if [ "$DESCRIPTOR" = "soap" ]; then
    if [ -n "$SOAP_R_CUT" ]; then
        CMD="$CMD --soap-r-cut $SOAP_R_CUT"
    fi
    CMD="$CMD --soap-n-max $SOAP_N_MAX --soap-l-max $SOAP_L_MAX"
    CMD="$CMD --pca-variance $PCA_VARIANCE"
fi

CMD="$CMD --dpi $DPI"
if [ -n "$TITLE" ]; then
    CMD="$CMD --title \"$TITLE\""
fi
if [ "$NO_DENSITY" = "true" ]; then
    CMD="$CMD --no-density"
fi
if [ -n "$FONT_FAMILY" ]; then
    CMD="$CMD --font-family \"$FONT_FAMILY\""
fi
if [ -n "$FONT_SIZE" ]; then
    CMD="$CMD --font-size $FONT_SIZE"
fi
if [ -n "$WORKERS" ]; then
    CMD="$CMD --workers $WORKERS"
fi
if [ -n "$CIFGEN_INP" ]; then
    CMD="$CMD --cifgen-inp $CIFGEN_INP"
fi
if [ -n "$WC_CMAP" ]; then
    CMD="$CMD --wc-cmap $WC_CMAP"
fi

echo "Running: $CMD"
eval $CMD

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo "============================================================"
echo "Visualization completed!"
echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"
echo "============================================================"
