#!/bin/bash

#SBATCH --job-name=run_selection
#SBATCH --output=run_selection.out
#SBATCH --error=run_selection.err
#SBATCH --time=2400:00:00
#SBATCH --cpus-per-task=64
#SBATCH -p <partition>
#SBATCH -N 1

# ============================================================
# GEWUM Unified Structure Selection Script
# Selects diverse structures for all dimensions (0D-3D)
# ============================================================

# Configuration - modify these parameters as needed
DIM="3d"               # Options: 3d, 2d, 1d, 0d
METHOD="medoid"        # Options: random, kmeans, medoid, maxmin, hdbscan
DESCRIPTOR="simple"    # Options: simple, soap (1D: soap only; 0D: +coulomb)
TARGET_COUNT=40        # Number of structures to keep (ignored for hdbscan)
WORKERS=64             # Number of parallel workers
INPUT_DIR="."          # Base directory

# HDBSCAN parameters (used when METHOD="hdbscan")
# Default: Moderate filtering mode (auto-adaptive)
MIN_CLUSTER_SIZE=3           # Minimum cluster size (5=moderate, 15+=strict, 3=relaxed)
MIN_SAMPLES=3                # Min samples for core points (3=moderate, 10+=strict)
CLUSTER_EPSILON=0.0          # Distance threshold (0.0=auto, 0.1-0.5=manual merge)
ALPHA=1.0                    # Density decay parameter (1.0=moderate, 2.0+=strict)
KEEP_NOISE="true"            # Keep noise structures: "true"=keep (default), "false"=discard

# SOAP parameters (used when DESCRIPTOR="soap")
SOAP_R_CUT=6.0         # SOAP cutoff radius (auto-adjusted per dimension)
SOAP_N_MAX=4           # SOAP radial basis
SOAP_L_MAX=4           # SOAP angular basis

# PCA parameters (auto-enabled for SOAP/Coulomb descriptors)
# PCA reduces high-dim features for better clustering
USE_PCA="auto"         # Options: auto, yes, no (auto=enable for soap/coulomb)
PCA_VARIANCE=0.9      # Variance ratio to preserve (default: 95%)
PCA_COMPONENTS=        # Fixed components (leave empty to use variance ratio)

# 0D Coulomb Matrix parameters (used when DESCRIPTOR="coulomb")
COULOMB_N_ATOMS_MAX=100

# 1D specific parameters
PERIODIC_DIR="z"       # Periodic direction for 1D: x, y, or z

# ============================================================
# Dimension-specific defaults
# ============================================================
# 3D: descriptor=simple/soap, r_cut=6.0
# 2D: descriptor=simple/soap, r_cut=5.0 (smaller due to vacuum)
# 1D: descriptor=soap only,   r_cut=4.0 (constrained by vacuum)
# 0D: descriptor=simple/soap/coulomb, r_cut=5.0 (non-periodic)
# ============================================================

start_time=$(date +%s)

echo "============================================================"
echo "GEWUM Unified Structure Selection"
echo "============================================================"
echo "Dimension: $DIM"
echo "Method: $METHOD"
echo "Descriptor: $DESCRIPTOR"
if [ "$METHOD" = "hdbscan" ]; then
    echo "HDBSCAN min_cluster_size: $MIN_CLUSTER_SIZE"
    echo "HDBSCAN min_samples: $MIN_SAMPLES"
    if [ "$CLUSTER_EPSILON" != "0.0" ]; then
        echo "HDBSCAN epsilon: $CLUSTER_EPSILON"
    fi
    if [ "$ALPHA" != "1.0" ]; then
        echo "HDBSCAN alpha: $ALPHA"
    fi
    echo "HDBSCAN keep_noise: $KEEP_NOISE"
else
    echo "Target count: $TARGET_COUNT"
fi
echo "Workers: $WORKERS"
echo "Input directory: $INPUT_DIR"
if [ "$DESCRIPTOR" = "soap" ]; then
    echo "SOAP params: r_cut=$SOAP_R_CUT, n_max=$SOAP_N_MAX, l_max=$SOAP_L_MAX"
    if [ "$USE_PCA" != "no" ]; then
        echo "PCA: enabled (variance=$PCA_VARIANCE)"
    fi
elif [ "$DESCRIPTOR" = "coulomb" ]; then
    echo "Coulomb Matrix: n_atoms_max=$COULOMB_N_ATOMS_MAX"
    if [ "$USE_PCA" != "no" ]; then
        echo "PCA: enabled (variance=$PCA_VARIANCE)"
    fi
fi
if [ "$DIM" = "1d" ]; then
    echo "Periodic direction: $PERIODIC_DIR"
fi
echo "============================================================"

# ============================================================
# Build selection command as a function for per-formula calling
# ============================================================
run_selection() {
    local FORMULA_DIR=$1

    CMD="python -m gewum.src.common.selection.structure_select \
        --dim $DIM \
        --method $METHOD \
        --target $TARGET_COUNT \
        --input-dir $FORMULA_DIR \
        --workers $WORKERS"

    if [ "$METHOD" = "hdbscan" ]; then
        CMD="$CMD --min-cluster-size $MIN_CLUSTER_SIZE"
        CMD="$CMD --min-samples $MIN_SAMPLES"
        if [ "$CLUSTER_EPSILON" != "0.0" ]; then
            CMD="$CMD --cluster-selection-epsilon $CLUSTER_EPSILON"
        fi
        if [ "$ALPHA" != "1.0" ]; then
            CMD="$CMD --alpha $ALPHA"
        fi
        if [ "$KEEP_NOISE" = "false" ]; then
            CMD="$CMD --no-keep-noise"
        fi
    fi

    if [ "$DESCRIPTOR" = "soap" ]; then
        CMD="$CMD --descriptor soap \
            --soap-r-cut $SOAP_R_CUT \
            --soap-n-max $SOAP_N_MAX \
            --soap-l-max $SOAP_L_MAX"
    elif [ "$DESCRIPTOR" = "coulomb" ]; then
        CMD="$CMD --descriptor coulomb \
            --coulomb-n-atoms-max $COULOMB_N_ATOMS_MAX"
    else
        CMD="$CMD --descriptor simple"
    fi

    if [ "$USE_PCA" = "no" ]; then
        CMD="$CMD --no-pca"
    elif [ -n "$PCA_COMPONENTS" ]; then
        CMD="$CMD --pca-components $PCA_COMPONENTS"
    else
        CMD="$CMD --pca-variance $PCA_VARIANCE"
    fi

    if [ "$DIM" = "1d" ]; then
        CMD="$CMD --periodic-dir $PERIODIC_DIR"
    fi

    eval $CMD
}

# ============================================================
# Per-formula processing loop: select directly from DB
# (No pack/unpack needed with SQLite storage)
# ============================================================

top_dir=$(pwd)
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
    echo "Selecting structures in: $formula_name"
    echo "============================================================"

    # Run selection on this formula directory (reads from structures.db if present)
    run_selection "$formula_dir"
done

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

# ============================================================
# Final summary across all formula directories
# Reads stage counts directly from each structures.db
# (sqlite3 CLI preferred, falls back to a tiny Python query).
# Stage semantics in DB mode (single label per row):
#   selected = stage != 'removed'   (kept by selection)
#   removed  = stage  = 'removed'   (filtered out)
#   total    = all rows
# ============================================================
db_count_stage() {
    local db=$1
    local where=$2
    local n
    if [ -n "$where" ]; then
        n=$(sqlite3 "$db" "SELECT COUNT(*) FROM structures WHERE $where;" 2>/dev/null)
    else
        n=$(sqlite3 "$db" "SELECT COUNT(*) FROM structures;" 2>/dev/null)
    fi
    if [ -z "$n" ]; then
        n=$(python -c "
import sqlite3, sys
try:
    c = sqlite3.connect(sys.argv[1])
    q = 'SELECT COUNT(*) FROM structures' + ((' WHERE ' + sys.argv[2]) if sys.argv[2] else '')
    print(c.execute(q).fetchone()[0])
except Exception:
    print(0)
" "$db" "$where" 2>/dev/null)
    fi
    echo "${n:-0}"
}

sum_total=0
sum_selected=0
sum_removed=0
formula_count=0

echo ""
echo "============================================================"
echo "Selection summary (per formula)"
echo "============================================================"
printf "%-30s %10s %10s %10s\n" "Formula" "Total" "Selected" "Removed"
printf "%-30s %10s %10s %10s\n" "------------------------------" "----------" "----------" "----------"
for formula_dir in "${formula_dirs[@]}"; do
    formula_name=$(basename "$formula_dir")
    db_path="${formula_dir%/}/structures.db"
    [ ! -f "$db_path" ] && continue
    n_total=$(db_count_stage    "$db_path" "")
    n_selected=$(db_count_stage "$db_path" "stage != 'removed'")
    n_removed=$(db_count_stage  "$db_path" "stage = 'removed'")
    printf "%-30s %10s %10s %10s\n" "$formula_name" "$n_total" "$n_selected" "$n_removed"
    sum_total=$((sum_total + n_total))
    sum_selected=$((sum_selected + n_selected))
    sum_removed=$((sum_removed + n_removed))
    formula_count=$((formula_count + 1))
done
printf "%-30s %10s %10s %10s\n" "------------------------------" "----------" "----------" "----------"
printf "%-30s %10s %10s %10s\n" "TOTAL ($formula_count formulas)" "$sum_total" "$sum_selected" "$sum_removed"

# Persist scalar totals to small text files (consistent with TOT_relaxed convention)
echo "$sum_total"    >> TOT_cif
echo "$sum_selected" >> TOT_selected
echo "$sum_removed"  >> TOT_removed

echo "============================================================"
echo "Selection completed!"
echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"
echo "============================================================"

# (Pack is handled per-formula in the loop above)
