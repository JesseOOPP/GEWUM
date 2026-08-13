#!/bin/bash

#SBATCH --job-name=run_srss
#SBATCH --output=run_srss.out
#SBATCH --error=run_srss.err
#SBATCH --time=2400:00:00
#SBATCH --cpus-per-task=64
#SBATCH -p <partition>
#SBATCH -N 1
module purge
module load cmake/<version>
module load gcc/<version>
module load intel/<version>
module load mpi/<version>
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate <env_name>
export PATH=/path/to/your/tools/bin:$PATH   # adapt to your cluster

# ============================================================
# GEWUM SRSS Automated Workflow (Single-Shot Full Pipeline)
# Usage: sbatch run_srss.sh
#
# IMPORTANT: Edit the individual scripts BEFORE submitting this:
#   cifgen.sh         
#   run_selection.sh  .
#   relax_umlip.sh    
#   post_relax.sh    
# ============================================================

# --- Ehull parameters ---
MP_DATA="/path/to/MPtrj_2022.9_full.json"   # Path to offline MP JSON file
USE_COMPAT="false"                            # Use MP2020 compatibility corrections (true/false)
SELF_HULL="false"                             # Self-hull mode: build hull from input structures only (true/false)
EHULL_API_KEY=""                              # MP API key for online mode (leave empty to use --mp-data)

echo "============================================"
echo "GEWUM SRSS Automated Workflow"
echo "============================================"
echo "  cifgen:  (see cifgen.sh)"
echo "  select:  (see run_selection.sh)"
echo "  relax:   (see relax_umlip.sh)"
echo "  post:    (see post_relax.sh)"
echo "  Ehull:   mp_data=$MP_DATA, compat=$USE_COMPAT, self-hull=$SELF_HULL"
echo "Start time: $(date)"
echo "============================================"

for script in cifgen.sh run_selection.sh relax_umlip.sh post_relax.sh; do
    if [ ! -f "$script" ]; then
        echo "[ERROR] Required script not found: $script"
        echo "  Please run the following commands first:"
        echo "    gewum RD --mode cifgen"
        echo "    gewum RD --mode select"
        echo "    gewum RD --mode relax"
        echo "    gewum RD --mode post"
        exit 1
    fi
done

if [ -z "$MP_DATA" ] && [ -z "$EHULL_API_KEY" ]; then
    echo "[ERROR] Neither MP_DATA nor EHULL_API_KEY is set."
    echo "  Please edit run_srss.sh and set either:"
    echo "    MP_DATA=\"/path/to/MPtrj_2022.9_full.json\""
    echo "    or"
    echo "    EHULL_API_KEY=\"your-mp-api-key\""
    exit 1
fi

echo "[Step 1/5] Submitting cifgen.sh..."
JOB_CIFGEN=$(sbatch --parsable cifgen.sh)
echo "  cifgen.sh submitted: Job ${JOB_CIFGEN}"

echo "[Step 2/5] Submitting run_selection.sh (after ${JOB_CIFGEN})..."
JOB_SELECT=$(sbatch --parsable --dependency=afterok:${JOB_CIFGEN} run_selection.sh)
echo "  run_selection.sh submitted: Job ${JOB_SELECT}"

echo "[Step 3/5] Submitting relax_umlip.sh (after ${JOB_SELECT})..."
JOB_RELAX=$(sbatch --parsable --dependency=afterok:${JOB_SELECT} relax_umlip.sh)
echo "  relax_umlip.sh submitted: Job ${JOB_RELAX}"

echo "[Step 4/5] Submitting post_relax.sh (after ${JOB_RELAX})..."
JOB_POST=$(sbatch --parsable --dependency=afterok:${JOB_RELAX} post_relax.sh)
echo "  post_relax.sh submitted: Job ${JOB_POST}"

echo "[Step 5/5] Generating Ehull job script..."

EHULL_CMD="gewum RD --mode Ehull -r"
if [ -n "$MP_DATA" ]; then
    EHULL_CMD="$EHULL_CMD --mp-data ${MP_DATA}"
else
    EHULL_CMD="$EHULL_CMD --api-key ${EHULL_API_KEY}"
fi
if [ "$USE_COMPAT" = "true" ]; then
    EHULL_CMD="$EHULL_CMD -cor"
fi
if [ "$SELF_HULL" = "true" ]; then
    EHULL_CMD="$EHULL_CMD --self-hull"
fi
echo "  Ehull command: $EHULL_CMD"

cat > _run_ehull_srss.sh << EHULLEOF
#!/bin/bash
#SBATCH --job-name=SRSS_Ehull
#SBATCH --output=SRSS_Ehull.out
#SBATCH --error=SRSS_Ehull.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=64
#SBATCH -p <partition>
#SBATCH -N 1

module purge
module load cmake/<version>
module load gcc/<version>
module load intel/<version>
module load mpi/<version>
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate <env_name>
export PATH=/path/to/your/tools/bin:$PATH   # adapt to your cluster

echo "============================================"
echo "SRSS Ehull Calculation"
echo "Start time: \$(date)"
echo "============================================"

echo "Running: $EHULL_CMD"
eval $EHULL_CMD
EHULL_STATUS=\$?

if [ \$EHULL_STATUS -ne 0 ]; then
    echo "[ERROR] Ehull calculation failed with exit code \$EHULL_STATUS."
    exit 1
fi

echo "============================================"
echo "SRSS Workflow Complete!"
echo "End time: \$(date)"
echo "Output: Hull_result.csv"
echo "============================================"
EHULLEOF

chmod +x _run_ehull_srss.sh

JOB_EHULL=$(sbatch --parsable --dependency=afterok:${JOB_POST} _run_ehull_srss.sh)
echo "  Ehull submitted: Job ${JOB_EHULL} (after ${JOB_POST})"

echo ""
echo "============================================"
echo "All pipeline jobs submitted successfully!"
echo "  cifgen:       Job ${JOB_CIFGEN}"
echo "  select:       Job ${JOB_SELECT} (after ${JOB_CIFGEN})"
echo "  relax:        Job ${JOB_RELAX} (after ${JOB_SELECT})"
echo "  post:         Job ${JOB_POST} (after ${JOB_RELAX})"
echo "  Ehull:        Job ${JOB_EHULL} (after ${JOB_POST})"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "============================================"
echo "End time: $(date)"
