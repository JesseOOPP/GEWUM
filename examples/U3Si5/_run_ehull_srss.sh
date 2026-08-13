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
echo "Start time: $(date)"
echo "============================================"

echo "Running: gewum RD --mode Ehull -r --mp-data /path/to/MPtrj_2022.9_full.json"
eval gewum RD --mode Ehull -r --mp-data /path/to/MPtrj_2022.9_full.json
EHULL_STATUS=$?

if [ $EHULL_STATUS -ne 0 ]; then
    echo "[ERROR] Ehull calculation failed with exit code $EHULL_STATUS."
    exit 1
fi

echo "============================================"
echo "SRSS Workflow Complete!"
echo "End time: $(date)"
echo "Output: Hull_result.csv"
echo "============================================"
