#!/bin/bash

#SBATCH --job-name=tc_post
#SBATCH --output=tc_post_%j.out
#SBATCH --error=tc_post_%j.err
#SBATCH --time=2400:00:00
#SBATCH --cpus-per-task=64
#SBATCH -p <partition>
#SBATCH -N 1

# ============================================================
# GEWUM Thermal Conductivity - Post-processing
# Computes thermal conductivity from FC3 using phono3py
# ============================================================

unset DISPLAY
start_time=$(date +%s)

# Configurable parameters
MESH=${1:-"5 5 5"}
TMIN=${2:-100}
TMAX=${3:-2000}
TSTEP=${4:-20}

echo "============================================================"
echo "GEWUM Thermal Conductivity - Post-processing"
echo "============================================================"
echo "  Mesh: $MESH"
echo "  Temperature range: $TMIN - $TMAX K (step: $TSTEP)"
echo "============================================================"

module purge
module load cmake/<version>
module load gcc/<version>
module load intel/<version>
module load mpi/<version>
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate <env_name>
export PATH=/path/to/parallel/build/bin:$PATH

export PARALLEL="--will-cite --no-notice"
export OMP_NUM_THREADS=64
export SLURM_CPUS_ON_NODE=64

get_real_cores() {
    if [ -n "$SLURM_CPUS_PER_TASK" ]; then
        echo $SLURM_CPUS_PER_TASK
    else
        echo 64
    fi
}

TOTAL_CORES=$(get_real_cores)
RESERVED_CORES=4
USABLE_CORES=$((TOTAL_CORES - RESERVED_CORES))

DIR_COUNT=$(find . -mindepth 1 -maxdepth 1 -type d -not -exec test -f {}/K.dat \; -print | wc -l)

if [ $DIR_COUNT -eq 0 ]; then
    echo "All directories already processed!"
    exit 0
elif [ $DIR_COUNT -le $((USABLE_CORES/2)) ]; then
    CORES_PER_TASK=$((USABLE_CORES / DIR_COUNT))
else
    CORES_PER_TASK=1
fi

MAX_JOBS=$((USABLE_CORES / CORES_PER_TASK))

echo "Core Allocation:"
echo "  Total cores: $TOTAL_CORES"
echo "  Directories to process: $DIR_COUNT"
echo "  Cores per task: $CORES_PER_TASK"
echo "  Parallel jobs: $MAX_JOBS"
echo "============================================================"

find . -mindepth 1 -maxdepth 1 -type d -not -exec test -f {}/K.dat \; -print0 | \
parallel -0 -j $MAX_JOBS --jobs $MAX_JOBS --joblog parallel_joblog.txt '
    echo "Processing directory: {}"
    
    cd {} || exit 1
    
    # Step 1: Force constant symmetrization
    OMP_NUM_THREADS=$CORES_PER_TASK MKL_NUM_THREADS=$CORES_PER_TASK \
    phono3py --fc-symmetry 2>&1 | tee -a phono3py.log
    
    # Step 2: Compute thermal conductivity
    OMP_NUM_THREADS=$CORES_PER_TASK MKL_NUM_THREADS=$CORES_PER_TASK \
    phono3py --fc3 --fc2 --dim="1 1 1" --mesh="'"$MESH"'" --br \
        --tmin '"$TMIN"' --tmax '"$TMAX"' --tstep '"$TSTEP"' > K.dat 2>&1

    # Cleanup HDF5 files after successful processing
    if [ -f K.dat ]; then
        echo "Cleaning up HDF5 files in: {}"
        find . -maxdepth 1 -type f -name "*.hdf5" -delete
    fi

    echo "Completed: {}"
'

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo "============================================================"
echo "Thermal conductivity calculation completed!"
echo "Results saved in: <structure>/K.dat"
echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"
echo "============================================================" | tee Time_tc_post.log
