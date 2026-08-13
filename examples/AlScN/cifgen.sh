#!/bin/bash

#SBATCH --job-name=cifgen
#SBATCH --output=cifgen.out
#SBATCH --error=cifgen.err
#SBATCH --time=2400:00:00
#SBATCH --cpus-per-task=64
#SBATCH -p <partition>
#SBATCH -N 1

# example:    sbatch cifgen.sh <dim> [max-atoms] [max-attempts]
#   sbatch cifgen.sh 0 36 30      # 0D, max-atoms=36,  max-attempts=30
#   sbatch cifgen.sh 3 24 150     # 3D, max-atoms=24,  max-attempts=150

DIM=${1:-3}
MAX_ATOMS=${2:-64}
MAX_ATTEMPTS=${3:-300}

if [[ ! "$DIM" =~ ^[0123]$ ]]; then
    echo "Error: dim must be 0, 1, 2, or 3"
    exit 1
fi

start_time=$(date +%s)

module purge
module load cmake/<version>
module load gcc/<version>
module load intel/<version>
module load mpi/<version>
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate <env_name>
export PATH=/path/to/your/tools/bin:$PATH   # adapt to your cluster

echo "Generating ${DIM}D crystals with max-atoms=${MAX_ATOMS}, max-attempts=${MAX_ATTEMPTS}"

python -m gewum.src.RDworkflows.cif_generate --dim $DIM --max-atoms $MAX_ATOMS --max-attempts $MAX_ATTEMPTS

# Count CIF files from SQLite databases
python -c "import sqlite3,glob; total=0
for db in glob.glob('*/structures.db'):
    conn=sqlite3.connect(db)
    total+=conn.execute('SELECT COUNT(*) FROM structures').fetchone()[0]
    conn.close()
print(total)" >> TOT_cif

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Total runtime: $hours hours, $minutes minutes, $seconds seconds" > Time_cif_generate.log
