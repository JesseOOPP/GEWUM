#!/bin/bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

# example:    sbatch cifgen.sh <dim> [max-atoms] [max-attempts]
#   sbatch cifgen.sh 0 36 30      # 0D, max-atoms=36,  max-attempts=30
#   sbatch cifgen.sh 3 24 150     # 3D, max-atoms=24,  max-attempts=150

DIM=${1:-3}
MAX_ATOMS=${2:-24}
MAX_ATTEMPTS=${3:-30}

if [[ ! "$DIM" =~ ^[0123]$ ]]; then
    echo "Error: dim must be 0, 1, 2, or 3"
    exit 1
fi

start_time=$(date +%s)

{{ENV_SETUP}}

echo "Generating ${DIM}D crystals with max-atoms=${MAX_ATOMS}, max-attempts=${MAX_ATTEMPTS}"

python -m gewum.src.RDworkflows.cif_generate --dim $DIM --max-atoms $MAX_ATOMS --max-attempts $MAX_ATTEMPTS

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
