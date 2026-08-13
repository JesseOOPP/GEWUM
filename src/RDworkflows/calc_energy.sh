#!/bin/bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

# ============================================================
# GEWUM Parallel Energy Calculation Script
# Calculates single-point energies for all CIF files
# Output format matches 0_final_result_tot.txt
# ============================================================

start_time=$(date +%s)

DEVICE=${1:-cpu}         
OUTPUT_FILE=${2:-"0_final_result_tot.txt"}

echo "============================================================"
echo "GEWUM Parallel Energy Calculation"
echo "============================================================"
echo "  DEVICE: $DEVICE"
echo "  OUTPUT_FILE: $OUTPUT_FILE"
echo "============================================================"

{{ENV_SETUP}}

top_dir=$(pwd)
echo "Working directory: $top_dir"

TEMP_DIR="$top_dir/.energy_temp"
mkdir -p "$TEMP_DIR"

run_energy_calc() {
    CIF_FILE=$1
    TEMP_DIR=$2
    CALC_DEVICE=$3
    
    if [ -z "$CIF_FILE" ] || [ -z "$TEMP_DIR" ]; then
        echo "Error: CIF_FILE or TEMP_DIR is empty. Skipping."
        return 1
    fi
    
    base_name=$(basename "$CIF_FILE" .cif)
    output_csv="$TEMP_DIR/${base_name}.csv"
    
    python -m gewum.src.RDworkflows.calc_energy_single "$CIF_FILE" "$output_csv" --device "$CALC_DEVICE" &>> calc_energy.out
}

export -f run_energy_calc

tasks_file="calc_tasks.txt"
> "$tasks_file"

for cif_file in $(ls "$top_dir" | grep '\.cif$'); do
    cif_file_path="$top_dir/$cif_file"
    
    if [ -f "$cif_file_path" ]; then
        echo "$cif_file_path $TEMP_DIR" >> "$tasks_file"
    fi
done

task_count=$(wc -l < "$tasks_file")
echo "Found $task_count CIF files to process"

if [ ! -s "$tasks_file" ]; then
    echo "No CIF files found in current directory."
    rm -rf "$TEMP_DIR"
    exit 1
fi

TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-64}
CORES_PER_TASK=1
NUM_TASKS=$((TOTAL_CPUS / CORES_PER_TASK))

export OMP_NUM_THREADS=$CORES_PER_TASK
export DEVICE

parallel -j $NUM_TASKS -a "$tasks_file" --colsep ' ' run_energy_calc {1} {2} $DEVICE

echo "Collecting results..."

python3 -c "
import os
import csv
import glob

temp_dir = '$TEMP_DIR'
output_file = '$OUTPUT_FILE'

all_data = []
for csv_file in glob.glob(os.path.join(temp_dir, '*.csv')):
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) == 6:
                all_data.append({
                    'formula': row[0],
                    'name': row[1],
                    'total': float(row[2]),
                    'per_atom': float(row[3]),
                    'path': row[4],
                    'sg': row[5]
                })

# Sort by energy per atom
all_data.sort(key=lambda x: x['per_atom'])

# Write to output file
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow([
        'Chemical_Formula',
        'CIF_Base_Name',
        'Total_Energy_eV',
        'Energy_per_Atom_eV',
        'Relaxed_CIF_Path',
        'SG_ori'
    ])
    for d in all_data:
        writer.writerow([
            d['formula'],
            d['name'],
            f\"{d['total']:.6f}\",
            f\"{d['per_atom']:.6f}\",
            d['path'],
            d['sg']
        ])

print(f'Total structures processed: {len(all_data)}')
"

success_count=$(ls "$TEMP_DIR"/*.csv 2>/dev/null | wc -l)

rm -rf "$TEMP_DIR"
rm -f "$tasks_file"

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo "============================================================"
echo "Energy calculation completed!"
echo "  Successful: $success_count / $task_count"
echo "  Results saved to: $OUTPUT_FILE"
echo "  Total runtime: ${hours}h ${minutes}m ${seconds}s"
echo "============================================================" | tee Time_calc_energy.log
