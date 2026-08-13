#!/bin/bash

#SBATCH --job-name=cal_qha
#SBATCH --output=cal_qha.out
#SBATCH --error=cal_qha.err
#SBATCH --time=2400:00:00
#SBATCH --cpus-per-task=64
#SBATCH -p <partition>
#SBATCH -N 1

unset DISPLAY
start_time=$(date +%s)

module purge
module load cmake/<version>
module load gcc/<version>
module load intel/<version>
module load mpi/<version>
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate <env_name>
export PATH=/path/to/parallel/build/bin:$PATH

main_dir=$(pwd)
scales=(0.95 0.96 0.97 0.98 0.99 1.00 1.01 1.02 1.03 1.04 1.05)

echo "============================================================"
echo "GEWUM QHA - Quasi-Harmonic Approximation Calculation"
echo "============================================================"

process_qha() {
    local subdir="$1"
    local main_dir="$2"
    
    echo ""
    echo "Processing: $subdir"
    echo "------------------------------------------------------------"
    
    cd "$subdir" || return 1
    
    python -m gewum.src.QHAworkflows.cif_pri
    mv primitive_POSCAR POSCAR
    echo -e "401\n1\n2 2 2" | vaspkit
    rm POSCAR
    mv SC222.vasp POSCAR
    
    for scale in "${scales[@]}"; do
        mkdir -p $scale
        cp POSCAR $scale/
    done
    
    echo "  Starting volume-scaling relaxations..."
    for scale in "${scales[@]}"; do
        (
            cd $scale
            python -m gewum.src.QHAworkflows.relax_volume $scale
        ) &
    done
    wait
    
    mkdir -p CP
    python -m gewum.src.QHAworkflows.collect_energy
    
    echo "  Starting phonon calculations..."
    for scale in "${scales[@]}"; do
        (
            cd $scale
            python -m gewum.src.QHAworkflows.run_phonon
        ) &
    done
    wait
    
    cd CP
    phonopy-qha v-e.dat thermal_properties-{1..11}.yaml --tmax 2002 --pressure 0 > thermo.out
    
    echo "  Completed: $subdir"
    cd "$main_dir" || return 1
}

for subdir in $(find . -maxdepth 1 -type d ! -path .); do
    process_qha "$subdir" "$main_dir"
done

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))

echo ""
echo "============================================================"
echo "QHA calculation completed!"
echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"
echo "============================================================"
