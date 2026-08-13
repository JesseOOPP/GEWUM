#!/bin/bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

unset DISPLAY
start_time=$(date +%s)

{{ENV_SETUP}}

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
    python -c "from ase.io import read,write; a=read('primitive_POSCAR'); write('POSCAR',a.repeat((2,2,2)),format='vasp')"
    rm primitive_POSCAR
    
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
