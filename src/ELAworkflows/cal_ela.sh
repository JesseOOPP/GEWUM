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

process_structure_strains() {
    local subdir_path="$1"

    cd "$subdir_path" || exit 1
    echo "Processing strains in: $(pwd) on $(hostname)"

    # Load the MatterSim model once and relax every strain_* structure in this directory
    python -m gewum.src.ELAworkflows.uMLIP_ela_sq --batch . > umlip_ela.log 2>&1

    echo "Completed: $subdir_path"
}

export -f process_structure_strains

for subdir in $(find . -maxdepth 1 -type d ! -path .); do
    if [ "$subdir" = "." ]; then
        continue
    fi
    
    cd "$subdir" || exit
    
    cp "${main_dir}/VPKIT.in1" .
    cp "${main_dir}/VPKIT.in2" .
    python -m gewum.src.ELAworkflows.cif_pri
    mv primitive_POSCAR POSCAR
    cp VPKIT.in1 VPKIT.in
    vaspkit -task 201
    
    cd "$main_dir" || exit
done

echo "Starting parallel processing of structures ..."

find . -maxdepth 1 -type d ! -path . | parallel -j 64  --joblog parallel_joblog.txt \
    --retries 2 --delay 1 "process_structure_strains {}"

if [ $? -ne 0 ]; then
    echo "Warning: Some parallel tasks may have failed. Check parallel_joblog.txt"
fi

for subdir in $(find . -maxdepth 1 -type d ! -path .); do
    if [ "$subdir" = "." ]; then
        continue
    fi
    
    cd "$subdir" || exit
    
    root_path=$(pwd)
    cd "$root_path" || exit 
    rm -f VPKIT.in
    cp VPKIT.in2 VPKIT.in
    vaspkit -task 201 > ela.dat 2>&1

    cd "$main_dir" || exit
done

echo "All done"

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Total runtime: $hours hours, $minutes minutes, $seconds seconds" > Time_ela.log
