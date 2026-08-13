#!/bin/bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --output={{JOB_NAME}}.out
#SBATCH --error={{JOB_NAME}}.err
#SBATCH --time={{SLURM_TIME}}
#SBATCH --cpus-per-task={{SLURM_CPUS}}
#SBATCH -p {{SLURM_PARTITION}}
#SBATCH -N {{SLURM_NODES}}

start_time=$(date +%s)

# ============ Stochastic Quenching Parameters ============
N_CONFIGS=${1:-20}          # Random configs per worker
SUPERCELL_X=${2:-2}
SUPERCELL_Y=${3:-2}
SUPERCELL_Z=${4:-2}
FMAX=${5:-0.05}
D_MIN_SCALE=${6:-0.7}
VARIABLE_CELL=${7:-1}       # 1 = variable cell (zero pressure), 0 = fixed cell

# ============ Ensemble Parameters ============
N_WORKERS=${8:-8}           # Parallel SQ workers (independent seed ranges)
N_SELECT=${10:-100}         # Final amorphous count kept after maxmin selection
SAVE_INITIAL=1              # 1 = also save pre-relaxation initial configs to
                            #     each worker's init_frames/ (paired by filename
                            #     with am_frames/); 0 = only keep relaxed frames

# ============ Perturbation-route Parameters (edit here) ============
# Two-route ensemble: the first N_PERTURB_WORKERS workers build their initial
# configurations by PERTURBING the crystal template (reusing the PT-module
# operators: rattle + optional lattice strain + optional atom rotation); the
# remaining workers use pure RANDOM packing. All frames merge into one pool.
# Set N_PERTURB_WORKERS=0 (default) for a pure-random ensemble.
N_PERTURB_WORKERS=${9:-0}   # Number of workers (of N_WORKERS) using perturb route
RATTLE_STDEV_MIN=0.4         # Min rattle displacement stdev (A); larger = more amorphous
RATTLE_STDEV_MAX=1.0         # Max rattle displacement stdev (A)
STRAIN_PROB=0.5              # Probability of applying random lattice strain
MAX_STRAIN=0.05              # Max lattice strain magnitude (0-1)
ROTATION_PROB=0.3            # Probability of applying atom rotation
ROTATION_PER_ATOM=0.15       # Per-atom rotation probability
MAX_ROTATION_ANGLE=60.0      # Max rotation angle (degrees)
PERTURB_MIN_DIST=1.0         # Min pairwise distance after perturbation (A)

# ==== Advanced pipeline (edit here) ====
# The SQ ensemble runs a six-stage pipeline: coarse relax (cheap) -> merge ->
# crystallinity gate (Q6 + SOAP) -> maxmin diversity selection -> final LBFGS
# refinement of the selected subset -> post-refinement crystallinity QC. This
# spends expensive relaxation only on confirmed-amorphous, diverse structures.
COARSE_FIRE=60               # FIRE steps for the cheap first-pass relaxation
FINAL_FIRE=50                # FIRE steps before the final LBFGS refinement
FINAL_LBFGS=200              # LBFGS steps for the final refinement
SELECT_METHOD=maxmin         # Diversity selection: maxmin (farthest-point)
SELECT_DESCRIPTOR=soap       # Structure descriptor for selection: soap
SOAP_R_CUT=6.0               # SOAP cutoff radius (A)
SOAP_N_MAX=4                 # SOAP radial basis size
SOAP_L_MAX=4                 # SOAP angular basis size
SELECT_WORKERS=$N_WORKERS    # Parallel workers for featurization/selection
REFINE_WORKERS=$N_WORKERS    # Parallel shards for the final refinement (Stage 5)
# NOTE: keep REFINE_WORKERS fixed across resumes. Changing the shard count after
#       some frames are already refined can drop their rows from the merged
#       am_ml_energy.csv (the frames themselves are preserved on disk).
Q6_REF_FRAC=0.5              # Crystalline if Q6 >= frac*Q6_template
Q6_FRAC_THRESHOLD=0.5        # Crystalline if crystalline-atom fraction >= this
SOAP_SIM_THRESHOLD=0.05      # Crystalline if SOAP cosine distance to template <= this
Q6_CUTOFF=                   # Neighbour cutoff (A); empty = auto (first RDF min)

echo "=============================================="
echo "GEWUM SQ Ensemble: ${N_WORKERS}-worker Stochastic Quenching"
echo "=============================================="
echo "Configs/worker: $N_CONFIGS"
echo "Supercell:      ${SUPERCELL_X}x${SUPERCELL_Y}x${SUPERCELL_Z}"
echo "Fmax:           $FMAX eV/A"
echo "D_min scale:    $D_MIN_SCALE"
echo "Variable cell:  $VARIABLE_CELL"
echo "N workers:      $N_WORKERS"
echo "N select:       $N_SELECT (final amorphous count)"
echo "Perturb workers: $N_PERTURB_WORKERS / $N_WORKERS"
if [ "$N_PERTURB_WORKERS" -gt 0 ]; then
    echo "  rattle stdev:  [$RATTLE_STDEV_MIN, $RATTLE_STDEV_MAX] A"
    echo "  strain prob:   $STRAIN_PROB (max $MAX_STRAIN)"
    echo "  rotation prob: $ROTATION_PROB (per-atom $ROTATION_PER_ATOM, max ${MAX_ROTATION_ANGLE} deg)"
    echo "  perturb d_min: $PERTURB_MIN_DIST A"
fi
echo "=============================================="

{{ENV_SETUP}}

# Resolve fixed-cell flag once for all workers
CELL_FLAG=""
if [ "$VARIABLE_CELL" -eq 0 ]; then
    CELL_FLAG="--fixed-cell"
fi

# Resolve save-initial flag once for all workers
SAVE_INIT_FLAG=""
if [ "$SAVE_INITIAL" -eq 1 ]; then
    SAVE_INIT_FLAG="--save-initial"
fi

cif_files=(*.cif)

if [ ${#cif_files[@]} -eq 0 ] || [ "${cif_files[0]}" == "*.cif" ]; then
    echo "Error: No .cif files found in current directory!"
    exit 1
fi

echo "Found ${#cif_files[@]} CIF file(s)"

for cif_file in "${cif_files[@]}"; do
    base_name="${cif_file%.cif}"
    ensemble_dir="sq_ensemble_${base_name}"
    merged_frames="${ensemble_dir}/am_frames"
    merged_energy="${ensemble_dir}/am_ml_energy.csv"

    echo ""
    echo "=============================================="
    echo "Processing: $cif_file -> $ensemble_dir"
    echo "=============================================="

    # ================================================================
    # Stage 1: Run N_WORKERS independent SQ workers (parallel + resume)
    # ================================================================
    echo ""
    echo "===== Stage 1: ${N_WORKERS} independent SQ workers (parallel) ====="

    worker_success=0
    declare -A pids

    for worker_id in $(seq 1 $N_WORKERS); do
        worker_label=$(printf "run_%02d" $worker_id)
        output_dir="${ensemble_dir}/${worker_label}"

        # Resume detection: skip if this worker already completed
        if [ -f "${output_dir}/am_ml_energy.csv" ] && [ -d "${output_dir}/am_frames" ]; then
            n_existing=$(ls "${output_dir}/am_frames"/*.cif 2>/dev/null | wc -l)
            if [ "$n_existing" -gt 0 ]; then
                echo "  ${worker_label}: already completed (${n_existing} frames), skipping"
                worker_success=$((worker_success + 1))
                continue
            fi
        fi

        # Non-overlapping seed range per worker
        seed_base=$((worker_id * 100000))

        # Route assignment: the first N_PERTURB_WORKERS workers perturb the
        # template; the remaining workers use pure random packing.
        if [ "$worker_id" -le "$N_PERTURB_WORKERS" ]; then
            route_args="--init-mode perturb --rattle-stdev-min $RATTLE_STDEV_MIN --rattle-stdev-max $RATTLE_STDEV_MAX --strain-prob $STRAIN_PROB --max-strain $MAX_STRAIN --rotation-prob $ROTATION_PROB --rotation-per-atom $ROTATION_PER_ATOM --max-rotation-angle $MAX_ROTATION_ANGLE --perturb-min-dist $PERTURB_MIN_DIST"
            echo "  ${worker_label}: launching PERTURB route (seed_base=${seed_base})..."
        else
            route_args="--init-mode random"
            echo "  ${worker_label}: launching RANDOM route (seed_base=${seed_base})..."
        fi

        python -m gewum.src.MDworkflows.sq_am "$cif_file" \
            --output-dir "$output_dir" \
            --n-configs $N_CONFIGS \
            --supercell $SUPERCELL_X $SUPERCELL_Y $SUPERCELL_Z \
            --fmax $FMAX \
            --d-min-scale $D_MIN_SCALE \
            --seed $seed_base \
            --steps-fire $COARSE_FIRE --steps-final 0 \
            $CELL_FLAG $SAVE_INIT_FLAG $route_args &
        pids[$worker_id]=$!
    done

    # Wait for all background jobs and tally results
    for worker_id in "${!pids[@]}"; do
        worker_label=$(printf "run_%02d" $worker_id)
        wait ${pids[$worker_id]}
        if [ $? -eq 0 ]; then
            echo "  ${worker_label}: completed"
            worker_success=$((worker_success + 1))
        else
            echo "  ${worker_label}: FAILED"
        fi
    done

    echo ""
    echo "Workers completed: ${worker_success}/${N_WORKERS}"

    if [ "$worker_success" -eq 0 ]; then
        echo "Error: All SQ workers failed for $cif_file, skipping."
        continue
    fi

    # ================================================================
    # Stages 2-6: advanced amorphous pipeline
    #   merge -> gate -> maxmin select -> refine -> QC
    # Resume safety:
    #   * am_final/am_order.csv present  -> whole pipeline done, skip 2-6
    #   * merged_frames/selection_log.txt present -> a selection was already
    #     committed (maxmin is non-deterministic), so DO NOT re-merge / re-gate
    #     / re-select; the selected subset already sits at the top level. Jump
    #     straight to refine (which itself resumes per-frame).
    # ================================================================
    am_final_dir="${ensemble_dir}/am_final"
    am_final_frames="${am_final_dir}/am_frames"
    selection_marker="${merged_frames}/selection_log.txt"

    Q6_CUTOFF_FLAG=""
    if [ -n "$Q6_CUTOFF" ]; then
        Q6_CUTOFF_FLAG="--q6-cutoff $Q6_CUTOFF"
    fi

    total_frames="(resumed)"

    if [ -f "${am_final_dir}/am_order.csv" ]; then
        echo ""
        echo "Advanced pipeline already completed for ${base_name} (am_final/am_order.csv exists), skipping stages 2-6."
    else
        if [ ! -f "$selection_marker" ]; then
            # ------------------------------------------------------------
            # Stage 2: Merge frames and energy CSVs
            # ------------------------------------------------------------
            echo ""
            echo "===== Stage 2: Merging ${worker_success} workers ====="

            mkdir -p "$merged_frames"

            total_frames=0
            first_csv=1

            for worker_id in $(seq 1 $N_WORKERS); do
                worker_label=$(printf "run_%02d" $worker_id)
                run_dir="${ensemble_dir}/${worker_label}"
                src_frames="${run_dir}/am_frames"
                src_energy="${run_dir}/am_ml_energy.csv"

                if [ -d "$src_frames" ]; then
                    for frame_cif in "$src_frames"/*.cif; do
                        if [ -f "$frame_cif" ]; then
                            frame_name=$(basename "$frame_cif")
                            cp "$frame_cif" "${merged_frames}/${worker_label}_${frame_name}"
                            total_frames=$((total_frames + 1))
                        fi
                    done
                fi

                if [ -f "$src_energy" ]; then
                    if [ "$first_csv" -eq 1 ]; then
                        head -n 1 "$src_energy" > "$merged_energy"
                        first_csv=0
                    fi
                    tail -n +2 "$src_energy" | awk -F, -v prefix="${worker_label}_" '{
                        $1 = prefix $1
                        print
                    }' OFS=, >> "$merged_energy"
                fi
            done

            echo "Total merged frames: $total_frames"
            echo "Merged energy CSV: $merged_energy"

            if [ "$total_frames" -eq 0 ]; then
                echo "Warning: No frames collected."
                continue
            fi

            # ------------------------------------------------------------
            # Stage 3: Crystallinity gate (pre-selection) - Q6 + SOAP
            # ------------------------------------------------------------
            echo ""
            echo "===== Stage 3: Crystallinity gate (prefilter) ====="
            python -m gewum.src.MDworkflows.am_order "$merged_frames" \
                --template "$cif_file" \
                --supercell $SUPERCELL_X $SUPERCELL_Y $SUPERCELL_Z \
                --q6-ref-frac $Q6_REF_FRAC --q6-frac-threshold $Q6_FRAC_THRESHOLD $Q6_CUTOFF_FLAG \
                --soap-threshold $SOAP_SIM_THRESHOLD \
                --soap-r-cut $SOAP_R_CUT --soap-n-max $SOAP_N_MAX --soap-l-max $SOAP_L_MAX \
                --mode prefilter --move-subdir crystalline

            n_clean=$(ls "${merged_frames}"/*.cif 2>/dev/null | wc -l)
            echo "Amorphous pool after gate: $n_clean frames"
            if [ "$n_clean" -eq 0 ]; then
                echo "Warning: crystallinity gate removed all frames; skipping $cif_file."
                continue
            fi

            # ------------------------------------------------------------
            # Stage 4: maxmin diversity selection (reduce to N_SELECT)
            # ------------------------------------------------------------
            echo ""
            echo "===== Stage 4: ${SELECT_METHOD} selection -> ${N_SELECT} ====="
            python -m gewum.src.common.selection.structure_select \
                --dim 3d --method $SELECT_METHOD --target $N_SELECT \
                --single-dir "$merged_frames" \
                --descriptor $SELECT_DESCRIPTOR \
                --soap-r-cut $SOAP_R_CUT --soap-n-max $SOAP_N_MAX --soap-l-max $SOAP_L_MAX \
                --workers $SELECT_WORKERS

            n_selected=$(ls "${merged_frames}"/*.cif 2>/dev/null | wc -l)
            echo "Selected amorphous subset: $n_selected frames"
        else
            echo ""
            echo "Resume: selection already committed (selection_log.txt found); skipping merge/gate/select for ${base_name}."
        fi

        # ------------------------------------------------------------
        # Stage 5: Final LBFGS refinement of the selected subset
        # Sharded across REFINE_WORKERS processes; all shards write frames into
        # the same am_final/am_frames (distinct names) but separate CSVs, which
        # are merged into am_final/am_ml_energy.csv afterwards.
        # ------------------------------------------------------------
        echo ""
        echo "===== Stage 5: Final refinement (${REFINE_WORKERS} shards) ====="
        declare -A refine_pids
        for refine_id in $(seq 1 $REFINE_WORKERS); do
            shard_idx=$((refine_id - 1))
            shard_csv=$(printf "am_ml_energy_w%02d.csv" $refine_id)
            refine_log=$(printf "${ensemble_dir}/refine_w%02d.log" $refine_id)
            python -m gewum.src.MDworkflows.sq_refine "$merged_frames" \
                --output-dir "$am_final_dir" \
                --fmax $FMAX \
                --steps-fire $FINAL_FIRE --steps-final $FINAL_LBFGS \
                --shard-index $shard_idx --shard-count $REFINE_WORKERS \
                --energy-csv-name "$shard_csv" \
                $CELL_FLAG > "$refine_log" 2>&1 &
            refine_pids[$refine_id]=$!
        done

        refine_ok=0
        for refine_id in "${!refine_pids[@]}"; do
            wait ${refine_pids[$refine_id]}
            if [ $? -eq 0 ]; then
                refine_ok=$((refine_ok + 1))
            else
                echo "  refine shard ${refine_id}: FAILED (see refine_w*.log)"
            fi
        done
        echo "Refinement shards completed: ${refine_ok}/${REFINE_WORKERS}"

        if [ ! -d "$am_final_frames" ]; then
            echo "Warning: refinement produced no frames; skipping QC for $cif_file."
            continue
        fi

        # Merge per-shard energy CSVs into the canonical am_ml_energy.csv
        merged_final_csv="${am_final_dir}/am_ml_energy.csv"
        first_final=1
        for refine_id in $(seq 1 $REFINE_WORKERS); do
            shard_csv=$(printf "${am_final_dir}/am_ml_energy_w%02d.csv" $refine_id)
            [ -f "$shard_csv" ] || continue
            if [ "$first_final" -eq 1 ]; then
                head -n 1 "$shard_csv" > "$merged_final_csv"
                first_final=0
            fi
            tail -n +2 "$shard_csv" >> "$merged_final_csv"
        done

        # Sanity check (detects a changed REFINE_WORKERS across a resume): every
        # refined cif, including any QC-moved into am_frames/crystalline/, should
        # have exactly one energy row.
        if [ -f "$merged_final_csv" ]; then
            refined_count=$(find "$am_final_frames" -name '*.cif' 2>/dev/null | wc -l)
            csv_rows=$(( $(wc -l < "$merged_final_csv") - 1 ))
            if [ "$refined_count" -ne "$csv_rows" ]; then
                echo "Warning: refined frames ($refined_count) != energy rows ($csv_rows) in am_ml_energy.csv."
                echo "         Energy table may be incomplete (did REFINE_WORKERS change across a resume?)."
            fi
        fi

        # ------------------------------------------------------------
        # Stage 6: Crystallinity QC (post-refinement)
        # ------------------------------------------------------------
        echo ""
        echo "===== Stage 6: Crystallinity QC (post-refinement) ====="
        python -m gewum.src.MDworkflows.am_order "$am_final_frames" \
            --template "$cif_file" \
            --supercell $SUPERCELL_X $SUPERCELL_Y $SUPERCELL_Z \
            --q6-ref-frac $Q6_REF_FRAC --q6-frac-threshold $Q6_FRAC_THRESHOLD $Q6_CUTOFF_FLAG \
            --soap-threshold $SOAP_SIM_THRESHOLD \
            --soap-r-cut $SOAP_R_CUT --soap-n-max $SOAP_N_MAX --soap-l-max $SOAP_L_MAX \
            --mode qc --move-subdir crystalline
    fi

    n_final=$(ls "${am_final_frames}"/*.cif 2>/dev/null | wc -l)

    echo ""
    echo "=============================================="
    echo "SQ ensemble completed: $ensemble_dir"
    echo "  Workers:        ${worker_success}/${N_WORKERS}"
    echo "  Merged frames:  $total_frames"
    echo "  Final amorphous: $n_final -> $am_final_frames"
    echo "  Final energy:   ${am_final_dir}/am_ml_energy.csv"
    echo "  Order table:    ${am_final_dir}/am_order.csv"
    echo "=============================================="

done

echo ""
echo "=============================================="
echo "SQ Ensemble workflow completed"
echo "=============================================="

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
hours=$((elapsed_time / 3600))
minutes=$(((elapsed_time % 3600) / 60))
seconds=$((elapsed_time % 60))
echo "Total runtime: $hours hours, $minutes minutes, $seconds seconds" | tee Time_sq_ensemble.log
