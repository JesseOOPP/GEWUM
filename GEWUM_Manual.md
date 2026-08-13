# GEWUM Manual

**General Exploration Workflow for the Utopia of Materials**

Version 1.0 | Contact: songjx@szlab.ac.cn | License: MIT | Citation: arXiv:2604.21401 (2026)

---

## 1. Quick Start

```bash
gewum <command> --mode <mode> [--dest <dir>]   # scripts are copied to <dir> (default: current directory)

# Example: copy the PT element-substitution template
gewum PT --mode sub

# Modes without copied scripts run directly with -r
gewum RD --mode sym -r
```

## 2. Command Overview

| Command | Purpose | Modes |
| --- | --- | --- |
| `cifgen` | CIF generation input | all, oxidation, substitute, mutate, dp |
| `RD` | Random structure search (SRSS) | cifgen, select, relax, refine, calc, post, reorder, dedup, Ehull, sym, sym2d, chiral, ph, phpost, relaxhp, posthp, viz, viz2, viz3, auto |
| `PT` | Perturbation-based search | supercell, sub, mutate, dp, relax, post, dedup, Ehull, sym, sym2d, chiral, ph, phpost, viz |
| `ELA` | Elastic constants (energy-strain) | pre, cal, post |
| `QHA` | Quasi-harmonic approximation | pre, cal |
| `TC` | Thermal conductivity (RTA) | fc3, post |
| `MD` | Molecular dynamics (NVT) & amorphous ensemble | nvt, post, sq |
| `DB` | Read-only structures.db inspection/export | info, list, stats, show, export |

Use `gewum <command> -h` for per-command options.

---

## 3. RD Workflow (Random Design)

### 3.1 Modes

| Mode | Function | Copied script(s) |
| --- | --- | --- |
| `cifgen` | Generate 0D-3D random crystals | cifgen.sh |
| `select` | Diversity selection (0D-3D) | run_selection.sh |
| `relax` | Structure relaxation with uMLIP | relax_umlip.sh |
| `refine` | Second-stage lattice refinement | refine_umlip.sh |
| `calc` | Single-point energy (no relaxation) | calc_energy.sh |
| `post` | Bond check + energy ranking | post_relax.sh |
| `reorder` | Reorder energies by composition | (direct: `-r`) |
| `dedup` | Remove duplicates (composition + structure) | (direct: `-r`) |
| `Ehull` | Convex-hull stability screening | ehull/*.py (or direct with --api-key/--mp-data) |
| `sym` | Space-group symmetry rename | (direct: `-r`) |
| `sym2d` | 2D layer-group classification | (direct: `-r`) |
| `chiral` | Extract chiral structures | (direct: `-r`) |
| `ph` | Phonon calculation | ph_cal.sh |
| `phpost` | Phonon post-processing (stable/unstable filter) | (direct: `-r`) |
| `relaxhp` | High-pressure relaxation | relax_umlip_hp.sh |
| `posthp` | High-pressure post-processing (enthalpy sort) | post_relax_hp.sh |
| `viz` | UMAP/t-SNE structure visualization | run_visualize.sh |
| `viz2` | Advanced analysis (RDF, SG curves) | run_viz2.sh |
| `viz3` | Formation-energy phase diagram | (direct: `-r`) |
| `auto` | Fully automated pipeline (SLURM-only) | cifgen.sh, relax_umlip.sh, post_relax.sh, run_srss.sh, run_selection.sh |

### 3.2 Typical Workflow

```bash
# 1. Generate candidate structures (edit cifgen.sh: conda env, attempts)
gewum cifgen --mode all
gewum RD --mode cifgen
sbatch cifgen.sh

# 2. Diversity selection (SRSS)
gewum RD --mode select
sbatch run_selection.sh

# 3. Relaxation with uMLIP (edit relax_umlip.sh: mode 1/2, fmax, steps)
gewum RD --mode relax
sbatch relax_umlip.sh

# 4. Post-processing: bond check + energy ranking
gewum RD --mode post
sbatch post_relax.sh

# 5. Single-point energy (optional)
gewum RD --mode calc
sbatch calc_energy.sh

# 6. Ehull screening (online / offline / filter)
gewum RD --mode Ehull --api-key YOUR_MP_KEY          # online
gewum RD --mode Ehull -r --mp-data /path/to/MPtrj.json   # offline
gewum RD --mode Ehull --post -t 0.2                   # extract below threshold

# 7. Lattice refinement of survivors (optional)
cd 0_cif
gewum RD --mode refine
sbatch refine_umlip.sh

# 8. Symmetry classification / dedup
gewum RD --mode sym -r            # 3D space groups
gewum RD --mode sym2d -r          # 2D layer groups
gewum RD --mode chiral -r         # chiral structures
gewum RD --mode dedup -r          # remove duplicates (--rdf merges similar motifs)

# 9. Phonon stability check (edit ph_cal.sh: supercell, band path)
gewum RD --mode ph
sbatch ph_cal.sh
gewum RD --mode phpost -r         # stable -> 0_final/, unstable -> 0_unstable/

# 10. Visualization (optional)
gewum RD --mode viz -r
gewum RD --mode viz2 -r
gewum RD --mode viz3 -r --mp-data /path/to/MP.json

# High-pressure variant
gewum RD --mode relaxhp           # edit relax_umlip_hp.sh (pressure, fmax, steps)
sbatch relax_umlip_hp.sh
gewum RD --mode posthp

# Fully automated pipeline: cifgen -> select -> relax -> post -> Ehull
# (chains jobs via sbatch --dependency; edit run_srss.sh first: MP_DATA or EHULL_API_KEY required)
gewum RD --mode auto
sbatch run_srss.sh
```

---

## 4. PT Workflow (Perturbation)

Generate perturbed / substituted / doped variants from known structures.

| Mode | Function | Copied file(s) |
| --- | --- | --- |
| `supercell` | Build supercells | (direct: `-r`) |
| `sub` | Element-substitution template | replacements.yaml |
| `mutate` | Perturbation parameter file | INPUT |
| `dp` | Doping template | doping.yaml |
| `relax` | Relaxation | relax_umlip_pt.sh |
| `post` | Energy collection | post_relax_pt.sh |
| `dedup` / `sym` / `sym2d` / `chiral` / `phpost` | Same downstream steps as RD | (direct: `-r`) |
| `Ehull` | Convex-hull stability screening | ehull/*.py (or direct with --api-key/--mp-data) |
| `ph` | Phonon calculation | ph_cal.sh |
| `viz` | PT structure visualization | run_viz_pt.sh |

```bash
# Optional: supercells from a template CIF
gewum PT --mode supercell -r -i input.cif --matrix 2 2 2

# Perturbation / substitution / doping: copy template, edit, generate
gewum PT --mode mutate
gewum cifgen --mode mutate
gewum PT --mode sub
gewum cifgen --mode substitute
gewum PT --mode dp
gewum cifgen --mode dp

# Relaxation & post-processing (edit relax_umlip_pt.sh as needed)
gewum PT --mode relax
sbatch relax_umlip_pt.sh
gewum PT --mode post

# Downstream: same as RD (Ehull, sym, ...)
gewum PT --mode viz -r
```

---

## 5. TC Workflow (Thermal Conductivity)

```bash
# FC3 calculation (edit tc_cal.sh: supercell)
gewum TC --mode fc3
sbatch tc_cal.sh

# Thermal conductivity (edit tc_post.sh: q-mesh, temperature range)
gewum TC --mode post
sbatch tc_post.sh
# Results: K.dat in each subdirectory
```

---

## 6. ELA Workflow (Elastic Constants)

```bash
# One working directory per CIF
gewum ELA --mode pre -r

# Calculation (copies cal_ela.sh + VPKIT.in1/in2; edit as needed)
gewum ELA --mode cal
sbatch cal_ela.sh

# Collect results (ela.dat / ela_tot.dat)
gewum ELA --mode post -r
```

---

## 7. QHA Workflow (Quasi-Harmonic Approximation)

```bash
# One working directory per CIF
gewum QHA --mode pre -r

# Volume sweep + phonon + QHA thermal properties (edit cal_qha.sh as needed)
gewum QHA --mode cal
sbatch cal_qha.sh
```

---

## 8. MD Workflow (Molecular Dynamics)

| Mode | Function | Copied script(s) |
| --- | --- | --- |
| `nvt` | NVT MD simulation | run_md_nvt.sh |
| `post` | Post-processing (energy-time plots) | post_md.sh |
| `sq` | Amorphous ensemble (stochastic quenching) | run_sq_ensemble.sh |

```bash
# NVT simulation (edit run_md_nvt.sh: T, steps, supercell)
gewum MD --mode nvt
sbatch run_md_nvt.sh

# Post-processing
gewum MD --mode post
sbatch post_md.sh

# Amorphous ensemble (optional)
gewum MD --mode sq
sbatch run_sq_ensemble.sh
```

---

## 9. DB Command (Read-Only Database Tools)

Inspect and export `structures.db` files. The DB command never writes to them.

| Subcommand | Function |
| --- | --- |
| `info` | Per-DB metadata: row counts by stage, energy range, schema |
| `list` | Filtered rows (table/csv/tsv/json) |
| `stats` | Group-by aggregation (1-2 dimensions: formula/stage/sg) |
| `show` | Dump one row's CIF content |
| `export` | Batch export CIFs to a directory or ZIP |

Common filters: `--formula`, `--sg`, `--name`, `--stage`, `--energy-min`, `--epa-min/max`, `--order-by`, `--limit`.

```bash
gewum DB info --root <work_dir>
gewum DB list --root <work_dir> --stage relaxed --limit 50
gewum DB stats --root <work_dir> --by formula,stage
gewum DB show --db <path>/structures.db --sg 225 --name xtal_1.cif
gewum DB export --root <work_dir> --stage relaxed --out filtered_cifs/
gewum DB export --root <work_dir> --zip all_relaxed.zip
```

---

## 10. Configuration

Place `slurm_config.yaml` in the target or working directory to customize the `#SBATCH` header and environment setup injected into every copied `.sh` script. Search order: target directory -> current directory -> GEWUM installation directory.

```yaml
slurm:
  time: "2400:00:00"        # wall time (HH:MM:SS)
  cpus_per_task: 64
  partition: "<partition>"
  nodes: 1

environment:
  module_purge: true
  modules:                  # modules to load, in order
    - "cmake/<version>"
    - "gcc/<version>"
  conda_path: "/path/to/anaconda3"
  conda_env: "gewum"

# optional: GNU Parallel for parallel batch jobs
parallel:
  path: "/path/to/parallel/bin"
```

Placeholders replaced at copy time:

| Placeholder | Meaning |
| --- | --- |
| `{{JOB_NAME}}` | Job name (derived from the script filename) |
| `{{SLURM_TIME}}` | Wall time |
| `{{SLURM_CPUS}}` | CPUs per task |
| `{{SLURM_PARTITION}}` | Partition |
| `{{SLURM_NODES}}` | Node count |
| `{{CONDA_PATH}}` | Conda installation path |
| `{{CONDA_ENV}}` | Conda environment name |
| `{{SLURM_HEADER}}` | Full generated #SBATCH block |
| `{{ENV_SETUP}}` | Full environment setup block |

---

## 11. Key Input/Output Files

### cifgen.inp (crystal generation input)

```python
# [elements], [stoichiometry], count-per-group, [optional: max atoms]
["Te", "O"], [1, 2], 5
["Cu", "Sn", "S"], [2, 1, 3], 3, 48
```

### replacements.yaml (PT element substitution)

```yaml
- old: A
  new: B
- old: C
  new: D
```

### doping.yaml (PT doping)

```yaml
dopant: Li          # dopant element
host: Na            # element to replace
concentration: 0.1  # doping fraction
num_structures: 20  # number of structures to generate
```

### MP offline data (Ehull --mp-data)

Materials Project trajectory (MPtrj) JSON, downloadable from the Materials Project website. Format: `{mp_id: {entry_id: {energy_per_atom, corrected_total_energy, uncorrected_total_energy, structure}}}`. First use builds a `<name>_index.json` index next to the file (takes a few minutes; needs `pip install ijson`); later runs load the index instantly.

### Energy results

`energy_final.txt` (per composition) and `0_final_result_tot.txt` (global) are CSV files with: `Chemical_Formula, CIF_Base_Name, Total_Energy_eV, Energy_per_Atom_eV, Relaxed_CIF_Path, SG_ori`.

High-pressure output (`posthp`) adds: `Final_Pressure_GPa, Enthalpy_per_Atom_eV, Corrected_Enthalpy_per_Atom_eV`.

---

## 12. FAQ

**Q: How do I get a Materials Project API key?**
A: Register at https://materialsproject.org/ and copy the key from your profile settings.

**Q: Which relaxation mode should I use?**
A: Mode 1 optimizes atomic positions only (fast, high-throughput screening); mode 2 also optimizes the lattice (final refinement). High-pressure relaxation uses the `relaxhp` mode instead.

**Q: What does Comp_CompatibilityError mean?**
A: The pyxtal Wyckoff positions cannot be balanced for a given composition/space-group combination; the script skips such assignments automatically.

**Q: How do the PT modes supercell / sub / mutate / dp differ?**
A: `supercell` expands the cell; `sub` performs element substitution; `mutate` applies random strain/atomic displacements; `dp` enumerates doping configurations by symmetry (no interstitial doping).

---

## 13. Contact

- Developer: Jiexi Song (songjx@szlab.ac.cn)
- Contributors: Diwei Shi, Yanqing Qin, Aixian She, Zhenyu Liu
- Acknowledgments: Fengyuan Xuan, Chongde Cao
