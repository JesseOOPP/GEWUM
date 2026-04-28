[README.md](https://github.com/user-attachments/files/27154806/README.md)
# GEWUM: General Exploration Workflow for the Utopia of Materials

**Version:** 1.0.0  
**Authors:** Jiexi Song, Aixian She, Changpeng Song, Diwei Shi, Fengyuan Xuan, Chongde Cao  
**Contact:** songjiexi@mail.nwpu.edu.cn  
**arXiv:** [2604.21401](https://doi.org/10.48550/arXiv.2604.21401) [cond-mat.mtrl-sci]  
**GitHub:** [https://github.com/JesseOOPP/GEWUM](https://github.com/JesseOOPP/GEWUM)
**License:** MIT

GEWUM is an integrated computational workflow system for materials science that automates crystal structure generation, structure relaxation via universal machine-learning interatomic potentials (uMLIP), thermodynamic property calculations, and high-throughput structure screening. It orchestrates the complete pipeline from random/perturbation-based structure exploration to phonon, elastic, quasi-harmonic, thermal conductivity, and molecular dynamics analyses.

---

## Directory Structure

```
gewum-3-31-3/
├── main.py                     # CLI entry point
├── config.py                   # Centralized path and configuration management
├── pyproject.toml              # Package metadata, dependencies, build config
├── setup.py                    # Legacy setuptools support
├── requirements.txt            # Dependency listing
├── slurm_config.yaml           # SLURM job scheduler configuration template
├── MANIFEST.in                 # Package data inclusion rules
│
├── commands/                   # CLI command modules (one per workflow)
│   ├── __init__.py
│   ├── base.py                 # Abstract base class for workflow commands
│   ├── template_utils.py       # SLURM template rendering utilities
│   ├── cifgen.py               # CIF generation command
│   ├── RD.py                   # Random Design workflow command
│   ├── PT.py                   # Perturbation workflow command
│   ├── ELA.py                  # Elastic constants workflow command
│   ├── QHA.py                  # Quasi-Harmonic Approximation command
│   ├── TC.py                   # Thermal Conductivity command
│   ├── MD.py                   # Molecular Dynamics command
│   └── FT.py                   # Fine-Tuning / active learning command
│
├── src/                        # Core computational modules
│   ├── cifgen_input/           # CIF generation engines
│   │   ├── chemi.py            # Chemical composition enumeration
│   │   ├── conv_all_cifgen.py  # Full composition CIF conversion
│   │   ├── conv_oxidation_cifgen.py  # Oxidation-state-filtered CIF conversion
│   │   ├── oxidation.py        # Oxidation state assignment
│   │   ├── substitute.py       # Element substitution generator
│   │   ├── mutate.py           # Structure perturbation from template CIF
│   │   └── doping.py           # Random doping structure generator
│   │
│   ├── RDworkflows/            # Random Design workflow scripts
│   │   ├── cif_generate.py     # Random crystal structure generation (via PyXtaL)
│   │   ├── cifgen.sh           # SLURM script for CIF generation
│   │   ├── relax_umlip.sh      # SLURM script for uMLIP relaxation
│   │   ├── refine_umlip.sh     # SLURM script for cell refinement
│   │   ├── post_relax.sh       # Post-relaxation processing script
│   │   ├── relax_umlip_hp.sh   # High-pressure relaxation script
│   │   ├── post_relax_hp.sh    # High-pressure post-processing script
│   │   ├── energy_post.py      # Energy extraction and ranking
│   │   ├── energy_post_hp.py   # High-pressure energy post-processing
│   │   ├── reorder_energy.py   # Reorder energy data by composition
│   │   └── cif_dedup.py        # Structure deduplication by matching
│   │
│   ├── PTworkflows/            # Perturbation workflow scripts
│   │   ├── supercell.py        # Supercell generation
│   │   ├── energy_post_pt.py   # Perturbation energy post-processing
│   │   ├── relax_umlip_pt.sh   # Perturbation relaxation script
│   │   ├── post_relax_pt.sh    # Perturbation post-processing script
│   │   ├── replacements.yaml   # Element substitution mapping
│   │   ├── doping.yaml         # Doping configuration
│   │   └── INPUT               # Mutation parameter file
│   │
│   ├── ELAworkflows/           # Elastic constants workflow
│   │   ├── ela_dir.py          # Directory setup for ELA calculations
│   │   ├── uMLIP_ela_sq.py     # uMLIP-based elastic tensor calculation
│   │   ├── cif_pri.py          # Primitive cell conversion
│   │   ├── cal_ela.sh          # SLURM script for ELA calculation
│   │   ├── post_ela_uMLIP.py   # Elastic property extraction
│   │   ├── VPKIT.in1           # VPKIT elastic input (part 1)
│   │   └── VPKIT.in2           # VPKIT elastic input (part 2)
│   │
│   ├── QHAworkflows/           # Quasi-Harmonic Approximation workflow
│   │   ├── qha_dir.py          # Directory setup for QHA
│   │   ├── relax_volume.py     # Volume-strain relaxation
│   │   ├── run_phonon.py       # Phonon calculation driver
│   │   ├── collect_energy.py   # Energy collection across volumes
│   │   ├── cif_pri.py          # Primitive cell conversion
│   │   └── cal_qha.sh          # SLURM script for QHA
│   │
│   ├── MDworkflows/            # Molecular Dynamics workflow
│   │   ├── md_nvt.py           # NVT MD simulation (Langevin thermostat)
│   │   ├── md_post.py          # MD trajectory post-processing
│   │   ├── run_md_nvt.sh       # SLURM script for NVT MD
│   │   └── post_md.sh          # SLURM script for MD post-processing
│   │
│   ├── FTworkflows/            # Fine-Tuning / Active Learning
│   │   ├── ft_uncertainty.py   # Uncertainty estimation for active learning
│   │   └── run_ft.sh           # SLURM script for FT data preparation
│   │
│   ├── common/                 # Shared modules across workflows
│   │   ├── ehull/              # Energy above hull (Ehull) calculation
│   │   │   ├── Ehull_compatibility.py    # Ehull with MP2020 corrections
│   │   │   ├── Ehull_no_compatibility.py # Ehull without corrections
│   │   │   ├── Ehull_post.py             # Ehull post-processing & filtering
│   │   │   └── mp_offline_loader.py      # Offline Materials Project data loader
│   │   ├── phonon/             # Phonon calculation utilities
│   │   │   ├── ph_uMLIP.py     # uMLIP-based force calculation for phonons
│   │   │   ├── ph_post.py      # Phonon post-processing
│   │   │   └── ph_cal.sh       # SLURM script for phonon calculations
│   │   ├── postprocess/        # Structure post-processing
│   │   │   ├── sym_rename.py   # Symmetrize and rename by space group (3D)
│   │   │   ├── layergroup_sym.py  # Layer group classification (2D)
│   │   │   ├── chiral_check.py # Chirality classification
│   │   │   └── bond_check.py   # Bond length validation
│   │   ├── relaxation/         # Structure relaxation
│   │   │   └── umlip_relax.py  # Universal MLIP relaxation engine
│   │   ├── selection/          # Structure diversity selection
│   │   │   ├── structure_select.py  # Main selection driver
│   │   │   ├── base_selector.py     # Base selector class
│   │   │   ├── selector_0d.py       # 0D structure selector
│   │   │   ├── selector_1d.py       # 1D structure selector
│   │   │   ├── selector_2d.py       # 2D structure selector
│   │   │   ├── selector_3d.py       # 3D structure selector
│   │   │   └── run_selection.sh     # SLURM script for selection
│   │   └── thermal/            # Thermal conductivity utilities
│   │       ├── fc3_workflow.py  # Third-order force constants workflow
│   │       ├── tc_fc3.py        # FC3 calculation via phono3py
│   │       ├── tc_cal.sh        # SLURM script for TC calculation
│   │       └── tc_post.sh       # SLURM script for TC post-processing
│   │
│   └── phonon_src/             # Low-level phonon calculation library
│       ├── ph/
│       │   ├── phonon.py        # Phonon dispersion & DOS calculation
│       │   └── ph3cal.py        # Third-order phonon calculation
│       └── utils/
│           └── phonon_helpers.py  # Phonon helper functions
```

---

## Installation

### Prerequisites

- **Python**: >= 3.8, <= 3.11
- **HPC environment**: SLURM job scheduler (for batch calculations)
- **Conda** (recommended for environment management)

### Basic Installation

```bash
# Clone the repository
git clone https://github.com/JesseOOPP/GEWUM.git
cd gewum

# Install with core dependencies
pip install .
```

### Install with ML Potential Support

To enable MatterSim universal machine-learning interatomic potential:

```bash
pip install ".[ml]"
```

### Development Installation

```bash
pip install -e ".[dev]"
```

### Verify Installation

```bash
gewum --version
gewum --help
```

---

## Configuration

### SLURM Configuration (`slurm_config.yaml`)

GEWUM uses a YAML configuration file to customize SLURM job parameters and environment settings. Place `slurm_config.yaml` in your working directory; GEWUM will automatically detect and apply it when generating scripts.

**Template:**

```yaml
# SLURM job parameters
slurm:
  time: "2400:00:00"        # Maximum wall time (HH:MM:SS)
  cpus_per_task: 64          # CPUs per task
  partition: "amd_256"       # SLURM partition name
  nodes: 1                   # Number of nodes

# Environment setup
environment:
  module_purge: true         # Purge loaded modules before setup
  modules:                   # Modules to load (in order)
    - "cmake/3.22.0"
    - "gcc/12.2"
    - "intel/2022.1"
    - "mpi/intel/2022.1"
  conda_path: "/path/to/anaconda3"
  conda_env: "mattersim"     # Conda environment name

# Optional: GNU Parallel path (for parallel batch jobs)
parallel:
  path: "/path/to/parallel/bin"
```

Users should customize `partition`, `conda_path`, `conda_env`, and `modules` to match their HPC system.

---

## Usage / Execution

### Command Format

```bash
gewum <command> --mode <mode> [options]
```

### Available Commands

| Command   | Description                                              |
|-----------|----------------------------------------------------------|
| `cifgen`  | Crystal structure (CIF) generation from chemical formulas |
| `RD`      | Random Design: random structure search workflow           |
| `PT`      | Perturbation: perturbation-based structure search workflow |
| `ELA`     | Elastic constants calculation                             |
| `QHA`     | Quasi-Harmonic Approximation for thermal properties       |
| `TC`      | Thermal conductivity via third-order force constants      |
| `MD`      | Molecular dynamics simulation (NVT ensemble)              |
| `FT`      | Fine-tuning data preparation via uncertainty estimation   |

Use `gewum <command> -h` for detailed help on each command.

### Command Details and Examples

#### 1. `cifgen` — CIF Generation

Generate candidate crystal structures from chemical compositions.

```bash
# Generate CIFs for all possible compositions
gewum cifgen --mode all

# Generate CIFs filtered by common oxidation states
gewum cifgen --mode oxidation

# Generate structures via element substitution
gewum cifgen --mode substitute

# Generate perturbed structures from a template CIF
gewum cifgen --mode mutate

# Generate doped structures with random atom replacements
gewum cifgen --mode dp
```

#### 2. `RD` — Random Design Workflow

End-to-end random crystal structure exploration.

```bash
# Step 1: Generate random structures (copy SLURM scripts)
gewum RD --mode cifgen

# Step 2: Relax structures with uMLIP
gewum RD --mode relax

# Step 3: Refine pre-relaxed structures
gewum RD --mode refine

# Step 4: Post-process relaxation results
gewum RD --mode post

# Step 5: Reorder energy data by composition
gewum RD --mode reorder -r

# Step 6: Deduplicate structures
gewum RD --mode dedup -r

# Step 7: Calculate energy above convex hull
gewum RD --mode Ehull --api-key mp-XXXXX          # Online
gewum RD --mode Ehull --mp-data /path/to/MP.json   # Offline
gewum RD --mode Ehull --post -t 0.2                # Filter by threshold

# Step 8: Symmetry classification
gewum RD --mode sym -r       # 3D space group classification
gewum RD --mode sym2d -r     # 2D layer group classification
gewum RD --mode chiral -r    # Chirality classification

# Step 9: Structure diversity selection
gewum RD --mode select

# Step 10: Phonon calculations
gewum RD --mode ph           # Copy phonon scripts
gewum RD --mode phpost -r    # Post-process phonon results
```

**High-pressure variant:**

```bash
gewum RD --mode relaxhp      # High-pressure relaxation
gewum RD --mode posthp       # High-pressure post-processing
```

#### 3. `PT` — Perturbation Workflow

Structure exploration via perturbation of known structures.

```bash
# Generate supercells from input CIF files
gewum PT --mode supercell -r -i input.cif --matrix 2 2 2
gewum PT --mode supercell -r --batch --matrix 3 3 1

# Copy substitution/mutation/doping config files
gewum PT --mode sub          # Copy replacements.yaml
gewum PT --mode mutate       # Copy INPUT parameter file
gewum PT --mode dp           # Copy doping.yaml

# Relaxation and post-processing
gewum PT --mode relax
gewum PT --mode post

# Ehull, symmetry, phonon (same interface as RD)
gewum PT --mode Ehull --api-key mp-XXXXX
gewum PT --mode sym -r
gewum PT --mode ph
```

#### 4. `ELA` — Elastic Constants

```bash
# Prepare directories (one per CIF)
gewum ELA --mode pre -r

# Copy calculation scripts (includes VPKIT input files)
gewum ELA --mode cal

# Extract elastic properties to ela_tot.dat
gewum ELA --mode post -r
```

#### 5. `QHA` — Quasi-Harmonic Approximation

```bash
# Prepare directories
gewum QHA --mode pre -r

# Copy QHA calculation scripts
gewum QHA --mode cal
```

#### 6. `TC` — Thermal Conductivity

```bash
# Copy FC3 calculation scripts
gewum TC --mode fc3

# Copy post-processing scripts
gewum TC --mode post
```

#### 7. `MD` — Molecular Dynamics

```bash
# Copy NVT simulation scripts
gewum MD --mode nvt
# Run: sbatch run_md_nvt.sh [TEMP] [STEPS] [TIMESTEP] [NX] [NY] [NZ]
# Example: sbatch run_md_nvt.sh 300 10000 1.0 2 2 1

# Copy post-processing scripts
gewum MD --mode post
```

#### 8. `FT` — Fine-Tuning Data Preparation

```bash
# Copy uncertainty estimation scripts
gewum FT --mode prepare
```

---

## Comprehensive Test Run Example

The following demonstrates a complete **Random Design (RD) workflow** for exploring novel crystal structures of a ternary system. This example covers all major stages from structure generation to stability screening.

### Prerequisites

1. GEWUM installed (`pip install .` or `pip install ".[ml]"`)
2. SLURM cluster access with `slurm_config.yaml` configured
3. (Optional) Materials Project API key for Ehull calculation

### Step 1: Prepare Working Directory

```bash
mkdir test_run && cd test_run
cp /path/to/gewum/slurm_config.yaml .
# Edit slurm_config.yaml to match your HPC environment
```

### Step 2: Generate Random Crystal Structures

```bash
# Generate CIF structures for all compositions of interest
gewum cifgen --mode all
```

**Input:** The `cifgen` module reads element combinations and uses PyXtaL to enumerate candidate structures.

**Output:** A set of `.cif` files representing candidate crystal structures.

### Step 3: Copy and Submit Relaxation Scripts

```bash
# Copy the random structure generation SLURM script
gewum RD --mode cifgen
sbatch cifgen.sh

# After structure generation completes, copy relaxation scripts
gewum RD --mode relax
sbatch relax_umlip.sh
```

The relaxation step uses the MatterSim universal machine-learning interatomic potential to optimize atomic positions and cell parameters.

### Step 4: Post-Process Relaxation Results

```bash
gewum RD --mode post
sbatch post_relax.sh
```

**Output:** `energy_final.txt` — a ranked list of structures with their total energies and energies per atom.

### Step 5: Reorder and Deduplicate

```bash
# Reorder by composition
gewum RD --mode reorder -r

# Remove duplicate structures
gewum RD --mode dedup -r
```

### Step 6: Thermodynamic Stability Screening (Ehull)

```bash
# Calculate energy above convex hull (online mode)
gewum RD --mode Ehull --api-key mp-XXXXX

# Or offline mode with pre-downloaded MP data
gewum RD --mode Ehull --mp-data /path/to/MP_entries.json

# Filter structures by Ehull threshold (e.g., 0.2 eV/atom)
gewum RD --mode Ehull --post -t 0.2
```

**Output:** `Ehull_results.csv` — Ehull values for all candidates. Structures below the threshold are extracted to `0_cif/`.

### Step 7: Symmetry Classification

```bash
gewum RD --mode sym -r
```

**Output:** Structures are renamed and organized by space group number (e.g., `Fm-3m_225_*.cif`).

### Step 8: Phonon Stability Check

```bash
gewum RD --mode ph
sbatch ph_cal.sh

# After phonon calculations complete:
gewum RD --mode phpost -r
```

**Output:** Phonon dispersion plots and `ph.log` files. Structures with no imaginary frequencies are dynamically stable.

### Expected Results Summary

| Stage               | Key Output                          | Description                                    |
|---------------------|--------------------------------------|------------------------------------------------|
| CIF generation      | `*.cif` files                       | Random candidate crystal structures            |
| Relaxation          | Relaxed `*.cif` files               | Energy-minimized structures                    |
| Post-processing     | `energy_final.txt`                  | Ranked energy list (eV/atom)                   |
| Deduplication       | Reduced structure set               | Unique structures by composition               |
| Ehull filtering     | `Ehull_results.csv`, `0_cif/`       | Thermodynamically viable candidates            |
| Symmetry            | Renamed CIF files by space group    | Classified by crystallographic symmetry        |
| Phonon              | Dispersion plots, `ph.log`          | Dynamic stability verification                 |

---

---

## Dependencies

### Core Dependencies

| Package       | Version   | Purpose                                    |
|---------------|-----------|--------------------------------------------|
| numpy         | >= 1.20.0 | Numerical computation                      |
| pandas        | >= 1.3.0  | Data manipulation and analysis             |
| matplotlib    | >= 3.4.0  | Plotting and visualization                 |
| scipy         | >= 1.7.0  | Scientific computing utilities             |
| ase           | >= 3.22.0 | Atomic Simulation Environment              |
| pymatgen      | >= 2022.0 | Materials analysis and structure handling   |
| pyxtal        | >= 0.6.0  | Random crystal structure generation        |
| spglib        | >= 2.0.0  | Space group detection and symmetry analysis|
| phonopy       | >= 2.20.0 | Phonon calculation framework               |
| phono3py      | >= 2.0.0  | Third-order phonon / thermal conductivity  |
| mp-api        | >= 0.30.0 | Materials Project API client               |
| scikit-learn  | >= 1.0.0  | Machine learning utilities for selection   |
| pyyaml        | >= 5.4.0  | YAML configuration parsing                 |
| tqdm          | (any)     | Progress bars                              |

### Optional Dependencies

| Package    | Install Command        | Purpose                                   |
|------------|------------------------|-------------------------------------------|
| mattersim  | `pip install ".[ml]"`  | MatterSim universal ML interatomic potential |

### Development Dependencies

| Package | Install Command         | Purpose           |
|---------|-------------------------|-------------------|
| pytest  | `pip install ".[dev]"`  | Testing framework |
| black   | `pip install ".[dev]"`  | Code formatting   |
| flake8  | `pip install ".[dev]"`  | Linting           |

---

## License

This project is distributed under the **MIT License**. See `LICENSE` for details.

## Citation

If you use GEWUM in your research, please cite:

> Jiexi Song, Aixian She, Changpeng Song, Diwei Shi, Fengyuan Xuan, Chongde Cao. *GEWUM: General Exploration Workflow for the Utopia of Materials: A Unified Platform for Automated Structure Generation, Selection, and Validation.* arXiv:2604.21401 [cond-mat.mtrl-sci], April 2026. DOI: [10.48550/arXiv.2604.21401](https://doi.org/10.48550/arXiv.2604.21401)

**BibTeX:**

```bibtex
@article{song2026gewum,
  title={GEWUM: General Exploration Workflow for the Utopia of Materials: A Unified Platform for Automated Structure Generation, Selection, and Validation},
  author={Song, Jiexi and She, Aixian and Song, Changpeng and Shi, Diwei and Xuan, Fengyuan and Cao, Chongde},
  journal={arXiv preprint arXiv:2604.21401},
  year={2026}
}
```
