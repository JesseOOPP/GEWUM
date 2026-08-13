# GEWUM: General Exploration Workflow for the Utopia of Materials

Version: 1.0.0
Authors: Jiexi Song, Aixian She, Changpeng Song, Diwei Shi, Fengyuan Xuan, Chongde Cao
Contact: songjx@szlab.ac.cn
GitHub: https://github.com/JesseOOPP/GEWUM
License: MIT

GEWUM is an integrated computational workflow platform for materials science that automates crystal structure generation, structure relaxation via universal machine-learning interatomic potentials (uMLIP), thermodynamic and dynamical stability screening, and property calculations (elastic constants, thermal conductivity, quasi-harmonic approximation, molecular dynamics). It orchestrates the complete pipeline from random/perturbation-based structure exploration to phonon, convex-hull, and thermophysical analyses, with native SLURM-based high-throughput execution.

---

## Directory Structure

```
gewum/
|-- main.py                     # CLI entry point (command dispatcher)
|-- config.py                   # Centralized path and configuration management
|-- pyproject.toml              # Package metadata, dependencies, build config
|-- setup.py                    # Legacy setuptools support
|-- requirements.txt            # Dependency listing
|-- environment.yml             # Tested conda environment (pinned versions)
|-- requirements-lock.txt       # pip-managed lockfile (pinned versions)
|-- slurm_config.yaml           # SLURM job scheduler configuration template
|-- MANIFEST.in                 # Package data inclusion rules
|
|-- commands/                   # CLI command modules (one per workflow)
|   |-- __init__.py
|   |-- base.py                 # BaseWorkflowCommand abstract base class
|   |-- template_utils.py       # SLURM template rendering utilities
|   |-- cifgen.py               # CIF generation command
|   |-- RD.py                   # Random Design workflow command
|   |-- PT.py                   # Perturbation workflow command
|   |-- ELA.py                  # Elastic constants workflow command
|   |-- QHA.py                  # Quasi-Harmonic Approximation command
|   |-- TC.py                   # Thermal Conductivity command
|   |-- MD.py                   # Molecular Dynamics command
|   `-- DB.py                   # Read-only structures.db inspection/export
|
|-- src/                        # Core computational modules
|   |-- cifgen_input/           # CIF generation engines (chemi, oxidation,
|   |                           #   substitute, mutate, doping)
|   |-- RDworkflows/            # Random Design workflow scripts
|   |-- PTworkflows/            # Perturbation workflow scripts
|   |-- ELAworkflows/           # Elastic constants workflow
|   |-- QHAworkflows/           # Quasi-Harmonic Approximation workflow
|   |-- MDworkflows/            # Molecular Dynamics workflow
|   |-- common/                 # Shared modules across workflows
|   |   |-- db_admin/           # structures.db management (query/formatter/stats/exporter)
|   |   |-- ehull/              # Energy above hull (Ehull) calculation
|   |   |-- phonon/             # Phonon calculation utilities
|   |   |-- postprocess/        # Symmetry rename, bond check, visualization
|   |   |-- relaxation/         # uMLIP relaxation engine
|   |   |-- selection/          # Structure diversity selection (0D-3D)
|   |   `-- thermal/            # Third-order force constants / thermal conductivity
|   `-- phonon_src/             # Low-level phonon calculation library
|
`-- tests/                      # pytest unit & CLI smoke tests
```

## Installation

### Prerequisites

- Python: >= 3.8
- Conda (recommended for environment management)
- SLURM job scheduler (only for batch execution; all workflow scripts can also be run locally)

### Basic Installation

```
git clone https://github.com/JesseOOPP/GEWUM.git
cd gewum

# Install with core dependencies
pip install .
```

### Tested Environment (Recommended)

The versions below were verified on the development cluster (Linux x86_64,
Python 3.10.20) where all GEWUM workflows were run. To recreate the exact
tested environment:

```
conda env create -f environment.yml
conda activate gewum
pip install .                    # or: pip install ".[ml]" for MatterSim support
```

[environment.yml](environment.yml) pins the conda-managed packages (numpy,
matplotlib, phonopy, phono3py, pyyaml, spglib);
[requirements-lock.txt](requirements-lock.txt) pins the pip-managed packages.
`pyproject.toml` additionally declares minimum version requirements so that
a plain `pip install .` works with any compatible combination.

### Install with ML Potential Support

To enable MatterSim universal machine-learning interatomic potential:

```
pip install ".[ml]"
```

### Development Installation

```
pip install -e ".[dev]"
```

### Verify Installation

```
gewum --version
gewum --help
```

## Configuration

### SLURM Configuration (slurm_config.yaml)

GEWUM uses a single YAML configuration file to customize SLURM job parameters and environment settings. Place `slurm_config.yaml` in your working directory; GEWUM will automatically detect and apply it when generating scripts (search order: target directory -> current directory -> installation directory).

```yaml
# SLURM job parameters
slurm:
  time: "2400:00:00"        # Maximum wall time (HH:MM:SS)
  cpus_per_task: 64          # CPUs per task
  partition: "<partition>"   # SLURM partition name
  nodes: 1                   # Number of nodes

# Environment setup
environment:
  module_purge: true         # Purge loaded modules before setup
  modules:                   # Modules to load (in order)
    - "cmake/<version>"
    - "gcc/<version>"
    - "intel/<version>"
    - "mpi/<version>"
  conda_path: "/path/to/anaconda3"
  conda_env: "<env_name>"    # Conda environment name

# Optional: GNU Parallel path (for parallel batch jobs)
parallel:
  path: "/path/to/parallel/bin"
```

## Usage / Execution

### Command Format

```
gewum <command> --mode <mode> [options]
```

### Available Commands

| Command | Description |
| --- | --- |
| cifgen | Generate CIF generation input from chemical formulas |
| RD | Selective Random Structure Search (SRSS) workflow |
| PT | Perturbation-based structure search workflow |
| ELA | Elastic constants calculation |
| QHA | Quasi-Harmonic Approximation for thermal properties |
| TC | Thermal conductivity via third-order force constants |
| MD | Molecular dynamics simulation (NVT ensemble) |
| DB | Read-only inspection and batch CIF export of structures.db |

Use `gewum <command> -h` for detailed help on each command.

### Command Details

#### 1. cifgen - CIF Generation Input

```
# Generate input for all possible compositions
gewum cifgen --mode all

# Generate input filtered by common oxidation states
gewum cifgen --mode oxidation

# Generate structures via element substitution
gewum cifgen --mode substitute

# Generate perturbed structures from a template CIF
gewum cifgen --mode mutate

# Generate doped structures with random atom replacements
gewum cifgen --mode dp
```

#### 2. RD - Random Design Workflow (SRSS)

End-to-end random crystal structure exploration. Each mode copies the corresponding workflow scripts into the working directory.

```
# Step 1: Generate random structures (copy cifgen.sh, then sbatch)
gewum RD --mode cifgen

# Step 2: Relax structures with uMLIP
gewum RD --mode relax

# Step 3: Cell refinement
gewum RD --mode refine

# Step 4: Post-process relaxation results
gewum RD --mode post

# Step 5: Single-point energy calculation
gewum RD --mode calc

# Step 6: Reorder energy data by composition
gewum RD --mode reorder -r

# Step 7: Deduplicate structures
gewum RD --mode dedup -r

# Step 8: Calculate energy above convex hull
gewum RD --mode Ehull --api-key mp-XXXXX          # Online
gewum RD --mode Ehull --mp-data /path/to/MP.json   # Offline
gewum RD --mode Ehull --post -t 0.2                # Filter by threshold

# Step 9: Symmetry classification
gewum RD --mode sym -r       # 3D space group classification
gewum RD --mode sym2d -r     # 2D layer group classification
gewum RD --mode chiral -r    # Chirality classification

# Step 10: Structure diversity selection
gewum RD --mode select

# Step 11: Phonon calculations
gewum RD --mode ph           # Copy phonon scripts
gewum RD --mode phpost -r    # Post-process phonon results

# Step 12: Visualization
gewum RD --mode viz          # Structure/energy visualization
gewum RD --mode viz2         # Advanced analysis (RDF, SG curves)
gewum RD --mode viz3         # Extended visualization
```

High-pressure variant:

```
gewum RD --mode relaxhp      # High-pressure relaxation
gewum RD --mode posthp       # High-pressure post-processing
```

Automated pipeline (all-in-one script set):

```
gewum RD --mode auto         # Copies cifgen.sh + relax_umlip.sh + post_relax.sh + run_srss.sh
```

#### 3. PT - Perturbation Workflow

Structure exploration via perturbation of known structures.

```
# Generate supercells from input CIF files
gewum PT --mode supercell -r -i input.cif --matrix 2 2 2

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
gewum PT --mode viz          # PT visualization
```

#### 4. ELA - Elastic Constants

```
# Prepare directories (one per CIF)
gewum ELA --mode pre -r

# Copy calculation scripts (includes VPKIT input files)
gewum ELA --mode cal

# Extract elastic properties to ela_tot.dat
gewum ELA --mode post -r
```

#### 5. QHA - Quasi-Harmonic Approximation

```
# Prepare directories
gewum QHA --mode pre -r

# Copy QHA calculation scripts
gewum QHA --mode cal
```

#### 6. TC - Thermal Conductivity

```
# Copy FC3 calculation scripts
gewum TC --mode fc3

# Copy post-processing scripts
gewum TC --mode post
```

#### 7. MD - Molecular Dynamics

```
# Copy NVT simulation scripts
gewum MD --mode nvt
# Run: sbatch run_md_nvt.sh [TEMP] [STEPS] [TIMESTEP] [NX] [NY] [NZ]
# Example: sbatch run_md_nvt.sh 300 10000 1.0 2 2 1

# Copy post-processing scripts
gewum MD --mode post

# Structural-quantity ensemble workflow
gewum MD --mode sq
```

#### 8. DB - structures.db Inspection & Export

Read-only inspection and batch CIF export of `structures.db` files (never writes to them).

```
# Per-DB metadata: row counts by stage, energy range, schema
gewum DB info --root <work_dir>

# List rows matching a filter (table/csv/json/tsv)
gewum DB list --root <work_dir> --formula 'Na*Cl*' --stage relaxed --format csv

# GROUP BY pivot (1 or 2 dimensions: formula / stage / sg)
gewum DB stats --root <work_dir> --by formula,stage

# Dump the CIF of a single row
gewum DB show --db <path>/structures.db --sg 225 --name xtal_001

# Batch export filtered CIFs to a directory or zip archive
gewum DB export --root <work_dir> --out exported_cifs/ --layout by-formula
gewum DB export --root <work_dir> --zip archive.zip --dry-run
```

## Comprehensive Test Run Example

The following demonstrates a complete Random Design (RD) workflow for exploring novel crystal structures of a ternary system. It covers all major stages from structure generation to stability screening.

### Prerequisites

1. GEWUM installed (`pip install .` or `pip install ".[ml]"`)
2. SLURM cluster access with `slurm_config.yaml` configured (or run scripts directly for small systems)
3. (Optional) Materials Project API key for Ehull calculation

### Step 1: Prepare Working Directory

```
mkdir test_run && cd test_run
cp /path/to/gewum/slurm_config.yaml .
# Edit slurm_config.yaml to match your HPC environment
```

### Local (non-SLURM) Execution

All single-stage workflow scripts (e.g., `cifgen.sh`, `relax_umlip.sh`, `post_relax.sh`) can be executed directly with `bash` on a laptop or single workstation - SLURM environment variables are used only when present (e.g., `TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-64}` falls back to a safe default). Only the fully automated pipeline `run_srss.sh` requires a SLURM cluster because it chains jobs via `sbatch --dependency`.

### Pipeline Overview

The table below summarizes the full RD pipeline used in the walkthrough; each command copies the corresponding workflow scripts into the working directory (see [Command Details](#command-details) for per-mode options).

| Step | Command(s) | Key Output |
| --- | --- | --- |
| 1. Generate structures | `gewum cifgen --mode all`; `gewum RD --mode cifgen`; `sbatch cifgen.sh` | `*.cif` candidates (0D-3D) |
| 2. Relax with uMLIP | `gewum RD --mode relax`; `sbatch relax_umlip.sh` | Relaxed CIFs (MatterSim) |
| 3. Post-process | `gewum RD --mode post`; `sbatch post_relax.sh` | `energy_final.txt` - ranked list (eV/atom) |
| 4. Reorder & deduplicate | `gewum RD --mode reorder -r`; `gewum RD --mode dedup -r` | Unique structures by composition |
| 5. Ehull screening | `gewum RD --mode Ehull --api-key mp-XXXX` (online) / `--mp-data` (offline); `--post -t 0.2` to filter | `Ehull_results.csv`, `0_cif/` |
| 6. Symmetry classification | `gewum RD --mode sym -r` | CIFs renamed by space group (e.g., `Fm-3m_225_*.cif`) |
| 7. Phonon stability check | `gewum RD --mode ph`; `sbatch ph_cal.sh`; `gewum RD --mode phpost -r` | Dispersion plots, `ph.log` |
| 8. Visualization (optional) | `gewum RD --mode viz -r` / `viz2 -r` / `viz3 -r` | UMAP/t-SNE maps, phase diagrams |

Structures that survive Step 7 are both thermodynamically (convex hull) and dynamically (no imaginary phonon modes) stable.

## Dependencies

### Core Dependencies

| Package | Purpose |
| --- | --- |
| numpy | Numerical computation |
| pandas | Data manipulation and analysis |
| matplotlib | Plotting and visualization |
| scipy | Scientific computing utilities |
| ase | Atomic Simulation Environment |
| pymatgen | Materials analysis and structure handling |
| pyxtal | Random crystal structure generation |
| spglib | Space group detection and symmetry analysis |
| phonopy | Phonon calculation framework |
| phono3py | Third-order phonon / thermal conductivity |
| mp-api | Materials Project API client |
| scikit-learn | Machine learning utilities for selection |
| dscribe | Structure feature descriptors |
| hdbscan | Density-based clustering |
| umap-learn | Dimensionality reduction for visualization |
| pyyaml | YAML configuration parsing |
| ijson | Incremental JSON parsing |
| tqdm | Progress bars |

### Optional Dependencies

| Package | Install Command | Purpose |
| --- | --- | --- |
| mattersim | `pip install ".[ml]"` | MatterSim universal ML interatomic potential |

### Development Dependencies

| Package | Install Command | Purpose |
| --- | --- | --- |
| pytest | `pip install ".[dev]"` | Testing framework |
| black | `pip install ".[dev]"` | Code formatting |
| flake8 | `pip install ".[dev]"` | Linting |

## Testing

```
pip install -e ".[dev]"
pytest tests/ -v
```

The test suite covers SLURM template rendering, configuration management, and CLI smoke tests. It does not require a SLURM cluster or any uMLIP backend.

## License

This project is distributed under the MIT License. See [LICENSE](LICENSE) for details.

## Citation

If you use GEWUM in your research, please cite:

> Jiexi Song, Aixian She, Changpeng Song, Diwei Shi, Fengyuan Xuan, Chongde Cao. GEWUM: General Exploration Workflow for the Utopia of Materials: A Unified Platform for Automated Structure Generation, Selection, and Validation. arXiv:2604.21401 [cond-mat.mtrl-sci], April 2026. DOI: 10.48550/arXiv.2604.21401

BibTeX:

```bibtex
@article{song2026gewum,
  title={GEWUM: General Exploration Workflow for the Utopia of Materials: A Unified Platform for Automated Structure Generation, Selection, and Validation},
  author={Song, Jiexi and She, Aixian and Song, Changpeng and Shi, Diwei and Xuan, Fengyuan and Cao, Chongde},
  journal={arXiv preprint arXiv:2604.21401},
  year={2026}
}
```
