# GEWUM Laptop Minimal Working Example

End-to-end run of the GEWUM **RD workflow** (structure generation ->
database -> diversity selection -> uMLIP relaxation -> energy ranking ->
convex-hull stability screening) on a **single workstation, without SLURM**.

Structures come from the real cifgen pipeline: this folder's `cifgen.inp`
(the Mg-Sn compositions you list, 20 structures per space group each) is
fed to `cif_generate`, which walks **all space groups (2-230)** exactly
like a production run. The selection step then keeps 2 of the 20
structures per space group, and only those are relaxed.

## 1. Installation (fresh environment)

```bash
# Python >= 3.8  #install a uMLIP (for example, pip install mattersim)
git clone https://github.com/JesseOOPP/GEWUM.git
cd GEWUM
pip install -e ".[ml]"        # GEWUM + mattersim (uMLIP) + all core deps
```

Verify the install:

```bash
gewum -h          # command-line help
gewum RD -h       # per-workflow help (modes and options)
```

## 2. Run

```bash
cd examples/laptop_mwe
mkdir mwe && cd mwe
python ../run_laptop_mwe.py         
```

No network access is needed. Elemental references for the hull come
from single-element rows in `cifgen.inp` - relaxed exactly like any
other structure, so nothing is preset and the hull is fully
uMLIP-consistent. The bundled `cifgen.inp` currently lists only the
MgSn compound; add rows such as `['Mg'], [1], 2` and `['Sn'], [1], 2`
to include elemental references in the hull. The structures themselves
are pyxtal random crystals (like any cifgen run), so re-runs differ;
the pipeline steps and their order are fixed.

## 3. What it runs

| Step | Command (what the script calls) | HPC counterpart | Output |
| --- | --- | --- | --- |
| 1. Generation | `cifgen.inp` + `python -m gewum.src.RDworkflows.cif_generate --dim 3 --max-atoms 12 --max-attempts 20` | `gewum RD --mode cifgen` (cifgen.sh) | `<formula>/structures.db` (20 CIF rows per SG, stage=initial) |
| 2. Selection | `python -m gewum.src.common.selection.structure_select --target 2` | `gewum RD --mode select` (run_selection.sh) | 2 per SG kept; rest marked `removed` in DB |
| 3. Relaxation | `python -m gewum.src.common.relaxation.relax_db_io worker/commit` | `gewum RD --mode relax` (relax_umlip.sh) | stage=relaxed rows + energy in DB, energy_results.csv |
| 4. Post | `python -m gewum.src.common.postprocess.bond_check` + `python -m gewum.src.RDworkflows.energy_post` | `gewum RD --mode post` (post_relax.sh) | energy_final.txt, 0_final_results.txt, 0_cif_final/, 0_final_result_tot.txt |
| 5. Ehull | `python -m gewum.src.common.ehull.Ehull_compatibility --self-hull` | `gewum RD --mode Ehull` | Hull_result.csv (formation energy, e_above_hull) |

## 4. Expected outputs

```
mwe/
|-- MgSn1/                         # one dir per composition in cifgen.inp
|   |-- structures.db              # sqlite: initial -> relaxed -> kept/removed
|   |-- energy_results.csv         # per-structure relaxed energies
|   |-- energy_final.txt           # all relaxed entries, ranked
|   |-- 0_final_results.txt        # selected subset (gap filter)
|   `-- 0_cif_final/               # selected CIFs
|-- cifgen.inp                     # generation input (copied by the script)
|-- cifgen.out                     # cifgen log
|-- 0_final_result_tot.txt         # global ranking (CSV)
`-- Hull_result.csv                # formation energy + e_above_hull
```
## 5. HPC counterpart (the same pipeline with SLURM)

On a cluster, the identical workflow is driven by the `gewum` CLI, which
copies the templated shell scripts (SLURM header + environment injected
from `slurm_config.yaml`) and submits them:

```bash
# 1. Random crystal generation (0D-3D)
gewum cifgen --mode all            # edit cifgen.inp
gewum RD --mode cifgen
sbatch cifgen.sh 3 24 150

# 2. Diversity selection
gewum RD --mode select             # edit run_selection.sh (dim, method, target)
sbatch run_selection.sh

# 3. uMLIP relaxation
gewum RD --mode relax              # edit relax_umlip.sh (mode 1/2, fmax, steps)
sbatch relax_umlip.sh

# 4. Post-processing: bond check + energy ranking
gewum RD --mode post
sbatch post_relax.sh

# 5. Convex-hull screening (online or offline MP data)
gewum RD --mode Ehull --api-key YOUR_MP_KEY
gewum RD --mode Ehull --mp-data /path/to/MPtrj.json

# ...or the fully automated pipeline (chains the 5 jobs):
gewum RD --mode auto               # edit run_srss.sh first
sbatch run_srss.sh
```
