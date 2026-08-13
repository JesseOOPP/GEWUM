"""
GEWUM MD (Molecular Dynamics) Workflow Command
Provides scripts for NVT molecular dynamics simulation using MatterSim
"""
import os
import shutil
from ..config import MD_REPOSITORY
from .template_utils import get_config_info, copy_with_template

MD_FILE_MAP = {
    "nvt": ["run_md_nvt.sh"],
    "post": ["post_md.sh"],
    "sq": ["run_sq_ensemble.sh"],
}

MODE_DESCRIPTIONS = """
  nvt: NVT molecular dynamics simulation (Langevin thermostat)
  post: MD post-processing (energy-time plots, statistics)
  sq: MD-free amorphous ensemble via stochastic quenching
      (coarse relax -> gate -> maxmin select -> refine -> gate; amorphous ensemble)
"""


def setup_args(parser):
    all_modes = list(MD_FILE_MAP.keys())
    
    parser.add_argument(
        "--mode", 
        required=True,
        choices=all_modes,
        help=f"GEWUM MD mode to copy:\n{MODE_DESCRIPTIONS}"
    )
    parser.add_argument(
        "--dest", 
        default=".",
        help="Target directory (default: current directory)"
    )


def execute(args, remaining_args=None):
    mode = args.mode
    dest = args.dest
    
    if mode not in MD_FILE_MAP:
        print(f"Error: Unknown mode '{mode}'")
        return
    
    files = MD_FILE_MAP[mode]
    if not files:
        print(f"No files configured for mode: {mode}")
        return
    
    os.makedirs(dest, exist_ok=True)
    
    config_path, config = get_config_info(dest)
    
    print(f"GEWUM MD [{mode}] - Copying files to: {dest}")
    if config_path:
        print(f"Using SLURM config: {config_path}")
    else:
        print("No slurm_config.yaml found, using default values in scripts")
    print("-" * 50)
    
    copied_count = 0
    for filename in files:
        src_path = os.path.join(MD_REPOSITORY, filename)
        dst_path = os.path.join(dest, filename)
        
        if os.path.exists(src_path):
            if filename.endswith('.sh') and config is not None:
                copy_with_template(src_path, dst_path, config)
            else:
                shutil.copy2(src_path, dst_path)
            print(f"  [OK] {filename}")
            copied_count += 1
        else:
            print(f"  [FAIL] {filename} (not found)")
    
    print("-" * 50)
    print(f"Copied {copied_count}/{len(files)} files")
    
    if mode == "nvt":
        print("\nUsage:")
        print("  sbatch run_md_nvt.sh [TEMP] [STEPS] [TIMESTEP] [NX] [NY] [NZ]")
        print("  Example: sbatch run_md_nvt.sh 300 10000 1.0 2 2 1")
    elif mode == "post":
        print("\nUsage:")
        print("  sbatch post_md.sh")
    elif mode == "sq":
        print("\nUsage:")
        print("  sbatch run_sq_ensemble.sh [N_CONFIGS] [NX] [NY] [NZ] [FMAX] [D_MIN_SCALE] [VARIABLE_CELL] [N_WORKERS] [N_PERTURB_WORKERS] [N_SELECT]")
        print("  Example: sbatch run_sq_ensemble.sh 20 2 2 2 0.05 0.7 1 8 0 100")
        print("  Note: keep N_WORKERS*N_CONFIGS comfortably larger than N_SELECT")
        print("        so the gate + maxmin selection keep a diverse margin.")
