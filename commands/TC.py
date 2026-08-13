"""
GEWUM TC (Thermal Conductivity) Workflow Command
Provides scripts for thermal conductivity calculation using phono3py
"""
import os
import shutil
from ..config import GEWUM_SOURCE_DIR
from .template_utils import get_config_info, copy_with_template

TC_REPOSITORY = os.path.join(GEWUM_SOURCE_DIR, "common", "thermal")

TC_FILE_MAP = {
    "fc3": ["tc_cal.sh"],
    "post": ["tc_post.sh"],
}

MODE_DESCRIPTIONS = """
  fc3: Calculate third-order force constants (FC3) using MatterSim
  post: Post-processing to compute thermal conductivity (kappa)
"""


def setup_args(parser):
    """Configure command arguments"""
    parser.add_argument(
        "--mode",
        required=True,
        choices=list(TC_FILE_MAP.keys()),
        help=f"GEWUM TC mode to copy:\n{MODE_DESCRIPTIONS}"
    )
    parser.add_argument(
        "--dest",
        default=".",
        help="Target directory (default: current directory)"
    )


def execute(args, remaining_args=None):
    """Copy GEWUM TC scripts to target directory"""
    copied = []
    missing = []
    
    os.makedirs(args.dest, exist_ok=True)
    
    config_path, config = get_config_info(args.dest)
    
    print(f" Preparing GEWUM TC {args.mode} scripts...")
    print(f" Source repository: {TC_REPOSITORY}")
    
    if config_path:
        print(f" Using SLURM config: {config_path}")
    else:
        print(" No slurm_config.yaml found, using default values in scripts")
    
    if not os.path.exists(TC_REPOSITORY):
        print(f"[ERROR] Source repository not found: {TC_REPOSITORY}")
        return
    
    for filename in TC_FILE_MAP[args.mode]:
        src_path = os.path.join(TC_REPOSITORY, filename)
        dst_path = os.path.join(args.dest, filename)
        
        if not os.path.exists(src_path):
            missing.append(filename)
            continue
        
        if filename.endswith('.sh') and config is not None:
            copy_with_template(src_path, dst_path, config)
        else:
            shutil.copy2(src_path, dst_path)
        copied.append(filename)
    
    print("\n" + "=" * 60)
    print(f" GEWUM TC {args.mode} scripts copy results")
    print("-" * 60)
    
    if copied:
        print("Copied files:")
        for f in copied:
            print(f"  - {f}")
    
    if missing:
        print("\n Missing files in GEWUM repository:")
        for f in missing:
            print(f"  - {f}")
    
    print(f"\n Target location: {os.path.abspath(args.dest)}")
    print("=" * 60)
