import os
import shutil
import subprocess
import sys
from ..config import CIFGEN_MODULE, ELA_REPOSITORY
from .template_utils import get_config_info, copy_with_template

GEWUM_FILE_MAP = {
    "pre": [],
    "cal": ["cal_ela.sh", "VPKIT.in1", "VPKIT.in2"],
    "post": [],
}

DIRECT_EXEC_MODULES = {
    "pre": "gewum.src.ELAworkflows.ela_dir",
    "post": "gewum.src.ELAworkflows.post_ela_uMLIP",
}

MODE_DESCRIPTIONS = """
  pre: Create a directory for each CIF file (-r)
  cal: ELA calculation scripts
  post: Extract elastic properties to ela_tot.dat (-r)
"""

def setup_args(parser):
    parser.add_argument("--mode", 
                        required=True,
                        choices=list(GEWUM_FILE_MAP.keys()),
                        help=f"GEWUM ELA mode:\n{MODE_DESCRIPTIONS}"
    )
    parser.add_argument("--dest", 
                        default=".",
                        help="Target directory (default: current directory)")
    parser.add_argument(
        "-r",
        action="store_true",
        dest="run",
        help="Directly execute (pre, post)"
    )

def execute(args, remaining_args=None):
    """Execute the ELA workflow command"""
    if hasattr(args, 'run') and args.run:
        if args.mode in DIRECT_EXEC_MODULES:
            module_name = DIRECT_EXEC_MODULES[args.mode]
            print(f"Executing {args.mode} in current directory...")
            try:
                result = subprocess.run([sys.executable, "-m", module_name], check=True)
                return
            except subprocess.CalledProcessError as e:
                print(f"Error executing {args.mode}: {e}")
                sys.exit(1)
        else:
            print(f"Warning: Mode '{args.mode}' does not support direct execution (-r).")
            print(f"Supported modes: {', '.join(DIRECT_EXEC_MODULES.keys())}")
            print("Falling back to copy mode...")
    
    copied = []
    missing = []
    
    os.makedirs(args.dest, exist_ok=True)
    
    config_path, config = get_config_info(args.dest)
    
    print(f" Preparing GEWUM ELA {args.mode} scripts...")
    print(f" Source repository: {ELA_REPOSITORY}")
    
    if config_path:
        print(f" Using SLURM config: {config_path}")
    else:
        print(" No slurm_config.yaml found, using default values in scripts")
    
    if not os.path.exists(ELA_REPOSITORY):
        print(f"[ERROR] Source repository not found: {ELA_REPOSITORY}")
        return
    
    for filename in GEWUM_FILE_MAP[args.mode]:
        src_path = os.path.join(ELA_REPOSITORY, filename)
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
    print(f" GEWUM ELA {args.mode} scripts copy results")
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
