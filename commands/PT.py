"""GEWUM PT (Perturbation) Workflow Command
Provides scripts for structure perturbation and processing
"""
import os
import shutil
import subprocess
import sys
from ..config import PT_REPOSITORY, CIFGEN_MODULE

PT_FILE_MAP = {
    "supercell": [],  
    "sub": ["replacements.yaml"],  
    "mutate": ["INPUT"],  
    "dp": ["doping.yaml"],  
    "relax": ["relax_umlip_pt.sh"],  
    "post": ["post_relax_pt.sh"],   
    "dedup": [],
    "Ehull": [],  
    "sym": [],    
    "sym2d": [],
    "chiral": [],
    "ph": [],     
    "phpost": [],
    "viz": ["run_viz_pt.sh"],
}

DIRECT_EXEC_MODULES = {
    "supercell": "gewum.src.PTworkflows.supercell",
    "dedup": "gewum.src.RDworkflows.cif_dedup",
    "sym": "gewum.src.common.postprocess.sym_rename",
    "sym2d": "gewum.src.common.postprocess.layergroup_sym",
    "chiral": "gewum.src.common.postprocess.chiral_check",
    "phpost": "gewum.src.common.phonon.ph_post",
    "viz": "gewum.src.PTworkflows.viz_pt",
}

EHULL_MODES = ["Ehull"]

CIFGEN_FILES = {}

MODE_DESCRIPTIONS = """
  supercell: Generate supercells from input CIF files (-r)
  sub: replacements.yaml
  mutate: INPUT
  dp: doping.yaml
  relax: Structure relaxation 
  post: Post-processing
  dedup: Deduplicate CIF structures by composition-aware structure matching (-r, --rdf --sim-t 0.2)
  Ehull: Calculate Ehull (-r, -cor, --mp-data,--api-key, --post -t 0.2)
  sym: Symmetrize by space group (-r)
  sym2d: Classify 2D layer group number (-r)
  chiral: Classify 3D chiral structures (-r)
  ph: Phonon calculation
  phpost: Phonon post-processing with ph.log (-r)
  viz: PT structure visualization
"""


def setup_args(parser):
    """Configure command arguments"""
    all_modes = list(PT_FILE_MAP.keys())
    
    parser.add_argument(
        "--mode", 
        required=True,
        choices=all_modes,
        help=f"GEWUM PT mode to copy:\n{MODE_DESCRIPTIONS}"
    )
    parser.add_argument(
        "--dest", 
        default=".",
        help="Target directory (default: ./)"
    )
    parser.add_argument(
        "-r",
        action="store_true",
        dest="run",
        help="Directly execute (supercell, sym, sym2d, chiral, phpost, Ehull)"
    )
    parser.add_argument(
        "-i",
        metavar="FILE",
        dest="input",
        help="Input"
    )
    parser.add_argument(
        "--matrix",
        nargs=3,
        type=int,
        default=[2, 2, 2],
        help="Scaling matrix(3 integers, default: 2 2 2)"
    )
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Batch process all CIF files in current directory"
    )
    parser.add_argument(
        "--api-key",
        help="Materials Project API key (online)"
    )
    parser.add_argument(
        "--mp-data",
        help="Path to offline MP JSON file (offline)"
    )
    parser.add_argument(
        "-cor",
        action="store_true",
        help="Use MP2020 energy compatibility corrections"
    )
    parser.add_argument(
        "--post",
        action="store_true",
        help="Run Ehull post-processing"
    )
    parser.add_argument(
        "-t",
        type=float,
        default=0.2,
        dest="threshold",
        help="Ehull threshold in eV/atom for filtering (default: 0.2)"
    )
    parser.add_argument(
        "-N",
        action="store_true",
        dest="no_hull",
        help="Skip hull calculation, directly extract all CIF files"
    )
    parser.add_argument(
        "--cif-dir",
        default="0_cif",
        help="Output directory(default: 0_cif)"
    )
    parser.add_argument(
        "-o",
        metavar="FILE",
        dest="output",
        help="Output file"
    )


def execute(args, remaining_args=None):
    """Execute the PT workflow command"""
    if args.mode == "Ehull":
        if args.api_key or args.mp_data or getattr(args, 'no_hull', False) or args.post:
            _execute_ehull(args)
            return

    if args.mode == "viz" and hasattr(args, 'run') and args.run:
        _execute_viz(args, remaining_args=remaining_args)
        return
    
    if hasattr(args, 'run') and args.run:
        if args.mode in DIRECT_EXEC_MODULES:
            module_name = DIRECT_EXEC_MODULES[args.mode]
            print(f"Executing {args.mode} in current directory...")
            
            cmd_args = [sys.executable, "-m", module_name]
            if args.mode == "supercell":
                if args.batch:
                    cmd_args.append("--batch")
                elif args.input:
                    cmd_args.extend(["-i", args.input])
                cmd_args.extend(["--matrix"] + [str(x) for x in args.matrix])
            if remaining_args:
                cmd_args.extend(remaining_args)
            
            try:
                result = subprocess.run(cmd_args, check=True)
                return
            except subprocess.CalledProcessError as e:
                print(f"Error executing {args.mode}: {e}")
                sys.exit(1)
        else:
            print(f"Warning: Mode '{args.mode}' does not support direct execution (-r).")
            print(f"Supported modes: {', '.join(DIRECT_EXEC_MODULES.keys())}")
            print("Falling back to copy mode...")
    
    from .base import BaseWorkflowCommand
    
    class PTCommand(BaseWorkflowCommand):
        @property
        def repository_path(self):
            return PT_REPOSITORY
        
        @property
        def file_map(self):
            return PT_FILE_MAP
        
        @property
        def workflow_name(self):
            return "PT"

        def get_common_files(self, mode: str) -> list:
            """PT viz uses its own script, skip common viz files."""
            if mode == "viz":
                return []
            return super().get_common_files(mode)
        
        def execute(self, args):
            """Override to handle CIFGEN_MODULE files for mutate mode"""
            super().execute(args)
            
            cifgen_files = CIFGEN_FILES.get(args.mode, [])
            if cifgen_files:
                print(f" Copying files from cifgen module...")
                for filename in cifgen_files:
                    src = os.path.join(CIFGEN_MODULE, filename)
                    dst = os.path.join(args.dest, filename)
                    if os.path.exists(src):
                        shutil.copy2(src, dst)
                        print(f"   [OK] {filename}")
                    else:
                        print(f"   [FAIL] {filename} (not found)")
    
    cmd = PTCommand()
    cmd.execute(args)


def _execute_ehull(args):
    """
    Execute Ehull calculation or post-processing directly.
    
    Usage:
        gewum PT --mode Ehull --api-key mp-xxxxx           # Online mode (no correction)
        gewum PT --mode Ehull --api-key mp-xxxxx -cor      # Online mode (with correction)
        gewum PT --mode Ehull --mp-data /path/to/MP.json   # Offline mode (no correction)
        gewum PT --mode Ehull --mp-data /path/to/MP.json -cor  # Offline mode (with correction)
        gewum PT --mode Ehull --post -t 0.2                # Post-processing
        gewum PT --mode Ehull -N                           # Skip hull, extract all
    """
    if hasattr(args, 'no_hull') and args.no_hull:
        post_module = "gewum.src.common.ehull.Ehull_post"
        cmd_args = [sys.executable, "-m", post_module, "-N"]
        if args.cif_dir:
            cmd_args.extend(["--cif-dir", args.cif_dir])
        if args.input:
            cmd_args.extend(["--raw-input", args.input])
        print(f"Extracting all CIF files (skipping hull calculation)...")
        try:
            subprocess.run(cmd_args, check=True)
            return
        except subprocess.CalledProcessError as e:
            print(f"Error: {e}")
            sys.exit(1)
    
    if hasattr(args, 'post') and args.post:
        post_module = "gewum.src.common.ehull.Ehull_post"
        cmd_args = [sys.executable, "-m", post_module]
        cmd_args.extend(["-t", str(args.threshold)])
        if args.cif_dir:
            cmd_args.extend(["--cif-dir", args.cif_dir])
        if args.input:
            cmd_args.extend(["-i", args.input])
        if args.output:
            cmd_args.extend(["-o", args.output])
        print(f"Running Ehull post-processing (threshold: {args.threshold} eV/atom)...")
        try:
            subprocess.run(cmd_args, check=True)
            return
        except subprocess.CalledProcessError as e:
            print(f"Error: {e}")
            sys.exit(1)
    
    if not args.api_key and not args.mp_data:
        print("Error: Either --api-key or --mp-data is required for Ehull calculation.")
        print("Usage examples:")
        print("  gewum PT --mode Ehull --api-key mp-xxxxx           # Online mode (no correction)")
        print("  gewum PT --mode Ehull --api-key mp-xxxxx -cor      # Online mode (with correction)")
        print("  gewum PT --mode Ehull --mp-data /path/to/MP.json   # Offline mode (no correction)")
        print("  gewum PT --mode Ehull --mp-data /path/to/MP.json -cor  # Offline mode (with correction)")
        print("  gewum PT --mode Ehull --post -t 0.2                # Post-process results")
        print("  gewum PT --mode Ehull -N                           # Skip hull, extract all CIFs")
        sys.exit(1)
    
    module_name = "gewum.src.common.ehull.Ehull_compatibility"
    if hasattr(args, 'cor') and args.cor:
        print("Using Ehull with MP2020 compatibility corrections...")
    else:
        print("Using Ehull without compatibility corrections...")
    
    cmd_args = [sys.executable, "-m", module_name]
    if args.api_key:
        cmd_args.extend(["--api-key", args.api_key])
        mode_info = "online"
    elif args.mp_data:
        cmd_args.extend(["--mp-data", args.mp_data])
        mode_info = "offline"
    if hasattr(args, 'cor') and args.cor:
        cmd_args.append("-cor")
    if args.input:
        cmd_args.extend(["-i", args.input])
    if args.output:
        cmd_args.extend(["-o", args.output])
    
    display_args = ["python", "-m", module_name]
    if args.api_key:
        display_args.extend(["--api-key", "******"])
    if args.mp_data:
        display_args.extend(["--mp-data", args.mp_data])
    print(f"Running ({mode_info} mode): {' '.join(display_args)}")
    
    try:
        subprocess.run(cmd_args, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        sys.exit(1)


def _execute_viz(args, remaining_args=None):
    """
    Execute PT visualization mode by forwarding args to viz_pt module.

    Usage:
        gewum PT --mode viz -r -- --cif-dir . --plot all
        gewum PT --mode viz -r -- --cif-dir . --plot sankey --ignore-p1
        gewum PT --mode viz -r -- --cif-dir . --plot violin --violin-top-n 15
    """
    module_name = DIRECT_EXEC_MODULES["viz"]
    cmd_args = [sys.executable, "-m", module_name]

    if remaining_args:
        cmd_args.extend(remaining_args)

    print(f"Running viz: {' '.join(cmd_args)}")

    try:
        subprocess.run(cmd_args, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing viz: {e}")
        sys.exit(1)
