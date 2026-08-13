"""GEWUM RD (Random Design) Workflow Command
Provides scripts for random crystal structure generation and processing
"""
import os
import shutil
import subprocess
import sys
from ..config import RD_REPOSITORY

RD_FILE_MAP = {
    "cifgen": ["cifgen.sh"],
    "relax": ["relax_umlip.sh"],
    "refine": ["refine_umlip.sh"],
    "calc": ["calc_energy.sh"],
    "post": ["post_relax.sh"],
    "relaxhp": ["relax_umlip_hp.sh"],
    "posthp": ["post_relax_hp.sh"],
    "reorder": [],
    "dedup": [],
    "Ehull": [],
    "sym": [],
    "sym2d": [],
    "chiral": [],
    "select": [],
    "ph": [],
    "phpost": [],
    "viz": [],
    "viz2": [],
    "viz3": [],
    "auto": ["cifgen.sh", "relax_umlip.sh", "post_relax.sh", "run_srss.sh"],
}

DIRECT_EXEC_MODULES = {
    "sym": "gewum.src.common.postprocess.sym_rename",
    "sym2d": "gewum.src.common.postprocess.layergroup_sym",
    "chiral": "gewum.src.common.postprocess.chiral_check",
    "phpost": "gewum.src.common.phonon.ph_post",
    "reorder": "gewum.src.RDworkflows.reorder_energy",
    "dedup": "gewum.src.RDworkflows.cif_dedup",
    "viz": "gewum.src.common.postprocess.visualization",
    "viz2": "gewum.src.common.postprocess.viz2_analysis",
    "viz3": "gewum.src.common.postprocess.phase_diagram",
}

EHULL_MODES = ["Ehull"]

MODE_DESCRIPTIONS = """
  cifgen: Random crystal structure generation 
  relax: Structure relaxation 
  refine: Cell relaxation for pre-relaxed structures
  calc: Single-point energy calculation
  post: Post-processing 
  relaxhp: High pressure structure relaxation 
  posthp: High pressure post-processing 
  reorder: Reorder energy_final.txt for specific composition (-r)
  dedup: Deduplicate CIF structures by composition-aware structure matching (-r, --rdf --sim-t 0.2)
  Ehull: Calculate Ehull (-r, -cor, --mp-data,--api-key, --post -t 0.2)
  sym: Symmetrize and classify 3D structures by space group (-r)
  sym2d: Classify 2D structures by layer group symmetry (-r)
  chiral: Classify 3D chiral structures (-r) 
  select: Unified structure diversity selection
  ph: Phonon dispersion calculation
  phpost: Phonon post-processing with ph.log (-r)
  viz: Structure visualization with UMAP/t-SNE. 
  viz2: Advanced structure analysis & visualization. 
  viz3: Formation energy phase diagram. 
  auto: Automated full pipeline (cifgen-select-relax-post-Ehull). 
"""


def setup_args(parser):
    """Configure command arguments"""
    all_modes = list(RD_FILE_MAP.keys())
    
    parser.add_argument(
        "--mode", 
        required=True,
        choices=all_modes,
        help=f"GEWUM RD mode to copy:\n{MODE_DESCRIPTIONS}"
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
        help="Directly execute (sym, sym2d, chiral, Ehull, phpost, reorder, dedup, viz3)"
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
        help="Ehull post-processing "
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
        "--self-hull",
        action="store_true",
        dest="self_hull",
        help="Self-hull mode: build convex hull from input structures only "
             "(no MP competing phases). MP data used only for elemental references."
    )
    parser.add_argument(
        "--element-dir",
        dest="element_dir",
        default=None,
        help="Directory with per-element CIFs (<Symbol>.cif) for uMLIP-computed "
             "elemental references in self-hull mode. Elements not found here "
             "fall back to MP offline data."
    )
    parser.add_argument(
        "--cif-dir",
        nargs='+',
        default=None,
        help="CIF directories. For Ehull: single dir (default: 0_cif)."
    )
    parser.add_argument(
        "-i",
        metavar="FILE",
        dest="input",
        help="Input file"
    )
    parser.add_argument(
        "-o",
        metavar="FILE",
        dest="output",
        help="Output file"
    )



def execute(args, remaining_args=None):
    """Execute the RD workflow command"""
    if args.mode == "Ehull":
        if args.api_key or args.mp_data or getattr(args, 'no_hull', False) or args.post:
            _execute_ehull(args)
            return

    if args.mode in ("viz", "viz2", "viz3") and hasattr(args, 'run') and args.run:
        _execute_viz(args, remaining_args=remaining_args)
        return

    if hasattr(args, 'run') and args.run:
        if args.mode in DIRECT_EXEC_MODULES:
            module_name = DIRECT_EXEC_MODULES[args.mode]
            print(f"Executing {args.mode} in current directory...")
            cmd_args = [sys.executable, "-m", module_name]
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
    
    class RDCommand(BaseWorkflowCommand):
        @property
        def repository_path(self):
            return RD_REPOSITORY
        
        @property
        def file_map(self):
            return RD_FILE_MAP
        
        @property
        def workflow_name(self):
            return "RD"
    
    cmd = RDCommand()
    cmd.execute(args)


def _execute_viz(args, remaining_args=None):
    """Execute visualization mode by forwarding args to visualization module.

    Usage:
        gewum RD --mode viz -r -- --dim 3d --total-dir ./total_cifs
        gewum RD --mode viz -r -- --dim 2d --descriptor soap --total-dir ./total --select-dir ./selected
        gewum RD --mode viz -r -- --dim 3d --load-features ./features.npz --reduction tsne
        gewum RD --mode viz3 -r --mp-data /path/to/MP.json
    """
    module_name = DIRECT_EXEC_MODULES[args.mode]
    cmd_args = [sys.executable, "-m", module_name]

    # For viz3, forward recognized RD-level args directly (no -- needed)
    if args.mode == "viz3":
        if getattr(args, 'mp_data', None):
            cmd_args.extend(["--mp-data", args.mp_data])
        if getattr(args, 'input', None):
            cmd_args.extend(["-i", args.input])
        if getattr(args, 'output', None):
            cmd_args.extend(["--output", args.output])

    if remaining_args:
        cmd_args.extend(remaining_args)

    print(f"Running {args.mode}: {' '.join(cmd_args)}")

    try:
        subprocess.run(cmd_args, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing {args.mode}: {e}")
        sys.exit(1)


def _execute_ehull(args):
    """
    Execute Ehull calculation or post-processing directly.
    
    Usage:
        gewum RD --mode Ehull --api-key mp-xxxxx           # Online mode (no correction)
        gewum RD --mode Ehull --api-key mp-xxxxx -cor      # Online mode (with correction)
        gewum RD --mode Ehull --mp-data /path/to/MP.json   # Offline mode (no correction)
        gewum RD --mode Ehull --mp-data /path/to/MP.json -cor  # Offline mode (with correction)
        gewum RD --mode Ehull --post -t 0.2                # Post-processing
        gewum RD --mode Ehull -N                           # Skip hull, extract all
    """
    if hasattr(args, 'no_hull') and args.no_hull:
        post_module = "gewum.src.common.ehull.Ehull_post"
        cmd_args = [sys.executable, "-m", post_module, "-N"]
        cif_dir_val = args.cif_dir[0] if args.cif_dir else "0_cif"
        cmd_args.extend(["--cif-dir", cif_dir_val])
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
            cmd_args.extend(["--cif-dir", args.cif_dir[0]])
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
        print("  gewum RD --mode Ehull --api-key mp-xxxxx           # Online mode (no correction)")
        print("  gewum RD --mode Ehull --api-key mp-xxxxx -cor      # Online mode (with correction)")
        print("  gewum RD --mode Ehull --mp-data /path/to/MP.json   # Offline mode (no correction)")
        print("  gewum RD --mode Ehull --mp-data /path/to/MP.json -cor  # Offline mode (with correction)")
        print("  gewum RD --mode Ehull --post -t 0.2                # Post-process results")
        print("  gewum RD --mode Ehull -N                           # Skip hull, extract all CIFs")
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
    if hasattr(args, 'self_hull') and args.self_hull:
        cmd_args.append("--self-hull")
    if getattr(args, 'element_dir', None):
        cmd_args.extend(["--element-dir", args.element_dir])
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
