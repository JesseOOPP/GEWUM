"""
GEWUM Phonon Calculation Module using MatterSim
Calculate phonon dispersion and check for imaginary modes
"""
import os
import numpy as np
from ase.io import read
from mattersim.forcefield.potential import MatterSimCalculator
from gewum.src.phonon_src.ph.phonon import PhononWorkflow
import warnings
import glob
import argparse
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.bandstructure import HighSymmKpath

warnings.filterwarnings("ignore", category=DeprecationWarning)


def get_auto_band_path(atoms):
    """
    Automatically generate high-symmetry band path from crystal structure.
    
    Args:
        atoms: ASE Atoms object
    
    Returns:
        List of (label, coordinates) tuples for band path
    """
    try:
        structure = AseAtomsAdaptor.get_structure(atoms)
        kpath = HighSymmKpath(structure)
        
        path_labels = kpath.kpath['path']
        kpoints = kpath.kpath['kpoints']
        
        band_path = []
        for segment in path_labels:
            for label in segment:
                coords = kpoints[label]
                display_label = "G" if label == "\\Gamma" or label == "GAMMA" else label
                band_path.append((display_label, list(coords)))
        
        cleaned_path = [band_path[0]]
        for i in range(1, len(band_path)):
            if band_path[i][0] != band_path[i-1][0]:
                cleaned_path.append(band_path[i])
        
        return cleaned_path
        
    except Exception as e:
        print(f"Warning: Auto band path failed ({e}), using default cubic path")
        return [
            ("G", [0, 0, 0]),
            ("X", [0.5, 0, 0]),
            ("M", [0.5, 0.5, 0]),
            ("R", [0.5, 0.5, 0.5]),
            ("G", [0, 0, 0])
        ]


def calculate_phonon(cif_file, work_directory="./tmp", model_path="MatterSim-v1.0.0-5M.pth",
                     device="cpu", supercell_matrix=None, band_path=None, 
                     band_npoints=201, amplitude=0.02, cleanup_yaml=True, auto_path=True,
                     band_color='purple', calc=None):
    """
    Calculate phonon dispersion for a CIF structure.
    
    Args:
        cif_file: Path to input CIF file
        work_directory: Working directory for intermediate files
        model_path: Path to MatterSim model
        device: Computation device ('cpu' or 'cuda')
        supercell_matrix: Supercell matrix for phonon calculation
        band_path: Custom band path for dispersion (overrides auto_path)
        band_npoints: Number of points along band path
        amplitude: Displacement amplitude for force constants
        cleanup_yaml: Whether to delete intermediate YAML files
        auto_path: Automatically determine band path from structure symmetry
        band_color: Color of the band dispersion line (default: purple)
        calc: Pre-loaded MatterSim calculator to reuse across CIFs. If None,
              a new calculator is created from model_path/device.
    
    Returns:
        Tuple of (has_imaginary_modes, phonon_frequencies)
    """
    if supercell_matrix is None:
        supercell_matrix = np.diag([2, 2, 2])
    
    struc = read(cif_file)
    if calc is None:
        calc = MatterSimCalculator(load_path=model_path, device=device)
    struc.calc = calc
    
    if band_path is None:
        if auto_path:
            band_path = get_auto_band_path(struc)
            print(f"  Auto band path: {' -> '.join([p[0] for p in band_path])}")
        else:
            band_path = [
                ("G", [0, 0, 0]),
                ("X", [0.5, 0, 0]),
                ("M", [0.5, 0.5, 0]),
                ("R", [0.5, 0.5, 0.5]),
                ("G", [0, 0, 0])
            ]
    
    base_name = os.path.splitext(os.path.basename(cif_file))[0]
    ph_work_dir = os.path.join(work_directory, base_name)
    
    ph = PhononWorkflow(
        atoms=struc,
        find_prim=True,
        work_dir=ph_work_dir,
        amplitude=amplitude,
        supercell_matrix=supercell_matrix,
        band_path=band_path,
        band_npoints=band_npoints,
        band_color=band_color
    )
    
    has_imag, phonons = ph.run()
    
    if cleanup_yaml:
        yaml_files = glob.glob(os.path.join(ph_work_dir, "*.yaml"))
        for yaml_file in yaml_files:
            try:
                os.remove(yaml_file)
            except OSError as e:
                print(f"Error deleting {yaml_file}: {e}")
    
    return has_imag, phonons


def main():
    """Main function for batch phonon calculation"""
    parser = argparse.ArgumentParser(description='Calculate phonon dispersion using MatterSim')
    parser.add_argument('--work-dir', '-w', default='./tmp',
                        help='Working directory for intermediate files')
    parser.add_argument('--model', '-m', default='MatterSim-v1.0.0-5M.pth',
                        help='Path to MatterSim model')
    parser.add_argument('--device', '-d', default='cpu', choices=['cpu', 'cuda'],
                        help='Computation device')
    parser.add_argument('--supercell', '-s', type=int, nargs=3, default=[2, 2, 2],
                        help='Supercell dimensions (default: 2 2 2)')
    parser.add_argument('--no-auto-path', action='store_true',
                        help='Disable auto band path detection, use default cubic path')
    parser.add_argument('--band-path', type=str, default=None,
                        help='Custom band path, format: "G:0,0,0|X:0.5,0,0|M:0.5,0.5,0"')
    parser.add_argument('--cif-dir', default='.',
                        help='Directory containing CIF files')
    parser.add_argument('--band-color', default='purple',
                        help='Color of the band dispersion line (default: purple)')
    args = parser.parse_args()
    
    supercell_matrix = np.diag(args.supercell)
    
    custom_band_path = None
    if args.band_path:
        try:
            custom_band_path = []
            for point in args.band_path.split('|'):
                label, coords = point.split(':')
                coords = [float(x) for x in coords.split(',')]
                custom_band_path.append((label.strip(), coords))
            print(f"Using custom band path: {' -> '.join([p[0] for p in custom_band_path])}")
        except Exception as e:
            print(f"Error parsing band path: {e}")
            print('Format: "G:0,0,0|X:0.5,0,0|M:0.5,0.5,0"')
            return
    
    auto_path = not args.no_auto_path
    
    cif_files = [f for f in os.listdir(args.cif_dir) if f.endswith('.cif')]
    
    if not cif_files:
        print(f"No CIF files found in {args.cif_dir}")
        return
    
    print(f"Found {len(cif_files)} CIF files to process")
    
    # Load the MatterSim model once and reuse it across all CIF files
    calc = MatterSimCalculator(load_path=args.model, device=args.device)
    
    results = []
    for cif_file in cif_files:
        cif_path = os.path.join(args.cif_dir, cif_file)
        print(f"Processing: {cif_file}")
        
        try:
            has_imag, phonons = calculate_phonon(
                cif_path,
                work_directory=args.work_dir,
                model_path=args.model,
                device=args.device,
                supercell_matrix=supercell_matrix,
                band_path=custom_band_path,
                auto_path=auto_path,
                band_color=args.band_color,
                calc=calc
            )
            
            print(f"  Has imaginary phonon: {has_imag}")
            print(f"  Phonon frequencies: {phonons}")
            results.append((cif_file, has_imag, phonons))
            
        except Exception as e:
            print(f"  Error: {e}")
            results.append((cif_file, None, str(e)))
    
    print("\n" + "=" * 60)
    print("PHONON CALCULATION SUMMARY")
    print("=" * 60)
    stable = sum(1 for r in results if r[1] == False)
    unstable = sum(1 for r in results if r[1] == True)
    failed = sum(1 for r in results if r[1] is None)
    print(f"Stable (no imaginary modes): {stable}")
    print(f"Unstable (has imaginary modes): {unstable}")
    print(f"Failed: {failed}")


if __name__ == "__main__":
    main()
