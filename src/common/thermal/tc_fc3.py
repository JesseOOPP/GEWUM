"""
GEWUM Thermal Conductivity Calculation Entry Point
Computes third-order force constants (FC3) using MatterSim + phono3py
"""
import os
import glob
import argparse
import numpy as np
from ase.io import read
from mattersim.forcefield.potential import MatterSimCalculator
from gewum.src.common.thermal.fc3_workflow import FC3CalculationWorkflow


def main():
    parser = argparse.ArgumentParser(
        description='Calculate FC3 for thermal conductivity using MatterSim'
    )
    parser.add_argument(
        '--cif', '-c',
        default=None,
        help='Input CIF file (default: first .cif in current directory)'
    )
    parser.add_argument(
        '--supercell', '-s',
        type=int, nargs=3,
        default=[2, 2, 2],
        help='Supercell matrix diagonal (default: 2 2 2)'
    )
    parser.add_argument(
        '--device', '-d',
        default='cpu',
        choices=['cpu', 'cuda'],
        help='Device for MatterSim (default: cpu)'
    )
    parser.add_argument(
        '--work-dir', '-w',
        default='./',
        help='Working directory for output (default: ./)'
    )
    parser.add_argument(
        '--find-prim',
        action='store_true',
        help='Find primitive cell before calculation'
    )
    args = parser.parse_args()
    
    if args.cif:
        cif_file = args.cif
    else:
        cif_files = glob.glob("*.cif")
        if not cif_files:
            raise FileNotFoundError("No .cif file found in current directory")
        cif_file = cif_files[0]
    
    print("=" * 60)
    print("GEWUM Thermal Conductivity - FC3 Calculation")
    print("=" * 60)
    print(f"  CIF file: {cif_file}")
    print(f"  Supercell: {args.supercell}")
    print(f"  Device: {args.device}")
    print("=" * 60)
    
    atoms = read(cif_file)
    
    model_path = "MatterSim-v1.0.0-1M.pth"
    calculator = MatterSimCalculator(load_path=model_path, device=args.device)
    atoms.calc = calculator
    
    workflow = FC3CalculationWorkflow(
        atoms=atoms,
        work_dir=args.work_dir,
        supercell_matrix=np.diag(args.supercell),
        find_prim=args.find_prim
    )
    
    workflow.run()
    
    print("\n" + "=" * 60)
    print("FC3 calculation completed!")
    print(f"Output files: phono3py_disp.yaml, FORCES_FC3")
    print("Next step: Run tc_post.sh to compute thermal conductivity")
    print("=" * 60)


if __name__ == "__main__":
    main()
