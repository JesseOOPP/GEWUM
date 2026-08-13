import sys
import numpy as np
from ase.io import read, write
from ase.optimize import BFGS
from ase.constraints import ExpCellFilter
from mattersim.forcefield import MatterSimCalculator

def relax(atoms, calc, pressure=0, maxstep=0.1, fmax=0.005, steps=1000):
    atoms.calc = calc
    mask = [False, False, False, False, False, False]  
    ucf = ExpCellFilter(atoms, scalar_pressure=pressure, mask=mask)
    dyn = BFGS(ucf, maxstep=maxstep)
    dyn.run(fmax=fmax, steps=steps)
    return atoms

def main():
    import argparse
    parser = argparse.ArgumentParser(description='QHA Volume-scaled Relaxation')
    parser.add_argument('scale', type=float, help='Volume scale factor')
    parser.add_argument('--input', '-i', default='POSCAR', help='Input POSCAR file')
    parser.add_argument('--output', '-o', default='CONTCAR', help='Output CONTCAR file')
    parser.add_argument('--fmax', type=float, default=0.005, help='Force convergence threshold')
    parser.add_argument('--device', default='cpu', help='Device for calculation')
    args = parser.parse_args()
    
    atoms = read(args.input)
    atoms.set_cell(atoms.cell * args.scale, scale_atoms=True)
    
    calc = MatterSimCalculator(load_path="MatterSim-v1.0.0-5M.pth", device=args.device)
    relaxed = relax(atoms, calc, pressure=0.0, fmax=args.fmax)
    
    write(args.output, relaxed, format="vasp")
    energy = relaxed.get_potential_energy()
    volume = relaxed.get_volume()
    
    with open("energy.dat", "w") as f:
        f.write(f"{volume} {energy}\n")
    
    print(f"Scale {args.scale:.2f}: Volume = {volume:.4f} A3, Energy = {energy:.6f} eV")


if __name__ == "__main__":
    main()
