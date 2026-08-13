"""
GEWUM MD NVT Module
Molecular dynamics simulation with NVT ensemble using MatterSim
Supports supercell expansion and Langevin thermostat
"""
import os
import sys
import logging
import argparse
import numpy as np
from ase.io import read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase import units
from ase.build import make_supercell
from ase.io.trajectory import Trajectory
from mattersim.forcefield.potential import MatterSimCalculator

logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s'
)


class MDLogger:
    """Custom MD logger for energy and temperature tracking"""
    
    def __init__(self, atoms, logfile, header=True):
        self.atoms = atoms
        self.logfile = logfile
        self.dyn = None
        self._file = open(logfile, 'w')
        if header:
            self._file.write("Step,Time_ps,Temperature_K,Potential_Energy_eV,"
                             "Kinetic_Energy_eV,Total_Energy_eV\n")
    
    def __call__(self):
        epot = self.atoms.get_potential_energy()
        ekin = self.atoms.get_kinetic_energy()
        temp = ekin / (1.5 * len(self.atoms) * units.kB)
        time_ps = self.dyn.nsteps * (self.dyn.dt / units.fs) / 1000.0
        self._file.write(f"{self.dyn.nsteps},{time_ps:.4f},{temp:.2f},"
                         f"{epot:.6f},{ekin:.6f},{epot + ekin:.6f}\n")
        self._file.flush()

    def close(self):
        """Close the log file handle."""
        if self._file and not self._file.closed:
            self._file.close()


def create_supercell(atoms, supercell_size):
    """
    Create supercell from unit cell
    
    Args:
        atoms: ASE atoms object
        supercell_size: Tuple of (nx, ny, nz) expansion factors
    
    Returns:
        Expanded supercell atoms object
    """
    P = np.diag(supercell_size)
    supercell = make_supercell(atoms, P)
    
    logging.info(f"Supercell: {supercell_size[0]}x{supercell_size[1]}x{supercell_size[2]}")
    logging.info(f"Atoms: {len(atoms)} -> {len(supercell)}")
    
    return supercell


def run_nvt_md(cif_file, output_dir, temperature=300, timestep=1.0, steps=10000,
               dump_interval=100, log_interval=1, friction=0.01,
               supercell_size=(2, 2, 1), model_path="mattersim-v1.0.0-5M.pth", device="cpu"):
    """
    Run NVT molecular dynamics simulation
    
    Args:
        cif_file: Input CIF file path
        output_dir: Output directory for results
        temperature: Temperature in Kelvin (default: 300)
        timestep: Time step in femtoseconds (default: 1.0)
        steps: Total number of MD steps (default: 10000)
        dump_interval: Trajectory output interval (default: 100)
        log_interval: Log output interval (default: 1)
        friction: Langevin friction coefficient (default: 0.01)
        supercell_size: Supercell expansion factors (nx, ny, nz)
        model_path: Path to MatterSim model (default: mattersim-v1.0.0-5M.pth)
        device: Computation device 'cpu' or 'cuda' (default: 'cpu')
    
    Returns:
        bool: True if simulation completed successfully
    """
    logger = None
    try:
        atoms = read(cif_file)
        original_atoms = len(atoms)
        base_name = os.path.basename(cif_file).replace('.cif', '')
        
        logging.info(f"Starting NVT MD for: {base_name}")
        logging.info(f"Original structure: {original_atoms} atoms")
        
        if any(s > 1 for s in supercell_size):
            atoms = create_supercell(atoms, supercell_size)
        
        if model_path:
            calc = MatterSimCalculator(load_path=model_path, device=device)
        else:
            calc = MatterSimCalculator(device=device)
        atoms.calc = calc
        
        MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)
        
        os.makedirs(output_dir, exist_ok=True)
        
        supercell_file = os.path.join(output_dir, f"{base_name}_supercell.cif")
        write(supercell_file, atoms)
        logging.info(f"Supercell saved: {supercell_file}")
        
        dyn = Langevin(atoms, timestep * units.fs, 
                      temperature_K=temperature, 
                      friction=friction)
        
        log_file = os.path.join(output_dir, f"{base_name}_md.log")
        logger = MDLogger(atoms, log_file)
        logger.dyn = dyn
        dyn.attach(logger, interval=log_interval)
        
        traj_file = os.path.join(output_dir, f"{base_name}_md.traj")
        traj_writer = Trajectory(traj_file, 'w', atoms)
        dyn.attach(traj_writer.write, interval=dump_interval)
        
        total_time_ps = steps * timestep / 1000.0
        logging.info(f"Temperature: {temperature} K")
        logging.info(f"Timestep: {timestep} fs, Steps: {steps}")
        logging.info(f"Total simulation time: {total_time_ps:.2f} ps")
        logging.info(f"Friction coefficient: {friction}")
        
        dyn.run(steps)
        
        final_file = os.path.join(output_dir, f"{base_name}_final.cif")
        write(final_file, atoms)
        
        traj_writer.close()
        
        logging.info(f"MD simulation completed: {base_name}")
        logging.info(f"Final structure: {final_file}")
        logging.info(f"Trajectory: {traj_file}")
        logging.info(f"Energy log: {log_file}")
        
        return True
        
    except Exception as e:
        logging.error(f"MD simulation failed for {cif_file}: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return False
    finally:
        if logger is not None:
            logger.close()


def main():
    """Command line interface for NVT MD simulation"""
    parser = argparse.ArgumentParser(
        description='GEWUM NVT Molecular Dynamics using MatterSim',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument('cif_file', help='Input CIF file path')
    parser.add_argument('output_dir', help='Output directory for results')
    
    parser.add_argument('--temperature', '-T', type=float, default=300,
                        help='Temperature in Kelvin (default: 300)')
    parser.add_argument('--timestep', '-dt', type=float, default=1.0,
                        help='Time step in femtoseconds (default: 1.0)')
    parser.add_argument('--steps', '-n', type=int, default=10000,
                        help='Number of MD steps (default: 10000)')
    parser.add_argument('--dump-interval', type=int, default=100,
                        help='Trajectory dump interval (default: 100)')
    parser.add_argument('--log-interval', type=int, default=1,
                        help='Log output interval (default: 1)')
    parser.add_argument('--friction', type=float, default=0.01,
                        help='Langevin friction coefficient (default: 0.01)')
    parser.add_argument('--supercell', type=int, nargs=3, default=[2, 2, 2],
                        metavar=('NX', 'NY', 'NZ'),
                        help='Supercell expansion factors (default: 2 2 2)')
    parser.add_argument('--model', type=str, default='mattersim-v1.0.0-5M.pth',
                        help='Path to MatterSim model file (default: mattersim-v1.0.0-5M.pth)')
    parser.add_argument('--device', type=str, default='cpu',
                        choices=['cpu', 'cuda'],
                        help='Computation device (default: cpu)')
    
    args = parser.parse_args()
    
    success = run_nvt_md(
        cif_file=args.cif_file,
        output_dir=args.output_dir,
        temperature=args.temperature,
        timestep=args.timestep,
        steps=args.steps,
        dump_interval=args.dump_interval,
        log_interval=args.log_interval,
        friction=args.friction,
        supercell_size=tuple(args.supercell),
        model_path=args.model,
        device=args.device
    )
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
