import os
import numpy as np
import yaml
from ase.io import read
from mattersim.forcefield.potential import MatterSimCalculator
from gewum.src.phonon_src.ph.phonon import PhononWorkflow


def main():
    import argparse
    parser = argparse.ArgumentParser(description='QHA Phonon Calculation')
    parser.add_argument('--input', '-i', default='CONTCAR', help='Input structure file')
    parser.add_argument('--device', default='cpu', help='Device for calculation')
    parser.add_argument('--amplitude', type=float, default=0.02, help='Displacement amplitude')
    parser.add_argument('--mesh', type=str, default='7,7,7', help='Mesh for phonon calculation')
    parser.add_argument('--t-max', type=int, default=2001, help='Maximum temperature')
    args = parser.parse_args()
    
    atoms = read(args.input)
    calc = MatterSimCalculator(load_path="MatterSim-v1.0.0-5M.pth", device=args.device)
    atoms.calc = calc
    
    ph = PhononWorkflow(
        atoms=atoms,
        find_prim=False,
        work_dir="phonon",
        amplitude=args.amplitude,
        supercell_matrix=np.diag([1, 1, 1]),
    )
    
    has_imag, phonon_obj = ph.run()
    
    mesh = [int(x) for x in args.mesh.split(',')]
    t_min = 0
    t_max = args.t_max
    t_step = 2
    
    phonon_obj.run_mesh(mesh, with_eigenvectors=False)
    phonon_obj.run_thermal_properties(t_min=t_min, t_max=t_max, t_step=t_step)
    
    tp_dict = phonon_obj.get_thermal_properties_dict()
    temperatures = tp_dict['temperatures']
    free_energy = tp_dict['free_energy']
    entropy = tp_dict['entropy']
    heat_capacity = tp_dict['heat_capacity']
    energy = tp_dict['internal_energy'] if 'internal_energy' in tp_dict else None
    
    natom = len(phonon_obj.primitive) if hasattr(phonon_obj, 'primitive') else len(atoms)
    cutoff_frequency = getattr(phonon_obj, 'cutoff_frequency', 0.0)
    num_modes = 3 * natom
    zero_point_energy = free_energy[0] if len(free_energy) > 0 else 0.0
    
    thermal_properties_data = {
        'unit': {
            'temperature': 'K',
            'free_energy': 'kJ/mol',
            'entropy': 'J/K/mol',
            'heat_capacity': 'J/K/mol'
        },
        'natom': natom,
        'cutoff_frequency': float(cutoff_frequency),
        'num_modes': num_modes,
        'zero_point_energy': float(zero_point_energy),
        'thermal_properties': []
    }
    
    for i, t in enumerate(temperatures):
        prop_data = {
            'temperature': float(t),
            'free_energy': float(free_energy[i]),
            'entropy': float(entropy[i]),
            'heat_capacity': float(heat_capacity[i])
        }
        
        if energy is not None and i < len(energy):
            prop_data['energy'] = float(energy[i])
        
        thermal_properties_data['thermal_properties'].append(prop_data)
    
    scale = os.path.basename(os.getcwd())
    scales = [0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05]
    idx = scales.index(float(scale)) + 1
    os.makedirs("../CP", exist_ok=True)
    
    output_file = f"../CP/thermal_properties-{idx}.yaml"
    with open(output_file, 'w') as f:
        yaml.dump(thermal_properties_data, f, default_flow_style=False, sort_keys=False)
    
    print(f"Phonon calculation for scale {scale} done. Thermal properties YAML saved to {output_file}.")


if __name__ == "__main__":
    main()