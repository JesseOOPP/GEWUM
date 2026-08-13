"""
GEWUM MD Post-processing Module
Processes MD simulation results and generates energy-time plots
"""
import os
import sys
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from pymatgen.core import Structure

plt.switch_backend('Agg')


def process_single_folder(folder, output_dir='./md_plots'):
    """
    Process a single MD output folder
    
    Args:
        folder: Path to MD output folder
        output_dir: Directory to save plots
    
    Returns:
        dict: Processing result with status and data
    """
    result = {
        'folder': folder,
        'success': False,
        'message': '',
        'data': None
    }
    
    log_files = glob.glob(os.path.join(folder, "*_md.log"))
    if not log_files:
        result['message'] = "No MD log file found"
        return result
    
    log_file = log_files[0]
    
    cif_files = glob.glob(os.path.join(folder, "*_final.cif"))
    if not cif_files:
        cif_files = glob.glob(os.path.join(folder, "*_supercell.cif"))
    
    if not cif_files:
        result['message'] = "No CIF file found for atom count"
        return result
    
    cif_file = cif_files[0]
    
    try:
        structure = Structure.from_file(cif_file)
        num_atoms = structure.num_sites
        
        df = pd.read_csv(log_file)
        
        required_cols = ['Time_ps', 'Potential_Energy_eV', 'Temperature_K']
        if not all(col in df.columns for col in required_cols):
            df = pd.read_csv(log_file, header=None,
                           names=['Step', 'Time_ps', 'Temperature_K',
                                 'Potential_Energy_eV', 'Kinetic_Energy_eV',
                                 'Total_Energy_eV'])
        
        df['Energy_per_Atom_eV'] = df['Potential_Energy_eV'] / num_atoms
        
        folder_name = os.path.basename(folder)
        
        fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        
        axes[0].plot(df['Time_ps'], df['Energy_per_Atom_eV'], 
                    linewidth=1.0, color='#1f77b4')
        axes[0].set_ylabel('Potential Energy (eV/atom)', fontsize=12)
        axes[0].grid(True, alpha=0.3)
        axes[0].set_title(f'MD Simulation: {folder_name}', fontsize=14)
        
        axes[1].plot(df['Time_ps'], df['Temperature_K'], 
                    linewidth=1.0, color='#d62728')
        axes[1].set_xlabel('Time (ps)', fontsize=12)
        axes[1].set_ylabel('Temperature (K)', fontsize=12)
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        os.makedirs(output_dir, exist_ok=True)
        png_path = os.path.join(output_dir, f"{folder_name}.png")
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        stats = {
            'num_atoms': num_atoms,
            'total_time_ps': df['Time_ps'].max(),
            'avg_temp': df['Temperature_K'].mean(),
            'std_temp': df['Temperature_K'].std(),
            'final_energy': df['Energy_per_Atom_eV'].iloc[-1],
            'avg_energy': df['Energy_per_Atom_eV'].mean()
        }
        
        result['success'] = True
        result['message'] = f"Plot saved: {png_path}"
        result['data'] = stats
        
    except Exception as e:
        result['message'] = f"Error: {str(e)}"
    
    return result


def process_all_folders(base_dir='.', output_dir='./md_plots', pattern='md_output_*'):
    """
    Process all MD output folders in base directory
    
    Args:
        base_dir: Base directory to search
        output_dir: Directory to save plots
        pattern: Glob pattern for MD output folders
    """
    print(f"Searching for MD output folders in: {base_dir}")
    print(f"Pattern: {pattern}")
    
    folders = glob.glob(os.path.join(base_dir, pattern))
    
    if not folders:
        all_dirs = [d for d in os.listdir(base_dir) 
                   if os.path.isdir(os.path.join(base_dir, d)) 
                   and not d.startswith('.')]
        folders = [os.path.join(base_dir, d) for d in all_dirs
                  if glob.glob(os.path.join(base_dir, d, '*_md.log'))]
    
    if not folders:
        print("No MD output folders found!")
        return
    
    print(f"Found {len(folders)} MD output folder(s)")
    
    results = []
    for folder in folders:
        print(f"\nProcessing: {os.path.basename(folder)}")
        result = process_single_folder(folder, output_dir)
        results.append(result)
        
        if result['success']:
            print(f"  [OK] {result['message']}")
            stats = result['data']
            print(f"    Atoms: {stats['num_atoms']}, Time: {stats['total_time_ps']:.2f} ps")
            print(f"    Avg T: {stats['avg_temp']:.1f} +/- {stats['std_temp']:.1f} K")
            print(f"    Final E: {stats['final_energy']:.6f} eV/atom")
        else:
            print(f"  [FAIL] {result['message']}")
    
    success_count = sum(1 for r in results if r['success'])
    print(f"\n{'='*50}")
    print(f"Processed: {success_count}/{len(folders)} folders successfully")
    print(f"Plots saved to: {output_dir}")
    
    if success_count > 0:
        summary_path = os.path.join(output_dir, 'md_summary.csv')
        with open(summary_path, 'w') as f:
            f.write("Folder,Atoms,Time_ps,Avg_Temp_K,Std_Temp_K,Final_Energy_eV_per_atom,Avg_Energy_eV_per_atom\n")
            for r in results:
                if r['success']:
                    s = r['data']
                    f.write(f"{r['folder']},{s['num_atoms']},{s['total_time_ps']:.2f},"
                           f"{s['avg_temp']:.2f},{s['std_temp']:.2f},"
                           f"{s['final_energy']:.6f},{s['avg_energy']:.6f}\n")
        print(f"Summary saved: {summary_path}")


def main():
    """Command line interface for MD post-processing"""
    parser = argparse.ArgumentParser(
        description='GEWUM MD Post-processing - Generate energy-time plots',
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument('--base-dir', '-d', type=str, default='.',
                        help='Base directory containing MD output folders (default: .)')
    parser.add_argument('--output-dir', '-o', type=str, default='./md_plots',
                        help='Output directory for plots (default: ./md_plots)')
    parser.add_argument('--pattern', '-p', type=str, default='md_output_*',
                        help='Glob pattern for MD output folders (default: md_output_*)')
    parser.add_argument('--single', '-s', type=str, default=None,
                        help='Process a single folder instead of all')
    
    args = parser.parse_args()
    
    if args.single:
        result = process_single_folder(args.single, args.output_dir)
        if result['success']:
            print(f"Success: {result['message']}")
        else:
            print(f"Failed: {result['message']}")
        sys.exit(0 if result['success'] else 1)
    else:
        process_all_folders(args.base_dir, args.output_dir, args.pattern)


if __name__ == "__main__":
    main()
