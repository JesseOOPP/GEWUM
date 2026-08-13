"""
GEWUM Phonon Post-processing Module
1. Filter stable structures (no imaginary frequencies) from ph.log
2. Copy phonon spectrum images (PNG) to the filtered directory
"""
import os
import re
import glob
import shutil
import argparse


def analyze_and_filter_stable(log_file, target_dir, unstable_dir=None):
    """
    Analyze phonon log file and copy stable/unstable CIF files to target directories.
    
    Args:
        log_file: Path to phonon calculation log (e.g., cif_ph.out)
        target_dir: Directory to copy stable structures to
        unstable_dir: Directory to copy unstable structures to (optional)
    
    Returns:
        Tuple of (stable_files, unstable_files) - lists of basenames
    """
    if not os.path.exists(log_file):
        print(f"Error: Log file {log_file} not found")
        return [], []
    
    os.makedirs(target_dir, exist_ok=True)
    
    with open(log_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    pattern = r'Processing: (.*?\.cif).*?Has imaginary phonon: (True|False)'
    matches = re.findall(pattern, content, re.DOTALL)
    
    stable_files = []
    unstable_files = []
    
    for cif_file, has_imaginary in matches:
        if has_imaginary == "False":
            stable_files.append(cif_file)
        else:
            unstable_files.append(cif_file)
    
    print(f"Stable structures (no imaginary): {len(stable_files)}")
    print(f"Unstable structures (imaginary): {len(unstable_files)}")
    
    copied_stable = []
    for cif_file in stable_files:
        if os.path.exists(cif_file):
            basename = os.path.basename(cif_file)
            dest_path = os.path.join(target_dir, basename)
            
            if os.path.exists(dest_path):
                base, ext = os.path.splitext(basename)
                counter = 1
                while os.path.exists(dest_path):
                    dest_path = os.path.join(target_dir, f"{base}_{counter}{ext}")
                    counter += 1
            
            shutil.copy2(cif_file, dest_path)
            copied_stable.append(basename)
            print(f"  Copied: {cif_file}")
        else:
            print(f"  Warning: {cif_file} not found")
    
    print(f"\nCopied {len(copied_stable)} stable CIF files to {target_dir}/")
    
    copied_unstable = []
    if unstable_dir:
        os.makedirs(unstable_dir, exist_ok=True)
        print(f"\n[Unstable Structures]")
        for cif_file in unstable_files:
            if os.path.exists(cif_file):
                basename = os.path.basename(cif_file)
                dest_path = os.path.join(unstable_dir, basename)
                
                if os.path.exists(dest_path):
                    base, ext = os.path.splitext(basename)
                    counter = 1
                    while os.path.exists(dest_path):
                        dest_path = os.path.join(unstable_dir, f"{base}_{counter}{ext}")
                        counter += 1
                
                shutil.copy2(cif_file, dest_path)
                copied_unstable.append(basename)
                print(f"  Copied: {cif_file}")
            else:
                print(f"  Warning: {cif_file} not found")
        
        print(f"\nCopied {len(copied_unstable)} unstable CIF files to {unstable_dir}/")
    
    return copied_stable, copied_unstable


def copy_phonon_images(target_dir, tmp_base="tmp"):
    """
    Copy phonon spectrum PNG images to the target directory.
    
    Args:
        target_dir: Directory containing filtered CIF files (e.g., 0_final)
        tmp_base: Base directory for tmp folders (relative to current working dir)
    """
    tmp_path = os.path.abspath(tmp_base)
    
    original_dir = os.getcwd()
    os.chdir(target_dir)
    
    copied_count = 0
    
    for cif_file in glob.glob('*.cif'):
        prefix = cif_file[:-4]  
        tmp_dir = os.path.join(tmp_path, prefix)
        
        if os.path.isdir(tmp_dir):
            png_files = glob.glob(os.path.join(tmp_dir, '*.png'))
            
            if png_files:
                source_png = png_files[0]
                target_png = f"{prefix}.png"
                shutil.copy2(source_png, target_png)
                print(f"  PNG: {prefix}.png")
                copied_count += 1
            else:
                print(f"  Warning: No PNG in {tmp_dir}")
        else:
            print(f"  Warning: tmp/{prefix} not found")
    
    os.chdir(original_dir)
    print(f"\nCopied {copied_count} PNG files to {target_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description='Phonon post-processing: filter stable structures and collect PNG images'
    )
    parser.add_argument(
        '--log', '-l',
        default='ph.log',
        help='Phonon calculation log file (default: ph.log)'
    )
    parser.add_argument(
        '--output', '-o',
        default='0_final',
        help='Output directory for stable structures (default: 0_final)'
    )
    parser.add_argument(
        '--unstable-output', '-u',
        default='0_unstable',
        help='Output directory for unstable structures (default: 0_unstable)'
    )
    parser.add_argument(
        '--tmp-dir', '-t',
        default='tmp',
        help='Temporary directory containing phonon results (default: tmp)'
    )
    parser.add_argument(
        '--skip-filter', 
        action='store_true',
        help='Skip filtering step, only copy PNG images'
    )
    parser.add_argument(
        '--skip-png',
        action='store_true',
        help='Skip PNG copying step'
    )
    parser.add_argument(
        '--no-unstable',
        action='store_true',
        help='Do not extract unstable structures'
    )
    args = parser.parse_args()
    
    print("=" * 60)
    print("GEWUM Phonon Post-processing")
    print("=" * 60)
    
    if not args.skip_filter:
        print("\n[Step 1] Filtering structures...")
        unstable_target = None if args.no_unstable else args.unstable_output
        analyze_and_filter_stable(args.log, args.output, unstable_target)
    else:
        print("\n[Step 1] Skipped (--skip-filter)")
    
    if not args.skip_png:
        print(f"\n[Step 2] Copying phonon spectrum images...")
        print(f"\n[Stable PNG]")
        if os.path.isdir(args.output):
            copy_phonon_images(args.output, tmp_base=args.tmp_dir)
        else:
            print(f"  Error: Output directory {args.output} not found")
        
        if not args.no_unstable and os.path.isdir(args.unstable_output):
            print(f"\n[Unstable PNG]")
            copy_phonon_images(args.unstable_output, tmp_base=args.tmp_dir)
    else:
        print("\n[Step 2] Skipped (--skip-png)")
    
    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
