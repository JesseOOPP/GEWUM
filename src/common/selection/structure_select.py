"""GEWUM Unified Structure Selection Script
Selects diverse crystal structures using clustering algorithms for all dimensions (0D-3D)

Usage:
    python structure_select.py --dim <dimension> --method <method> --target <count> [options]
    
Dimensions:
    3d - 3D periodic crystals (space groups 2-230)
    2d - 2D layered materials (layer groups 2-80)
    1d - 1D chains/nanowires (rod groups 2-75)
    0d - 0D molecules/clusters
    
Methods:
    random  - Random sampling
    kmeans  - K-means clustering with random selection from each cluster
    medoid  - K-means clustering with medoid selection (recommended)
    maxmin  - Maximum-minimum distance selection for maximum diversity

Descriptors:
    simple  - Handcrafted features (fast, available for 0D/2D/3D)
    soap    - SOAP descriptor (accurate, available for all dimensions)
    coulomb - Coulomb Matrix (fast, only for 0D)
"""
import os
import sys
import argparse
import shutil
import multiprocessing as mp
from functools import partial

def get_selector_class(dimension):
    """Get the appropriate selector class for the dimension."""
    try:
        from gewum.src.common.selection.selector_3d import StructureSelector3D
        from gewum.src.common.selection.selector_2d import StructureSelector2D
        from gewum.src.common.selection.selector_1d import StructureSelector1D
        from gewum.src.common.selection.selector_0d import StructureSelector0D
        selector_map = {
            '3d': StructureSelector3D,
            '2d': StructureSelector2D,
            '1d': StructureSelector1D,
            '0d': StructureSelector0D
        }
        if dimension in selector_map:
            return selector_map[dimension]
        raise ValueError(f"Unknown dimension: {dimension}")
    except ImportError:
        import sys
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        if dimension == '3d':
            from selector_3d import StructureSelector3D
            return StructureSelector3D
        elif dimension == '2d':
            from selector_2d import StructureSelector2D
            return StructureSelector2D
        elif dimension == '1d':
            from selector_1d import StructureSelector1D
            return StructureSelector1D
        elif dimension == '0d':
            from selector_0d import StructureSelector0D
            return StructureSelector0D
        else:
            raise ValueError(f"Unknown dimension: {dimension}")


def get_group_range(dimension):
    """Get the symmetry group range for the dimension."""
    ranges = {
        '3d': (2, 231),   # Space groups 2-230
        '2d': (2, 81),    # Layer groups 2-80
        '1d': (2, 76),    # Rod groups 2-75
        '0d': None        # No symmetry groups for 0D
    }
    return ranges.get(dimension)


def get_descriptor_choices(dimension):
    """Get available descriptor choices for the dimension."""
    choices = {
        '3d': ['simple', 'soap'],
        '2d': ['simple', 'soap'],
        '1d': ['soap'],  
        '0d': ['simple', 'soap', 'coulomb']
    }
    return choices.get(dimension, ['simple'])


def _build_log_section(dimension, method, feature_settings, selector, target_count,
                       hdbscan_params, cluster_count, noise_count,
                       all_cif_files, selected_files, unselected_files):
    """Build one selection_log.txt section (as a single string).

    Used for both DB and FS modes so the text content stays identical.
    For DB mode, several sections are concatenated into a single per-formula log.
    """
    lines = []
    if selector._use_db and selector._sg_number is not None:
        lines.append(f"===== SG {selector._sg_number} =====\n")
    lines.append(f"Dimension: {dimension.upper()}\n")
    lines.append(f"Selection method: {method}\n")
    lines.append(f"Descriptor: {feature_settings.get('descriptor_type', 'simple')}\n")
    if selector._use_db:
        lines.append(f"Storage: SQLite (SG={selector._sg_number})\n")
    if method == 'hdbscan':
        hdbscan_keep = hdbscan_params.get('keep_noise', True) if hdbscan_params else True
        lines.append(f"HDBSCAN clusters: {cluster_count}\n")
        lines.append(f"HDBSCAN noise points: {noise_count}\n")
        lines.append(f"HDBSCAN keep_noise: {hdbscan_keep}\n")
    else:
        lines.append(f"Target count: {target_count}\n")
    lines.append(f"Total CIF files: {len(all_cif_files)}\n")
    lines.append(f"Selected: {len(selected_files)}\n")
    lines.append(f"Removed: {len(unselected_files)}\n\n")
    lines.append("Selected files:\n")
    for file in selected_files:
        lines.append(f"  {file}\n")
    lines.append("\n")
    return "".join(lines)


def process_single_directory(input_dir, dimension, method, target_count, feature_settings, folder_name="", hdbscan_params=None):
    """Process a single directory containing CIF files.

    Returns:
        (success, n_selected, n_removed, log_info)
        log_info is either None (FS mode wrote its own log already) or a dict
        {'target_path': '<formula_dir>/selection_log.txt', 'text': '<section>'} 
        for DB mode, so the caller can aggregate sections per formula and write once.
    """
    print(f"\n{'='*50}")
    print(f"Processing: {folder_name or input_dir}")
    print(f"{'='*50}")

    SelectorClass = get_selector_class(dimension)
    selector = SelectorClass(input_dir, feature_settings)

    try:
        file_list, features = selector.load_and_featurize()

        if method != 'hdbscan' and len(file_list) <= target_count:
            print(f"Skipping: only {len(file_list)} files (target: {target_count})")
            return True, len(file_list), 0, None

        cluster_count = None
        noise_count = None
        if method == 'hdbscan':
            hdbscan_params = hdbscan_params or {}
            result = selector.select(method=method, **hdbscan_params)
            selected_files, cluster_count, noise_count = result
            print(f"HDBSCAN auto-selected {len(selected_files)} structures ({cluster_count} clusters)")
        else:
            selected_files = selector.select(method=method, target_count=target_count)
            print(f"Selected {len(selected_files)} structures using {method} method")

        all_cif_files = [f for f in os.listdir(input_dir) if f.endswith('.cif')] if not selector._use_db else file_list
        unselected_files = [f for f in all_cif_files if f not in selected_files]

        log_info = None

        if selector._use_db and selector._sg_number is not None:
            parent_dir = os.path.dirname(os.path.abspath(input_dir))
            from gewum.src.common.cif_db import CifDatabase
            with CifDatabase(parent_dir) as db:
                db.update_stage_batch(
                    [(selector._sg_number, file) for file in unselected_files],
                    'removed'
                )
            print(f"Marked {len(unselected_files)} structures as 'removed' in DB")

            section = _build_log_section(
                dimension, method, feature_settings, selector, target_count,
                hdbscan_params, cluster_count, noise_count,
                all_cif_files, selected_files, unselected_files,
            )
            log_info = {
                'target_path': os.path.join(parent_dir, 'selection_log.txt'),
                'text': section,
            }
        else:
            remove_dir = os.path.join(input_dir, 'remove')
            os.makedirs(remove_dir, exist_ok=True)
            for file in unselected_files:
                src_path = os.path.join(input_dir, file)
                dst_path = os.path.join(remove_dir, file)
                shutil.move(src_path, dst_path)
            print(f"Moved {len(unselected_files)} files to {remove_dir}")

            log_file = os.path.join(input_dir, 'selection_log.txt')
            with open(log_file, 'w') as f:
                f.write(_build_log_section(
                    dimension, method, feature_settings, selector, target_count,
                    hdbscan_params, cluster_count, noise_count,
                    all_cif_files, selected_files, unselected_files,
                ))

        return True, len(selected_files), len(unselected_files), log_info

    except Exception as e:
        print(f"Error processing {input_dir}: {e}")
        import traceback
        traceback.print_exc()
        return False, 0, 0, None


def process_directory_wrapper(args, dimension, method, target_count, feature_settings, hdbscan_params=None):
    """Wrapper for multiprocessing."""
    folder_path, folder_name = args
    return process_single_directory(folder_path, dimension, method, target_count, feature_settings, folder_name, hdbscan_params)


def _flush_db_logs(results):
    """Aggregate DB-mode log sections per formula directory and write once."""
    grouped = {}
    for r in results:
        if len(r) < 4 or r[3] is None:
            continue
        info = r[3]
        grouped.setdefault(info['target_path'], []).append(info['text'])
    for target_path, sections in grouped.items():
        try:
            with open(target_path, 'w') as f:
                for sec in sections:
                    f.write(sec)
        except OSError as e:
            print(f"Warning: failed to write {target_path}: {e}")


def find_directories_to_process(base_dir, dimension):
    """Find all directories with CIF files based on dimension.

    Supports two input levels:
      1. base_dir = top working directory containing formula subdirs
         e.g. base_dir/Na2Cl2/1/*.cif
      2. base_dir = a single formula directory directly
         e.g. base_dir/1/*.cif  (auto-detected when group subdirs exist)

    Also supports SQLite-based storage via structures.db.
    """
    directories = []
    group_range = get_group_range(dimension)

    db_path = os.path.join(base_dir, 'structures.db')
    use_db = os.path.isfile(db_path)

    if use_db:
        import sqlite3
        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row
        try:
            rows = conn.execute(
                "SELECT DISTINCT sg_number FROM structures WHERE stage = 'initial' ORDER BY sg_number"
            ).fetchall()
        finally:
            conn.close()

        folder_name = os.path.basename(os.path.normpath(base_dir))
        for row in rows:
            sg = row['sg_number']
            subdir_path = os.path.join(base_dir, str(sg))
            directories.append((subdir_path, f"{folder_name}/{sg}"))
        if directories:
            return directories

    if group_range:
        has_group_subdirs = any(
            os.path.isdir(os.path.join(base_dir, str(num)))
            for num in range(group_range[0], group_range[1])
        )
        if has_group_subdirs:
            folder_name = os.path.basename(os.path.normpath(base_dir))
            for num in range(group_range[0], group_range[1]):
                subdir_path = os.path.join(base_dir, str(num))
                if os.path.exists(subdir_path) and os.path.isdir(subdir_path):
                    cif_files = [f for f in os.listdir(subdir_path) if f.endswith('.cif')]
                    if cif_files:
                        directories.append((subdir_path, f"{folder_name}/{num}"))
            if directories:
                return directories

    for folder in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder)

        if not os.path.isdir(folder_path):
            continue

        formula_db = os.path.join(folder_path, 'structures.db')
        if os.path.isfile(formula_db) and group_range:
            import sqlite3
            conn = sqlite3.connect(formula_db)
            conn.row_factory = sqlite3.Row
            try:
                rows = conn.execute(
                    "SELECT DISTINCT sg_number FROM structures WHERE stage = 'initial' ORDER BY sg_number"
                ).fetchall()
            finally:
                conn.close()
            prev_count = len(directories)
            for row in rows:
                sg = row['sg_number']
                subdir_path = os.path.join(folder_path, str(sg))
                directories.append((subdir_path, f"{folder}/{sg}"))
            if len(directories) > prev_count:
                continue  

        if group_range:
            for num in range(group_range[0], group_range[1]):
                subdir_path = os.path.join(folder_path, str(num))
                if os.path.exists(subdir_path) and os.path.isdir(subdir_path):
                    cif_files = [f for f in os.listdir(subdir_path) if f.endswith('.cif')]
                    if cif_files:
                        directories.append((subdir_path, f"{folder}/{num}"))
        else:
            cif_files = [f for f in os.listdir(folder_path) if f.endswith('.cif')]
            if cif_files:
                directories.append((folder_path, folder))
            else:
                for subdir in os.listdir(folder_path):
                    subdir_path = os.path.join(folder_path, subdir)
                    if os.path.isdir(subdir_path):
                        cif_files = [f for f in os.listdir(subdir_path) if f.endswith('.cif')]
                        if cif_files:
                            directories.append((subdir_path, f"{folder}/{subdir}"))

    return directories


def main():
    parser = argparse.ArgumentParser(
        description="GEWUM Unified Structure Selection - Select diverse structures for all dimensions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 3D structures (space groups)
  python structure_select.py --dim 3d --method medoid --target 30
  
  # 2D materials (layer groups) with SOAP
  python structure_select.py --dim 2d --method medoid --target 20 --descriptor soap
  
  # 1D chains (SOAP required)
  python structure_select.py --dim 1d --method medoid --target 15
  
  # 0D molecules with Coulomb Matrix
  python structure_select.py --dim 0d --method medoid --target 20 --descriptor coulomb
  
  # HDBSCAN auto-selection (automatic cluster count)
  python structure_select.py --dim 3d --method hdbscan --min-cluster-size 5
  
  # Single directory mode
  python structure_select.py --dim 3d --method maxmin --target 50 --single-dir ./cifs
        """
    )
    
    parser.add_argument('--dim', '-D',
                        required=True,
                        choices=['3d', '2d', '1d', '0d'],
                        help='Structure dimension (3d/2d/1d/0d)')
    parser.add_argument('--method', '-m', 
                        choices=['random', 'kmeans', 'medoid', 'maxmin', 'hdbscan'],
                        default='medoid',
                        help='Selection method (default: medoid)')
    parser.add_argument('--target', '-t', 
                        type=int, default=30,
                        help='Target number of structures per directory (default: 30, ignored for hdbscan)')
    parser.add_argument('--input-dir', '-i',
                        default='.',
                        help='Base input directory (default: current directory)')
    parser.add_argument('--single-dir', '-s',
                        help='Process a single directory instead of searching for subdirectories')
    parser.add_argument('--workers', '-w',
                        type=int, default=1,
                        help='Number of parallel workers (default: 1)')
    parser.add_argument('--descriptor', '-d',
                        default=None,
                        help='Feature descriptor type (default depends on dimension)')
    parser.add_argument('--soap-r-cut',
                        type=float, default=None,
                        help='SOAP cutoff radius')
    parser.add_argument('--soap-n-max',
                        type=int, default=4,
                        help='SOAP radial basis (default: 4)')
    parser.add_argument('--soap-l-max',
                        type=int, default=4,
                        help='SOAP angular basis (default: 4)')
    parser.add_argument('--coulomb-n-atoms-max',
                        type=int, default=100,
                        help='Max atoms for Coulomb Matrix (default: 100, 0D only)')
    parser.add_argument('--periodic-dir', '-p',
                        choices=['x', 'y', 'z'],
                        default='z',
                        help='Periodic direction for 1D structure (default: z)')
    parser.add_argument('--min-cluster-size',
                        type=int, default=5,
                        help='HDBSCAN minimum cluster size (default: 5, moderate filtering)')
    parser.add_argument('--min-samples',
                        type=int, default=3,
                        help='HDBSCAN min samples for core points (default: 3)')
    parser.add_argument('--cluster-selection-epsilon',
                        type=float, default=0.0,
                        help='HDBSCAN distance threshold for merging clusters (default: 0.0, auto)')
    parser.add_argument('--alpha',
                        type=float, default=1.0,
                        help='HDBSCAN density decay parameter (default: 1.0, moderate)')
    parser.add_argument('--no-keep-noise',
                        dest='keep_noise',
                        action='store_false',
                        default=True,
                        help='Discard HDBSCAN noise points (default: keep noise points)')
    parser.add_argument('--no-pca',
                        action='store_true',
                        help='Disable auto PCA for high-dim features (SOAP/Coulomb)')
    parser.add_argument('--pca-variance',
                        type=float, default=0.95,
                        help='PCA variance ratio to preserve (default: 0.95 = 95%%)')
    parser.add_argument('--pca-components',
                        type=int, default=None,
                        help='Fixed PCA components (overrides pca-variance)')
    
    args = parser.parse_args()
    
    valid_descriptors = get_descriptor_choices(args.dim)
    if args.descriptor is None:
        args.descriptor = valid_descriptors[0]
    elif args.descriptor not in valid_descriptors:
        print(f"Error: Descriptor '{args.descriptor}' not available for {args.dim.upper()}")
        print(f"Available descriptors: {', '.join(valid_descriptors)}")
        sys.exit(1)
    
    if args.soap_r_cut is None:
        default_r_cuts = {'3d': 6.0, '2d': 5.0, '1d': 4.0, '0d': 5.0}
        args.soap_r_cut = default_r_cuts[args.dim]
    
    feature_settings = {
        'descriptor_type': args.descriptor,
        'soap_r_cut': args.soap_r_cut,
        'soap_n_max': args.soap_n_max,
        'soap_l_max': args.soap_l_max,
        'coulomb_n_atoms_max': args.coulomb_n_atoms_max,
        'periodic_direction': args.periodic_dir,
        'include_lattice_params': True,
        'include_position_stats': True,
        'include_basic_features': True,
        'use_pca': None if not args.no_pca else False,
        'pca_variance': args.pca_variance,
        'pca_n_components': args.pca_components
    }
    
    hdbscan_params = {
        'min_cluster_size': args.min_cluster_size,
        'min_samples': args.min_samples,
        'cluster_selection_epsilon': args.cluster_selection_epsilon,
        'alpha': args.alpha,
        'keep_noise': args.keep_noise
    }
    
    dim_names = {'3d': '3D Crystal', '2d': '2D Material', '1d': '1D Chain', '0d': '0D Cluster'}
    
    print(f"GEWUM {dim_names[args.dim]} Structure Selection")
    print(f"{'='*50}")
    print(f"Dimension: {args.dim.upper()}")
    print(f"Method: {args.method}")
    print(f"Descriptor: {args.descriptor}")
    if args.method == 'hdbscan':
        print(f"HDBSCAN min_cluster_size: {args.min_cluster_size}")
        if args.min_samples:
            print(f"HDBSCAN min_samples: {args.min_samples}")
        if args.cluster_selection_epsilon > 0:
            print(f"HDBSCAN epsilon: {args.cluster_selection_epsilon}")
        if args.alpha != 1.0:
            print(f"HDBSCAN alpha: {args.alpha}")
        print(f"HDBSCAN keep_noise: {args.keep_noise}")
    else:
        print(f"Target count: {args.target}")
    print(f"Input directory: {os.path.abspath(args.input_dir)}")
    if args.descriptor == 'soap':
        print(f"SOAP params: r_cut={args.soap_r_cut}, n_max={args.soap_n_max}, l_max={args.soap_l_max}")
    elif args.descriptor == 'coulomb':
        print(f"Coulomb Matrix: n_atoms_max={args.coulomb_n_atoms_max}")
    if args.dim == '1d':
        print(f"Periodic direction: {args.periodic_dir}")
    print(f"{'='*50}")
    
    if args.single_dir:
        success, selected, removed, log_info = process_single_directory(
            args.single_dir, args.dim, args.method, args.target, feature_settings,
            hdbscan_params=hdbscan_params if args.method == 'hdbscan' else None
        )
        if log_info is not None:
            _flush_db_logs([(success, selected, removed, log_info)])
        print(f"\nResult: {'Success' if success else 'Failed'}")
        print(f"Selected: {selected}, Removed: {removed}")
        return
    
    directories = find_directories_to_process(args.input_dir, args.dim)
    
    if not directories:
        print("No directories with CIF files found!")
        group_range = get_group_range(args.dim)
        if group_range:
            print(f"Expected structure: input_dir/formula/group_number/*.cif")
        else:
            print(f"Expected structure: input_dir/formula/*.cif")
        return
    
    print(f"Found {len(directories)} directories to process")
    
    if args.workers > 1:
        print(f"Using {args.workers} parallel workers")
        with mp.Pool(processes=args.workers) as pool:
            worker_func = partial(
                process_directory_wrapper,
                dimension=args.dim,
                method=args.method,
                target_count=args.target,
                feature_settings=feature_settings,
                hdbscan_params=hdbscan_params if args.method == 'hdbscan' else None
            )
            results = pool.map(worker_func, directories)
    else:
        results = []
        for dir_path, dir_name in directories:
            result = process_single_directory(
                dir_path, args.dim, args.method, args.target, feature_settings, dir_name,
                hdbscan_params=hdbscan_params if args.method == 'hdbscan' else None
            )
            results.append(result)
    
    _flush_db_logs(results)
    
    success_count = sum(1 for r in results if r[0])
    total_selected = sum(r[1] for r in results)
    total_removed = sum(r[2] for r in results)
    
    print(f"\n{'='*50}")
    print(f"{args.dim.upper()} Selection Summary")
    print(f"{'='*50}")
    print(f"Directories processed: {len(directories)}")
    print(f"Successful: {success_count}")
    print(f"Total structures selected: {total_selected}")
    print(f"Total structures removed: {total_removed}")


if __name__ == "__main__":
    main()
