"""
GEWUM PT Visualization Module
Generates space group Sankey diagrams, energy violin plots, and RDF comparison
plots for the PT (Perturbation) workflow.

PT directory structure:
    work_dir/
    |-- *.cif                   # Initial/perturbed CIF files (root level)
    |-- relaxed/                # Relaxed structures
    |   |-- *_relaxed.cif
    |   `-- bond_mis/           # Structures with short bonds (excluded)
    |-- energy_results.csv      # 4-col energy data
    `-- 0_final_result_tot.txt  # 6-col final result

Usage:
    python -m gewum.src.PTworkflows.viz_pt --cif-dir . --plot all
"""
import os
import sys
import csv
import argparse
import datetime
import numpy as np
from pathlib import Path
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from tqdm import tqdm


def _find_initial_cifs(directory):
    """
    Find initial CIF files in PT work directory (root level only).
    Excludes relaxed/, bond_mis/, supercell_structures/, etc.

    Args:
        directory: PT work directory.

    Returns:
        Sorted list of CIF file paths.
    """
    excluded_dirs = {'relaxed', 'bond_mis', 'supercell_structures',
                     'relaxed_symmetry', 'error_str', '0_cif'}
    cif_files = []
    for f in os.listdir(directory):
        if f.endswith('.cif'):
            cif_files.append(os.path.join(directory, f))
    for sub in os.listdir(directory):
        sub_path = os.path.join(directory, sub)
        if os.path.isdir(sub_path) and sub not in excluded_dirs:
            for f in os.listdir(sub_path):
                if f.endswith('.cif'):
                    cif_files.append(os.path.join(sub_path, f))
    return sorted(cif_files)


def _find_relaxed_cifs(directory):
    """
    Find relaxed CIF files in PT work directory.
    Looks in relaxed/ subdirectory, excluding bond_mis/.

    Args:
        directory: PT work directory.

    Returns:
        Sorted list of relaxed CIF file paths.
    """
    relaxed_dir = os.path.join(directory, 'relaxed')
    if not os.path.isdir(relaxed_dir):
        print(f"[PT-viz] Warning: relaxed/ directory not found in {directory}")
        return []

    cif_files = []
    for f in os.listdir(relaxed_dir):
        fpath = os.path.join(relaxed_dir, f)
        if f.endswith('.cif') and os.path.isfile(fpath):
            cif_files.append(fpath)
    return sorted(cif_files)


def _collect_cifs(cif_dirs):
    """
    Collect initial and relaxed CIF files from one or more PT directories.

    Args:
        cif_dirs: list of PT work directory paths.

    Returns:
        dict: {'initial': [...], 'relaxed': [...]}
    """
    initial_all = []
    relaxed_all = []

    for cif_dir in cif_dirs:
        print(f"[PT-viz] Scanning directory: {cif_dir}")
        initial = _find_initial_cifs(cif_dir)
        relaxed = _find_relaxed_cifs(cif_dir)
        print(f"[PT-viz]   Initial CIFs: {len(initial)}")
        print(f"[PT-viz]   Relaxed CIFs: {len(relaxed)}")
        initial_all.extend(initial)
        relaxed_all.extend(relaxed)

    result = {}
    if initial_all:
        result['initial'] = initial_all
    if relaxed_all:
        result['relaxed'] = relaxed_all

    print(f"[PT-viz] --- Collection summary ---")
    for stage, paths in result.items():
        print(f"[PT-viz]   {stage}: {len(paths)} CIF files")
    return result

def _identify_sg_from_structure(cif_path, symprec=0.1):
    """
    Identify space group from a CIF file using spglib.

    Returns:
        (symbol, number): e.g. ('Fm-3m', 225)
    """
    from pymatgen.io.cif import CifParser
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    try:
        parser = CifParser(cif_path)
        structure = parser.parse_structures(primitive=False)[0]
    except Exception:
        from pymatgen.core import Structure
        structure = Structure.from_file(cif_path)

    sga = SpacegroupAnalyzer(structure, symprec=symprec)
    sg_symbol = sga.get_space_group_symbol()
    sg_number = sga.get_space_group_number()
    return sg_symbol, sg_number


def _get_initial_sg(cif_path, energy_data_map=None):
    """
    Get initial space group for a relaxed CIF.

    Strategy:
    1. If energy_data_map is available and has SG_ori, use it.
    2. Otherwise try to infer from the relaxed CIF basename -> initial CIF parent dir.
    3. Fall back to 'unknown'.

    Args:
        cif_path: path to relaxed CIF.
        energy_data_map: dict from _load_energy_data().

    Returns:
        str: space group number or label.
    """
    if energy_data_map:
        base = Path(cif_path).stem
        match_key = base.replace('_relaxed', '')
        for key, data in energy_data_map.items():
            if match_key in key or key in match_key:
                sg = data.get('sg_ori', '')
                if sg and sg != 'Unknown':
                    return sg

    work_dir = str(Path(cif_path).parent.parent)
    dir_name = os.path.basename(work_dir)
    if dir_name.isdigit():
        return dir_name
    return 'unknown'


def _analyse_sg_worker(args):
    """Worker: identify space group for a single CIF (picklable)."""
    cif_path, symprec = args
    try:
        rel_sym, rel_num = _identify_sg_from_structure(cif_path, symprec=symprec)
        return cif_path, rel_sym, rel_num
    except Exception as e:
        print(f"[PT-viz] Warning: SG analysis failed for {os.path.basename(cif_path)}: {e}")
        return cif_path, 'P1', 1

def _load_energy_data(cif_dirs):
    """
    Load energy data from PT work directories.

    Searches for:
    1. 0_final_result_tot.txt (6-col with Chemical_Formula, SG_ori)
    2. energy_results.csv (4-col fallback)

    Returns:
        dict: {cif_basename: {'formula': str, 'total_energy': float,
               'energy_per_atom': float, 'cif_path': str, 'sg_ori': str}}
    """
    energy_data = {}

    for cif_dir in cif_dirs:
        # Try 6-col format first
        result_file = os.path.join(cif_dir, '0_final_result_tot.txt')
        if os.path.isfile(result_file):
            print(f"[PT-viz] Loading energy data from: {result_file}")
            try:
                with open(result_file, 'r', encoding='utf-8') as f:
                    reader = csv.reader(f)
                    header = next(reader, None)
                    for row in reader:
                        if len(row) < 6:
                            continue
                        formula = row[0].strip()
                        base_name = row[1].strip()
                        total_energy = float(row[2].strip())
                        energy_per_atom = float(row[3].strip())
                        cif_path = row[4].strip()
                        sg_ori = row[5].strip()
                        energy_data[base_name] = {
                            'formula': formula,
                            'total_energy': total_energy,
                            'energy_per_atom': energy_per_atom,
                            'cif_path': cif_path,
                            'sg_ori': sg_ori,
                        }
            except Exception as e:
                print(f"[PT-viz] Warning: failed to read {result_file}: {e}")
            continue

        csv_file = os.path.join(cif_dir, 'energy_results.csv')
        if os.path.isfile(csv_file):
            print(f"[PT-viz] Loading energy data from: {csv_file}")
            try:
                with open(csv_file, 'r', encoding='utf-8') as f:
                    reader = csv.reader(f)
                    header = next(reader, None)
                    for row in reader:
                        if len(row) < 4:
                            continue
                        base_name = row[0].strip()
                        total_energy = float(row[1].strip())
                        energy_per_atom = float(row[2].strip())
                        cif_path = row[3].strip()
                        energy_data[base_name] = {
                            'formula': '',
                            'total_energy': total_energy,
                            'energy_per_atom': energy_per_atom,
                            'cif_path': cif_path,
                            'sg_ori': os.path.basename(cif_dir),
                        }
            except Exception as e:
                print(f"[PT-viz] Warning: failed to read {csv_file}: {e}")

    print(f"[PT-viz] Energy data loaded: {len(energy_data)} entries")
    return energy_data


def _match_energy_for_cif(cif_path, energy_data):
    """
    Match energy data entry for a given relaxed CIF path.

    Returns:
        dict with energy info, or None if not matched.
    """
    p = Path(cif_path)
    cif_stem = p.stem 
    match_key = cif_stem.replace('_relaxed', '')

    if match_key in energy_data:
        return energy_data[match_key]

    for key, data in energy_data.items():
        if key in cif_stem or cif_stem in key or match_key in key or key in match_key:
            return data

    cif_resolved = str(p.resolve())
    for key, data in energy_data.items():
        try:
            if str(Path(data['cif_path']).resolve()) == cif_resolved:
                return data
        except Exception:
            pass

    return None

def _compute_rdf(structure, r_max=10.0, bin_size=0.1):
    """
    Compute radial distribution function for a structure.

    Returns:
        (r_values, g_r): arrays of bin centers and RDF values.
    """
    n_bins = int(r_max / bin_size)
    r_values = np.linspace(bin_size / 2, r_max - bin_size / 2, n_bins)
    counts = np.zeros(n_bins)

    n_atoms = len(structure)
    volume = structure.volume
    rho = n_atoms / volume

    for i, site in enumerate(structure):
        neighbors = structure.get_neighbors(site, r=r_max)
        for neighbor in neighbors:
            dist = neighbor.nn_distance
            bin_idx = int(dist / bin_size)
            if 0 <= bin_idx < n_bins:
                counts[bin_idx] += 1

    g_r = np.zeros(n_bins)
    for i in range(n_bins):
        r_inner = i * bin_size
        r_outer = (i + 1) * bin_size
        shell_vol = (4 / 3) * np.pi * (r_outer ** 3 - r_inner ** 3)
        if shell_vol > 0 and rho > 0 and n_atoms > 0:
            g_r[i] = counts[i] / (n_atoms * rho * shell_vol)

    return r_values, g_r


def _compute_rdf_worker(args):
    """Worker: compute RDF for a single CIF file."""
    path, r_max, bin_size = args
    try:
        from pymatgen.core import Structure
        struct = Structure.from_file(path)
        r, gr = _compute_rdf(struct, r_max=r_max, bin_size=bin_size)
        return r, gr, None
    except Exception as e:
        return None, None, f"{os.path.basename(path)}: {e}"

def _cubic_bezier_y(x, y0, y1):
    """Compute smooth Bezier interpolation."""
    t = x
    return y0 * (1 - t) ** 3 + y0 * 3 * t * (1 - t) ** 2 + y1 * 3 * t ** 2 * (1 - t) + y1 * t ** 3


def _merge_others(column_totals, column_nodes, min_fraction=0.01, min_keep=20):
    """Merge tiny nodes into 'Others'."""
    if len(column_nodes) <= 30:
        return column_nodes, set()
    total = sum(column_totals[n] for n in column_nodes)
    sorted_nodes = sorted(column_nodes, key=lambda n: column_totals[n], reverse=True)
    tiny_set = set()
    keep = []
    for i, n in enumerate(sorted_nodes):
        if i < min_keep:
            keep.append(n)
        elif column_totals[n] / max(total, 1) < min_fraction:
            tiny_set.add(n)
        else:
            keep.append(n)
    return keep, tiny_set


def plot_sankey(cif_paths_relaxed, output_path, energy_data=None,
               symprec=0.1, dpi=300, ignore_p1=False, workers=1,
               cmap_left=None, cmap_right=None):
    """
    Draw a space-group Sankey diagram for PT relaxed structures.

    Shows flow from initial SG (from energy data SG_ori or directory name)
    to relaxed SG (identified via spglib).

    Args:
        cif_paths_relaxed: list of relaxed CIF file paths.
        output_path: Output PNG path.
        energy_data: dict from _load_energy_data().
        symprec: Symmetry precision.
        dpi: Figure DPI.
        ignore_p1: Exclude P1 from diagram.
        workers: Number of parallel workers.
        cmap_left: Colormap name for left column (initial SG). None = default pink-blue.
        cmap_right: Colormap name for right column (relaxed SG). None = default blue-pink.
    """
    import matplotlib.cm as cm

    if not cif_paths_relaxed:
        print("[PT-viz] Warning: no relaxed CIF files for Sankey diagram")
        return

    flows = []
    if workers > 1 and len(cif_paths_relaxed) > 1:
        tasks = [(p, symprec) for p in cif_paths_relaxed]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(tqdm(
                pool.map(_analyse_sg_worker, tasks,
                         chunksize=max(1, len(tasks) // workers // 4)),
                total=len(tasks), desc="[PT-viz] Analysing space groups"))
        for cif_path, rel_sym, rel_num in results:
            init_sg = _get_initial_sg(cif_path, energy_data)
            flows.append((f"SG {init_sg}", f"SG {rel_num}"))
    else:
        for cif_path in tqdm(cif_paths_relaxed, desc="[PT-viz] Analysing space groups"):
            rel_sym, rel_num = _identify_sg_from_structure(cif_path, symprec=symprec)
            init_sg = _get_initial_sg(cif_path, energy_data)
            flows.append((f"SG {init_sg}", f"SG {rel_num}"))

    if ignore_p1:
        flows = [(l, r) for l, r in flows if r != "SG 1"]
        if not flows:
            print("[PT-viz] Warning: all flows lead to P1, nothing to display")
            return

    if not flows:
        print("[PT-viz] Warning: no valid flows for Sankey diagram")
        return

    flow_counts = Counter(flows)

    left_totals = defaultdict(int)
    right_totals = defaultdict(int)
    for (l, r), cnt in flow_counts.items():
        left_totals[l] += cnt
        right_totals[r] += cnt

    others_label = "Others"
    left_nodes_sorted = sorted(left_totals.keys(), key=lambda k: -left_totals[k])
    right_nodes_sorted = sorted(right_totals.keys(), key=lambda k: -right_totals[k])
    left_keep, left_tiny = _merge_others(left_totals, left_nodes_sorted)
    right_keep, right_tiny = _merge_others(right_totals, right_nodes_sorted)

    new_flows = defaultdict(int)
    for (lf, rf), cnt in flow_counts.items():
        lf_mapped = others_label if lf in left_tiny else lf
        rf_mapped = others_label if rf in right_tiny else rf
        new_flows[(lf_mapped, rf_mapped)] += cnt
    top_flows = list(new_flows.items())

    left_totals = defaultdict(int)
    right_totals = defaultdict(int)
    for (lf, rf), cnt in top_flows:
        left_totals[lf] += cnt
        right_totals[rf] += cnt

    left_nodes = sorted([n for n in left_totals if n != others_label],
                        key=lambda k: -left_totals[k])
    if others_label in left_totals:
        left_nodes.append(others_label)

    right_nodes = sorted([n for n in right_totals if n != others_label],
                         key=lambda k: -right_totals[k])
    if others_label in right_totals:
        right_nodes.append(others_label)

    total_left = sum(left_totals[n] for n in left_nodes)
    total_right = sum(right_totals[n] for n in right_nodes)
    total_max = max(total_left, total_right)

    gap = 0.008
    bar_w = 0.06
    x_left = 0.0
    x_right = 0.6

    def _compute_positions(nodes, totals, total_val, gap_val):
        positions = {}
        cursor = 0.0
        for node in nodes:
            frac = totals[node] / total_val * (1 - gap_val * (len(nodes) - 1))
            positions[node] = (cursor, cursor + frac)
            cursor += frac + gap_val
        return positions

    left_pos = _compute_positions(left_nodes, left_totals, total_max, gap)
    right_pos = _compute_positions(right_nodes, right_totals, total_max, gap)

    if cmap_left:
        _left_cmap = plt.get_cmap(cmap_left)
    else:
        _left_cmap = LinearSegmentedColormap.from_list('left_pink_blue', ['#F3C1D3', '#A8D8EA'])
    if cmap_right:
        _right_cmap = plt.get_cmap(cmap_right)
    else:
        _right_cmap = LinearSegmentedColormap.from_list('right_blue_pink', ['#A8D8EA', '#F3C1D3'])

    non_others_left = [n for n in left_nodes if n != others_label]
    left_colors = {}
    for i, n in enumerate(non_others_left):
        left_colors[n] = _left_cmap(i / max(len(non_others_left) - 1, 1))
    if others_label in left_nodes:
        left_colors[others_label] = (0.85, 0.85, 0.85, 1.0)

    non_others_right = [n for n in right_nodes if n != others_label]
    right_colors = {}
    for i, n in enumerate(non_others_right):
        right_colors[n] = _right_cmap(i / max(len(non_others_right) - 1, 1))
    if others_label in right_nodes:
        right_colors[others_label] = (0.85, 0.85, 0.85, 1.0)

    n_max_nodes = max(len(left_nodes), len(right_nodes))
    fig_h = max(8, n_max_nodes * 0.35)
    fig, ax = plt.subplots(figsize=(10, fig_h))

    fontsize = max(6, min(9, 200 // max(n_max_nodes, 1)))

    for node in left_nodes:
        y0, y1 = left_pos[node]
        ax.barh(y=(y0 + y1) / 2, width=bar_w, height=y1 - y0,
                left=x_left, color=left_colors[node],
                edgecolor='white', linewidth=0.5)
        count = left_totals[node]
        name = node[3:] if node.startswith('SG ') else node
        ax.text(x_left - 0.01, (y0 + y1) / 2,
                f"{name}:{count}",
                ha='right', va='center', fontsize=fontsize, fontweight='bold')

    for node in right_nodes:
        y0, y1 = right_pos[node]
        ax.barh(y=(y0 + y1) / 2, width=bar_w, height=y1 - y0,
                left=x_right, color=right_colors[node],
                edgecolor='white', linewidth=0.5)
        count = right_totals[node]
        name = node[3:] if node.startswith('SG ') else node
        ax.text(x_right + bar_w + 0.01, (y0 + y1) / 2,
                f"{name}:{count}",
                ha='left', va='center', fontsize=fontsize, fontweight='bold')

    left_cursor = {n: left_pos[n][0] for n in left_nodes}
    right_cursor = {n: right_pos[n][0] for n in right_nodes}
    x_vals = np.linspace(0, 1, 200)

    for (l_node, r_node), cnt in sorted(top_flows, key=lambda x: -x[1]):
        if l_node not in left_pos or r_node not in right_pos:
            continue
        frac = cnt / total_max * (1 - gap * (n_max_nodes - 1))

        l_y0 = left_cursor[l_node]
        l_y1 = l_y0 + frac
        left_cursor[l_node] = l_y1

        r_y0 = right_cursor[r_node]
        r_y1 = r_y0 + frac
        right_cursor[r_node] = r_y1

        t = x_vals
        s = 3 * t ** 2 - 2 * t ** 3
        top_curve = l_y1 + s * (r_y1 - l_y1)
        bot_curve = l_y0 + s * (r_y0 - l_y0)

        x_draw = (x_left + bar_w) + t * (x_right - x_left - bar_w)
        ax.fill_between(x_draw, bot_curve, top_curve,
                        alpha=0.35, color=left_colors[l_node], edgecolor='none')

    ax.set_xlim(-0.15, x_right + bar_w + 0.15)
    y_max = max(
        max((p[1] for p in left_pos.values()), default=1),
        max((p[1] for p in right_pos.values()), default=1),
    )
    ax.set_ylim(-0.05, y_max + 0.05)
    ax.invert_yaxis()
    ax.axis('off')

    ax.text(x_left + bar_w / 2, -0.03, "Initial",
            ha='center', va='bottom', fontsize=12, fontweight='bold')
    ax.text(x_right + bar_w / 2, -0.03, "Relaxed",
            ha='center', va='bottom', fontsize=12, fontweight='bold')

    fig.set_facecolor('white')
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f"[PT-viz] Saved: {output_path}")

def plot_violin_energy(cif_paths_relaxed, output_path, energy_data=None,
                       symprec=0.1, dpi=300, ignore_p1=False,
                       workers=1, top_n=20):
    """
    Draw a violin plot showing energy_per_atom distribution grouped by relaxed space group.

    Args:
        cif_paths_relaxed: list of relaxed CIF file paths.
        output_path: Output PNG path.
        energy_data: dict from _load_energy_data().
        symprec: Symmetry precision.
        dpi: Figure DPI.
        ignore_p1: Exclude P1 from plot.
        workers: Number of parallel workers.
        top_n: Show top-N most populated space groups.
    """
    if not cif_paths_relaxed:
        print("[PT-viz] Warning: no relaxed CIF files for violin plot")
        return

    if not energy_data:
        print("[PT-viz] Warning: no energy data found, skipping violin plot")
        return

    cif_energy_map = {}
    for cif_path in cif_paths_relaxed:
        matched = _match_energy_for_cif(cif_path, energy_data)
        if matched and matched.get('energy_per_atom') is not None:
            cif_energy_map[cif_path] = matched['energy_per_atom']

    if not cif_energy_map:
        print("[PT-viz] Warning: no energy data matched to relaxed CIFs, skipping violin")
        return

    paths_with_energy = list(cif_energy_map.keys())
    print(f"[PT-viz] Matched energy for {len(paths_with_energy)} / {len(cif_paths_relaxed)} relaxed CIFs")

    sg_map = {}
    if workers > 1 and len(paths_with_energy) > 1:
        tasks = [(p, symprec) for p in paths_with_energy]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(tqdm(
                pool.map(_analyse_sg_worker, tasks,
                         chunksize=max(1, len(tasks) // workers // 4)),
                total=len(tasks), desc="[PT-viz] Analysing space groups (violin)"))
        for cif_path, rel_sym, rel_num in results:
            sg_map[cif_path] = rel_num
    else:
        for cif_path in tqdm(paths_with_energy, desc="[PT-viz] Analysing space groups (violin)"):
            _, rel_num = _identify_sg_from_structure(cif_path, symprec=symprec)
            sg_map[cif_path] = rel_num

    sg_energies = defaultdict(list)
    for cif_path in paths_with_energy:
        if cif_path not in sg_map:
            continue
        rel_num = sg_map[cif_path]
        sg_energies[rel_num].append(cif_energy_map[cif_path])

    if ignore_p1 and 1 in sg_energies:
        del sg_energies[1]

    if not sg_energies:
        print("[PT-viz] Warning: no data for violin plot")
        return

    sorted_by_count = sorted(sg_energies.items(), key=lambda x: len(x[1]), reverse=True)
    selected = sorted_by_count[:top_n]
    selected.sort(key=lambda x: x[0])

    sg_labels = [f"SG {sg}" for sg, _ in selected]
    sg_values = [vals for _, vals in selected]
    sg_counts = [len(vals) for vals in sg_values]

    print(f"[PT-viz] Violin plot: {len(selected)} space groups, "
          f"{sum(sg_counts)} structures total")

    n_sg = len(selected)
    fig_w = max(10, n_sg * 0.6)
    fig, ax = plt.subplots(figsize=(fig_w, 6))

    positions = list(range(n_sg))
    cmap = plt.cm.get_cmap('Pastel2', n_sg)

    all_values = [v for vals in sg_values for v in vals]
    overall_median = float(np.median(all_values))
    ax.axhline(y=overall_median, color='#888888', linestyle='--', linewidth=0.8,
               label=f'Overall median ({overall_median:.3f})', zorder=1)

    rng = np.random.default_rng(42)

    for i, ((_sg_num, vals), count) in enumerate(zip(selected, sg_counts)):
        if len(vals) < 2:
            ax.scatter([i], vals, s=20, color=cmap(i), edgecolors='gray',
                       linewidths=0.5, zorder=4)
        else:
            parts = ax.violinplot([vals], positions=[i],
                                  showmeans=True, showmedians=True,
                                  widths=0.7)
            for pc in parts['bodies']:
                pc.set_facecolor(cmap(i))
                pc.set_edgecolor('gray')
                pc.set_alpha(0.7)
            for key in ('cmeans', 'cmedians', 'cbars', 'cmins', 'cmaxes'):
                if key in parts:
                    parts[key].set_edgecolor('gray')
                    parts[key].set_linewidth(0.8)

        jitter = rng.uniform(-0.15, 0.15, len(vals))
        ax.scatter(np.array([i] * len(vals)) + jitter, vals,
                   s=3, alpha=0.3, color='gray', zorder=3)

        y_top = max(vals)
        ax.text(i, y_top, f"n={count}", ha='center', va='bottom',
                fontsize=7, color='#555555')

    ax.set_xticks(positions)
    ax.set_xticklabels(sg_labels, rotation=45, ha='right', fontweight='bold')
    ax.set_ylabel('Energy (eV/atom)', fontweight='bold')
    ax.set_xlabel('Space Group', fontweight='bold')
    ax.legend(loc='upper right', fontsize=8)

    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"[PT-viz] Saved: {output_path}")

def plot_rdf_comparison(cif_dict, output_path, r_max=10.0, bin_size=0.1,
                        max_structures=200, dpi=300, workers=1):
    """
    Plot average RDF comparison between initial and relaxed structures.

    Args:
        cif_dict: {'initial': [paths], 'relaxed': [paths]}
        output_path: Output PNG path.
        r_max: Maximum radius.
        bin_size: Bin width.
        max_structures: Cap per stage.
        dpi: Figure DPI.
        workers: Number of parallel workers.
    """
    stage_colors = {
        'initial': '#AAAAAA',
        'relaxed': '#4ECDC4',
    }
    stage_order = ['initial', 'relaxed']

    fig, ax = plt.subplots(figsize=(10, 6))
    has_data = False

    for stage in stage_order:
        paths = cif_dict.get(stage, [])
        if not paths:
            continue

        if len(paths) > max_structures:
            rng = np.random.default_rng(42)
            paths = list(rng.choice(paths, size=max_structures, replace=False))

        all_gr = []
        r_vals = None

        if workers > 1 and len(paths) > 1:
            tasks = [(p, r_max, bin_size) for p in paths]
            with ProcessPoolExecutor(max_workers=workers) as pool:
                rdf_results = list(tqdm(
                    pool.map(_compute_rdf_worker, tasks,
                             chunksize=max(1, len(tasks) // workers // 4)),
                    total=len(tasks), desc=f"[PT-viz] RDF {stage}"))
            for r, gr, err in rdf_results:
                if err is not None:
                    print(f"[PT-viz] Warning: skipping {err}")
                else:
                    all_gr.append(gr)
                    if r_vals is None:
                        r_vals = r
        else:
            from pymatgen.core import Structure
            for p in tqdm(paths, desc=f"[PT-viz] RDF {stage}"):
                try:
                    struct = Structure.from_file(p)
                    r, gr = _compute_rdf(struct, r_max=r_max, bin_size=bin_size)
                    all_gr.append(gr)
                    if r_vals is None:
                        r_vals = r
                except Exception as e:
                    print(f"[PT-viz] Warning: skipping {os.path.basename(p)}: {e}")

        if not all_gr:
            print(f"[PT-viz] Warning: no valid RDF data for stage '{stage}'")
            continue

        gr_array = np.array(all_gr)
        mean_gr = gr_array.mean(axis=0)
        std_gr = gr_array.std(axis=0)

        color = stage_colors.get(stage, '#999999')
        label = f"{stage.capitalize()} (n={len(all_gr)})"
        ax.plot(r_vals, mean_gr, color=color, linewidth=1.8, label=label)
        ax.fill_between(r_vals, mean_gr - std_gr, mean_gr + std_gr,
                        color=color, alpha=0.2)

        if len(mean_gr) > 2:
            peak_idx = np.argmax(mean_gr[1:]) + 1
            ax.annotate(f'{r_vals[peak_idx]:.2f} \u00c5',
                        xy=(r_vals[peak_idx], mean_gr[peak_idx]),
                        fontsize=9, ha='center', va='bottom',
                        color=color, fontweight='bold')
        has_data = True

    if not has_data:
        print("[PT-viz] Warning: no RDF data to plot")
        plt.close(fig)
        return

    ax.set_xlabel(r'$r$ ($\mathrm{\AA}$)', fontsize=13)
    ax.set_ylabel(r'$g(r)$', fontsize=13)
    ax.set_xlim(left=0)
    ax.legend(fontsize=11, framealpha=0.9)

    ax.tick_params(labelsize=11)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.set_facecolor('white')
    fig.set_facecolor('white')

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f"[PT-viz] Saved: {output_path}")

def generate_report(cif_dict, energy_data, output_path):
    """Write a summary report."""
    lines = [
        "=" * 60,
        "GEWUM PT Visualization Report",
        "=" * 60,
        f"Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "--- Dataset Statistics ---",
    ]
    for stage, paths in cif_dict.items():
        lines.append(f"  {stage}: {len(paths)} CIF files")

    lines.append(f"  Energy entries: {len(energy_data)}")
    lines.append("")
    lines.append("--- Output Files ---")

    output_dir = os.path.dirname(output_path)
    if os.path.isdir(output_dir):
        for fname in sorted(os.listdir(output_dir)):
            fpath = os.path.join(output_dir, fname)
            if os.path.isfile(fpath):
                size_kb = os.path.getsize(fpath) / 1024
                lines.append(f"  {fname}  ({size_kb:.1f} KB)")

    lines.append("=" * 60)

    with open(output_path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f"[PT-viz] Report saved: {output_path}")

def _build_parser():
    parser = argparse.ArgumentParser(
        description='GEWUM PT Visualization - Space group Sankey, energy violin, RDF comparison',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python -m gewum.src.PTworkflows.viz_pt --cif-dir .
  python -m gewum.src.PTworkflows.viz_pt --cif-dir . --plot sankey --ignore-p1
  python -m gewum.src.PTworkflows.viz_pt --cif-dir . --plot violin --violin-top-n 15
  python -m gewum.src.PTworkflows.viz_pt --cif-dir ./dir1 ./dir2 --plot rdf
""",
    )
    parser.add_argument('--cif-dir', nargs='+', required=True,
                        help='One or more PT work directories')
    parser.add_argument('--plot', choices=['all', 'sankey', 'violin', 'rdf'],
                        default='all',
                        help='Which plot(s) to generate (default: all)')
    parser.add_argument('--symprec', type=float, default=0.1,
                        help='Symmetry precision for spglib (default: 0.1)')
    parser.add_argument('--rdf-rmax', type=float, default=10.0,
                        help='RDF maximum radius in Angstrom (default: 10.0)')
    parser.add_argument('--rdf-bin', type=float, default=0.1,
                        help='RDF bin size in Angstrom (default: 0.1)')
    parser.add_argument('--max-structures', type=int, default=200,
                        help='Max structures per stage for RDF (default: 200)')
    parser.add_argument('--ignore-p1', action='store_true', default=False,
                        help='Exclude P1 (SG 1) from Sankey and violin plots')
    parser.add_argument('--violin-top-n', type=int, default=20,
                        help='Number of top space groups in violin plot (default: 20)')
    parser.add_argument('--sankey-cmap-left', type=str, default=None,
                        help='Colormap for left column (initial SG) in sankey_spacegroup (e.g. viridis)')
    parser.add_argument('--sankey-cmap-right', type=str, default=None,
                        help='Colormap for right column (relaxed SG) in sankey_spacegroup (e.g. viridis_r)')
    parser.add_argument('-o', '--output', default='./viz_pt_output',
                        help='Output directory (default: ./viz_pt_output)')
    parser.add_argument('--dpi', type=int, default=300,
                        help='Figure DPI (default: 300)')
    parser.add_argument('--font-family', type=str, default='Liberation Sans',
                        help='Font family for figures (default: Liberation Sans)')
    parser.add_argument('--font-size', type=int, default=None,
                        help='Global font size for figures (default: matplotlib default)')
    parser.add_argument('--workers', type=int, default=min(os.cpu_count() or 4, 4),
                        help='Number of parallel workers (default: min(cpu_count, 4))')
    return parser


def main():
    parser = _build_parser()
    args = parser.parse_args()

    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = [args.font_family] + matplotlib.rcParams['font.sans-serif']
    if args.font_size is not None:
        matplotlib.rcParams['font.size'] = args.font_size

    os.makedirs(args.output, exist_ok=True)

    cif_dict = _collect_cifs(args.cif_dir)
    if not cif_dict:
        print("[PT-viz] Error: no CIF files found in any provided directory.")
        sys.exit(1)

    energy_data = _load_energy_data(args.cif_dir)

    plot_type = args.plot
    workers = args.workers
    print(f"[PT-viz] Workers: {workers}")

    if plot_type in ('all', 'sankey'):
        if 'relaxed' in cif_dict and cif_dict['relaxed']:
            sankey_path = os.path.join(args.output, 'sankey_spacegroup.png')
            plot_sankey(cif_dict['relaxed'], sankey_path,
                        energy_data=energy_data,
                        symprec=args.symprec, dpi=args.dpi,
                        ignore_p1=args.ignore_p1, workers=workers,
                        cmap_left=args.sankey_cmap_left,
                        cmap_right=args.sankey_cmap_right)
        else:
            print("[PT-viz] Skipping Sankey: no relaxed structures found")

    if plot_type in ('all', 'violin'):
        if 'relaxed' in cif_dict and cif_dict['relaxed']:
            violin_path = os.path.join(args.output, 'violin_energy.png')
            plot_violin_energy(cif_dict['relaxed'], violin_path,
                               energy_data=energy_data,
                               symprec=args.symprec, dpi=args.dpi,
                               ignore_p1=args.ignore_p1, workers=workers,
                               top_n=args.violin_top_n)
        else:
            print("[PT-viz] Skipping Violin: no relaxed structures found")

    if plot_type in ('all', 'rdf'):
        rdf_path = os.path.join(args.output, 'rdf_comparison.png')
        plot_rdf_comparison(cif_dict, rdf_path,
                            r_max=args.rdf_rmax, bin_size=args.rdf_bin,
                            max_structures=args.max_structures,
                            dpi=args.dpi, workers=workers)

    report_path = os.path.join(args.output, 'viz_pt_report.txt')
    generate_report(cif_dict, energy_data, report_path)

    print(f"\n[PT-viz] All outputs written to {os.path.abspath(args.output)}")

    import gc
    plt.close('all')
    gc.collect()


if __name__ == '__main__':
    main()
