"""
Standalone formation energy phase diagram plotting module.

Supports binary (2-element) convex hull diagrams, ternary (3-element) triangle diagrams,
and quaternary (4-element) tetrahedral diagrams.

Usage:
    python -m gewum.src.common.postprocess.phase_diagram --dir . --mp-data /path/to/MP.json -o phase_diagram.png
"""

import argparse
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize as mplNormalize
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd


# Font configuration
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'DejaVu Sans', 'Arial', 'sans-serif']
plt.rcParams['font.size'] = 14


def _get_reduced_formula(comp_str: str) -> str:
    """Get reduced formula using pymatgen (Be8Pt8 -> BePt)."""
    try:
        from pymatgen.core import Composition
        return Composition(comp_str).reduced_formula
    except Exception:
        return comp_str


def _format_composition_label(comp: str) -> str:
    """Convert composition string to matplotlib LaTeX subscript format (upright).

    Uses reduced formula first (Be8Pt8 -> BePt), then applies subscript.
    Example: 'Li6Al2C3' -> '$\\mathrm{Li_{6}Al_{2}C_{3}}$', 'BePt' -> '$\\mathrm{BePt}$'
    """
    reduced = _get_reduced_formula(comp)
    formatted = re.sub(r'(\d+)', r'_{\1}', reduced)
    return f'$\\mathrm{{{formatted}}}$'


def _ternary_to_cartesian(a, b, c):
    """Convert ternary coordinates (a, b, c) to 2D Cartesian (x, y).

    a = bottom-left vertex, b = bottom-right vertex, c = top vertex.
    """
    total = a + b + c
    x = 0.5 * (2 * b + c) / total
    y = (np.sqrt(3) / 2) * c / total
    return x, y


def _quaternary_to_cartesian3d(fractions):
    """Convert quaternary mole fractions to 3D Cartesian coordinates on a regular tetrahedron.

    Vertices of the tetrahedron:
        element[0] -> (0, 0, 0)
        element[1] -> (1, 0, 0)
        element[2] -> (0.5, sqrt(3)/2, 0)
        element[3] -> (0.5, sqrt(3)/6, sqrt(6)/3)

    Args:
        fractions: (N, 4) array of mole fractions summing to 1.

    Returns:
        x, y, z: 1D arrays of shape (N,).
    """
    f = np.asarray(fractions)
    x = f[:, 1] + f[:, 2] * 0.5 + f[:, 3] * 0.5
    y = f[:, 2] * np.sqrt(3) / 2 + f[:, 3] * np.sqrt(3) / 6
    z = f[:, 3] * np.sqrt(6) / 3
    return x, y, z


def lower_convex_hull(points_sorted_by_x):
    """Compute lower convex hull of points sorted by x using Andrew's monotone chain."""
    hull = []
    for p in points_sorted_by_x:
        while len(hull) >= 2:
            o, a = hull[-2], hull[-1]
            cross = (a[0] - o[0]) * (p[1] - o[1]) - (a[1] - o[1]) * (p[0] - o[0])
            if cross <= 0:
                hull.pop()
            else:
                break
        hull.append(p)
    return np.array(hull)


def main():
    parser = argparse.ArgumentParser(
        description='Formation energy phase diagram plotting (binary / ternary / quaternary).'
    )
    parser.add_argument('-i', '--input', default='./Hull_result.csv',
                        help='Path to input CSV file (default: ./Hull_result.csv)')
    parser.add_argument('--mp-data', default=None, help='Path to offline MP JSON data file (optional)')
    parser.add_argument('--output', '-o', default='phase_diagram.png', help='Output file name (default: phase_diagram.png)')
    args = parser.parse_args()

    csv_path = args.input
    if not os.path.isfile(csv_path):
        print(f"Error: {csv_path} not found.")
        sys.exit(1)

    df = pd.read_csv(csv_path)
    if 'Chemical_Formula' not in df.columns or 'formation_energy_per_atom' not in df.columns:
        print("Error: Hull_result.csv must contain 'Chemical_Formula' and 'formation_energy_per_atom' columns.")
        sys.exit(1)

    df = df[df['formation_energy_per_atom'].notna()].copy()
    if len(df) == 0:
        print("Error: No valid formation energy data found.")
        sys.exit(1)

    from pymatgen.core import Composition

    all_elements = set()
    for comp_str in df['Chemical_Formula']:
        comp = Composition(str(comp_str))
        all_elements.update(comp.get_el_amt_dict().keys())
    elements = sorted(all_elements)

    if len(elements) < 2 or len(elements) > 4:
        print(f"Phase diagram requires 2-4 elements, found {len(elements)}: {elements}")
        print("Phase diagram generation is only supported for binary, ternary, and quaternary systems.")
        sys.exit(0)

    fractions = []
    for comp_str in df['Chemical_Formula']:
        comp = Composition(str(comp_str))
        el_dict = comp.fractional_composition.get_el_amt_dict()
        frac = [el_dict.get(el, 0.0) for el in elements]
        fractions.append(frac)
    fractions = np.array(fractions)

    mp_stable_data = []
    if args.mp_data and os.path.isfile(args.mp_data):
        try:
            from gewum.src.common.ehull.mp_offline_loader import MPOfflineLoader
            from pymatgen.analysis.phase_diagram import PhaseDiagram

            loader = MPOfflineLoader(args.mp_data)
            # Load as PDEntry - formation energies are self-consistent.
            # Note: In Ehull calculation, MP entries are always used raw (uncorrected)
            # regardless of -cor flag. The -cor only corrects USER entries.
            # So both -cor and non-cor Hull_result.csv share the same MP reference.
            mp_entries = loader.get_entries_in_chemsys(elements, as_pd_entry=True)

            if mp_entries:
                pd_mp = PhaseDiagram(mp_entries)
                for entry in pd_mp.stable_entries:
                    formula = entry.composition.reduced_formula
                    form_energy = pd_mp.get_form_energy_per_atom(entry)
                    comp_obj = entry.composition
                    el_dict = comp_obj.fractional_composition.get_el_amt_dict()
                    frac = [el_dict.get(el, 0.0) for el in elements]
                    # Skip pure elements
                    if sum(1 for f in frac if f > 1e-6) <= 1:
                        continue
                    mp_stable_data.append({
                        'composition': formula,
                        'formation_energy': form_energy,
                        'fractions': frac,
                    })
        except Exception as e:
            print(f"Warning: Failed to load MP data ({e}). Proceeding without MP reference.")

    output_path = args.output
    if not os.path.isabs(output_path):
        output_path = os.path.join(os.path.dirname(csv_path) or '.', output_path)

    user_y = df['formation_energy_per_atom'].values

    if len(elements) == 2:
        _plot_binary(elements, fractions, user_y, df, mp_stable_data, output_path)
    elif len(elements) == 3:
        _plot_ternary(elements, fractions, user_y, df, mp_stable_data, output_path)
    else:
        _plot_quaternary(elements, fractions, user_y, df, mp_stable_data, output_path)

    print(f"Phase diagram saved to: {output_path}")


def _is_on_hull(px, py, hull_pts, tol=1e-4):
    """Check if point (px, py) lies on the lower convex hull within tolerance."""
    for i in range(len(hull_pts) - 1):
        x0, y0 = hull_pts[i]
        x1, y1 = hull_pts[i + 1]
        if x0 - tol <= px <= x1 + tol:
            if abs(x1 - x0) < 1e-12:
                hull_y_at_px = min(y0, y1)
            else:
                t = (px - x0) / (x1 - x0)
                hull_y_at_px = y0 + t * (y1 - y0)
            if py <= hull_y_at_px + tol:
                return True
    return False


def _dedup_by_composition(df, fractions, user_y):
    """Deduplicate: keep only the lowest energy entry per composition."""
    comp_best = {}
    for i, comp_str in enumerate(df['Chemical_Formula']):
        comp_key = str(comp_str)
        if comp_key not in comp_best or user_y[i] < user_y[comp_best[comp_key]]:
            comp_best[comp_key] = i
    best_indices = sorted(comp_best.values())
    return best_indices


def _plot_binary(elements, fractions, user_y, df, mp_stable_data, output_path):
    """Plot binary formation energy convex hull diagram."""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_facecolor('#FAFAFA')

    best_indices = _dedup_by_composition(df, fractions, user_y)
    user_x_all = fractions[:, 1]
    user_x = fractions[best_indices, 1]
    user_y_dedup = user_y[best_indices]
    df_dedup = df.iloc[best_indices].reset_index(drop=True)

    mp_x = np.array([d['fractions'][1] for d in mp_stable_data]) if mp_stable_data else np.array([])
    mp_y = np.array([d['formation_energy'] for d in mp_stable_data]) if mp_stable_data else np.array([])

    hull_x = np.concatenate([[0.0], mp_x, [1.0]])
    hull_y = np.concatenate([[0.0], mp_y, [0.0]])
    sorted_idx = np.argsort(hull_x)
    sx = hull_x[sorted_idx]
    sy = hull_y[sorted_idx]
    pts_sorted = np.column_stack([sx, sy])
    lower_hull_pts = lower_convex_hull(pts_sorted)

    ax.plot(lower_hull_pts[:, 0], lower_hull_pts[:, 1],
            color='#AAAAAA', linewidth=1.5, linestyle='--', zorder=3)

    combined_hull_x = np.concatenate([[0.0], mp_x, user_x, [1.0]])
    combined_hull_y = np.concatenate([[0.0], mp_y, user_y_dedup, [0.0]])
    combined_sorted_idx = np.argsort(combined_hull_x)
    combined_pts = np.column_stack([combined_hull_x[combined_sorted_idx],
                                    combined_hull_y[combined_sorted_idx]])
    combined_lower_hull = lower_convex_hull(combined_pts)

    ax.plot(combined_lower_hull[:, 0], combined_lower_hull[:, 1],
            color='#222222', linewidth=2.0, linestyle='-', zorder=4)

    user_on_combined_hull = np.array([_is_on_hull(user_x[i], user_y_dedup[i], combined_lower_hull)
                                      for i in range(len(user_x))])
    stable_mask = user_on_combined_hull & (user_y_dedup < -0.01)
    unstable_mask = ~stable_mask
    mp_formula_set = set(_get_reduced_formula(d['composition']) for d in mp_stable_data)

    if len(mp_x) > 0:
        ax.scatter(mp_x, mp_y, c='#89CFF0', s=90, marker='s',
                   edgecolors='#5B9BD5', linewidths=1.0, zorder=5,
                   label='Stable (MP)')
        for d in mp_stable_data:
            if d['formation_energy'] < -0.01:
                label = _format_composition_label(d['composition'])
                ax.annotate(label,
                            xy=(d['fractions'][1], d['formation_energy']),
                            xytext=(5, 6), textcoords='offset points',
                            fontsize=10, color='#3A7CA5',
                            bbox=dict(boxstyle='round,pad=0.2',
                                      facecolor='#E8F4FD',
                                      edgecolor='#89CFF0', alpha=0.85))

    if unstable_mask.any():
        ax.scatter(user_x[unstable_mask], user_y_dedup[unstable_mask],
                   c='#FFB6C1', s=40, marker='o', alpha=0.4,
                   edgecolors='none', zorder=3,
                   label='Calculated')

    if stable_mask.any():
        stable_indices = np.where(stable_mask)[0]
        confirmed_idx = []
        novel_idx = []
        for idx in stable_indices:
            comp_name = str(df_dedup.iloc[idx]['Chemical_Formula'])
            reduced = _get_reduced_formula(comp_name)
            if reduced in mp_formula_set:
                confirmed_idx.append(idx)
            else:
                novel_idx.append(idx)

        if confirmed_idx:
            cidx = np.array(confirmed_idx)
            ax.scatter(user_x[cidx], user_y_dedup[cidx],
                       c='#9B59B6', s=110, marker='*',
                       edgecolors='#6C3483', linewidths=0.8, zorder=7,
                       label='Confirmed (MP + Predicted)')
            for idx in confirmed_idx:
                comp_name = str(df_dedup.iloc[idx]['Chemical_Formula'])
                label = _format_composition_label(comp_name)
                ax.annotate(label, xy=(user_x[idx], user_y_dedup[idx]),
                            xytext=(5, -14), textcoords='offset points',
                            fontsize=10, fontweight='bold', color='#6C3483',
                            bbox=dict(boxstyle='round,pad=0.2',
                                      facecolor='#F4ECF7',
                                      edgecolor='#9B59B6', alpha=0.85))

        if novel_idx:
            nidx = np.array(novel_idx)
            ax.scatter(user_x[nidx], user_y_dedup[nidx],
                       c='#FF69B4', s=80, marker='o',
                       edgecolors='#C71585', linewidths=1.0, zorder=6,
                       label='Predicted (stable)')
            for idx in novel_idx:
                comp_name = str(df_dedup.iloc[idx]['Chemical_Formula'])
                label = _format_composition_label(comp_name)
                ax.annotate(label, xy=(user_x[idx], user_y_dedup[idx]),
                            xytext=(5, -14), textcoords='offset points',
                            fontsize=10, fontweight='bold', color='#C71585',
                            bbox=dict(boxstyle='round,pad=0.2',
                                      facecolor='#FFF0F5',
                                      edgecolor='#FF69B4', alpha=0.85))

    ax.axhline(y=0, color='#AAAAAA', linestyle='--', linewidth=0.8, alpha=0.6)

    ax.scatter([0, 1], [0, 0], c='#333333', s=100, marker='D', zorder=7)
    ax.text(0, 0.01, elements[0], ha='center', va='bottom',
            fontsize=13, fontweight='bold', color='#333333')
    ax.text(1, 0.01, elements[1], ha='center', va='bottom',
            fontsize=13, fontweight='bold', color='#333333')

    ax.set_xlabel(f'Mole fraction of {elements[1]}', fontsize=14)
    ax.set_ylabel('Formation Energy (eV/atom)', fontsize=14)
    ax.set_xlim(-0.05, 1.05)

    # Dynamic y-axis limits: focus on physically meaningful energy range
    # Collect all plotted energies (user + MP)
    all_plotted_y = list(user_y_dedup)
    if len(mp_y) > 0:
        all_plotted_y.extend(mp_y.tolist())
    all_plotted_y = np.array(all_plotted_y)
    y_min_data = np.min(all_plotted_y)
    # Upper limit: cap at max(0.5, 30% of |y_min|) but never exceed 2.0 eV/atom
    # This keeps focus on the hull region while showing some unstable context
    y_upper = min(max(0.5, -0.3 * y_min_data), 2.0)
    # Lower limit: 10% padding below the minimum
    y_lower = y_min_data * 1.15 if y_min_data < 0 else y_min_data - 0.1
    ax.set_ylim(y_lower, y_upper)

    ax.legend(loc='best', framealpha=0.9, fontsize=11)
    ax.grid(True, alpha=0.15, linestyle='-')

    plt.tight_layout()
    fig.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close(fig)


def _draw_ternary_hull_lines(ax, pd_obj, elements):
    """Draw convex hull tie-lines on ternary diagram from PhaseDiagram facets."""
    from pymatgen.core import Composition

    if hasattr(pd_obj, 'qhull_entries'):
        qhull_entries = pd_obj.qhull_entries
    elif hasattr(pd_obj, '_qhull_entries'):
        qhull_entries = pd_obj._qhull_entries
    else:
        qhull_entries = sorted(pd_obj.stable_entries,
                               key=lambda e: e.composition.reduced_formula)

    entry_coords = {}
    for idx, entry in enumerate(qhull_entries):
        el_dict = entry.composition.fractional_composition.get_el_amt_dict()
        a = el_dict.get(elements[0], 0.0)
        b = el_dict.get(elements[1], 0.0)
        c = el_dict.get(elements[2], 0.0)
        x, y = _ternary_to_cartesian(np.array([a]), np.array([b]), np.array([c]))
        entry_coords[idx] = (float(x[0]), float(y[0]))

    drawn_edges = set()
    for facet in pd_obj.facets:
        for i in range(len(facet)):
            for j in range(i + 1, len(facet)):
                edge_key = tuple(sorted([facet[i], facet[j]]))
                if edge_key in drawn_edges:
                    continue
                drawn_edges.add(edge_key)
                if facet[i] in entry_coords and facet[j] in entry_coords:
                    x1, y1 = entry_coords[facet[i]]
                    x2, y2 = entry_coords[facet[j]]
                    ax.plot([x1, x2], [y1, y2], color='#555555',
                            linewidth=1.2, linestyle='-', alpha=0.7, zorder=2)


def _plot_ternary(elements, fractions, user_y, df, mp_stable_data, output_path):
    """Plot ternary formation energy triangle diagram."""
    fig, ax = plt.subplots(figsize=(9, 8))
    fig.set_facecolor('white')
    ax.set_facecolor('white')

    best_indices = _dedup_by_composition(df, fractions, user_y)
    fractions_dedup = fractions[best_indices]
    user_y_dedup = user_y[best_indices]
    df_dedup = df.iloc[best_indices].reset_index(drop=True)

    triangle = plt.Polygon([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]],
                           fill=False, edgecolor='#555555', linewidth=1.5)
    ax.add_patch(triangle)

    user_cart_x, user_cart_y = _ternary_to_cartesian(
        fractions_dedup[:, 0], fractions_dedup[:, 1], fractions_dedup[:, 2])
    user_fe = user_y_dedup

    user_stable_mask = np.zeros(len(user_fe), dtype=bool)
    pd_full = None
    try:
        from pymatgen.core import Composition
        from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry

        user_entries = []
        for i in range(len(df_dedup)):
            comp = Composition(str(df_dedup.iloc[i]['Chemical_Formula']))
            energy = user_fe[i] * comp.num_atoms
            entry = PDEntry(comp, energy, name=f'user_{i}')
            user_entries.append(entry)

        all_entries = user_entries[:]
        if mp_stable_data:
            for d in mp_stable_data:
                comp = Composition(d['composition'])
                energy = d['formation_energy'] * comp.num_atoms
                all_entries.append(PDEntry(comp, energy, name='mp'))
        for el in elements:
            all_entries.append(PDEntry(Composition(el), 0.0, name=f'pure_{el}'))

        pd_full = PhaseDiagram(all_entries)
        for i, entry in enumerate(user_entries):
            e_hull = pd_full.get_e_above_hull(entry)
            if e_hull < 1e-4 and user_fe[i] < -0.01:
                user_stable_mask[i] = True

        stable_entries_for_hull = []
        for i in range(len(df_dedup)):
            if user_stable_mask[i]:
                comp = Composition(str(df_dedup.iloc[i]['Chemical_Formula']))
                energy = user_fe[i] * comp.num_atoms
                stable_entries_for_hull.append(PDEntry(comp, energy, name=f'user_{i}'))
        if mp_stable_data:
            for d in mp_stable_data:
                comp = Composition(d['composition'])
                energy = d['formation_energy'] * comp.num_atoms
                stable_entries_for_hull.append(PDEntry(comp, energy, name='mp'))
        for el in elements:
            stable_entries_for_hull.append(PDEntry(Composition(el), 0.0, name=f'pure_{el}'))

        if len(stable_entries_for_hull) >= 4:  
            pd_hull = PhaseDiagram(stable_entries_for_hull)
            _draw_ternary_hull_lines(ax, pd_hull, elements)

    except Exception as e:
        import traceback
        traceback.print_exc()
        top_k = min(5, len(user_fe))
        top_indices = np.argsort(user_fe)[:top_k]
        user_stable_mask[top_indices] = True

    unstable_mask = ~user_stable_mask

    mp_formula_set = set(_get_reduced_formula(d['composition']) for d in mp_stable_data)

    if mp_stable_data:
        mp_fracs = np.array([d['fractions'] for d in mp_stable_data])
        mp_cart_x, mp_cart_y = _ternary_to_cartesian(
            mp_fracs[:, 0], mp_fracs[:, 1], mp_fracs[:, 2])

        ax.scatter(mp_cart_x, mp_cart_y, c='#89CFF0', s=90, marker='s',
                   edgecolors='#5B9BD5', linewidths=1.0, zorder=4,
                   label='Stable (MP)')
        for i, d in enumerate(mp_stable_data):
            if d['formation_energy'] < -0.01:
                label = _format_composition_label(d['composition'])
                ax.annotate(label, xy=(mp_cart_x[i], mp_cart_y[i]),
                            xytext=(5, 5), textcoords='offset points',
                            fontsize=10, color='#3A7CA5',
                            bbox=dict(boxstyle='round,pad=0.2',
                                      facecolor='#E8F4FD',
                                      edgecolor='#89CFF0', alpha=0.85))

    if unstable_mask.any():
        ax.scatter(user_cart_x[unstable_mask], user_cart_y[unstable_mask],
                   c='#FFB6C1', s=40, marker='o', alpha=0.4,
                   edgecolors='none', zorder=3,
                   label='Calculated')

    if user_stable_mask.any():
        stable_indices = np.where(user_stable_mask)[0]
        confirmed_idx = []
        novel_idx = []
        for idx in stable_indices:
            comp_name = str(df_dedup.iloc[idx]['Chemical_Formula'])
            reduced = _get_reduced_formula(comp_name)
            if reduced in mp_formula_set:
                confirmed_idx.append(idx)
            else:
                novel_idx.append(idx)

        if confirmed_idx:
            cidx = np.array(confirmed_idx)
            ax.scatter(user_cart_x[cidx], user_cart_y[cidx],
                       c='#9B59B6', s=110, marker='*',
                       edgecolors='#6C3483', linewidths=0.8, zorder=7,
                       label='Confirmed (MP + Predicted)')
            for idx in confirmed_idx:
                comp_name = str(df_dedup.iloc[idx]['Chemical_Formula'])
                label = _format_composition_label(comp_name)
                ax.annotate(label, xy=(user_cart_x[idx], user_cart_y[idx]),
                            xytext=(5, -10), textcoords='offset points',
                            fontsize=10, fontweight='bold', color='#6C3483',
                            bbox=dict(boxstyle='round,pad=0.2',
                                      facecolor='#F4ECF7',
                                      edgecolor='#9B59B6', alpha=0.85))

        if novel_idx:
            nidx = np.array(novel_idx)
            ax.scatter(user_cart_x[nidx], user_cart_y[nidx],
                       c='#FF69B4', s=80, marker='o',
                       edgecolors='#C71585', linewidths=1.0, zorder=6,
                       label='Predicted (stable)')
            for idx in novel_idx:
                comp_name = str(df_dedup.iloc[idx]['Chemical_Formula'])
                label = _format_composition_label(comp_name)
                ax.annotate(label, xy=(user_cart_x[idx], user_cart_y[idx]),
                            xytext=(5, -10), textcoords='offset points',
                            fontsize=10, fontweight='bold', color='#C71585',
                            bbox=dict(boxstyle='round,pad=0.2',
                                      facecolor='#FFF0F5',
                                      edgecolor='#FF69B4', alpha=0.85))

    offset = 0.03
    ax.text(0 - offset, 0 - offset, elements[0], fontsize=14,
            fontweight='bold', ha='center', color='#333333')
    ax.text(1 + offset, 0 - offset, elements[1], fontsize=14,
            fontweight='bold', ha='center', color='#333333')
    ax.text(0.5, np.sqrt(3) / 2 + offset, elements[2], fontsize=14,
            fontweight='bold', ha='center', color='#333333')

    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, np.sqrt(3) / 2 + 0.15)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.patch.set_visible(False)
    ax.legend(loc='upper right', framealpha=0.9, fontsize=11)

    plt.tight_layout()
    fig.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close(fig)


def _plot_quaternary(elements, fractions, user_y, df, mp_stable_data, output_path):
    """Plot quaternary formation energy tetrahedron diagram (3D)."""
    from itertools import combinations
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    fig = plt.figure(figsize=(10, 9))
    ax = fig.add_subplot(111, projection='3d')

    # --- Deduplicate ---
    best_indices = _dedup_by_composition(df, fractions, user_y)
    fractions_dedup = fractions[best_indices]
    user_y_dedup = user_y[best_indices]
    df_dedup = df.iloc[best_indices].reset_index(drop=True)

    # --- Tetrahedron wireframe ---
    # Vertices: element[0]=(0,0,0), element[1]=(1,0,0),
    #           element[2]=(0.5, sqrt3/2, 0), element[3]=(0.5, sqrt3/6, sqrt6/3)
    tet_vertices = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, np.sqrt(3) / 2, 0.0],
        [0.5, np.sqrt(3) / 6, np.sqrt(6) / 3],
    ])
    for i, j in combinations(range(4), 2):
        ax.plot([tet_vertices[i, 0], tet_vertices[j, 0]],
                [tet_vertices[i, 1], tet_vertices[j, 1]],
                [tet_vertices[i, 2], tet_vertices[j, 2]],
                color='#888888', linewidth=1.0, alpha=0.6)

    # Vertex labels
    offsets = [(-0.06, -0.04, -0.02), (0.04, -0.04, -0.02),
               (0.0, 0.04, -0.02), (0.0, 0.0, 0.04)]
    for i, el in enumerate(elements):
        ax.text(tet_vertices[i, 0] + offsets[i][0],
                tet_vertices[i, 1] + offsets[i][1],
                tet_vertices[i, 2] + offsets[i][2],
                el, fontsize=13, fontweight='bold', color='#333333')

    # --- Convert user data to 3D coordinates ---
    user_x3, user_y3, user_z3 = _quaternary_to_cartesian3d(fractions_dedup)

    # --- Stability determination via pymatgen PhaseDiagram ---
    user_stable_mask = np.zeros(len(user_y_dedup), dtype=bool)
    try:
        from pymatgen.core import Composition
        from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry

        user_entries = []
        for i in range(len(df_dedup)):
            comp = Composition(str(df_dedup.iloc[i]['Chemical_Formula']))
            energy = user_y_dedup[i] * comp.num_atoms
            entry = PDEntry(comp, energy, name=f'user_{i}')
            user_entries.append(entry)

        all_entries = user_entries[:]
        if mp_stable_data:
            for d in mp_stable_data:
                comp = Composition(d['composition'])
                energy = d['formation_energy'] * comp.num_atoms
                all_entries.append(PDEntry(comp, energy, name='mp'))
        for el in elements:
            all_entries.append(PDEntry(Composition(el), 0.0, name=f'pure_{el}'))

        pd_full = PhaseDiagram(all_entries)
        for i, entry in enumerate(user_entries):
            e_hull = pd_full.get_e_above_hull(entry)
            if e_hull < 1e-4 and user_y_dedup[i] < -0.01:
                user_stable_mask[i] = True

    except Exception as e:
        import traceback
        traceback.print_exc()
        # Fallback: mark top-5 lowest energy as stable
        top_k = min(5, len(user_y_dedup))
        top_indices = np.argsort(user_y_dedup)[:top_k]
        user_stable_mask[top_indices] = True

    mp_formula_set = set(_get_reduced_formula(d['composition']) for d in mp_stable_data)
    unstable_mask = ~user_stable_mask

    # --- Plot MP stable phases ---
    if mp_stable_data:
        mp_fracs = np.array([d['fractions'] for d in mp_stable_data])
        mp_x3, mp_y3, mp_z3 = _quaternary_to_cartesian3d(mp_fracs)
        ax.scatter(mp_x3, mp_y3, mp_z3, c='#89CFF0', s=90, marker='s',
                   edgecolors='#5B9BD5', linewidths=1.0, depthshade=False,
                   label='Stable (MP)', zorder=5)
        for i, d in enumerate(mp_stable_data):
            if d['formation_energy'] < -0.01:
                label = _format_composition_label(d['composition'])
                ax.text(mp_x3[i], mp_y3[i], mp_z3[i], f'  {label}',
                        fontsize=8, color='#3A7CA5', zorder=6)

    # --- Plot unstable user data ---
    if unstable_mask.any():
        ax.scatter(user_x3[unstable_mask], user_y3[unstable_mask], user_z3[unstable_mask],
                   c='#FFB6C1', s=30, marker='o', alpha=0.4,
                   edgecolors='none', depthshade=False,
                   label='Calculated', zorder=3)

    # --- Plot stable user data (confirmed vs novel) ---
    if user_stable_mask.any():
        stable_indices = np.where(user_stable_mask)[0]
        confirmed_idx = []
        novel_idx = []
        for idx in stable_indices:
            comp_name = str(df_dedup.iloc[idx]['Chemical_Formula'])
            reduced = _get_reduced_formula(comp_name)
            if reduced in mp_formula_set:
                confirmed_idx.append(idx)
            else:
                novel_idx.append(idx)

        if confirmed_idx:
            cidx = np.array(confirmed_idx)
            ax.scatter(user_x3[cidx], user_y3[cidx], user_z3[cidx],
                       c='#9B59B6', s=110, marker='*',
                       edgecolors='#6C3483', linewidths=0.8, depthshade=False,
                       label='Confirmed (MP + Predicted)', zorder=7)
            for idx in confirmed_idx:
                comp_name = str(df_dedup.iloc[idx]['Chemical_Formula'])
                label = _format_composition_label(comp_name)
                ax.text(user_x3[idx], user_y3[idx], user_z3[idx],
                        f'  {label}', fontsize=9, fontweight='bold',
                        color='#6C3483', zorder=8)

        if novel_idx:
            nidx = np.array(novel_idx)
            ax.scatter(user_x3[nidx], user_y3[nidx], user_z3[nidx],
                       c='#FF69B4', s=80, marker='o',
                       edgecolors='#C71585', linewidths=1.0, depthshade=False,
                       label='Predicted (stable)', zorder=6)
            for idx in novel_idx:
                comp_name = str(df_dedup.iloc[idx]['Chemical_Formula'])
                label = _format_composition_label(comp_name)
                ax.text(user_x3[idx], user_y3[idx], user_z3[idx],
                        f'  {label}', fontsize=9, fontweight='bold',
                        color='#C71585', zorder=8)

    # --- Axis settings ---
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_zlabel('')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('none')
    ax.yaxis.pane.set_edgecolor('none')
    ax.zaxis.pane.set_edgecolor('none')
    ax.xaxis.line.set_color('none')
    ax.yaxis.line.set_color('none')
    ax.zaxis.line.set_color('none')
    ax.grid(False)
    ax.legend(loc='upper left', framealpha=0.9, fontsize=10)

    # Set view angle for good visualization
    ax.view_init(elev=20, azim=30)

    plt.tight_layout()
    fig.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
