"""
GEWUM viz2 - Structure Analysis Visualization
Generates space group Sankey diagrams, RDF comparison plots, and structure funnel charts.
"""
import os
import sys
import argparse
import datetime
import matplotlib
import numpy as np
from pathlib import Path
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

from gewum.src.common.cif_archive import (
    find_cifs_db_aware, detect_stages_db_aware,
    load_structure, entry_to_path, entry_basename,
)

def _find_cif_files(directory, mode=None):
    """
    Recursively find CIF files with mode-based filtering.
    Transparently reads from structures.zip if present (no unpack needed).

    Args:
        directory: Root directory to search.
        mode: Collection strategy
            - None or 'all': all CIF files recursively
            - 'total': exclude relaxed/ subdirectories
            - 'selected': exclude remove/ and relaxed/ subdirectories
            - 'relaxed': only from relaxed/ subdirectories
    Returns:
        Sorted list of CIF entries (file paths or (zip_path, member) tuples).
    """
    mode_info = {
        'total': 'excluding relaxed/ subdirectories',
        'selected': 'excluding remove/ and relaxed/ subdirectories',
        'relaxed': 'only from relaxed/ subdirectories (excluding bond_mis/)',
    }
    effective_mode = mode or 'all'
    print(f"[viz2] Searching for CIF files in: {directory} (mode: {effective_mode})")
    if effective_mode in mode_info:
        print(f"[viz2] Filter: {mode_info[effective_mode]}")

    entries = find_cifs_db_aware(directory, mode=mode)
    print(f"[viz2] Found {len(entries)} CIF files")
    return entries


def _auto_collect_cifs(cif_dirs):
    """
    Auto-detect workflow stages for each directory and collect CIF paths.

    Checks for ``remove/`` and ``relaxed/`` sub-directories to infer stages.

    Args:
        cif_dirs: list of directory paths.

    Returns:
        dict  ``{'total': [...], 'selected': [...], 'relaxed': [...]}``
        (only stages that are detected are included)
    """
    dir_info = {}
    has_remove_any = False
    has_relaxed_any = False

    for cif_dir in cif_dirs:
        info = detect_stages_db_aware(cif_dir)
        dir_info[cif_dir] = info
        has_remove_any = has_remove_any or info['remove']
        has_relaxed_any = has_relaxed_any or info['relaxed']

        print(f"[viz2] Auto-detecting workflow stages in: {cif_dir}")
        if info['remove']:
            print("[viz2]   Found 'remove/' \u2192 extracting selected stage")
        if info['relaxed']:
            print("[viz2]   Found 'relaxed/' \u2192 extracting relaxed stage")

    stages = ['total']
    if has_remove_any:
        stages.append('selected')
    if has_relaxed_any:
        stages.append('relaxed')

    stage_paths = {s: [] for s in stages}

    for cif_dir in cif_dirs:
        info = dir_info[cif_dir]
        for stage in stages:
            if stage == 'selected' and not info['remove']:
                continue
            if stage == 'relaxed' and not info['relaxed']:
                continue
            if not has_remove_any and not has_relaxed_any:
                m = None
            else:
                m = stage
            paths = _find_cif_files(cif_dir, mode=m)
            if paths:
                stage_paths[stage].extend(paths)

    print("[viz2] --- Collection summary ---")
    for stage in stages:
        n = len(stage_paths[stage])
        print(f"[viz2]   {stage}: {n} CIF files")

    return {k: v for k, v in stage_paths.items() if v}

def _extract_sg_from_path(cif_path):
    """
    Extract the *initial* space-group number from the CIF path (cifgen dir name).

    Strategy:
    - If 'relaxed' is in the path parts, take the part immediately before it.
    - If 'remove' is in the path parts, take the part immediately before it.
    - Otherwise use the parent directory name.
    - The chosen part must be a pure-digit string (space-group number).

    Returns:
        str - space-group number as string, or 'unknown'.
    """
    vpath = entry_to_path(cif_path)
    parts = Path(vpath).parts

    for anchor in ('relaxed', 'remove'):
        if anchor in parts:
            idx = parts.index(anchor)
            if idx > 0 and parts[idx - 1].isdigit():
                return parts[idx - 1]

    parent = Path(vpath).parent.name
    if parent.isdigit():
        return parent

    for part in reversed(parts):
        if part.isdigit():
            return part

    return 'unknown'


def _identify_sg_from_structure(cif_path, symprec=0.1):
    """
    Identify the space group of a CIF structure using spglib via pymatgen.

    Args:
        cif_path: Path/entry to a CIF file (str or (zip_path, member) tuple).
        symprec: Symmetry precision for SpacegroupAnalyzer.

    Returns:
        (sg_symbol, sg_number) - e.g. ('P2_1/c', 14).
        On failure returns ('P1', 1).
    """
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    try:
        structure = load_structure(cif_path)
        sga = SpacegroupAnalyzer(structure, symprec=symprec)
        sg_symbol = sga.get_space_group_symbol()
        sg_number = sga.get_space_group_number()
        return (sg_symbol, sg_number)
    except Exception as e:
        print(f"[viz2] Warning: spglib failed for {entry_basename(cif_path)}: {e}")
        return ('P1', 1)


def _analyse_sg_worker(args):
    """Worker: identify space group for a single CIF (picklable, module-level)."""
    cif_path, symprec = args
    try:
        init_sg = _extract_sg_from_path(cif_path)
        rel_sym, rel_num = _identify_sg_from_structure(cif_path, symprec=symprec)
        return cif_path, init_sg, rel_sym, rel_num
    except Exception as e:
        print(f"[viz2] Warning: SG worker failed for {entry_basename(cif_path)}: {e}")
        return cif_path, 'unknown', 'P1', 1


def _compute_rdf_worker(args):
    """Worker: compute RDF for a single CIF file (picklable, module-level)."""
    path, r_max, bin_size = args
    try:
        struct = load_structure(path)
        r, gr = _compute_rdf(struct, r_max=r_max, bin_size=bin_size)
        return r, gr, None
    except Exception as e:
        return None, None, f"{entry_basename(path)}: {e}"


def _compute_ehull_for_system(args):
    """
    Worker: compute Ehull for all structures in one chemical system.

    Each worker creates its own MP data loader (not shared across processes).
    Returns list of (cif_path, ehull_value) tuples and failure count.
    """
    import re
    import time
    from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

    (chemsys, entries_in_system, mp_data_path, mp_api_key, use_compat) = args

    results = []
    failed = 0

    mp_loader = None
    mp_rester = None
    if mp_data_path:
        from gewum.src.common.ehull.mp_offline_loader import MPOfflineLoader
        mp_loader = MPOfflineLoader(mp_data_path)
    elif mp_api_key:
        from mp_api.client import MPRester
        mp_rester = MPRester(mp_api_key)

    compatibility = None
    if use_compat:
        from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
        from pymatgen.entries.computed_entries import ComputedEntry
        from pymatgen.core.composition import Composition
        from gewum.src.common.ehull.Ehull_compatibility import get_mp2020_potcar_symbols
        compatibility = MaterialsProject2020Compatibility()

    try:
        if mp_loader is not None:
            mp_entries = mp_loader.get_entries_in_chemsys(
                chemsys.split('-'), as_pd_entry=(not use_compat)
            )
        else:
            mp_entries = mp_rester.get_entries_in_chemsys(chemsys.split('-'))
            time.sleep(1.5)
    except Exception as e:
        print(f"[viz2] Warning: failed to get MP entries for {chemsys}: {e}")
        return [], len(entries_in_system)

    for cif_path, formula, total_energy in entries_in_system:
        try:
            if use_compat:
                comp = Composition(formula)
                potcar_symbols = [get_mp2020_potcar_symbols(el.symbol) for el in comp.elements]
                current_entry = ComputedEntry(
                    composition=comp,
                    energy=total_energy,
                    entry_id=f"VIZ2_{formula}",
                    parameters={
                        "run_type": "GGA",
                        "is_hubbard": False,
                        "potcar_symbols": potcar_symbols,
                        "hubbards": {},
                        "potcar_spec": [{"titel": f"PAW_PBE {sym} 06Sep2000"} for sym in potcar_symbols]
                    }
                )
                try:
                    processed = compatibility.process_entries([current_entry])
                    if processed:
                        current_entry = processed[0]
                    else:
                        raise ValueError("Compatibility processing returned empty")
                except Exception as e:
                    print(f"[viz2] Warning: Compatibility error for {formula}: {e}")
                    failed += 1
                    continue
                all_entries = list(mp_entries) + [current_entry]
                pd_diag = PhaseDiagram(all_entries)
                ehull = pd_diag.get_e_above_hull(current_entry)
            else:
                user_entry = PDEntry(composition=formula, energy=total_energy)
                all_entries = list(mp_entries) + [user_entry]
                pd_diag = PhaseDiagram(all_entries)
                ehull = pd_diag.get_e_above_hull(user_entry)

            results.append((cif_path, ehull))
        except Exception as e:
            print(f"[viz2] Warning: Ehull computation failed for {entry_basename(cif_path)}: {e}")
            failed += 1

    return results, failed

# ---------------------------------------------------------------------------
# Module-level helpers for _collect_energy_data (must be picklable for Pool)
# ---------------------------------------------------------------------------

ENERGY_FILES_SG = ['energy_results.csv']
ENERGY_FILES_FORMULA = ['energy_final.txt', '0_final_results.txt']
ENERGY_FILES_CWD = ['energy_final.txt', '0_final_result_tot.txt', '0_final_results.txt']


def _find_energy_file_for_cif(cif_path):
    """Find the energy data file for a given CIF path (module-level, picklable)."""
    vpath = entry_to_path(cif_path)
    parts = Path(vpath).parts
    sg_dir = None
    formula_dir = None

    if 'relaxed' in parts:
        idx = parts.index('relaxed')
        sg_dir = Path(*parts[:idx]) if idx > 0 else None
        formula_dir = Path(*parts[:idx - 1]) if idx > 1 else None
    else:
        sg_dir = Path(vpath).parent.parent
        formula_dir = Path(vpath).parent.parent.parent

    if sg_dir:
        for fname in ENERGY_FILES_SG:
            candidate = sg_dir / fname
            if candidate.exists():
                return str(candidate)
    if formula_dir:
        for fname in ENERGY_FILES_FORMULA:
            candidate = formula_dir / fname
            if candidate.exists():
                return str(candidate)
    for fname in ENERGY_FILES_CWD:
        candidate = Path.cwd() / fname
        if candidate.exists():
            return str(candidate)
    return None


def _collect_energy_batch_worker(args):
    """Worker: read one energy CSV and match a batch of CIFs (picklable)."""
    import csv
    csv_path_str, cifs_info = args
    # cifs_info = [(cif_path, cif_base, cif_resolved), ...]

    rows = []
    fmt = '4col'
    try:
        with open(csv_path_str, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if header is not None:
                header_clean = [h.strip() for h in header]
                if len(header_clean) >= 6 or (header_clean and header_clean[0] == 'Chemical_Formula'):
                    fmt = '6col'
                rows = [row for row in reader if len(row) >= 3]
    except Exception as e:
        print(f"[viz2] Warning: failed to read {csv_path_str}: {e}")

    if fmt == '6col':
        idx_formula = 0
        idx_cif_base = 1
        idx_total_energy = 2
        idx_energy_per_atom = 3
        idx_relaxed_path = 4
    else:
        idx_formula = None
        idx_cif_base = 0
        idx_total_energy = 1
        idx_energy_per_atom = 2
        idx_relaxed_path = 3

    results = []
    for cif_path, cif_base, cif_resolved in cifs_info:
        if not rows:
            results.append((cif_path, None))
            continue

        matched_row = None

        if idx_relaxed_path is not None:
            for row in rows:
                if len(row) > idx_relaxed_path:
                    try:
                        row_path_resolved = str(Path(row[idx_relaxed_path].strip()).resolve())
                        if row_path_resolved == cif_resolved:
                            matched_row = row
                            break
                    except Exception:
                        pass

        if matched_row is None:
            for row in rows:
                if len(row) > idx_cif_base and row[idx_cif_base].strip() == cif_base:
                    matched_row = row
                    break

        if matched_row is None:
            for row in rows:
                if len(row) > idx_cif_base:
                    row_base = row[idx_cif_base].strip()
                    if row_base and (row_base in cif_base or cif_base in row_base):
                        matched_row = row
                        break

        if matched_row is None:
            results.append((cif_path, None))
            continue

        try:
            total_energy = float(matched_row[idx_total_energy].strip())
            energy_per_atom = float(matched_row[idx_energy_per_atom].strip())
        except (ValueError, IndexError):
            results.append((cif_path, None))
            continue

        formula = None
        if idx_formula is not None and len(matched_row) > idx_formula:
            formula = matched_row[idx_formula].strip()

        if not formula:
            try:
                structure = load_structure(cif_path)
                formula = structure.composition.formula
            except Exception:
                results.append((cif_path, None))
                continue

        results.append((cif_path, {
            'formula': formula,
            'total_energy': total_energy,
            'energy_per_atom': energy_per_atom,
        }))

    return results


def _collect_energy_data(cif_paths_relaxed, workers=1):
    """
    Collect energy data for relaxed CIF files.

    Searches for energy data files in the following order:
    1. energy_results.csv in the SG directory (4-column format)
    2. energy_final.txt in the formula directory (6-column format)
    3. 0_final_results.txt in the formula directory (6-column format, selected only)
    4. CWD fallback: energy_final.txt, 0_final_result_tot.txt, 0_final_results.txt

    Args:
        cif_paths_relaxed: list of relaxed CIF file paths.
        workers: Number of parallel workers for CSV matching (1 = serial).

    Returns:
        dict: {cif_path: {'formula': str, 'total_energy': float, 'energy_per_atom': float}}
    """
    from tqdm import tqdm

    if cif_paths_relaxed:
        print(f"[viz2] First CIF path sample: {entry_to_path(cif_paths_relaxed[0])}")
        p0 = Path(entry_to_path(cif_paths_relaxed[0]))
        print(f"[viz2] Path parts: {p0.parts}")

    # ------------------------------------------------------------------
    # Group CIFs by energy file (one unique CSV -> one worker batch)
    # ------------------------------------------------------------------
    csv_groups = defaultdict(list)  # {csv_path: [(cif_path, cif_base, cif_resolved), ...]}
    no_file_count = 0

    for cif_path in cif_paths_relaxed:
        csv_path_str = _find_energy_file_for_cif(cif_path)
        if csv_path_str is None:
            no_file_count += 1
            continue
        vpath = entry_to_path(cif_path)
        p = Path(vpath)
        cif_resolved = str(p.resolve())
        cif_base = p.stem
        csv_groups[csv_path_str].append((cif_path, cif_base, cif_resolved))

    n_csv_files = len(csv_groups)
    print(f"[viz2] Energy files found: {n_csv_files} unique CSVs, "
          f"{no_file_count} CIFs without energy file")

    if n_csv_files == 0:
        print(f"[viz2] Energy data: 0 matched, {len(cif_paths_relaxed)} skipped (no energy files found)")
        return {}

    # ------------------------------------------------------------------
    # Match CIFs to energy rows (parallel when workers > 1)
    # ------------------------------------------------------------------
    energy_data = {}
    matched = 0
    skipped = no_file_count

    if workers > 1 and n_csv_files > 1:
        tasks = list(csv_groups.items())
        with ProcessPoolExecutor(max_workers=workers) as pool:
            batch_results = list(tqdm(
                pool.map(_collect_energy_batch_worker, tasks,
                         chunksize=max(1, len(tasks) // workers // 4)),
                total=len(tasks), desc="[viz2] Collecting energy data"))
        for batch in batch_results:
            for cif_path, info in batch:
                if info is None:
                    skipped += 1
                else:
                    energy_data[cif_path] = info
                    matched += 1
    else:
        # Serial fallback: process each CSV file in a single thread
        for csv_path_str, cifs_info in tqdm(csv_groups.items(),
                                            desc="[viz2] Collecting energy data"):
            batch = _collect_energy_batch_worker((csv_path_str, cifs_info))
            for cif_path, info in batch:
                if info is None:
                    skipped += 1
                else:
                    energy_data[cif_path] = info
                    matched += 1

    print(f"[viz2] Energy data: {matched} matched, {skipped} skipped out of {len(cif_paths_relaxed)} CIFs")
    if matched == 0:
        print("[viz2] Hint: No energy data matched. Possible causes:")
        print("[viz2]   1. Energy files (energy_final.txt / 0_final_result_tot.txt) not found")
        print("[viz2]   2. CIF file names don't match entries in the energy file")
        print("[viz2]   3. Energy file is at an unexpected directory level")
    return energy_data


def _compute_ehull_batch(energy_data, mp_api_key=None, mp_data_path=None, use_compat=False, workers=1):
    """
    Compute energy above hull (Ehull) for a batch of structures.

    Reuses GEWUM's existing MP data infrastructure.

    Args:
        energy_data: dict from _collect_energy_data(),
                     {cif_path: {'formula': str, 'total_energy': float, ...}}
        mp_api_key: Materials Project API key (online mode).
        mp_data_path: Path to offline MP JSON file.
        use_compat: If True, apply MP2020 compatibility corrections
                    (ComputedEntry mode). If False, use raw PDEntry mode.
        workers: Number of parallel workers (1 = serial).

    Returns:
        dict: {cif_path: ehull_value} (eV/atom).
              Entries where Ehull computation failed are excluded.
    """
    import re
    from tqdm import tqdm

    if not mp_api_key and not mp_data_path:
        raise ValueError("[viz2] Either mp_api_key or mp_data_path must be provided for Ehull calculation.")

    if use_compat:
        print("[viz2] Ehull mode: MP2020 compatibility corrections enabled")
    else:
        print("[viz2] Ehull mode: raw PDEntry (no compatibility corrections)")

    system_groups = defaultdict(list)  
    for cif_path, info in energy_data.items():
        formula = info['formula']
        elements = sorted(set(re.findall(r'[A-Z][a-z]*', formula)))
        chemsys = '-'.join(elements)
        system_groups[chemsys].append((cif_path, formula, info['total_energy']))

    if mp_data_path:
        print(f"[viz2] Using offline MP data from: {mp_data_path}")
    else:
        print("[viz2] Using online Materials Project API")

    ehull_results = {}
    failed = 0

    if workers > 1 and len(system_groups) > 1:
        tasks = [
            (chemsys, entries, mp_data_path, mp_api_key, use_compat)
            for chemsys, entries in system_groups.items()
        ]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(_compute_ehull_for_system, t): t[0]
                for t in tasks
            }
            from concurrent.futures import as_completed
            for future in tqdm(as_completed(futures), total=len(futures),
                               desc="[viz2] Computing Ehull by chemical system"):
                chemsys_name = futures[future]
                try:
                    sys_results, sys_failed = future.result(timeout=600)
                    for cif_path, ehull in sys_results:
                        ehull_results[cif_path] = ehull
                    failed += sys_failed
                except TimeoutError:
                    print(f"[viz2] WARNING: Ehull worker timed out for {chemsys_name}")
                    failed += 1
                except Exception as e:
                    print(f"[viz2] WARNING: Ehull worker failed for {chemsys_name}: {e}")
                    failed += 1
    else:
        import time
        from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram

        compatibility = None
        if use_compat:
            from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
            from pymatgen.entries.computed_entries import ComputedEntry
            from pymatgen.core.composition import Composition
            from gewum.src.common.ehull.Ehull_compatibility import POTCAR_MAP, get_mp2020_potcar_symbols
            compatibility = MaterialsProject2020Compatibility()

        mp_loader = None
        mp_rester = None
        if mp_data_path:
            from gewum.src.common.ehull.mp_offline_loader import MPOfflineLoader
            mp_loader = MPOfflineLoader(mp_data_path)
        else:
            from mp_api.client import MPRester
            mp_rester = MPRester(mp_api_key)

        for chemsys in tqdm(system_groups, desc="[viz2] Computing Ehull by chemical system"):
            entries_in_system = system_groups[chemsys]

            try:
                if mp_loader is not None:
                    mp_entries = mp_loader.get_entries_in_chemsys(
                        chemsys.split('-'), as_pd_entry=(not use_compat)
                    )
                else:
                    mp_entries = mp_rester.get_entries_in_chemsys(chemsys.split('-'))
                    time.sleep(1.5)
            except Exception as e:
                print(f"[viz2] Warning: failed to get MP entries for {chemsys}: {e}")
                failed += len(entries_in_system)
                continue

            for cif_path, formula, total_energy in entries_in_system:
                try:
                    if use_compat:
                        comp = Composition(formula)
                        potcar_symbols = [get_mp2020_potcar_symbols(el.symbol) for el in comp.elements]
                        current_entry = ComputedEntry(
                            composition=comp,
                            energy=total_energy,
                            entry_id=f"VIZ2_{formula}",
                            parameters={
                                "run_type": "GGA",
                                "is_hubbard": False,
                                "potcar_symbols": potcar_symbols,
                                "hubbards": {},
                                "potcar_spec": [{"titel": f"PAW_PBE {sym} 06Sep2000"} for sym in potcar_symbols]
                            }
                        )
                        try:
                            processed = compatibility.process_entries([current_entry])
                            if processed:
                                current_entry = processed[0]
                            else:
                                raise ValueError("Compatibility processing returned empty")
                        except Exception as e:
                            print(f"[viz2] Warning: Compatibility error for {formula}: {e}")
                            failed += 1
                            continue

                        all_entries = list(mp_entries) + [current_entry]
                        pd_diag = PhaseDiagram(all_entries)
                        ehull = pd_diag.get_e_above_hull(current_entry)
                    else:
                        user_entry = PDEntry(composition=formula, energy=total_energy)
                        all_entries = list(mp_entries) + [user_entry]
                        pd_diag = PhaseDiagram(all_entries)
                        ehull = pd_diag.get_e_above_hull(user_entry)

                    ehull_results[cif_path] = ehull
                except Exception as e:
                    print(f"[viz2] Warning: Ehull computation failed for {entry_basename(cif_path)}: {e}")
                    failed += 1
                    continue

    print(f"[viz2] Ehull computed: {len(ehull_results)} succeeded, {failed} failed")
    return ehull_results

def _fmt_label(node, count):
    """Format node label: 'SG 4' + 192 -> '4:192'; non-SG labels kept as-is."""
    name = node[3:] if node.startswith('SG ') else node
    return f"{name}:{count}"


def _cubic_bezier_y(x, y0, y1):
    """Compute cubic Bezier interpolation between y0 and y1 over x in [0,1]."""
    t = x
    return y0 * (1 - t) ** 3 + y0 * 3 * t * (1 - t) ** 2 + y1 * 3 * t ** 2 * (1 - t) + y1 * t ** 3


def _merge_others(column_totals, column_nodes, min_fraction=0.01, min_keep=20):
    """Merge tiny nodes into 'Others' when node count > 30 and fraction < min_fraction.

    The top *min_keep* nodes (by flow) are always retained regardless of their
    fraction.  Only nodes ranked after *min_keep* **and** below *min_fraction*
    are merged into the Others group.

    Args:
        column_totals: dict mapping node -> count.
        column_nodes: list of node names (order preserved for keep list).
        min_fraction: threshold below which a node is considered tiny.
        min_keep: minimum number of independent nodes to keep (default 20).

    Returns:
        (keep_nodes, tiny_set): list of nodes to keep, set of nodes merged into Others.
    """
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


def plot_sankey(cif_paths_relaxed, output_path, symprec=0.1, top_n=15, dpi=300, ignore_p1=False, workers=1, cmap_left=None, cmap_right=None):
    """
    Draw a space-group alluvial / flow diagram for relaxed structures.

    For each relaxed CIF the initial SG is read from the path and the
    relaxed SG is identified with spglib.  The flow between (initial, relaxed)
    pairs is visualised as a Sankey-like diagram.

    Args:
        cif_paths_relaxed: list of relaxed CIF file paths.
        output_path: Output PNG path.
        symprec: Symmetry precision.
        top_n: Show only the top-N flows.
        dpi: Figure DPI.
        ignore_p1: Exclude flows where relaxed SG is P1.
        workers: Number of parallel workers (1 = serial).
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from tqdm import tqdm

    flows = []
    if workers > 1 and len(cif_paths_relaxed) > 1:
        tasks = [(p, symprec) for p in cif_paths_relaxed]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(tqdm(
                pool.map(_analyse_sg_worker, tasks,
                         chunksize=max(1, len(tasks) // workers // 4)),
                total=len(tasks), desc="[viz2] Analysing space groups"))
        for cif_path, init_sg, rel_sym, rel_num in results:
            flows.append((f"SG {init_sg}", f"SG {rel_num}"))
    else:
        for cif_path in tqdm(cif_paths_relaxed, desc="[viz2] Analysing space groups"):
            init_sg = _extract_sg_from_path(cif_path)
            rel_sym, rel_num = _identify_sg_from_structure(cif_path, symprec=symprec)
            init_label = f"SG {init_sg}"
            rel_label = f"SG {rel_num}"
            flows.append((init_label, rel_label))

    if ignore_p1:
        flows = [(l, r) for l, r in flows if r != "SG 1"]
        if not flows:
            print("[viz2] Warning: all flows lead to P1, nothing to display with --ignore-p1")
            return

    if not flows:
        print("[viz2] Warning: no valid flows for Sankey diagram")
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
    total_max = max(total_left, total_right, 1)

    gap = 0.01
    bar_w = 0.08

    def _assign_positions(nodes, totals, total_max, gap):
        positions = {}
        y = 0
        for node in nodes:
            h = totals[node] / total_max * (1 - gap * (len(nodes) - 1))
            positions[node] = (y, y + h)
            y += h + gap
        return positions

    left_pos = _assign_positions(left_nodes, left_totals, total_max, gap)
    right_pos = _assign_positions(right_nodes, right_totals, total_max, gap)

    _others_color = (0.831, 0.929, 0.855, 1.0)  # #D4EDDA
    if cmap_left:
        _left_cmap = plt.get_cmap(cmap_left)
        _left_others_color = _left_cmap(1.0)
        left_colors = {}
        for i, n in enumerate(left_nodes):
            if n == others_label:
                left_colors[n] = _left_others_color
            else:
                left_colors[n] = _left_cmap(i / max(len(left_nodes) - 1, 1))
    else:
        _pastel2_l = cm.get_cmap('Pastel2')
        left_colors = {n: (_others_color if n == others_label
                           else _pastel2_l(i % _pastel2_l.N / _pastel2_l.N))
                       for i, n in enumerate(left_nodes)}
    if cmap_right:
        _right_cmap = plt.get_cmap(cmap_right)
        _right_others_color = _right_cmap(1.0)
        right_colors = {}
        for i, n in enumerate(right_nodes):
            if n == others_label:
                right_colors[n] = _right_others_color
            else:
                right_colors[n] = _right_cmap(i / max(len(right_nodes) - 1, 1))
    else:
        _pastel2_r = cm.get_cmap('Pastel2')
        right_colors = {n: (_others_color if n == others_label
                            else _pastel2_r(i % _pastel2_r.N / _pastel2_r.N))
                        for i, n in enumerate(right_nodes)}

    n_left = len(left_nodes)
    n_right = len(right_nodes)
    max_nodes = max(n_left, n_right)

    if max_nodes <= 10:
        fontsize = 10
    elif max_nodes <= 25:
        fontsize = 8
    elif max_nodes <= 50:
        fontsize = 6
    else:
        fontsize = 5

    fig_height = max(8, max_nodes * 0.35)
    fig_width = 14
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    for node in left_nodes:
        y0, y1 = left_pos[node]
        ax.barh(y=(y0 + y1) / 2, width=bar_w, height=y1 - y0,
                left=0, color=left_colors[node], edgecolor='white', linewidth=0.5)
        label_text = _fmt_label(node, left_totals[node])
        ax.text(-0.02, (y0 + y1) / 2, label_text,
                ha='right', va='center', fontsize=fontsize, fontweight='bold')

    for node in right_nodes:
        y0, y1 = right_pos[node]
        ax.barh(y=(y0 + y1) / 2, width=bar_w, height=y1 - y0,
                left=1 - bar_w, color=right_colors[node], edgecolor='white',
                linewidth=0.5, alpha=0.9)
        label_text = _fmt_label(node, right_totals[node])
        ax.text(1 + 0.02, (y0 + y1) / 2, label_text,
                ha='left', va='center', fontsize=fontsize, fontweight='bold')

    left_cursor = {n: left_pos[n][0] for n in left_nodes}
    right_cursor = {n: right_pos[n][0] for n in right_nodes}

    x_vals = np.linspace(0, 1, 200)

    for (l_node, r_node), cnt in top_flows:
        frac = cnt / total_max * (1 - gap * (max(len(left_nodes), len(right_nodes)) - 1))

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

        x_draw = bar_w + t * (1 - 2 * bar_w)

        ax.fill_between(x_draw, bot_curve, top_curve,
                        alpha=0.35, color=left_colors[l_node], edgecolor='none')

    ax.set_xlim(-0.35, 1.35)
    y_max = max(
        max((p[1] for p in left_pos.values()), default=1),
        max((p[1] for p in right_pos.values()), default=1),
    )
    ax.set_ylim(-0.05, y_max + 0.05)
    ax.invert_yaxis()
    ax.axis('off')

    fig.set_facecolor('white')
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f"[viz2] Saved: {output_path}")


# ---------------------------------------------------------------------------
# 3c. Three-column Sankey with Ehull
# ---------------------------------------------------------------------------

def plot_sankey_ehull(cif_paths_relaxed, output_path,
                      mp_api_key=None, mp_data_path=None,
                      symprec=0.1, top_n=15, dpi=300,
                      ignore_p1=False, use_compat=False, workers=1,
                      cmap=None, energy_data_cache=None):
    """
    Draw a three-column Sankey diagram: Initial SG -> Relaxed SG -> Ehull bin.

    Args:
        cif_paths_relaxed: list of relaxed CIF file paths.
        output_path: Output PNG path.
        mp_api_key: Materials Project API key (online mode).
        mp_data_path: Path to offline MP JSON file.
        symprec: Symmetry precision for spglib.
        top_n: Top-N flows to show per segment.
        dpi: Figure DPI.
        ignore_p1: Exclude flows where relaxed SG is P1.
        use_compat: If True, apply MP2020 compatibility corrections.
        workers: Number of parallel workers (1 = serial).
        cmap: Matplotlib colormap name for left/right columns.
              Middle column uses the reversed colormap.
              Default: None (blue->pink 'A8D8EA->F3C1D3' scheme).
        energy_data_cache: Pre-collected energy data dict (from _collect_energy_data).
              If None, energy data will be collected inside this function.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from tqdm import tqdm

    # ------------------------------------------------------------------
    # Step 1: Collect energy + Ehull data
    # ------------------------------------------------------------------
    try:
        if energy_data_cache is not None:
            energy_data = energy_data_cache
        else:
            energy_data = _collect_energy_data(cif_paths_relaxed, workers=workers)
        if not energy_data:
            print("[viz2] WARNING: No energy data collected, skipping sankey_ehull")
            return

        ehull_results = _compute_ehull_batch(energy_data, mp_api_key, mp_data_path,
                                             use_compat=use_compat, workers=workers)
        if not ehull_results:
            print("[viz2] WARNING: No Ehull results computed, skipping sankey_ehull")
            return
    except Exception as e:
        print(f"[viz2] ERROR in Ehull calculation: {e}")
        return

    # ------------------------------------------------------------------
    # Step 2: Build per-structure records (init_sg, relaxed_sg, ehull_bin)
    # ------------------------------------------------------------------
    records = []  

    sg_map = {} 
    if workers > 1 and len(cif_paths_relaxed) > 1:
        tasks = [(p, symprec) for p in cif_paths_relaxed]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(tqdm(
                pool.map(_analyse_sg_worker, tasks,
                         chunksize=max(1, len(tasks) // workers // 4)),
                total=len(tasks), desc="[viz2] Analysing (ehull sankey)"))
        for cif_path, init_sg, rel_sym, rel_num in results:
            sg_map[cif_path] = (init_sg, rel_num)
    else:
        for cif_path in tqdm(cif_paths_relaxed, desc="[viz2] Analysing (ehull sankey)"):
            init_sg = _extract_sg_from_path(cif_path)
            rel_sym, rel_num = _identify_sg_from_structure(cif_path, symprec=symprec)
            sg_map[cif_path] = (init_sg, rel_num)

    for cif_path in cif_paths_relaxed:
        if cif_path not in sg_map or cif_path not in ehull_results:
            continue
        init_sg, rel_num = sg_map[cif_path]
        init_label = f"SG {init_sg}"
        rel_label = f"SG {rel_num}"
        ehull_val = ehull_results[cif_path]

        if ehull_val >= 0.2:
            ehull_label = ">0.2"
        else:
            bin_idx = int(ehull_val / 0.02)
            bin_idx = min(bin_idx, 9)
            low = bin_idx * 0.02
            high = low + 0.02
            low_s = "0" if low == 0 else f"{low:.2f}"
            high_s = f"{high:.2f}"
            ehull_label = f"{low_s}~{high_s}"

        records.append((init_label, rel_label, ehull_label))

    if ignore_p1:
        records = [r for r in records if r[1] != "SG 1"]

    if not records:
        print("[viz2] Warning: no valid records for Ehull Sankey")
        return

    # ------------------------------------------------------------------
    # Step 3: Count flows for each segment (no top_n truncation)
    # ------------------------------------------------------------------
    seg1_flows = Counter((r[0], r[1]) for r in records)
    seg2_flows = Counter((r[1], r[2]) for r in records)

    left_totals = defaultdict(int)
    mid_totals = defaultdict(int)
    right_totals = defaultdict(int)
    for init_sg, rel_sg, ehull_label in records:
        left_totals[init_sg] += 1
        mid_totals[rel_sg] += 1
        right_totals[ehull_label] += 1

    # ------------------------------------------------------------------
    # Step 4: Others merge (left & mid only; right uses fixed bins)
    # ------------------------------------------------------------------
    others_label = "Others"

    def _ehull_sort_key(label):
        if label.startswith(">"):
            return 9999
        try:
            return float(label.split('~')[0])
        except (ValueError, IndexError):
            return 9999

    ALL_EHULL_BINS = []
    for _idx in range(10):
        _low = _idx * 0.02
        _high = _low + 0.02
        _low_s = "0" if _low == 0 else f"{_low:.2f}"
        _high_s = f"{_high:.2f}"
        ALL_EHULL_BINS.append(f"{_low_s}~{_high_s}")
    ALL_EHULL_BINS.append(">0.2")

    left_nodes_sorted = sorted(left_totals.keys(), key=lambda k: -left_totals[k])
    left_keep, left_tiny = _merge_others(left_totals, left_nodes_sorted)

    mid_nodes_sorted = sorted(mid_totals.keys(), key=lambda k: -mid_totals[k])
    mid_keep, mid_tiny = _merge_others(mid_totals, mid_nodes_sorted)

    right_tiny = set()  

    new_seg1 = defaultdict(int)
    new_seg2 = defaultdict(int)
    for (l, m), cnt in seg1_flows.items():
        if l in left_tiny:
            l = others_label
        if m in mid_tiny:
            m = others_label
        new_seg1[(l, m)] += cnt

    for (m, r), cnt in seg2_flows.items():
        if m in mid_tiny:
            m = others_label
        new_seg2[(m, r)] += cnt

    seg1_flows = new_seg1
    seg2_flows = new_seg2

    left_totals = defaultdict(int)
    mid_totals = defaultdict(int)
    right_totals = defaultdict(int)
    for (l, m), cnt in seg1_flows.items():
        left_totals[l] += cnt
    for init_sg, rel_sg, ehull_label in records:
        rel_mapped = others_label if rel_sg in mid_tiny else rel_sg
        mid_totals[rel_mapped] += 1
    for (m, r), cnt in seg2_flows.items():
        right_totals[r] += cnt

    for _bin_label in ALL_EHULL_BINS:
        if _bin_label not in right_totals:
            right_totals[_bin_label] = 0

    left_nodes = sorted([n for n in left_totals if n != others_label],
                        key=lambda k: _sg_sort_key(k))
    if others_label in left_totals:
        left_nodes.append(others_label)

    mid_nodes = sorted([n for n in mid_totals if n != others_label],
                       key=lambda k: _sg_sort_key(k))
    if others_label in mid_totals:
        mid_nodes.append(others_label)

    right_nodes = list(ALL_EHULL_BINS)

    seg1_top = list(seg1_flows.items())
    seg2_top = list(seg2_flows.items())

    # ------------------------------------------------------------------
    # Step 5: Position calculation
    # ------------------------------------------------------------------
    total_max = max(
        sum(left_totals.values()),
        sum(mid_totals.values()),
        sum(right_totals.values()),
        1,
    )
    gap = 0.01
    bar_w = 0.06

    def _assign_pos(nodes, totals, total_max, gap, min_height=0.0):
        positions = {}
        y = 0
        for node in nodes:
            h = totals[node] / total_max * (1 - gap * (len(nodes) - 1))
            if h == 0 and min_height > 0:
                h = min_height
            positions[node] = (y, y + h)
            y += h + gap
        return positions

    _total_height = 1 - gap * (max(len(right_nodes) - 1, 0))
    _min_bin_height = _total_height * 0.005

    left_pos = _assign_pos(left_nodes, left_totals, total_max, gap)
    mid_pos = _assign_pos(mid_nodes, mid_totals, total_max, gap)
    right_pos = _assign_pos(right_nodes, right_totals, total_max, gap,
                            min_height=_min_bin_height)

    # ------------------------------------------------------------------
    # Step 6: Figure setup
    # ------------------------------------------------------------------
    max_nodes = max(len(left_nodes), len(mid_nodes), len(right_nodes))
    if max_nodes <= 10:
        fontsize = 10
    elif max_nodes <= 25:
        fontsize = 8
    elif max_nodes <= 50:
        fontsize = 6
    else:
        fontsize = 5

    fig_height = max(8, max_nodes * 0.35)
    fig_width = 20
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    x_left = 0.0
    x_mid = 0.45
    x_right = 0.90

    # ------------------------------------------------------------------
    # Step 7: Node colours
    # ------------------------------------------------------------------
    from matplotlib.colors import LinearSegmentedColormap

    # ---- Resolve colormaps ----
    if cmap:
        # User-provided colormap: left/right = cmap, mid = cmap_r
        _left_cmap = plt.get_cmap(cmap)
        _mid_cmap = plt.get_cmap(cmap + '_r')
        _ehull_cmap = plt.get_cmap(cmap)
        _others_color_left = _left_cmap(1.0)
        _others_color_mid  = _mid_cmap(1.0)
        _others_color_right = _left_cmap(1.0)
    else:
        # Default: blue -> pink gradient
        _others_color_left = (0.97, 0.85, 0.89, 1.0)
        _others_color_mid  = (0.85, 0.93, 0.97, 1.0)
        _others_color_right = (0.9, 0.9, 0.9, 1.0)

        _left_cmap = LinearSegmentedColormap.from_list(
            'lr_blue_pink', ['#A8D8EA', '#F3C1D3'],
        )
        _mid_cmap = LinearSegmentedColormap.from_list(
            'mid_pink_blue', ['#F3C1D3', '#A8D8EA'],
        )
        _ehull_cmap = LinearSegmentedColormap.from_list(
            'ehull_blue_pink', ['#A8D8EA', '#F3C1D3'],
        )

    non_others_left = [n for n in left_nodes if n != others_label]
    n_left = len(non_others_left)
    left_colors = {}
    for i, n in enumerate(non_others_left):
        left_colors[n] = _left_cmap(i / max(n_left - 1, 1))
    if others_label in left_nodes:
        left_colors[others_label] = _others_color_left

    non_others_mid = [n for n in mid_nodes if n != others_label]
    n_mid = len(non_others_mid)
    mid_colors = {}
    for i, n in enumerate(non_others_mid):
        mid_colors[n] = _mid_cmap(i / max(n_mid - 1, 1))
    if others_label in mid_nodes:
        mid_colors[others_label] = _others_color_mid

    _empty_bin_color = (0.9, 0.9, 0.9, 1.0)
    right_colors = {}
    n_all_bins = len(ALL_EHULL_BINS)
    for bin_idx, bin_label in enumerate(ALL_EHULL_BINS):
        if right_totals.get(bin_label, 0) == 0:
            right_colors[bin_label] = _empty_bin_color
        else:
            right_colors[bin_label] = _ehull_cmap(bin_idx / max(n_all_bins - 1, 1))
    if others_label in right_nodes:
        right_colors[others_label] = _others_color_right

    for node in left_nodes:
        y0, y1 = left_pos[node]
        ax.barh(y=(y0 + y1) / 2, width=bar_w, height=y1 - y0,
                left=x_left, color=left_colors[node],
                edgecolor='white', linewidth=0.5)
        ax.text(x_left - 0.01, (y0 + y1) / 2,
                _fmt_label(node, left_totals[node]),
                ha='right', va='center', fontsize=fontsize, fontweight='bold')

    for node in mid_nodes:
        y0, y1 = mid_pos[node]
        ax.barh(y=(y0 + y1) / 2, width=bar_w, height=y1 - y0,
                left=x_mid, color=mid_colors[node],
                edgecolor='white', linewidth=0.5)
        ax.text(x_mid + bar_w + 0.01, (y0 + y1) / 2,
                _fmt_label(node, mid_totals[node]),
                ha='left', va='center', fontsize=fontsize, fontweight='bold')

    for node in right_nodes:
        y0, y1 = right_pos[node]
        ax.barh(y=(y0 + y1) / 2, width=bar_w, height=y1 - y0,
                left=x_right, color=right_colors[node],
                edgecolor='white', linewidth=0.5)
        ax.text(x_right + bar_w + 0.01, (y0 + y1) / 2,
                _fmt_label(node, right_totals[node]),
                ha='left', va='center', fontsize=fontsize, fontweight='bold')

    left_cursor = {n: left_pos[n][0] for n in left_nodes}
    mid_cursor_left = {n: mid_pos[n][0] for n in mid_nodes}

    x_vals = np.linspace(0, 1, 200)

    for (l_node, m_node), cnt in seg1_top:
        if l_node not in left_pos or m_node not in mid_pos:
            continue
        frac = cnt / total_max * (1 - gap * (max(len(left_nodes), len(mid_nodes)) - 1))

        l_y0 = left_cursor[l_node]
        l_y1 = l_y0 + frac
        left_cursor[l_node] = l_y1

        m_y0 = mid_cursor_left[m_node]
        m_y1 = m_y0 + frac
        mid_cursor_left[m_node] = m_y1

        t = x_vals
        s = 3 * t ** 2 - 2 * t ** 3
        top_curve = l_y1 + s * (m_y1 - l_y1)
        bot_curve = l_y0 + s * (m_y0 - l_y0)

        x_draw = (x_left + bar_w) + t * (x_mid - x_left - bar_w)
        ax.fill_between(x_draw, bot_curve, top_curve,
                        alpha=0.35, color=left_colors[l_node], edgecolor='none')

    mid_cursor_right = {n: mid_pos[n][0] for n in mid_nodes}
    right_cursor = {n: right_pos[n][0] for n in right_nodes}

    for (m_node, r_node), cnt in seg2_top:
        if m_node not in mid_pos or r_node not in right_pos:
            continue
        frac = cnt / total_max * (1 - gap * (max(len(mid_nodes), len(right_nodes)) - 1))

        m_y0 = mid_cursor_right[m_node]
        m_y1 = m_y0 + frac
        mid_cursor_right[m_node] = m_y1

        r_y0 = right_cursor[r_node]
        r_y1 = r_y0 + frac
        right_cursor[r_node] = r_y1

        t = x_vals
        s = 3 * t ** 2 - 2 * t ** 3
        top_curve = m_y1 + s * (r_y1 - m_y1)
        bot_curve = m_y0 + s * (r_y0 - m_y0)

        x_draw = (x_mid + bar_w) + t * (x_right - x_mid - bar_w)
        ax.fill_between(x_draw, bot_curve, top_curve,
                        alpha=0.35, color=right_colors[r_node], edgecolor='none')

    ax.set_xlim(-0.25, x_right + bar_w + 0.25)
    y_max = max(
        max((p[1] for p in left_pos.values()), default=1),
        max((p[1] for p in mid_pos.values()), default=1),
        max((p[1] for p in right_pos.values()), default=1),
    )
    ax.set_ylim(-0.05, y_max + 0.05)
    ax.invert_yaxis()
    ax.axis('off')

    ax.text(x_left + bar_w / 2, -0.03, "Initial SG",
            ha='center', va='bottom', fontsize=12, fontweight='bold')
    ax.text(x_mid + bar_w / 2, -0.03, "Relaxed SG",
            ha='center', va='bottom', fontsize=12, fontweight='bold')
    ax.text(x_right + bar_w / 2, -0.03, "E$_{\mathrm{hull}}$ (eV/atom)",
            ha='center', va='bottom', fontsize=12, fontweight='bold')

    fig.set_facecolor('white')
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f"[viz2] Saved: {output_path}")


# ---------------------------------------------------------------------------
# 3b. Chord diagram
# ---------------------------------------------------------------------------

def _sg_sort_key(label):
    """Extract numeric space group number from 'SG xxx' label for sorting."""
    try:
        return int(label.replace("SG ", ""))
    except (ValueError, AttributeError):
        return 9999  

def plot_chord(cif_paths_relaxed, output_path, symprec=0.1, top_n=15, dpi=300, ignore_p1=False, workers=1):
    """
    Draw a chord diagram showing space group flow between initial and relaxed states.

    Nodes (space groups) are arranged on a circle. Ribbons connect nodes
    that have flow between them, with width proportional to count.

    Args:
        cif_paths_relaxed: list of relaxed CIF file paths.
        output_path: Output PNG path.
        symprec: Symmetry precision.
        top_n: Show only the top-N flows.
        dpi: Figure DPI.
        ignore_p1: Exclude flows where relaxed SG is P1.
        workers: Number of parallel workers (1 = serial).
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.path as mpath
    import matplotlib.cm as cm
    from tqdm import tqdm

    flows = []
    if workers > 1 and len(cif_paths_relaxed) > 1:
        tasks = [(p, symprec) for p in cif_paths_relaxed]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(tqdm(
                pool.map(_analyse_sg_worker, tasks,
                         chunksize=max(1, len(tasks) // workers // 4)),
                total=len(tasks), desc="[viz2] Analysing space groups (chord)"))
        for cif_path, init_sg, rel_sym, rel_num in results:
            flows.append((f"SG {init_sg}", f"SG {rel_num}"))
    else:
        for cif_path in tqdm(cif_paths_relaxed, desc="[viz2] Analysing space groups (chord)"):
            init_sg = _extract_sg_from_path(cif_path)
            rel_sym, rel_num = _identify_sg_from_structure(cif_path, symprec=symprec)
            init_label = f"SG {init_sg}"
            rel_label = f"SG {rel_num}"
            flows.append((init_label, rel_label))

    if ignore_p1:
        flows = [(l, r) for l, r in flows if r != "SG 1"]
        if not flows:
            print("[viz2] Warning: all flows lead to P1, nothing to display with --ignore-p1")
            return

    if not flows:
        print("[viz2] Warning: no valid flows for chord diagram")
        return

    flow_counts = Counter(flows)
    top_flows = flow_counts.most_common(top_n)

    node_flow = defaultdict(int)  
    for (l, r), cnt in top_flows:
        node_flow[l] += cnt
        node_flow[r] += cnt

    nodes = sorted(node_flow.keys(), key=lambda k: _sg_sort_key(k))
    n_nodes = len(nodes)

    if n_nodes == 0:
        print("[viz2] Warning: no nodes for chord diagram")
        return

    node_colors = {n: cm.tab20(i / max(n_nodes, 1)) for i, n in enumerate(nodes)}

    total_flow = sum(node_flow[n] for n in nodes)
    gap_deg = 1.5  
    total_gap = gap_deg * n_nodes
    available_deg = 360.0 - total_gap
    if available_deg < 10:
        gap_deg = 0.5
        total_gap = gap_deg * n_nodes
        available_deg = 360.0 - total_gap

    node_arc_start = {}
    node_arc_end = {}
    angle = 0.0
    for n in nodes:
        span = (node_flow[n] / total_flow) * available_deg
        node_arc_start[n] = np.radians(angle)
        node_arc_end[n] = np.radians(angle + span)
        angle += span + gap_deg

    node_cursor = {n: node_arc_start[n] for n in nodes}

    radius = 1.0

    def _angle_to_xy(theta, r=radius):
        return r * np.cos(theta), r * np.sin(theta)

    def _arc_points(theta0, theta1, r=radius, n_pts=50):
        angles = np.linspace(theta0, theta1, n_pts)
        return np.column_stack([r * np.cos(angles), r * np.sin(angles)])

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal')

    arc_width = 0.08
    for n in nodes:
        theta0 = node_arc_start[n]
        theta1 = node_arc_end[n]
        deg0 = np.degrees(theta0)
        deg1 = np.degrees(theta1)
        wedge = mpatches.Wedge(
            (0, 0), radius, deg0, deg1,
            width=arc_width,
            facecolor=node_colors[n],
            edgecolor='white',
            linewidth=0.5,
        )
        ax.add_patch(wedge)

    inner_r = radius - arc_width  

    for (l_node, r_node), cnt in top_flows:
        if l_node not in node_flow or r_node not in node_flow:
            continue

        l_span = (cnt / total_flow) * available_deg
        r_span = (cnt / total_flow) * available_deg
        l_theta0 = node_cursor[l_node]
        l_theta1 = l_theta0 + np.radians(l_span)
        node_cursor[l_node] = l_theta1

        r_theta0 = node_cursor[r_node]
        r_theta1 = r_theta0 + np.radians(r_span)
        node_cursor[r_node] = r_theta1


        src_arc = _arc_points(l_theta0, l_theta1, r=inner_r, n_pts=30)
        dst_arc = _arc_points(r_theta0, r_theta1, r=inner_r, n_pts=30)

        verts = []
        codes = []

        for i, pt in enumerate(src_arc):
            verts.append(pt)
            codes.append(mpath.Path.MOVETO if i == 0 else mpath.Path.LINETO)

        verts.append((0, 0))
        codes.append(mpath.Path.CURVE3)
        verts.append(dst_arc[0])
        codes.append(mpath.Path.CURVE3)

        for pt in dst_arc:
            verts.append(pt)
            codes.append(mpath.Path.LINETO)

        verts.append((0, 0))
        codes.append(mpath.Path.CURVE3)
        verts.append(src_arc[0])
        codes.append(mpath.Path.CURVE3)

        verts.append(src_arc[0])
        codes.append(mpath.Path.CLOSEPOLY)

        path = mpath.Path(verts, codes)
        patch = mpatches.PathPatch(
            path,
            facecolor=node_colors[l_node],
            edgecolor='none',
            alpha=0.5,
        )
        ax.add_patch(patch)

    label_r = radius + 0.06
    fontsize = 9 if n_nodes <= 15 else 7 if n_nodes <= 25 else 5
    for n in nodes:
        mid_angle = (node_arc_start[n] + node_arc_end[n]) / 2
        x, y = _angle_to_xy(mid_angle, r=label_r)
        deg = np.degrees(mid_angle) % 360

        if 90 < deg < 270:
            rotation = deg + 180
            ha = 'right'
        else:
            rotation = deg
            ha = 'left'

        ax.text(x, y, _fmt_label(n, node_flow[n]),
                ha=ha, va='center',
                fontsize=fontsize, fontweight='bold',
                rotation=rotation, rotation_mode='anchor')

    margin = 0.45
    lim = radius + margin
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axis('off')

    fig.set_facecolor('white')
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f"[viz2] Saved: {output_path}")


# ---------------------------------------------------------------------------
# 4. RDF comparison
# ---------------------------------------------------------------------------

def _compute_rdf(structure, r_max=10.0, bin_size=0.1):
    """
    Compute the radial distribution function g(r) for a single structure.

    Args:
        structure: pymatgen Structure object.
        r_max: Maximum radius in Angstrom.
        bin_size: Histogram bin width in Angstrom.

    Returns:
        (r_values, g_r) - 1-D numpy arrays.
    """
    bins = np.arange(0, r_max + bin_size, bin_size)
    r_values = (bins[:-1] + bins[1:]) / 2.0
    counts = np.zeros(len(r_values))

    all_neighbors = structure.get_all_neighbors(r_max)
    n_atoms = len(structure)

    for site_neighbors in all_neighbors:
        for neighbor in site_neighbors:
            d = neighbor.nn_distance
            idx = int(d / bin_size)
            if 0 <= idx < len(counts):
                counts[idx] += 1

    volume = structure.volume
    rho = n_atoms / volume  

    g_r = np.zeros_like(counts)
    for i, r in enumerate(r_values):
        r_inner = bins[i]
        r_outer = bins[i + 1]
        shell_vol = (4.0 / 3.0) * np.pi * (r_outer ** 3 - r_inner ** 3)
        if shell_vol > 0 and rho > 0 and n_atoms > 0:
            g_r[i] = counts[i] / (n_atoms * rho * shell_vol)

    return r_values, g_r


def plot_rdf_comparison(cif_dict, output_path, r_max=10.0, bin_size=0.1,
                        max_structures=200, dpi=300, workers=1):
    """
    Plot average RDF comparison across workflow stages.

    For each stage the per-structure RDF is computed and averaged;
    +/-1 std-dev is drawn as a shaded band.

    Args:
        cif_dict: ``{'total': [paths], 'selected': [paths], ...}``
        output_path: Output PNG path.
        r_max: Maximum r for RDF.
        bin_size: Bin width.
        max_structures: Cap per stage (random sample if exceeded).
        dpi: Figure DPI.
        workers: Number of parallel workers (1 = serial).
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from tqdm import tqdm

    stage_colors = {
        'total': '#AAAAAA',
        'selected': '#FF6B6B',
        'relaxed': '#4ECDC4',
    }
    stage_order = ['total', 'selected', 'relaxed']

    fig, ax = plt.subplots(figsize=(10, 6))
    has_data = False

    for stage in stage_order:
        paths = cif_dict.get(stage, [])
        if not paths:
            continue

        if len(paths) > max_structures:
            rng = np.random.default_rng(42)
            indices = rng.choice(len(paths), size=max_structures, replace=False)
            paths = [paths[i] for i in indices]

        all_gr = []
        r_vals = None

        if workers > 1 and len(paths) > 1:
            tasks = [(p, r_max, bin_size) for p in paths]
            with ProcessPoolExecutor(max_workers=workers) as pool:
                rdf_results = list(tqdm(
                    pool.map(_compute_rdf_worker, tasks,
                             chunksize=max(1, len(tasks) // workers // 4)),
                    total=len(tasks), desc=f"[viz2] RDF {stage}"))
            for r, gr, err in rdf_results:
                if err is not None:
                    print(f"[viz2] Warning: skipping {err}")
                else:
                    all_gr.append(gr)
                    if r_vals is None:
                        r_vals = r
        else:
            for p in tqdm(paths, desc=f"[viz2] RDF {stage}"):
                try:
                    struct = load_structure(p)
                    r, gr = _compute_rdf(struct, r_max=r_max, bin_size=bin_size)
                    all_gr.append(gr)
                    if r_vals is None:
                        r_vals = r
                except Exception as e:
                    print(f"[viz2] Warning: skipping {entry_basename(p)}: {e}")

        if not all_gr:
            print(f"[viz2] Warning: no valid RDF data for stage '{stage}'")
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
            peak_idx = np.argmax(mean_gr[1:]) + 1  # skip r~0
            ax.annotate(f'{r_vals[peak_idx]:.2f} \u00c5',
                        xy=(r_vals[peak_idx], mean_gr[peak_idx]),
                        fontsize=9, ha='center', va='bottom',
                        color=color, fontweight='bold')
        has_data = True

    if not has_data:
        print("[viz2] Warning: no RDF data to plot")
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
    print(f"[viz2] Saved: {output_path}")


# ---------------------------------------------------------------------------
# 5. Structure funnel chart
# ---------------------------------------------------------------------------

def plot_funnel(cif_dict, output_path, dpi=300):
    """
    Draw a funnel chart showing structure counts across workflow stages.

    Args:
        cif_dict: ``{'total': [paths], 'selected': [paths], ...}``
        output_path: Output PNG path.
        dpi: Figure DPI.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    stage_order = ['total', 'selected', 'relaxed']
    stage_colors = {
        'total': '#AAAAAA',
        'selected': '#FF6B6B',
        'relaxed': '#4ECDC4',
    }
    stage_labels = {
        'total': 'Total (generated)',
        'selected': 'Selected',
        'relaxed': 'Relaxed',
    }

    stages = []
    counts = []
    for s in stage_order:
        if s in cif_dict and cif_dict[s]:
            stages.append(s)
            counts.append(len(cif_dict[s]))

    if not stages:
        print("[viz2] Warning: no data for funnel chart")
        return

    max_count = max(counts)
    total_count = counts[0] if counts else 1

    fig, ax = plt.subplots(figsize=(8, 6))

    n_stages = len(stages)
    bar_height = 0.6
    y_gap = 1.0

    for i, (stage, count) in enumerate(zip(stages, counts)):
        y_center = -i * y_gap
        width = count / max_count 

        color = stage_colors.get(stage, '#CCCCCC')
        rect = mpatches.FancyBboxPatch(
            (-width / 2, y_center - bar_height / 2), width, bar_height,
            boxstyle="round,pad=0.02", facecolor=color, edgecolor='white',
            linewidth=2, alpha=0.9,
        )
        ax.add_patch(rect)

        retention = count / total_count * 100
        label = stage_labels.get(stage, stage.capitalize())
        ax.text(0, y_center, f"{label}\n{count} structures ({retention:.1f}%)",
                ha='center', va='center', fontsize=12, fontweight='bold',
                color='white' if count / max_count > 0.3 else 'black')

        if i < n_stages - 1:
            next_width = counts[i + 1] / max_count
            y_bot = y_center - bar_height / 2
            y_next_top = -(i + 1) * y_gap + bar_height / 2

            trap_x = [
                -width / 2, width / 2,
                next_width / 2, -next_width / 2,
            ]
            trap_y = [
                y_bot, y_bot,
                y_next_top, y_next_top,
            ]
            ax.fill(trap_x, trap_y, color=color, alpha=0.15, edgecolor='none')

    if len(counts) >= 2:
        overall = counts[-1] / total_count * 100
        ax.text(0, -(n_stages - 1) * y_gap - bar_height / 2 - 0.35,
                f"Overall retention: {overall:.1f}%",
                ha='center', va='top', fontsize=12, fontstyle='italic', color='#333333')

    margin = 0.15
    ax.set_xlim(-0.5 - margin, 0.5 + margin)
    ax.set_ylim(-(n_stages - 1) * y_gap - bar_height / 2 - 0.6,
                bar_height / 2 + 0.3)
    ax.axis('off')

    fig.set_facecolor('white')
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f"[viz2] Saved: {output_path}")

def generate_report(cif_dict, output_path, sg_mapping=None):
    """
    Generate a plain-text summary report ``viz2_report.txt``.

    Args:
        cif_dict: ``{'total': [paths], ...}``
        output_path: Path for the report file.
        sg_mapping: Optional dict mapping cif_path -> (init_sg, relaxed_sg).
    """
    lines = [
        "=" * 60,
        "GEWUM viz2 Report",
        "=" * 60,
        f"Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "--- Stage Statistics ---",
    ]

    stage_order = ['total', 'selected', 'relaxed']
    total_n = len(cif_dict.get('total', []))
    if total_n == 0:
        for s in stage_order:
            paths = cif_dict.get(s, [])
            if paths:
                total_n = len(paths)
                break

    for stage in stage_order:
        paths = cif_dict.get(stage, [])
        if not paths:
            continue
        n = len(paths)
        pct = (n / total_n * 100) if total_n > 0 else 0
        lines.append(f"  {stage:>10s}: {n:6d} structures  ({pct:5.1f}%)")

    if sg_mapping:
        lines.append("")
        lines.append("--- Space Group Changes ---")
        kept = sum(1 for (i, r) in sg_mapping.values() if str(i) == str(r))
        changed = len(sg_mapping) - kept
        lines.append(f"  Kept original SG : {kept}")
        lines.append(f"  Changed SG       : {changed}")
        if sg_mapping:
            pct_kept = kept / len(sg_mapping) * 100
            lines.append(f"  Retention rate   : {pct_kept:.1f}%")

    lines.append("")
    lines.append("=" * 60)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f"[viz2] Saved report: {output_path}")


# ---------------------------------------------------------------------------
# 6b. Violin plot: energy distribution by space group
# ---------------------------------------------------------------------------

def plot_violin_energy(cif_paths_relaxed, output_path,
                       symprec=0.1, dpi=300, ignore_p1=False,
                       workers=1, top_n=20, energy_data_cache=None):
    """
    Draw a violin plot showing energy_per_atom distribution grouped by space group.

    Args:
        cif_paths_relaxed: list of relaxed CIF file paths.
        output_path: Output PNG path.
        symprec: Symmetry precision for spglib.
        dpi: Figure DPI.
        ignore_p1: If True, exclude P1 (SG 1) from the plot.
        workers: Number of parallel workers (1 = serial).
        top_n: Only show the top-N most populated space groups (default: 20).
        energy_data_cache: Pre-collected energy data dict (from _collect_energy_data).
              If None, energy data will be collected inside this function.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from tqdm import tqdm

    energy_data = energy_data_cache if energy_data_cache is not None else _collect_energy_data(cif_paths_relaxed, workers=workers)
    if not energy_data:
        print("[viz2] Warning: no energy data found, skipping violin plot")
        return

    sg_map = {}  
    paths_with_energy = [p for p in cif_paths_relaxed if p in energy_data]
    if not paths_with_energy:
        print("[viz2] Warning: no CIF files with energy data, skipping violin plot")
        return

    if workers > 1 and len(paths_with_energy) > 1:
        tasks = [(p, symprec) for p in paths_with_energy]
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(tqdm(
                pool.map(_analyse_sg_worker, tasks,
                         chunksize=max(1, len(tasks) // workers // 4)),
                total=len(tasks), desc="[viz2] Analysing space groups (violin)"))
        for cif_path, _init_sg, rel_sym, rel_num in results:
            sg_map[cif_path] = (rel_sym, rel_num)
    else:
        for cif_path in tqdm(paths_with_energy, desc="[viz2] Analysing space groups (violin)"):
            rel_sym, rel_num = _identify_sg_from_structure(cif_path, symprec=symprec)
            sg_map[cif_path] = (rel_sym, rel_num)

    sg_energies = defaultdict(list)  
    for cif_path in paths_with_energy:
        if cif_path not in sg_map:
            continue
        _rel_sym, rel_num = sg_map[cif_path]
        epa = energy_data[cif_path].get('energy_per_atom')
        if epa is None:
            continue
        sg_energies[rel_num].append(epa)

    if ignore_p1 and 1 in sg_energies:
        del sg_energies[1]

    if not sg_energies:
        print("[viz2] Warning: no data remaining for violin plot")
        return

    sorted_by_count = sorted(sg_energies.items(), key=lambda x: len(x[1]), reverse=True)
    selected = sorted_by_count[:top_n]
    selected.sort(key=lambda x: x[0])  # sort by SG number

    sg_labels = [f"SG {sg}" for sg, _ in selected]
    sg_values = [vals for _, vals in selected]
    sg_counts = [len(vals) for vals in sg_values]

    print(f"[viz2] Violin plot: {len(selected)} space groups, "
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
    print(f"[viz2] Saved: {output_path}")


def _build_parser():
    parser = argparse.ArgumentParser(
        description='GEWUM viz2 - Structure Analysis Visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python -m gewum.src.common.postprocess.viz2_analysis --cif-dir ./Na2Cl2 --dim 3d
  python -m gewum.src.common.postprocess.viz2_analysis --cif-dir ./Na2Cl2 ./Na2Cl3 --plot sankey
  python -m gewum.src.common.postprocess.viz2_analysis --cif-dir ./Na2Cl2 --plot rdf --rdf-rmax 8.0
""",
    )
    parser.add_argument('--cif-dir', action='append', nargs='+', required=True,
                        help='One or more CIF root directories. Accepts both '
                             '"--cif-dir A B C" and repeated "--cif-dir A --cif-dir B".')
    parser.add_argument('--dim', choices=['0d', '1d', '2d', '3d'], default='3d',
                        help='Structure dimension (default: 3d)')
    parser.add_argument('--plot', choices=['all', 'sankey', 'chord', 'rdf', 'funnel', 'ehull', 'violin'], default='all',
                        help='Which plot(s) to generate (default: all)')
    parser.add_argument('--symprec', type=float, default=0.1,
                        help='Symmetry precision for spglib (default: 0.1)')
    parser.add_argument('--rdf-rmax', type=float, default=10.0,
                        help='RDF maximum radius in Angstrom (default: 10.0)')
    parser.add_argument('--rdf-bin', type=float, default=0.1,
                        help='RDF bin size in Angstrom (default: 0.1)')
    parser.add_argument('--top-n', type=int, default=15,
                        help='Top-N flows to show in Sankey (default: 15)')
    parser.add_argument('--max-structures', type=int, default=200,
                        help='Max structures per stage for RDF (default: 200)')
    parser.add_argument('-o', '--output', default='./viz2_output',
                        help='Output directory (default: ./viz2_output)')
    parser.add_argument('--dpi', type=int, default=300,
                        help='Figure DPI (default: 300)')
    parser.add_argument('--ignore-p1', action='store_true', default=False,
                        help='Exclude flows where relaxed space group is P1 from Sankey diagram')
    parser.add_argument('--api-key',
                        help='Materials Project API key (for Ehull calculation)')
    parser.add_argument('--mp-data',
                        help='Path to offline MP JSON file (for Ehull calculation)')
    parser.add_argument('--ehull-compat', action='store_true', default=False,
                        help='Apply MP2020 compatibility corrections for Ehull calculation')
    parser.add_argument('--ehull-cmap', default=None,
                        help='Colormap for sankey_ehull.png columns (e.g. viridis, plasma). '
                             'Left/right columns share the same colormap; '
                             'middle column uses the reversed variant. '
                             'Default: blue-pink gradient (#A8D8EA -> #F3C1D3)')
    parser.add_argument('--font-family', type=str, default='Liberation Sans',
                        help='Font family used for figures (default: Liberation Sans)')
    parser.add_argument('--font-size', type=int, default=None,
                        help='Global font size for figures (default: matplotlib default)')
    parser.add_argument('--violin-top-n', type=int, default=20,
                        help='Number of top space groups to display in violin plot (default: 20)')
    parser.add_argument('--workers', type=int, default=min(os.cpu_count() or 4, 4),
                        help='Number of parallel workers (default: min(cpu_count, 4))')
    parser.add_argument('--sankey-cmap-left', default=None,
                        help='Colormap for left column in sankey_spacegroup.png (e.g. viridis, Pastel2)')
    parser.add_argument('--sankey-cmap-right', default=None,
                        help='Colormap for right column in sankey_spacegroup.png (e.g. viridis_r, Pastel2)')
    return parser


def main():
    parser = _build_parser()
    args = parser.parse_args()

    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = [args.font_family] + matplotlib.rcParams['font.sans-serif']
    if args.font_size is not None:
        matplotlib.rcParams['font.size'] = args.font_size

    os.makedirs(args.output, exist_ok=True)

    cif_dirs = [d for sub in args.cif_dir for d in sub]
    cif_dict = _auto_collect_cifs(cif_dirs)

    if not cif_dict:
        print("[viz2] Error: no CIF files found in any provided directory.")
        sys.exit(1)

    plot_type = args.plot
    workers = args.workers
    print(f"[viz2] Workers: {workers}")

    # Collect energy data once for ehull + violin (avoid duplicate I/O on large CIF sets)
    energy_data_cache = None
    if plot_type in ('all', 'ehull', 'violin'):
        if 'relaxed' in cif_dict and cif_dict['relaxed']:
            energy_data_cache = _collect_energy_data(cif_dict['relaxed'], workers=workers)

    if plot_type in ('all', 'sankey'):
        if 'relaxed' in cif_dict and cif_dict['relaxed']:
            sankey_path = os.path.join(args.output, 'sankey_spacegroup.png')
            plot_sankey(cif_dict['relaxed'], sankey_path,
                        symprec=args.symprec, top_n=args.top_n,
                        dpi=args.dpi,
                        ignore_p1=args.ignore_p1,
                        workers=workers,
                        cmap_left=args.sankey_cmap_left,
                        cmap_right=args.sankey_cmap_right)
        else:
            print("[viz2] Skipping Sankey: no relaxed structures found")

    if plot_type in ('all', 'chord'):
        if 'relaxed' in cif_dict and cif_dict['relaxed']:
            chord_path = os.path.join(args.output, 'chord_spacegroup.png')
            plot_chord(cif_dict['relaxed'], chord_path,
                       symprec=args.symprec, top_n=args.top_n,
                       dpi=args.dpi,
                       ignore_p1=args.ignore_p1,
                       workers=workers)
        else:
            print("[viz2] Skipping Chord: no relaxed structures found")

    if plot_type in ('all', 'ehull'):
        if 'relaxed' in cif_dict and cif_dict['relaxed']:
            if not args.api_key and not args.mp_data:
                print("[viz2] Skipping Ehull Sankey: need --api-key or --mp-data")
            else:
                try:
                    ehull_path = os.path.join(args.output, 'sankey_ehull.png')
                    plot_sankey_ehull(cif_dict['relaxed'], ehull_path,
                                      mp_api_key=args.api_key, mp_data_path=args.mp_data,
                                      symprec=args.symprec, top_n=args.top_n,
                                      dpi=args.dpi, ignore_p1=args.ignore_p1,
                                      use_compat=args.ehull_compat,
                                      workers=workers,
                                      cmap=args.ehull_cmap,
                                      energy_data_cache=energy_data_cache)
                except Exception as e:
                    print(f"[viz2] ERROR: sankey_ehull failed: {e}")
        else:
            print("[viz2] Skipping Ehull Sankey: no relaxed structures found")

    if plot_type in ('all', 'rdf'):
        rdf_path = os.path.join(args.output, 'rdf_comparison.png')
        plot_rdf_comparison(cif_dict, rdf_path,
                            r_max=args.rdf_rmax, bin_size=args.rdf_bin,
                            max_structures=args.max_structures,
                            dpi=args.dpi,
                            workers=workers)

    if plot_type in ('all', 'funnel'):
        funnel_path = os.path.join(args.output, 'structure_funnel.png')
        plot_funnel(cif_dict, funnel_path, dpi=args.dpi)

    if plot_type in ('all', 'violin'):
        if 'relaxed' in cif_dict and cif_dict['relaxed']:
            violin_path = os.path.join(args.output, 'violin_energy.png')
            plot_violin_energy(cif_dict['relaxed'], violin_path,
                               symprec=args.symprec, dpi=args.dpi,
                               ignore_p1=args.ignore_p1, workers=workers,
                               top_n=args.violin_top_n,
                               energy_data_cache=energy_data_cache)
        else:
            print("[viz2] Skipping Violin: no relaxed structures found")

    report_path = os.path.join(args.output, 'viz2_report.txt')
    generate_report(cif_dict, report_path)

    print(f"\n[viz2] All outputs written to {os.path.abspath(args.output)}")

    import gc
    import matplotlib.pyplot as plt
    plt.close('all')
    gc.collect()


if __name__ == '__main__':
    main()
