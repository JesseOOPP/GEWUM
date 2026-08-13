"""
GEWUM CIF Deduplication Module
Composition- and space-group-aware structure deduplication:
1. Identify reduced formula and space group for each structure
2. Group by (reduced formula, space group): structures from different
   compositions are never compared; at least one structure per SG is kept
3. Within each group, use StructureMatcher to remove duplicates
4. (Optional --rdf) after exact dedup, merge structurally similar
   structures within each (composition, space group) using RDF fingerprint
   distance, so that near-identical conformers that StructureMatcher cannot
   superimpose are collapsed to one representative

Output filenames are prefixed with the reduced formula and SG number
(e.g. Fe2O3_sg167_xxx.cif), so structures from different compositions can
be told apart at a glance.

Usage:
    python -m gewum.src.RDworkflows.cif_dedup
    python -m gewum.src.RDworkflows.cif_dedup --input-dir ./relaxed
    python -m gewum.src.RDworkflows.cif_dedup --rdf --sim-t 0.2
"""
import os
import shutil
import argparse
from pathlib import Path
from collections import defaultdict

import numpy as np

from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher


def _read_sg_from_cif_header(cif_path):
    """Read the space group number from the CIF header.

    Looks for ``_symmetry_Int_Tables_number`` / ``_space_group_IT_number``,
    consistent with the sym module renaming convention. Returns None if the
    header is missing or unparseable.
    """
    try:
        with open(cif_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if '_symmetry_Int_Tables_number' in line:
                    parts = line.split()
                    if len(parts) >= 2:
                        val = parts[-1].strip("'\"")
                        if val.isdigit():
                            return int(val)
                    next_line = next(f, '').strip()
                    if next_line.isdigit():
                        return int(next_line)
                if '_space_group_IT_number' in line:
                    parts = line.split()
                    if len(parts) >= 2:
                        val = parts[-1].strip("'\"")
                        if val.isdigit():
                            return int(val)
                    next_line = next(f, '').strip()
                    if next_line.isdigit():
                        return int(next_line)
    except Exception:
        pass
    return None


def identify_spacegroup(cif_path, std_struct, sg_analyzed=None):
    """Determine the space group number used for grouping.

    Priority order:
    1. CIF header (``_symmetry_Int_Tables_number``), consistent with the sym
       module renaming convention.
    2. If the header is missing or is a placeholder P1 (typical of raw random
       CIFs), fall back to the actual symmetry analysis of the standardized
       structure, so that isomorphous structures with a placeholder header
       are still grouped together and deduplicated.

    Returns:
        int space group number (1 if unknown).
    """
    header_sg = _read_sg_from_cif_header(cif_path)
    if header_sg and header_sg != 1:
        return header_sg
    if sg_analyzed is not None:
        return sg_analyzed
    return header_sg or 1


def get_standard_structure(structure):
    """Get the standardized primitive structure for comparison.

    Returns:
        (std_struct, sg_analyzed): standardized primitive structure and its
        space group number from actual symmetry analysis (None on failure).
    """
    try:
        sga = SpacegroupAnalyzer(structure, symprec=0.01)
        return sga.get_primitive_standard_structure(), sga.get_space_group_number()
    except Exception as e:
        print(f"  Warning: Standardization failed, using primitive cell. Error: {e}")
        return structure.get_primitive_structure(), None


def get_unique_filename(base_name, used_names):
    """Generate a unique filename to avoid collisions."""
    stem = base_name.rsplit(".", 1)[0]
    suffix = ".cif"
    counter = 1
    new_name = base_name
    while new_name in used_names:
        new_name = f"{stem}_{counter}{suffix}"
        counter += 1
    used_names.add(new_name)
    return new_name


def collect_cif_files(root=".", recursive=True):
    """Collect CIF files from directory."""
    cif_files = []
    if recursive:
        for dirpath, _, filenames in os.walk(root):
            for name in sorted(filenames):
                if name.lower().endswith(".cif"):
                    cif_files.append(os.path.join(dirpath, name))
    else:
        root_path = Path(root)
        for f in sorted(root_path.iterdir()):
            if f.is_file() and f.suffix.lower() == '.cif':
                cif_files.append(str(f))
    return cif_files


def build_output_name(formula, sg_num, original_name):
    """Prefix the output filename with composition and SG for easy grouping.

    The original name is kept as the tail so the source file stays traceable.
    """
    if not original_name.startswith(f"{formula}_sg"):
        return f"{formula}_sg{sg_num}_{original_name}"
    return original_name


def _structure_to_rdf(structure, r_max=8.0, n_bins=80):
    """Convert a pymatgen Structure to a normalised RDF feature vector.

    Creates a supercell large enough to capture all neighbours within r_max
    under periodic boundary conditions, collects all pairwise distances from
    original-cell atoms, then bins into a histogram normalised to unit sum.

    Independent copy of scf_viz._structure_to_rdf (avoids cross-module import).

    Args:
        structure: pymatgen Structure
        r_max:     maximum pair distance (A)
        n_bins:    number of histogram bins

    Returns:
        1D numpy array of length n_bins, normalised to unit sum.
    """
    lattice = structure.lattice
    a_norm = float(np.linalg.norm(lattice.matrix[0]))
    b_norm = float(np.linalg.norm(lattice.matrix[1]))
    c_norm = float(np.linalg.norm(lattice.matrix[2]))

    na = max(1, min(int(np.ceil(2 * r_max / a_norm)), 4))
    nb = max(1, min(int(np.ceil(2 * r_max / b_norm)), 4))
    nc = max(1, min(int(np.ceil(2 * r_max / c_norm)), 4))

    supercell = structure * (na, nb, nc)
    n_orig = len(structure)
    orig_coords = supercell.cart_coords[:n_orig]
    all_coords = supercell.cart_coords

    distances = []
    for i in range(n_orig):
        diff = all_coords - orig_coords[i]
        dists = np.sqrt(np.sum(diff ** 2, axis=1))
        mask = (dists > 0.01) & (dists <= r_max)
        distances.extend(dists[mask].tolist())

    if not distances:
        return np.zeros(n_bins)

    hist, _ = np.histogram(distances, bins=n_bins, range=(0, r_max))
    total = float(np.sum(hist))
    if total > 0:
        hist = hist / total
    return hist.astype(float)


def merge_similar_by_fingerprint(entries, threshold):
    """Merge structurally similar structures within one (formula, SG) group.

    Greedy single-link clustering on RDF fingerprint distance: each entry is
    compared against the already-accepted representatives in order; if its
    distance to the nearest representative is below ``threshold`` it is
    merged (dropped), otherwise it becomes a new representative.

    Args:
        entries: list of (cif_path, std_struct, sg_num)
        threshold: RDF fingerprint distance threshold for merging

    Returns:
        (kept, merged): kept entries, and a list of
        (merged_index, rep_index, distance) for logging.
    """
    n = len(entries)
    if n < 2:
        return entries, []

    fps = [_structure_to_rdf(std) for _, std, _ in entries]
    kept_flags = [True] * n
    reps = [0]
    merged = []

    for i in range(1, n):
        dists = [float(np.linalg.norm(fps[i] - fps[j])) for j in reps]
        d_min = min(dists)
        if d_min < threshold:
            kept_flags[i] = False
            merged.append((i, reps[dists.index(d_min)], d_min))
        else:
            reps.append(i)

    kept = [e for e, keep in zip(entries, kept_flags) if keep]
    return kept, merged


def main():
    parser = argparse.ArgumentParser(
        description='GEWUM CIF Deduplication - composition & space-group-aware structure matching',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Strategy:
  1. Read space group number from CIF header (_symmetry_Int_Tables_number);
     if missing or placeholder P1, fall back to actual symmetry analysis
  2. Group structures by (reduced formula, space group)
  3. Keep at least one structure per space group
  4. Within each group, use StructureMatcher to remove duplicates
  5. Output filenames carry the formula/SG prefix (Fe2O3_sg167_xxx.cif)

--rdf: after exact StructureMatcher dedup, merge structurally
similar structures within each (composition, space group) using RDF
fingerprint distance (--sim-t, default 0.2). Near-identical conformers
that cannot be superimposed exactly are collapsed to one representative
per fingerprint cluster. The log lists every merged file with its distance
so the effect can be inspected.

Examples:
  python -m gewum.src.RDworkflows.cif_dedup
  python -m gewum.src.RDworkflows.cif_dedup --input-dir ./relaxed
  python -m gewum.src.RDworkflows.cif_dedup --ltol 0.03
  python -m gewum.src.RDworkflows.cif_dedup --rdf --sim-t 0.2
""",
    )
    parser.add_argument('--input-dir', '-d', default='.',
                        help='Directory containing CIF files (default: current directory)')
    parser.add_argument('--output-dir', '-o', default='0_s_cif',
                        help='Output directory for unique structures (default: 0_s_cif)')
    parser.add_argument('--ltol', type=float, default=0.1,
                        help='StructureMatcher length tolerance (default: 0.1)')
    parser.add_argument('--stol', type=float, default=0.1,
                        help='StructureMatcher site tolerance (default: 0.1)')
    parser.add_argument('--angle-tol', type=float, default=1,
                        help='StructureMatcher angle tolerance (default: 1)')
    parser.add_argument('--rdf', action='store_true',
                        help='After exact dedup, merge structurally similar structures within each '
                             'composition using RDF fingerprint distance (threshold: --sim-t)')
    parser.add_argument('--sim-t', type=float, default=0.2,
                        help='RDF fingerprint distance threshold for similar-structure merging '
                             '(default: 0.2)')
    parser.add_argument('--no-recursive', action='store_true',
                        help='Do not search subdirectories (flat mode)')
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    recursive = not args.no_recursive
    cif_files = collect_cif_files(str(input_dir), recursive=recursive)

    if not cif_files:
        print(f"[dedup] No CIF files found in {input_dir}")
        return

    print(f"[dedup] Found {len(cif_files)} CIF files in {input_dir}")
    print(f"[dedup] StructureMatcher: ltol={args.ltol}, stol={args.stol}, angle_tol={args.angle_tol}")
    print()

    # ------------------------------------------------------------------
    # Step 1: Parse all structures, identify composition + space group
    # ------------------------------------------------------------------
    # comp_groups: {formula: {sg_num: [(cif_path, std_struct, sg_num), ...]}}
    comp_groups = defaultdict(lambda: defaultdict(list))
    parse_errors = 0

    for idx, cif_path in enumerate(cif_files, 1):
        rel_path = os.path.relpath(cif_path)
        print(f"[{idx}/{len(cif_files)}] Parsing: {rel_path}", end="")

        try:
            cif_parser = CifParser(cif_path)
            structures = cif_parser.get_structures(primitive=False)
            if not structures:
                print(" -> No structure, skipped.")
                parse_errors += 1
                continue
            struct = structures[0]
            std_struct, sg_analyzed = get_standard_structure(struct)
            sg_num = identify_spacegroup(cif_path, std_struct, sg_analyzed)
            formula = std_struct.composition.reduced_formula
            comp_groups[formula][sg_num].append((cif_path, std_struct, sg_num))
            print(f" -> {formula} / SG {sg_num}")
        except Exception as e:
            print(f" -> Error: {e}")
            parse_errors += 1
            continue

    total_parsed = sum(len(v) for sg_d in comp_groups.values() for v in sg_d.values())
    print(f"\n[dedup] Parsed successfully: {total_parsed}")
    print(f"[dedup] Parse errors: {parse_errors}")
    print(f"[dedup] Unique compositions found: {len(comp_groups)}")
    unique_sgs = {sg for sg_d in comp_groups.values() for sg in sg_d}
    print(f"[dedup] Unique space groups found: {len(unique_sgs)}")
    print()

    # ------------------------------------------------------------------
    # Step 2: Within each (formula, SG) group, deduplicate with StructureMatcher
    # ------------------------------------------------------------------
    matcher = StructureMatcher(
        ltol=args.ltol,
        stol=args.stol,
        angle_tol=args.angle_tol,
        primitive_cell=True,
    )

    kept_by_formula = defaultdict(list)
    total_kept = 0
    total_duplicates = 0

    for formula in sorted(comp_groups.keys()):
        sg_dict = comp_groups[formula]
        n_comp = sum(len(v) for v in sg_dict.values())
        print(f"=== Composition {formula}: {n_comp} structures ===")

        for sg_key in sorted(sg_dict.keys()):
            entries = sg_dict[sg_key]
            sg_label = str(sg_key)
            print(f"--- {formula} / SG {sg_label}: {len(entries)} structures ---")

            unique_in_sg = []

            for cif_path, std_struct, sg_num in entries:
                is_duplicate = False
                for ref_std in unique_in_sg:
                    try:
                        if matcher.fit(ref_std, std_struct):
                            is_duplicate = True
                            break
                    except Exception:
                        continue

                if is_duplicate:
                    print(f"  Duplicate: {os.path.basename(cif_path)}")
                    total_duplicates += 1
                else:
                    unique_in_sg.append(std_struct)
                    kept_by_formula[formula].append((cif_path, std_struct, sg_num))
                    print(f"  Kept: {os.path.basename(cif_path)}")
                    total_kept += 1

            print(f"  -> Kept {len(unique_in_sg)}/{len(entries)}")
            print()

    # ------------------------------------------------------------------
    # Step 3 (optional): merge structurally similar structures by fingerprint
    # Within each (formula, SG) group: RDF fingerprints carry no symmetry
    # information, so merging across space groups could collapse genuinely
    # different conformers that happen to share similar pair-distance
    # histograms. Keeping the merge inside each SG preserves the
    # at-least-one-per-SG guarantee.
    # ------------------------------------------------------------------
    total_merged = 0
    if args.rdf:
        print(f"[dedup] --rdf: RDF fingerprint threshold = {args.sim_t}")
        for formula in sorted(kept_by_formula.keys()):
            sg_groups = defaultdict(list)
            for e in kept_by_formula[formula]:
                sg_groups[e[2]].append(e)
            kept_entries = []
            for sg_num in sorted(sg_groups.keys()):
                entries = sg_groups[sg_num]
                kept, merged = merge_similar_by_fingerprint(entries, args.sim_t)
                for m_idx, r_idx, dist in merged:
                    print(f"  Similar-merged: {os.path.basename(entries[m_idx][0])} "
                          f"(RDF distance {dist:.3f} vs {os.path.basename(entries[r_idx][0])})")
                total_merged += len(merged)
                kept_entries.extend(kept)
            kept_by_formula[formula] = kept_entries
        print()

    # ------------------------------------------------------------------
    # Step 4: copy the surviving structures to the output directory
    # ------------------------------------------------------------------
    saved_filenames = set()
    for formula in sorted(kept_by_formula.keys()):
        for cif_path, std_struct, sg_num in kept_by_formula[formula]:
            original_name = os.path.basename(cif_path)
            safe_name = get_unique_filename(
                build_output_name(formula, sg_num, original_name),
                saved_filenames,
            )
            shutil.copy2(cif_path, output_dir / safe_name)

    print("=" * 60)
    print("[dedup] Deduplication complete!")
    print(f"  Compositions: {len(comp_groups)}")
    print(f"  Total structures: {total_parsed}")
    print(f"  Kept (unique): {total_kept - total_merged}")
    print(f"  Removed (duplicates): {total_duplicates}")
    if args.rdf:
        print(f"  Merged (similar): {total_merged}")
    print(f"  Output directory: {output_dir.absolute()}")
    print("=" * 60)


if __name__ == "__main__":
    main()
