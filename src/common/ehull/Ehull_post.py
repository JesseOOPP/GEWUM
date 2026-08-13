"""
GEWUM Energy Hull Post-processing Module
Filter structures by energy above hull threshold and copy selected CIF files
"""
import pandas as pd
import os
import csv
import zipfile
from shutil import copy
import argparse


def filter_and_save_csv(input_file, output_file, threshold):
    """
    Filter structures by e_above_hull threshold.
    
    Args:
        input_file: Path to input CSV with e_above_hull column
        output_file: Path to save filtered results
        threshold: Maximum e_above_hull value (eV/atom)
    """
    df = pd.read_csv(input_file)
    filtered_df = df[df['e_above_hull'] < threshold]
    filtered_df.to_csv(output_file, index=False)
    print(f"Filtered {len(filtered_df)}/{len(df)} structures with e_above_hull < {threshold} eV")
    return filtered_df


def get_unique_filename(destination_path):
    """Generate unique filename to avoid overwriting"""
    base, ext = os.path.splitext(destination_path)
    counter = 1
    while os.path.exists(destination_path):
        destination_path = f"{base}_{counter}{ext}"
        counter += 1
    return destination_path


def _extract_from_zip(source_path):
    """
    Try to extract a CIF file from its formula directory's structures.db or structures.zip.

    Tries SQLite DB first, then falls back to ZIP.
    When CIF files are archived via cif_db or cif_archive, the original disk path
    no longer exists. This function locates the corresponding DB/ZIP and
    extracts the single file to a temporary location so it can be copied.

    Args:
        source_path: Absolute path like /work/Na2Cl2/34/relaxed/xtal_1_relaxed.cif

    Returns:
        Extracted file path if successful, None otherwise
    """
    import sqlite3
    parts = os.path.normpath(source_path).split(os.sep)
    for i in range(len(parts) - 1, 0, -1):
        candidate = os.sep.join(parts[:i])

        db_path = os.path.join(candidate, 'structures.db')
        if os.path.isfile(db_path):
            rel_parts = os.path.relpath(source_path, candidate).replace(os.sep, '/').split('/')
            try:
                sg_number = int(rel_parts[0])
            except (ValueError, IndexError):
                continue
            cif_name = rel_parts[-1]
            try:
                conn = sqlite3.connect(db_path)
                conn.row_factory = sqlite3.Row
                try:
                    row = conn.execute(
                        "SELECT cif_content FROM structures WHERE sg_number = ? AND cif_name = ?",
                        (sg_number, cif_name)
                    ).fetchone()
                    if row:
                        extracted = os.path.join(candidate, *rel_parts)
                        os.makedirs(os.path.dirname(extracted), exist_ok=True)
                        with open(extracted, 'w', encoding='utf-8') as f:
                            f.write(row['cif_content'])
                        if os.path.isfile(extracted):
                            return extracted
                finally:
                    conn.close()
            except (sqlite3.Error, KeyError):
                pass

        zip_path = os.path.join(candidate, 'structures.zip')
        if os.path.isfile(zip_path):
            rel_path = os.path.relpath(source_path, candidate)
            rel_path_zip = rel_path.replace(os.sep, '/')
            try:
                with zipfile.ZipFile(zip_path, 'r') as zf:
                    if rel_path_zip in zf.namelist():
                        zf.extract(rel_path_zip, candidate)
                        extracted = os.path.join(candidate, rel_path)
                        if os.path.isfile(extracted):
                            return extracted
            except (zipfile.BadZipFile, KeyError):
                pass
    return None


def _parse_db_uri(uri):
    """Parse 'db://<sg>/<cif_name>' into (sg:int, cif_name:str). Returns (None, None) on bad URI."""
    if not isinstance(uri, str) or not uri.startswith('db://'):
        return None, None
    rest = uri[len('db://'):]
    if '/' not in rest:
        return None, None
    sg_str, cif_name = rest.split('/', 1)
    try:
        return int(sg_str), cif_name
    except ValueError:
        return None, None


def _collect_db_dirs(base_dir):
    """Recursively find all directories under base_dir that contain a structures.db.

    Uses glob for robust nested-directory discovery, covering layouts like:
      - base_dir/structures.db
      - base_dir/<formula>/structures.db
      - base_dir/round_N/<formula>/structures.db
    """
    import glob as _glob

    base_dir = os.path.abspath(base_dir)
    found = set()
    try:
        matches = _glob.glob(os.path.join(base_dir, '**', 'structures.db'), recursive=True)
    except Exception:
        return list(found)
    for db_path in matches:
        d = os.path.dirname(db_path)
        if d not in found:
            found.add(d)
    return list(found)


def _read_db_uri(sg_number, cif_name, search_dirs, formula_name=None):
    """Resolve (sg, cif_name) by querying candidate structures.db files.

    When *formula_name* is provided, only directories whose basename
    matches the formula (normalised via pymatgen Composition) are
    searched.  This prevents cross-formula collisions where different
    compositions happen to share (sg_number, cif_name) tuples.

    Returns the cif_content text on hit, or None when no DB has the row.
    """
    import sqlite3

    # Filter by formula when known
    if formula_name:
        from pymatgen.core.composition import Composition as _Comp
        try:
            target_formula = _Comp(formula_name).reduced_formula
        except Exception:
            target_formula = formula_name
        search_dirs = [
            d for d in search_dirs
            if _Comp(os.path.basename(d)).reduced_formula == target_formula
        ]

    for d in search_dirs:
        db_path = os.path.join(d, 'structures.db')
        if not os.path.isfile(db_path):
            continue
        try:
            conn = sqlite3.connect(db_path)
            try:
                row = conn.execute(
                    "SELECT cif_content FROM structures WHERE sg_number=? AND cif_name=?",
                    (sg_number, cif_name),
                ).fetchone()
                if row and row[0]:
                    return row[0]
            finally:
                conn.close()
        except sqlite3.Error:
            continue
    return None


def process_compositions(input_file, output_dir):
    """
    Copy CIF files for filtered compositions.

    Supports three sources for Relaxed_CIF_Path:
      1) plain filesystem path (legacy disk layout)
      2) 'db://<sg>/<cif_name>' URI (stage-2 zero-disk layout)
      3) archived path inside structures.zip / structures.db (legacy archive)

    Args:
        input_file: Path to filtered CSV with Relaxed_CIF_Path column
        output_dir: Directory to copy/materialize CIF files to
    """
    os.makedirs(output_dir, exist_ok=True)

    base_dir = os.path.dirname(os.path.abspath(input_file)) or '.'
    db_dirs = _collect_db_dirs(base_dir)

    copied_count = 0
    extracted_files = []  
    with open(input_file, mode='r') as file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            source_path = row.get('Relaxed_CIF_Path', '').strip()
            if not source_path:
                print(f"Warning: Empty Relaxed_CIF_Path")
                continue

            if source_path.startswith('db://'):
                sg_number, cif_name = _parse_db_uri(source_path)
                if sg_number is None:
                    print(f"Warning: malformed db URI: {source_path}")
                    continue

                formula = (row.get('Chemical_Formula')
                           or row.get('chemical_formula') or '').strip()

                cif_text = _read_db_uri(sg_number, cif_name, db_dirs,
                                        formula_name=formula)
                if cif_text is None:
                    print(f"Warning: db URI not found in any structures.db: {source_path}")
                    continue

                dest_path = os.path.join(output_dir, cif_name)
                unique_dest = get_unique_filename(dest_path)
                try:
                    with open(unique_dest, 'w', encoding='utf-8') as fh:
                        fh.write(cif_text)
                    copied_count += 1
                    print(f"Resolved: {source_path} -> {unique_dest}")
                except OSError as e:
                    print(f"Error writing {unique_dest}: {e}")
                continue

            actual_path = source_path
            if not os.path.isfile(source_path):
                extracted = _extract_from_zip(source_path)
                if extracted:
                    actual_path = extracted
                    extracted_files.append(extracted)
                else:
                    print(f"Warning: Invalid or missing Relaxed_CIF_Path: {source_path}")
                    continue

            cif_filename = os.path.basename(actual_path)
            dest_path = os.path.join(output_dir, cif_filename)
            unique_dest = get_unique_filename(dest_path)

            try:
                copy(actual_path, unique_dest)
                copied_count += 1
                print(f"Copied: {source_path} -> {unique_dest}")
            except Exception as e:
                print(f"Error copying {source_path}: {e}")

    for f in extracted_files:
        try:
            os.remove(f)
        except OSError:
            pass

    print(f"\nSuccessfully copied {copied_count} CIF files to {output_dir}")
    return copied_count


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(description='Post-process energy hull results')
    parser.add_argument('--input', '-i', default='Hull_result.csv',
                        help='Input CSV file with hull results')
    parser.add_argument('--output', '-o', default='final_result_02.csv',
                        help='Output filtered CSV file')
    parser.add_argument('--threshold', '-t', type=float, default=0.2,
                        help='E_above_hull threshold in eV/atom (default: 0.2)')
    parser.add_argument('--cif-dir', '-d', default='0_cif',
                        help='Directory to copy selected CIF files')
    parser.add_argument('-N', '--no-hull', action='store_true',
                        help='Skip hull filtering, extract all CIFs from 0_final_result_tot.txt directly')
    parser.add_argument('--raw-input', default='0_final_result_tot.txt',
                        help='Raw input file for -N mode (default: 0_final_result_tot.txt)')
    args = parser.parse_args()
    
    if args.no_hull:
        print(f"Mode: No hull filtering. Extracting all CIFs from {args.raw_input}")
        process_compositions(args.raw_input, args.cif_dir)
    else:
        filter_and_save_csv(args.input, args.output, args.threshold)
        process_compositions(args.output, args.cif_dir)


if __name__ == "__main__":
    main()
