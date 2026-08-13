"""
GEWUM Energy Post-processing (DB-driven, stage 2).

Reads stage='relaxed' entries from the per-formula structures.db,
sorts by energy_per_atom, applies the energy-gap greedy filter and
the MAX_N cap, then materializes the selected structures' cif_content
into the 0_cif_final/ subdirectory.

Inputs:
  - cwd is a formula directory containing structures.db.

Outputs (under cwd):
  - energy_final.txt        : audit log of all relaxed entries
  - 0_final_results.txt     : selected subset (post-gap filtering)
  - 0_cif_final/*.cif       : the actual selected CIF files

Both output CSVs store Relaxed_CIF_Path as 'db://<sg>/<cif_name>' so
that downstream tools can always retrieve the source CIF from the DB.
"""

import os
import csv
import argparse
import tempfile

from gewum.src.common.cif_db import CifDatabase


def _formula_from_cif_text(cif_text):
    """Fallback formula parser for DB rows whose `formula` column is empty."""
    try:
        from pymatgen.core import Structure
    except Exception:
        return 'Unknown'
    fd, tmp_path = tempfile.mkstemp(suffix='.cif')
    os.close(fd)
    try:
        with open(tmp_path, 'w', encoding='utf-8') as fh:
            fh.write(cif_text)
        s = Structure.from_file(tmp_path)
        return s.composition.formula.replace(' ', '')
    except Exception:
        return 'Unknown'
    finally:
        try:
            os.remove(tmp_path)
        except OSError:
            pass


def _base_name(cif_name):
    return cif_name[:-4] if cif_name.lower().endswith('.cif') else cif_name


def main():
    parser = argparse.ArgumentParser(description='GEWUM Energy Post-processing (DB-driven)')
    parser.add_argument('--gap', '-g', type=float, default=0.005,
                        help='Minimum energy difference between selected structures in eV/atom (default: 0.005)')
    parser.add_argument('--max-n', '-n', type=int, default=5,
                        help='Maximum number of structures to select (default: 5)')
    args = parser.parse_args()

    current_dir = os.getcwd()
    chem_dir_name = os.path.basename(current_dir)
    output_dir = os.path.join(current_dir, '0_cif_final')
    os.makedirs(output_dir, exist_ok=True)

    db_path = os.path.join(current_dir, 'structures.db')
    if not os.path.isfile(db_path):
        print(f"[energy_post] no structures.db in {current_dir}, nothing to do")
        return

    with CifDatabase(current_dir) as db:
        rows = db.query_relaxed_for_post()

    if not rows:
        print("[energy_post] no relaxed entries with energy found.")
        return

    all_data = []
    for r in rows:
        formula = r.get('formula') or _formula_from_cif_text(r['cif_content'])
        sg = r['sg_number']
        name = r['cif_name']
        all_data.append({
            'chemical_formula': formula,
            'space_group': str(sg),
            'base_name': _base_name(name),
            'total_energy': r['energy'],
            'energy_per_atom': r['energy_per_atom'],
            'cif_content': r['cif_content'],
            'cif_name': name,
            'db_uri': f"db://{sg}/{name}",
        })

    with open(os.path.join(current_dir, 'energy_final.txt'), 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['Chemical_Formula', 'CIF_Base_Name', 'Total_Energy_eV',
                         'Energy_per_Atom_eV', 'Relaxed_CIF_Path', 'SG_ori'])
        for d in all_data:
            writer.writerow([
                d['chemical_formula'], d['base_name'],
                d['total_energy'], d['energy_per_atom'],
                d['db_uri'], d['space_group'],
            ])

    all_data.sort(key=lambda x: x['energy_per_atom'])

    print(f"\nTotal entries: {len(all_data)}")
    print("First 10 energies:")
    for i, d in enumerate(all_data[:10]):
        print(f"{i+1}. {d['energy_per_atom']:.6f} eV/atom")

    selected = []
    last_energy = None
    for entry in all_data:
        if len(selected) >= args.max_n:
            break
        if last_energy is None or (entry['energy_per_atom'] - last_energy) >= args.gap:
            selected.append(entry)
            last_energy = entry['energy_per_atom']

    print(f"\nSelected {len(selected)} structures with energy gap >= {args.gap} eV/atom")

    with open(os.path.join(current_dir, '0_final_results.txt'), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Chemical_Formula', 'CIF_Base_Name', 'Total_Energy_eV',
                         'Energy_per_Atom_eV', 'Relaxed_CIF_Path', 'SG_ori'])
        for d in selected:
            writer.writerow([
                d['chemical_formula'], d['base_name'],
                d['total_energy'], d['energy_per_atom'],
                d['db_uri'], d['space_group'],
            ])

    for idx, d in enumerate(selected, 1):
        out_name = (f"{chem_dir_name}_{d['chemical_formula']}_{d['space_group']}_"
                    f"{d['base_name']}_spaced_{idx}.cif")
        out_path = os.path.join(output_dir, out_name)
        with open(out_path, 'w', encoding='utf-8') as fh:
            fh.write(d['cif_content'])
        print(f"Wrote: {out_path}")

    print("\nFinal selected structures:")
    for i, d in enumerate(selected):
        diff = 0.0 if i == 0 else d['energy_per_atom'] - selected[i-1]['energy_per_atom']
        print(f"{i+1}. E = {d['energy_per_atom']:.6f} eV/atom (D = {diff:.6f}) "
              f"| {d['chemical_formula']} | SG: {d['space_group']}")


if __name__ == '__main__':
    main()
