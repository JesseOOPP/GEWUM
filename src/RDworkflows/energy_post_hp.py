"""
GEWUM High Pressure Energy Post-processing (DB-driven, stage 2).

Reads stage='relaxed' entries from structures.db, filters by pressure
tolerance against the target pressure, sorts by corrected enthalpy per
atom, applies the enthalpy-gap greedy filter and MAX_N cap, then writes
the selected CIFs to 0_cif_final/ from the DB cif_content.

Required DB columns: final_pressure, enthalpy_per_atom,
corrected_enthalpy_per_atom (created by CifDatabase._ensure_columns).
"""

import os
import csv
import argparse
import tempfile

from gewum.src.common.cif_db import CifDatabase


def _formula_from_cif_text(cif_text):
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
    parser = argparse.ArgumentParser(description='High pressure post-processing (DB-driven)')
    parser.add_argument('--pressure', type=float, default=100.0,
                        help='Target pressure in GPa (default: 100.0)')
    parser.add_argument('--tolerance', type=float, default=5.0,
                        help='Pressure tolerance in GPa (default: 5.0)')
    parser.add_argument('--gap', type=float, default=0.005,
                        help='Minimum enthalpy gap between selected structures in eV/atom (default: 0.005)')
    parser.add_argument('--max-n', type=int, default=5,
                        help='Maximum number of structures to select (default: 5)')
    args = parser.parse_args()

    current_dir = os.getcwd()
    chem_dir_name = os.path.basename(current_dir)
    output_dir = os.path.join(current_dir, '0_cif_final')
    os.makedirs(output_dir, exist_ok=True)

    db_path = os.path.join(current_dir, 'structures.db')
    if not os.path.isfile(db_path):
        print(f"[energy_post_hp] no structures.db in {current_dir}, nothing to do")
        return

    with CifDatabase(current_dir) as db:
        rows = db.query_relaxed_for_post()

    if not rows:
        print("[energy_post_hp] no relaxed entries found.")
        return

    all_data = []
    skipped_no_hp = 0
    skipped_pressure = 0
    for r in rows:
        if r.get('final_pressure') is None or r.get('corrected_enthalpy_per_atom') is None:
            skipped_no_hp += 1
            continue
        pressure_diff = abs(r['final_pressure'] - args.pressure)
        if pressure_diff > args.tolerance:
            skipped_pressure += 1
            continue
        formula = r.get('formula') or _formula_from_cif_text(r['cif_content'])
        sg = r['sg_number']
        name = r['cif_name']
        all_data.append({
            'chemical_formula': formula,
            'space_group': str(sg),
            'base_name': _base_name(name),
            'total_energy': r['energy'],
            'energy_per_atom': r['energy_per_atom'],
            'final_pressure': r['final_pressure'],
            'enthalpy_per_atom': r.get('enthalpy_per_atom'),
            'corrected_enthalpy': r['corrected_enthalpy_per_atom'],
            'cif_content': r['cif_content'],
            'cif_name': name,
            'db_uri': f"db://{sg}/{name}",
        })

    if skipped_no_hp:
        print(f"[energy_post_hp] skipped {skipped_no_hp} entries with missing HP fields")
    if skipped_pressure:
        print(f"[energy_post_hp] skipped {skipped_pressure} entries outside pressure tolerance")

    if not all_data:
        print("No valid data found within pressure tolerance.")
        return

    with open(os.path.join(current_dir, 'energy_final.txt'), 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['Chemical_Formula', 'CIF_Base_Name', 'Total_Energy_eV',
                         'Energy_per_Atom_eV', 'Final_Pressure_GPa',
                         'Enthalpy_per_Atom_eV', 'Corrected_Enthalpy_per_Atom_eV',
                         'Relaxed_CIF_Path', 'SG_ori'])
        for d in all_data:
            writer.writerow([
                d['chemical_formula'], d['base_name'],
                d['total_energy'], d['energy_per_atom'],
                d['final_pressure'], d['enthalpy_per_atom'], d['corrected_enthalpy'],
                d['db_uri'], d['space_group'],
            ])

    all_data.sort(key=lambda x: x['corrected_enthalpy'])

    print(f"\nTotal entries within pressure tolerance: {len(all_data)}")
    print("First 10 enthalpies:")
    for i, d in enumerate(all_data[:10]):
        print(f"{i+1}. H = {d['corrected_enthalpy']:.6f} eV/atom, P = {d['final_pressure']:.1f} GPa")

    selected = []
    last_h = None
    for entry in all_data:
        if len(selected) >= args.max_n:
            break
        if last_h is None or (entry['corrected_enthalpy'] - last_h) >= args.gap:
            selected.append(entry)
            last_h = entry['corrected_enthalpy']

    print(f"\nSelected {len(selected)} structures with enthalpy gap >= {args.gap} eV/atom")

    with open(os.path.join(current_dir, '0_final_results.txt'), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Chemical_Formula', 'CIF_Base_Name', 'Total_Energy_eV',
                         'Energy_per_Atom_eV', 'Final_Pressure_GPa',
                         'Enthalpy_per_Atom_eV', 'Corrected_Enthalpy_per_Atom_eV',
                         'Relaxed_CIF_Path', 'SG_ori'])
        for d in selected:
            writer.writerow([
                d['chemical_formula'], d['base_name'],
                d['total_energy'], d['energy_per_atom'],
                d['final_pressure'], d['enthalpy_per_atom'], d['corrected_enthalpy'],
                d['db_uri'], d['space_group'],
            ])

    for idx, d in enumerate(selected, 1):
        out_name = (f"{chem_dir_name}_{d['chemical_formula']}_{d['space_group']}_"
                    f"{d['base_name']}_hp_{idx}.cif")
        out_path = os.path.join(output_dir, out_name)
        with open(out_path, 'w', encoding='utf-8') as fh:
            fh.write(d['cif_content'])
        print(f"Wrote: {out_path}")

    print(f"\nFinal selected structures (Target P = {args.pressure} GPa):")
    for i, d in enumerate(selected):
        diff = 0.0 if i == 0 else d['corrected_enthalpy'] - selected[i-1]['corrected_enthalpy']
        print(f"{i+1}. H = {d['corrected_enthalpy']:.6f} eV/atom (D = {diff:.6f}) "
              f"| P = {d['final_pressure']:.1f} GPa | {d['chemical_formula']} "
              f"| SG: {d['space_group']}")


if __name__ == '__main__':
    main()
