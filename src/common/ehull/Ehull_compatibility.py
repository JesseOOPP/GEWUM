"""GEWUM Energy Hull Calculation Module (Unified)
Calculate formation energy and energy above hull using Materials Project data
Supports both online (MP API) and offline (local JSON) modes
Supports both with and without MP2020 compatibility corrections
"""
import logging
import re
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedEntry 
from pymatgen.core.composition import Composition 
import pandas as pd
import numpy as np
from collections import defaultdict
import warnings
import argparse
import time
import os
import json
import hashlib

warnings.filterwarnings("ignore", category=UserWarning, module='mp_api.client.core.client')

logging.basicConfig(filename='process.log', level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

# MP2020 POTCAR symbol mapping
POTCAR_MAP = {
    "Ac": "Ac", "Ag": "Ag", "Al": "Al", "Am": "Am", "Ar": "Ar", 
    "As": "As", "At": "At", "Au": "Au", "B": "B", "Ba": "Ba_sv", 
    "Be": "Be_sv", "Bi": "Bi", "Br": "Br", "C": "C", "Ca": "Ca_sv", 
    "Cd": "Cd", "Ce": "Ce", "Cf": "Cf", "Cl": "Cl", "Cm": "Cm", 
    "Co": "Co", "Cr": "Cr_pv", "Cs": "Cs_sv", "Cu": "Cu_pv", 
    "Dy": "Dy_3", "Er": "Er_3", "Eu": "Eu", "F": "F", "Fe": "Fe_pv", 
    "Fr": "Fr_sv", "Ga": "Ga_d", "Gd": "Gd", "Ge": "Ge_d", "H": "H", 
    "He": "He", "Hf": "Hf_pv", "Hg": "Hg", "Ho": "Ho_3", "I": "I", 
    "In": "In_d", "Ir": "Ir", "K": "K_sv", "Kr": "Kr", "La": "La", 
    "Li": "Li_sv", "Lu": "Lu_3", "Mg": "Mg_pv", "Mn": "Mn_pv", 
    "Mo": "Mo_pv", "N": "N", "Na": "Na_pv", "Nb": "Nb_pv", 
    "Nd": "Nd_3", "Ne": "Ne", "Ni": "Ni_pv", "Np": "Np", "O": "O", 
    "Os": "Os_pv", "P": "P", "Pa": "Pa", "Pb": "Pb_d", "Pd": "Pd", 
    "Pm": "Pm_3", "Po": "Po", "Pr": "Pr_3", "Pt": "Pt", "Pu": "Pu", 
    "Ra": "Ra_sv", "Rb": "Rb_sv", "Re": "Re_pv", "Rh": "Rh_pv", 
    "Rn": "Rn", "Ru": "Ru_pv", "S": "S", "Sb": "Sb", "Sc": "Sc_sv", 
    "Se": "Se", "Si": "Si", "Sm": "Sm_3", "Sn": "Sn_d", "Sr": "Sr_sv", 
    "Ta": "Ta_pv", "Tb": "Tb_3", "Tc": "Tc_pv", "Te": "Te", "Th": "Th", 
    "Ti": "Ti_pv", "Tl": "Tl_d", "Tm": "Tm_3", "U": "U", "V": "V_pv", 
    "W": "W_pv", "Xe": "Xe", "Y": "Y_sv", "Yb": "Yb_3", "Zn": "Zn", 
    "Zr": "Zr_sv"
}


def get_mp2020_potcar_symbols(element):
    """Get MP2020-compatible POTCAR symbol for an element"""
    return POTCAR_MAP.get(element, element)


def get_above_hull_and_formation_energy(file_path, mp_api_key=None, mp_data_path=None, use_compatibility=True):
    """
    Calculate energy above hull and formation energy for compounds.
    
    Args:
        file_path: Path to CSV file with formula and energy columns
        mp_api_key: Materials Project API key (for online mode)
        mp_data_path: Path to offline MP JSON file (for offline mode)
        use_compatibility: If True, use ComputedEntry + MP2020Compatibility corrections.
                          If False, use PDEntry with raw energies (no corrections).
        
    Note:
        Either mp_api_key or mp_data_path must be provided.
        If mp_data_path is provided, offline mode is used.
    
    Returns:
        DataFrame with added formation_energy_per_atom, e_above_hull columns
        (plus energy_correction and status columns when use_compatibility=True)
    """
    if not mp_api_key and not mp_data_path:
        raise ValueError("Either mp_api_key or mp_data_path must be provided")
    
    use_offline = mp_data_path is not None
    
    if use_compatibility:
        try:
            data = pd.read_csv(file_path, header=0)
            formulas = data['Chemical_Formula'].tolist()
            total_energies = data['Total_Energy_eV'].tolist()
            total_rows = len(data)
            logging.info(f"Loaded {total_rows} compounds from {file_path}")
        except Exception as e:
            logging.error(f"Loading error: {e}")
            raise
    else:
        try:
            data = pd.read_csv(file_path, sep=',')
            formulas = data.iloc[:, 0].tolist()
            total_energies = data.iloc[:, 2].tolist()
            total_rows = len(data)
            logging.info(f"Loaded {total_rows} compounds from {file_path} (no-compat mode)")
        except Exception as e:
            logging.error(f"Loading error: {e}")
            raise

    formation_energies = [np.nan] * total_rows
    e_above_hulls = [np.nan] * total_rows

    if use_compatibility:
        correction_values = [np.nan] * total_rows
        status_messages = [""] * total_rows

    unique_formulas = set(formulas)
    element_systems = defaultdict(list)
    for formula in unique_formulas:
        try:
            if use_compatibility:
                comp = Composition(formula)
                elements = sorted([el.symbol for el in comp.elements])
            else:
                elements = sorted(set(re.findall(r'[A-Z][a-z]*', formula)))
            element_system = '-'.join(elements)
            element_systems[element_system].append(formula)
        except Exception as e:
            logging.warning(f"Failed to parse unique formula {formula}: {e}")
            continue

    mp_entries_cache = {}

    if use_compatibility:
        compatibility = MaterialsProject2020Compatibility()
    
    if use_offline:
        logging.info("Using offline mode with local MP data")
        from .mp_offline_loader import MPOfflineLoader
        loader = MPOfflineLoader(mp_data_path)
        
        for element_system in element_systems.keys():
            logging.info(f"Loading entries for: {element_system}")
            try:
                if use_compatibility:
                    # Load with uncorrected energy so all MP entries pass through
                    # the SAME MP2020Compatibility pipeline as user entries
                    entries = loader.get_entries_in_chemsys(
                        element_system.split('-'), as_pd_entry=False, use_uncorrected=True
                    )
                    # Add potcar_symbols derived from composition elements
                    for entry in entries:
                        if isinstance(entry, ComputedEntry) and '_mp_elements' in entry.data:
                            elements = entry.data.pop('_mp_elements')
                            entry.parameters['run_type'] = 'GGA'
                            entry.parameters['is_hubbard'] = False
                            entry.parameters['hubbards'] = {}
                            entry.parameters['potcar_symbols'] = [
                                get_mp2020_potcar_symbols(el) for el in elements
                            ]
                            entry.parameters['potcar_spec'] = [
                                {"titel": f"PAW_PBE {get_mp2020_potcar_symbols(el)} 06Sep2000"}
                                for el in elements
                            ]
                    # Run MP entries through the same compatibility pipeline
                    processed_mp = compatibility.process_entries(entries)
                    if processed_mp:
                        entries = processed_mp
                    else:
                        logging.warning(
                            f"Compatibility processing returned empty for {element_system}, "
                            f"using uncorrected entries"
                        )
                else:
                    entries = loader.get_entries_in_chemsys(element_system.split('-'), as_pd_entry=True)
                mp_entries_cache[element_system] = entries
                logging.info(f"Found {len(entries)} entries for {element_system}")
            except Exception as e:
                logging.error(f"Error loading {element_system}: {e}")
                mp_entries_cache[element_system] = []
    else:
        logging.info("Using online mode with MP API")
        if use_compatibility:
            from mp_api.client import MPRester
            
            with MPRester(mp_api_key) as mpr:
                for element_system in element_systems.keys():
                    logging.info(f"Fetching: {element_system}")
                    try:
                        entries = mpr.get_entries_in_chemsys(element_system.split('-'))
                        processed_entries = []
                        for entry in entries:
                            if hasattr(entry, 'parameters') and 'potcar_symbols' in entry.parameters:
                                processed_entries.append(entry)
                            else:
                                try:
                                    new_entry = ComputedEntry(
                                        composition=entry.composition,
                                        energy=entry.energy,
                                        entry_id=entry.entry_id,
                                        parameters=getattr(entry, 'parameters', {})
                                    )
                                    processed_entries.append(new_entry)
                                except Exception:
                                    processed_entries.append(entry)
                        mp_entries_cache[element_system] = processed_entries
                        logging.info(f"Found {len(processed_entries)} entries for {element_system}")
                        time.sleep(1.5)
                    except Exception as e:
                        logging.error(f"Error fetching {element_system}: {e}")
                        mp_entries_cache[element_system] = []
        else:
            from mp_api.client import MPRester
            
            with MPRester(mp_api_key) as mpr:
                for element_system in element_systems.keys():
                    logging.info(f"Fetching MP entries for system: {element_system}")
                    try:
                        entries = mpr.get_entries_in_chemsys(element_system)
                        mp_entries_cache[element_system] = entries
                        logging.info(f"Found {len(entries)} MP entries for {element_system}")
                        time.sleep(1.5)
                    except Exception as e:
                        logging.error(f"Error fetching {element_system}: {e}")
                        mp_entries_cache[element_system] = []

    for idx, (formula, total_energy) in enumerate(zip(formulas, total_energies)):
        try:
            if use_compatibility:
                comp = Composition(formula)
                elements = sorted([el.symbol for el in comp.elements])
            else:
                elements = sorted(set(re.findall(r'[A-Z][a-z]*', formula)))
            element_system = '-'.join(elements)

            base_entries = mp_entries_cache.get(element_system, [])
            if not base_entries:
                if use_compatibility:
                    status_messages[idx] = f"No MP data for system: {element_system}"
                logging.warning(f"No MP data for system: {element_system}")
                continue

            if use_compatibility:
                potcar_symbols = [get_mp2020_potcar_symbols(el.symbol) for el in comp.elements]
                current_entry = ComputedEntry(
                    composition=comp,
                    energy=total_energy,
                    entry_id=f"CUSTOM_{idx}_{formula}",
                    parameters={
                        "run_type": "GGA",
                        "is_hubbard": False,
                        "potcar_symbols": potcar_symbols,
                        "hubbards": {},
                        "potcar_spec": [{"titel": f"PAW_PBE {sym} 06Sep2000"} for sym in potcar_symbols]
                    }
                )

                try:
                    adjustments = compatibility.get_adjustments(current_entry)
                    total_correction = sum(adj.value for adj in adjustments) if adjustments else 0.0
                    correction_values[idx] = total_correction

                    processed_entries = compatibility.process_entries([current_entry])
                    if processed_entries:
                        current_entry = processed_entries[0]
                    else:
                        raise ValueError("Compatibility processing returned empty list")
                except Exception as e:
                    status_messages[idx] = f"Compatibility error: {e}"
                    logging.warning(f"Compatibility error for {formula} (row {idx}): {e}")
                    continue
            else:
                current_entry = PDEntry(composition=formula, energy=total_energy)

            all_entries = base_entries + [current_entry]
            pd_calc = PhaseDiagram(all_entries)
            formation_energy = pd_calc.get_form_energy_per_atom(current_entry)
            e_above_hull = pd_calc.get_e_above_hull(current_entry)

            formation_energies[idx] = formation_energy
            e_above_hulls[idx] = e_above_hull

            if use_compatibility:
                status_messages[idx] = "Success"
                logging.info(f"Row {idx} {formula}: FE={formation_energy:.6f}, EH={e_above_hull:.6f}, corr={total_correction:.6f}")
            else:
                logging.info(
                    f"{formula}: Formation energy = {formation_energy:.4f} eV/atom, "
                    f"E_above_hull = {e_above_hull:.4f} eV/atom"
                )

        except Exception as e:
            if use_compatibility:
                status_messages[idx] = f"Processing error: {e}"
            logging.exception(f"Error at row {idx} ({formula}): {e}")

    data['formation_energy_per_atom'] = formation_energies
    data['e_above_hull'] = e_above_hulls

    if use_compatibility:
        data['energy_correction'] = correction_values
        data['status'] = status_messages

    return data


def load_element_references_from_dir(element_dir, needed_elements):
    """Relax user-provided elemental CIFs with uMLIP for hull references.

    For each element symbol in *needed_elements*, look for
    ``<element_dir>/<Symbol>.cif``. When present, relax it with the uMLIP
    (atoms + cell) and build a single-atom PDEntry from the resulting
    energy-per-atom. This lets a user override the default MP DFT elemental
    references so formation energies are computed on a fully uMLIP-consistent
    basis.

    Args:
        element_dir: directory holding per-element CIFs (e.g. ./element).
        needed_elements: iterable of element symbols required by the inputs.

    Returns:
        dict mapping element symbol -> PDEntry (only for CIFs found & relaxed).

    Caching:
        Relaxed energies are cached in ``<element_dir>/.element_energy_cache.json``
        keyed by element, guarded by the CIF content SHA-1 and the relax
        parameter signature. On subsequent rounds an unchanged CIF is served
        from cache (no re-relaxation). Editing/replacing the CIF invalidates it.
    """
    refs = {}
    if not element_dir or not os.path.isdir(element_dir):
        return refs

    try:
        from gewum.src.common.relaxation.umlip_relax import optimize_from_content
    except Exception as e:
        logging.warning(f"[self-hull] uMLIP relaxation unavailable, "
                        f"cannot load ./element references: {e}")
        return refs

    # Cross-round cache: avoid re-relaxing unchanged element CIFs. The signature
    # ties a cached energy to both the CIF content and the relax parameters, so
    # any change to either transparently forces a recompute.
    _RELAX_PARAMS_SIG = "mode2_fmax0.05_max200"
    cache_path = os.path.join(element_dir, ".element_energy_cache.json")
    cache = {}
    try:
        if os.path.isfile(cache_path):
            with open(cache_path, 'r', encoding='utf-8') as fh:
                loaded = json.load(fh)
            if isinstance(loaded, dict):
                cache = loaded
    except Exception as e:
        logging.warning(f"[self-hull] Ignoring unreadable element cache "
                        f"{cache_path}: {e}")
        cache = {}

    cache_dirty = False
    for element in sorted(set(needed_elements)):
        cif_path = os.path.join(element_dir, f"{element}.cif")
        if not os.path.isfile(cif_path):
            continue
        try:
            with open(cif_path, 'r', encoding='utf-8') as fh:
                cif_text = fh.read()
            cif_sha1 = hashlib.sha1(cif_text.encode('utf-8')).hexdigest()

            cached = cache.get(element)
            if (isinstance(cached, dict)
                    and cached.get('cif_sha1') == cif_sha1
                    and cached.get('params') == _RELAX_PARAMS_SIG
                    and 'energy_per_atom' in cached):
                energy_per_atom = float(cached['energy_per_atom'])
                refs[element] = PDEntry(
                    composition=Composition(element),
                    energy=energy_per_atom,
                    name=f"USER_ELEM_{element}",
                )
                logging.info(f"[self-hull] Element reference {element}: "
                             f"{energy_per_atom:.4f} eV/atom (cache hit)")
                continue

            _, _energy, energy_per_atom, _ = optimize_from_content(
                cif_text, f"element_{element}", mode=2, fmax=0.05, max_steps=200
            )
            energy_per_atom = float(energy_per_atom)
            refs[element] = PDEntry(
                composition=Composition(element),
                energy=energy_per_atom,
                name=f"USER_ELEM_{element}",
            )
            cache[element] = {
                'energy_per_atom': energy_per_atom,
                'cif_sha1': cif_sha1,
                'params': _RELAX_PARAMS_SIG,
            }
            cache_dirty = True
            logging.info(f"[self-hull] Element reference {element}: "
                         f"{energy_per_atom:.4f} eV/atom (user CIF + uMLIP)")
        except Exception as e:
            logging.warning(f"[self-hull] Failed to relax element CIF "
                            f"{cif_path}: {e}")

    if cache_dirty:
        try:
            tmp_path = cache_path + ".tmp"
            with open(tmp_path, 'w', encoding='utf-8') as fh:
                json.dump(cache, fh, indent=2, sort_keys=True)
            os.replace(tmp_path, cache_path)
        except Exception as e:
            logging.warning(f"[self-hull] Failed to write element cache "
                            f"{cache_path}: {e}")

    return refs


def get_self_hull(file_path, mp_data_path=None, use_compatibility=False,
                  element_dir=None):
    """
    Calculate formation energy and e_above_hull using ONLY self-computed structures.
    
    Unlike the standard function that mixes MP structures into the convex hull,
    this mode builds the hull exclusively from the input structures plus
    elemental reference energies (from MP offline data for proper formation
    energy definition).
    
    This is designed for COMP iterative workflows where mixing uMLIP-level
    structures with DFT-level MP data creates an unfair comparison.
    
    Args:
        file_path: Path to CSV with Chemical_Formula, Energy_per_Atom_eV, Total_Energy_eV columns
                   (or any 3-column format: formula, energy_per_atom, total_energy)
        mp_data_path: Path to offline MP JSON file (for elemental reference energies only).
                      If None, elemental references are derived from the input data itself.
        use_compatibility: If True, apply MP2020 corrections to user entries.
                          Elemental references from MP remain as PDEntry (already in MP framework).
        element_dir: Optional directory with per-element CIFs (<Symbol>.cif).
                     When present, each element's reference energy is computed
                     by relaxing the user CIF with uMLIP (highest priority);
                     elements not found here fall back to MP offline data.
    
    Returns:
        DataFrame with added formation_energy_per_atom and e_above_hull columns
    """
    try:
        data = pd.read_csv(file_path, sep=',')
        if 'Chemical_Formula' in data.columns and 'Total_Energy_eV' in data.columns:
            formulas = data['Chemical_Formula'].tolist()
            total_energies = data['Total_Energy_eV'].tolist()
        else:
            formulas = data.iloc[:, 0].tolist()
            total_energies = data.iloc[:, 2].tolist()
        total_rows = len(data)
        logging.info(f"[self-hull] Loaded {total_rows} compounds from {file_path}")
    except Exception as e:
        logging.error(f"[self-hull] Loading error: {e}")
        raise

    if use_compatibility:
        compatibility = MaterialsProject2020Compatibility()

    all_entries = []
    entry_indices = []  
    for idx, (formula, energy) in enumerate(zip(formulas, total_energies)):
        try:
            if use_compatibility:
                comp = Composition(formula)
                potcar_symbols = [get_mp2020_potcar_symbols(el.symbol) for el in comp.elements]
                entry = ComputedEntry(
                    composition=comp,
                    energy=float(energy),
                    entry_id=f"SELF_{idx}_{formula}",
                    parameters={
                        "run_type": "GGA",
                        "is_hubbard": False,
                        "potcar_symbols": potcar_symbols,
                        "hubbards": {},
                        "potcar_spec": [{"titel": f"PAW_PBE {sym} 06Sep2000"} for sym in potcar_symbols]
                    }
                )
                processed = compatibility.process_entries([entry])
                if processed:
                    entry = processed[0]
                else:
                    logging.warning(f"[self-hull] Compatibility processing failed for row {idx} ({formula}), using raw energy")
                    entry = PDEntry(composition=formula, energy=float(energy),
                                    name=f"SELF_{idx}_{formula}")
            else:
                entry = PDEntry(composition=formula, energy=float(energy),
                                name=f"SELF_{idx}_{formula}")
            all_entries.append(entry)
            entry_indices.append(idx)
        except Exception as e:
            logging.warning(f"[self-hull] Failed to create entry for row {idx} ({formula}): {e}")

    # Collect all elements present in the input compositions.
    all_elements = set()
    for formula in formulas:
        try:
            comp = Composition(formula)
            for el in comp.elements:
                all_elements.add(el.symbol)
        except Exception:
            continue

    # Layered elemental reference resolution:
    #   (1) user-provided ./element/<X>.cif relaxed with uMLIP (highest priority)
    #   (2) fall back to MP offline DFT elemental energy
    user_refs = load_element_references_from_dir(element_dir, all_elements)

    loader = None
    if mp_data_path is not None:
        from .mp_offline_loader import MPOfflineLoader
        loader = MPOfflineLoader(mp_data_path)

    for element in all_elements:
        if element in user_refs:
            all_entries.append(user_refs[element])
            continue
        if loader is not None:
            try:
                ref_entries = loader.get_entries_in_chemsys([element], as_pd_entry=True)
                for ref_entry in ref_entries:
                    ref_comp = Composition(ref_entry.composition)
                    if len(ref_comp.elements) == 1:
                        all_entries.append(ref_entry)
            except Exception as e:
                logging.warning(f"[self-hull] Failed to load reference for {element}: {e}")

    formation_energies = [np.nan] * total_rows
    e_above_hulls = [np.nan] * total_rows

    try:
        pd_calc = PhaseDiagram(all_entries)
        logging.info(f"[self-hull] Built PhaseDiagram with {len(all_entries)} entries "
                     f"(compatibility={'ON' if use_compatibility else 'OFF'})")

        for i, entry_idx in enumerate(entry_indices):
            entry = all_entries[i]
            try:
                formation_energies[entry_idx] = pd_calc.get_form_energy_per_atom(entry)
                e_above_hulls[entry_idx] = pd_calc.get_e_above_hull(entry)
            except Exception as e:
                logging.warning(
                    f"[self-hull] Error computing hull distance for row {entry_idx} "
                    f"({formulas[entry_idx]}): {e}"
                )
    except Exception as e:
        logging.error(f"[self-hull] Failed to build PhaseDiagram: {e}")
        raise

    data['formation_energy_per_atom'] = formation_energies
    data['e_above_hull'] = e_above_hulls

    return data


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(description='Calculate energy above hull (unified module)')
    parser.add_argument('--input', '-i', default='0_final_result_tot.txt',
                        help='Input CSV file path')
    parser.add_argument('--output', '-o', default='Hull_result.csv',
                        help='Output CSV file path')
    parser.add_argument('--api-key',
                        help='Materials Project API key (for online mode)')
    parser.add_argument('--mp-data',
                        help='Path to offline MP JSON file (for offline mode)')
    parser.add_argument('-cor', action='store_true', default=False,
                        help='Use MP2020 energy compatibility corrections')
    parser.add_argument('--self-hull', action='store_true', default=False,
                        help='Self-hull mode: build convex hull from input structures only '
                             '(no MP competing phases). MP data used only for elemental references.')
    parser.add_argument('--element-dir', default=None,
                        help='Directory with per-element CIFs (<Symbol>.cif) for '
                             'uMLIP-computed elemental references (self-hull mode only). '
                             'Elements not found fall back to MP offline data.')
    args = parser.parse_args()
    
    if not args.api_key and not args.mp_data:
        parser.error("Either --api-key or --mp-data must be provided")
    
    try:
        if args.self_hull:
            result_data = get_self_hull(
                file_path=args.input,
                mp_data_path=args.mp_data,
                use_compatibility=args.cor,
                element_dir=args.element_dir
            )
        else:
            result_data = get_above_hull_and_formation_energy(
                file_path=args.input, 
                mp_api_key=args.api_key,
                mp_data_path=args.mp_data,
                use_compatibility=args.cor
            )
        result_data.to_csv(args.output, index=False)
        
        if args.cor and 'status' in result_data.columns:
            success_count = sum(result_data['status'] == "Success")
            logging.info(f"Processing completed: {success_count}/{len(result_data)} compounds successful")
            print(f"All Done. Successfully processed {success_count}/{len(result_data)} compounds")
        else:
            logging.info("Processing completed successfully")
            print(f"All Done. Processed {len(result_data)} compounds")
    except Exception as e:
        logging.exception("Error")
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
