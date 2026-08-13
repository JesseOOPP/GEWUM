"""
GEWUM MP Offline Data Loader
Load Materials Project data from offline JSON file (MPtrj format)
Generates a compact index file for fast subsequent loading
"""
import json
import os
import logging
from collections import defaultdict
from typing import Dict, List, Optional, Tuple
from pathlib import Path

try:
    import ijson
    IJSON_AVAILABLE = True
except ImportError:
    IJSON_AVAILABLE = False

from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PDEntry


logger = logging.getLogger(__name__)


def get_index_path(mp_data_path: str) -> str:
    """Get the index file path for a given MP data file"""
    base_path = Path(mp_data_path)
    return str(base_path.parent / f"{base_path.stem}_index.json")


def extract_composition_from_structure(structure_dict: dict) -> str:
    """
    Extract chemical formula from structure dict.
    
    Args:
        structure_dict: Structure dictionary containing 'sites' field
        
    Returns:
        Chemical formula string (e.g., 'Fe2O3')
    """
    element_counts = defaultdict(int)
    sites = structure_dict.get('sites', [])
    
    for site in sites:
        species_list = site.get('species', [])
        for species in species_list:
            element = species.get('element', '')
            occu = float(species.get('occu', 1))
            if element:
                element_counts[element] += occu
    
    formula_parts = []
    for elem in sorted(element_counts.keys()):
        count = element_counts[elem]
        if count == int(count):
            count = int(count)
        if count == 1:
            formula_parts.append(elem)
        else:
            formula_parts.append(f"{elem}{count}")
    
    return ''.join(formula_parts)


def get_chemsys_from_composition(composition: str) -> str:
    """
    Get chemical system string from composition.
    
    Args:
        composition: Chemical formula string
        
    Returns:
        Sorted element system string (e.g., 'Fe-O')
    """
    try:
        comp = Composition(composition)
        elements = sorted([el.symbol for el in comp.elements])
        return '-'.join(elements)
    except Exception:
        return ''


def build_index_streaming(mp_data_path: str, index_path: str,
                         store_structures: bool = False) -> Dict:
    """
    Build index from large MP JSON file using streaming parser.
    Selects lowest energy configuration for each mp_id.
    
    Args:
        mp_data_path: Path to MPtrj JSON file
        index_path: Path to save index file
        store_structures: If True, embed structure dicts in the index
            (enables instant structure retrieval, but makes index larger)
        
    Returns:
        Index dictionary {mp_id: {composition, energy, chemsys, [structure]}}
    """
    if not IJSON_AVAILABLE:
        raise ImportError("ijson is required for streaming large JSON files. "
                         "Install with: pip install ijson")
    
    logger.info(f"Building index from {mp_data_path} (streaming mode)...")
    logger.info("This may take several minutes for large files...")
    
    index = {}
    processed_count = 0
    
    with open(mp_data_path, 'rb') as f:
        parser = ijson.kvitems(f, '')
        
        for mp_id, sub_entries in parser:
            if not isinstance(sub_entries, dict):
                continue
            
            best_entry = None
            best_energy = float('inf')
            
            for sub_id, entry_data in sub_entries.items():
                if not isinstance(entry_data, dict):
                    continue
                
                energy_per_atom = entry_data.get('energy_per_atom')
                if energy_per_atom is None:
                    continue
                
                try:
                    energy_val = float(energy_per_atom)
                except (TypeError, ValueError):
                    continue
                
                if energy_val < best_energy:
                    best_energy = energy_val
                    
                    structure = entry_data.get('structure', {})
                    composition = extract_composition_from_structure(structure)
                    corrected_energy = entry_data.get('corrected_total_energy')
                    uncorrected_energy = entry_data.get('uncorrected_total_energy')
                    
                    n_atoms = len(structure.get('sites', []))
                    
                    # Resolve elements from composition for potcar_symbols reconstruction
                    try:
                        elements = sorted([el.symbol for el in Composition(composition).elements])
                    except Exception:
                        elements = []
                    
                    # total_energy: prefer corrected for backward compatibility
                    total_energy = corrected_energy if corrected_energy is not None else uncorrected_energy
                    
                    best_entry = {
                        'composition': composition,
                        'energy_per_atom': energy_val,
                        'total_energy': float(total_energy) if total_energy else energy_val * n_atoms,
                        'uncorrected_total_energy': float(uncorrected_energy) if uncorrected_energy is not None else None,
                        'n_atoms': n_atoms,
                        'chemsys': get_chemsys_from_composition(composition),
                        'elements': elements,
                        'entry_id': mp_id
                    }
                    
                    if store_structures:
                        best_entry['structure'] = structure
            
            if best_entry and best_entry['composition']:
                index[mp_id] = best_entry
                processed_count += 1
                
                if processed_count % 10000 == 0:
                    logger.info(f"Processed {processed_count} materials...")
    
    logger.info(f"Index building complete. Total materials: {len(index)}")
    if store_structures:
        logger.info("Structures embedded in index (fast subsequent loads)")
    
    logger.info(f"Saving index to {index_path}...")
    with open(index_path, 'w', encoding='utf-8') as f:
        json.dump(index, f, ensure_ascii=False)
    logger.info("Index saved successfully.")
    
    return index


def load_index(index_path: str) -> Dict:
    """Load pre-built index from file"""
    logger.info(f"Loading index from {index_path}...")
    with open(index_path, 'r', encoding='utf-8') as f:
        index = json.load(f)
    logger.info(f"Loaded index with {len(index)} materials.")
    return index


class MPOfflineLoader:
    """
    Offline loader for Materials Project data.
    Manages index creation and provides query interface.
    """
    
    def __init__(self, mp_data_path: str, force_rebuild: bool = False,
                 store_structures: bool = False):
        """
        Initialize offline loader.
        
        Args:
            mp_data_path: Path to MPtrj JSON file
            force_rebuild: Force rebuild index even if exists
            store_structures: If True, embed structure dicts in index
                (one-time overhead, enables instant structure retrieval)
        """
        self.mp_data_path = mp_data_path
        self.index_path = get_index_path(mp_data_path)
        self.index = None
        self.chemsys_cache = None
        self._has_structures = store_structures
        
        self._load_or_build_index(force_rebuild, store_structures)
        self._build_chemsys_cache()
    
    def _load_or_build_index(self, force_rebuild: bool,
                             store_structures: bool = False):
        """Load existing index or build new one"""
        # Rebuild if forced OR if we want structures but existing index lacks them
        need_rebuild = force_rebuild
        if store_structures and not force_rebuild and os.path.exists(self.index_path):
            # Quick check: does the first entry have a structure key?
            try:
                with open(self.index_path, 'r', encoding='utf-8') as f:
                    peek = json.load(f)
                first_entry = next(iter(peek.values()), {}) if peek else {}
                if 'structure' not in first_entry:
                    logger.info("Existing index lacks structures; rebuilding...")
                    need_rebuild = True
                else:
                    self._has_structures = True
            except Exception:
                need_rebuild = True
        
        if need_rebuild or not os.path.exists(self.index_path):
            self.index = build_index_streaming(
                self.mp_data_path, self.index_path,
                store_structures=store_structures or self._has_structures
            )
        else:
            self.index = load_index(self.index_path)
    
    def _build_chemsys_cache(self):
        """Build chemsys -> mp_ids mapping for fast lookup"""
        logger.info("Building chemical system cache...")
        self.chemsys_cache = defaultdict(list)
        
        for mp_id, entry_data in self.index.items():
            chemsys = entry_data.get('chemsys', '')
            if chemsys:
                self.chemsys_cache[chemsys].append(mp_id)
        
        logger.info(f"Cache built with {len(self.chemsys_cache)} chemical systems.")
    
    def get_entries_in_chemsys(self, elements: List[str], 
                               as_pd_entry: bool = False,
                               use_uncorrected: bool = False) -> List:
        """
        Get all entries in a chemical system (including subsystems).
        Keeps the lowest energy_per_atom entry for each composition to ensure
        ground-state phases are used for convex hull construction.
        
        Args:
            elements: List of element symbols (e.g., ['Fe', 'O'])
            as_pd_entry: If True, return PDEntry objects; else ComputedEntry
            use_uncorrected: If True, use uncorrected_total_energy instead of
                             total_energy (corrected). When combined with
                             as_pd_entry=False, ComputedEntry are returned with
                             uncorrected energy ready for unified
                             MP2020Compatibility reprocessing.
            
        Returns:
            List of entry objects for phase diagram construction
        """
        from itertools import combinations
        
        elements = sorted(elements)
        target_systems = set()
        
        for i in range(1, len(elements) + 1):
            for combo in combinations(elements, i):
                target_systems.add('-'.join(sorted(combo)))
        
        # Collect the lowest energy_per_atom entry for each composition
        best_per_composition: Dict[str, Dict] = {}
        
        for chemsys in target_systems:
            mp_ids = self.chemsys_cache.get(chemsys, [])
            
            for mp_id in mp_ids:
                entry_data = self.index.get(mp_id, {})
                composition = entry_data.get('composition', '')
                if not composition:
                    continue
                
                epa = entry_data.get('energy_per_atom', float('inf'))
                
                if composition not in best_per_composition or epa < best_per_composition[composition].get('energy_per_atom', float('inf')):
                    best_per_composition[composition] = entry_data
        
        # Create entry objects from best entries
        entries = []
        for composition, entry_data in best_per_composition.items():
            total_energy = entry_data.get('total_energy', 0)
            uncorrected_energy = entry_data.get('uncorrected_total_energy')
            entry_id = entry_data.get('entry_id', '')
            entry_elements = entry_data.get('elements', [])
            
            # Choose which energy to use
            if use_uncorrected and uncorrected_energy is not None:
                energy = uncorrected_energy
            else:
                energy = total_energy
            
            try:
                if as_pd_entry:
                    entry = PDEntry(
                        composition=composition,
                        energy=energy,
                        name=entry_id
                    )
                else:
                    entry = ComputedEntry(
                        composition=composition,
                        energy=energy,
                        entry_id=entry_id
                    )
                    # Attach elements list for caller to reconstruct potcar_symbols
                    if use_uncorrected and entry_elements:
                        entry.data['_mp_elements'] = entry_elements
                entries.append(entry)
            except Exception as e:
                logger.warning(f"Failed to create entry for {entry_id}: {e}")
        
        logger.info(f"Found {len(entries)} entries for system {'-'.join(elements)}")
        return entries
    
    def get_entry_by_mp_id(self, mp_id: str) -> Optional[Dict]:
        """Get entry data by mp_id"""
        return self.index.get(mp_id)
    
    def get_all_chemsys(self) -> List[str]:
        """Get all available chemical systems"""
        return list(self.chemsys_cache.keys())
    
    def __len__(self):
        return len(self.index) if self.index else 0
