"""
GEWUM 3D Structure Selector
Feature extraction optimized for 3D periodic crystals
"""
import numpy as np
from scipy.spatial.distance import pdist
try:
    from .base_selector import BaseStructureSelector
except ImportError:
    from base_selector import BaseStructureSelector


class StructureSelector3D(BaseStructureSelector):
    """
    Structure selector for 3D periodic crystals.
    
    Features optimized for 3D structures:
    - Volume, density, atom count
    - Lattice parameters (a, b, c)
    - Fractional coordinate statistics
    
    Supported descriptors:
    - simple: Handcrafted features (fast)
    - soap: SOAP descriptor (accurate)
    """
    
    @property
    def dimension_name(self) -> str:
        return "3D"
    
    def get_descriptor_choices(self) -> list:
        return ['simple', 'soap']
    
    def extract_structure_features(self, structure) -> np.ndarray:
        """Extract features based on descriptor type."""
        if self.descriptor_type == 'soap':
            r_cut = self.feature_settings.get('soap_r_cut', 6.0)
            return self._extract_soap_features(structure, periodic=True, r_cut=r_cut)
        else:
            return self._extract_simple_features(structure)
    
    def _extract_simple_features(self, structure) -> np.ndarray:
        """Extract handcrafted features for 3D crystals."""
        features = []
        
        features.append(structure.volume)
        features.append(structure.density)
        features.append(len(structure))
        
        if self.feature_settings.get('include_lattice_params', True):
            lattice = structure.lattice
            features.extend([lattice.a, lattice.b, lattice.c])

        if self.feature_settings.get('include_position_stats', True):
            all_coords = np.array([site.frac_coords for site in structure])
            
            features.extend(np.mean(all_coords, axis=0))  
            features.extend(np.std(all_coords, axis=0))   
            
            if len(all_coords) > 1:
                max_dist = np.max(pdist(all_coords, metric='euclidean'))
                features.append(max_dist)
            else:
                features.append(0.0)

        return np.array(features)
