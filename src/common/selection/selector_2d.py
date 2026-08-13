"""
GEWUM 2D Structure Selector
Feature extraction optimized for 2D layered materials
"""
import numpy as np
from scipy.spatial.distance import pdist
try:
    from .base_selector import BaseStructureSelector
except ImportError:
    from base_selector import BaseStructureSelector


class StructureSelector2D(BaseStructureSelector):
    """
    Structure selector for 2D layered materials.
    
    Features optimized for 2D structures:
    - Area (a * b) instead of volume
    - 2D density (atoms per area)
    - Z-direction thickness
    - XY plane coordinate statistics
    
    Supported descriptors:
    - simple: Handcrafted 2D features (fast)
    - soap: SOAP descriptor with smaller r_cut (accurate)
    """
    
    @property
    def dimension_name(self) -> str:
        return "2D"
    
    def get_descriptor_choices(self) -> list:
        return ['simple', 'soap']
    
    def extract_structure_features(self, structure) -> np.ndarray:
        """Extract features based on descriptor type."""
        if self.descriptor_type == 'soap':
            r_cut = self.feature_settings.get('soap_r_cut', 5.0)  
            return self._extract_soap_features(structure, periodic=True, r_cut=r_cut)
        else:
            return self._extract_simple_features(structure)
    
    def _extract_simple_features(self, structure) -> np.ndarray:
        """Extract handcrafted features optimized for 2D structures."""
        features = []
        lattice = structure.lattice
        
        area = lattice.a * lattice.b
        features.append(area)
        
        num_atoms = len(structure)
        features.append(num_atoms / area)  
        features.append(num_atoms)
        
        if self.feature_settings.get('include_lattice_params', True):
            features.extend([lattice.a, lattice.b])
            features.append(lattice.a / lattice.b)  
        
        if self.feature_settings.get('include_position_stats', True):
            all_coords = np.array([site.frac_coords for site in structure])
            
            xy_coords = all_coords[:, :2]
            
            features.extend(np.mean(xy_coords, axis=0))  
            features.extend(np.std(xy_coords, axis=0))  
            
            if all_coords.shape[0] > 0:
                z_coords = all_coords[:, 2]
                z_thickness = (np.max(z_coords) - np.min(z_coords)) * lattice.c
                features.append(z_thickness)
            else:
                features.append(0.0)
                
            if len(xy_coords) > 1:
                cart_coords = structure.cart_coords[:, :2]  
                max_dist = np.max(pdist(cart_coords, metric='euclidean'))
                features.append(max_dist)
                
                avg_dist = np.mean(pdist(cart_coords, metric='euclidean'))
                features.append(avg_dist)
            else:
                features.extend([0.0, 0.0])

        return np.array(features)
