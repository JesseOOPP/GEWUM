"""
GEWUM 1D Structure Selector
Feature extraction optimized for 1D structures (chains, nanowires, nanotubes)
Note: 1D structures require SOAP descriptor due to vacuum layers in non-periodic directions
"""
import numpy as np
from scipy.spatial.distance import pdist
try:
    from .base_selector import BaseStructureSelector, DSCRIBE_AVAILABLE
except ImportError:
    from base_selector import BaseStructureSelector, DSCRIBE_AVAILABLE


class StructureSelector1D(BaseStructureSelector):
    """
    Structure selector for 1D structures (chains, nanowires, nanotubes).
    
    Features optimized for 1D structures:
    - SOAP descriptor with small r_cut to avoid vacuum interaction
    - Basic geometric features combined with SOAP
    
    Supported descriptors:
    - soap: SOAP descriptor (required for 1D, simple not supported)
    
    Note: Simple features alone are insufficient for 1D structures
          due to vacuum layers in non-periodic directions.
    """
    
    def __init__(self, cif_directory, feature_settings=None):
        """Initialize with SOAP as default descriptor."""
        feature_settings = feature_settings or {}
        feature_settings['descriptor_type'] = 'soap'
        
        if not DSCRIBE_AVAILABLE:
            raise ImportError(
                "1D structure selection requires SOAP descriptor. "
                "Install dscribe with: pip install dscribe"
            )
        
        super().__init__(cif_directory, feature_settings)
        
        self.periodic_direction = feature_settings.get('periodic_direction', 'z')
    
    @property
    def dimension_name(self) -> str:
        return "1D"
    
    def get_descriptor_choices(self) -> list:
        return ['soap']  # Only SOAP supported for 1D
    
    def _get_periodic_axis(self):
        """Get the index of the periodic direction."""
        direction_map = {'x': 0, 'y': 1, 'z': 2}
        return direction_map.get(self.periodic_direction, 2)
    
    def _get_non_periodic_axes(self):
        """Get indices of non-periodic directions (vacuum directions)."""
        axis = self._get_periodic_axis()
        return [i for i in range(3) if i != axis]
    
    def extract_structure_features(self, structure) -> np.ndarray:
        """Extract SOAP features with optional basic geometric features."""
        features = []
        
        if self.feature_settings.get('include_basic_features', True):
            basic_features = self._extract_basic_1d_features(structure)
            features.extend(basic_features)
        
        soap_features = self._extract_1d_soap_features(structure)
        features.extend(soap_features)
        
        return np.array(features)
    
    def _extract_basic_1d_features(self, structure) -> list:
        """Extract basic geometric features for 1D structures."""
        features = []
        lattice = structure.lattice
        
        p_axis = self._get_periodic_axis()
        np_axes = self._get_non_periodic_axes()
        
        abc = [lattice.a, lattice.b, lattice.c]
        
        chain_length = abc[p_axis]
        features.append(chain_length)
        
        num_atoms = len(structure)
        linear_density = num_atoms / chain_length
        features.append(linear_density)
        features.append(num_atoms)
        
        cross_section = abc[np_axes[0]] * abc[np_axes[1]]
        features.append(cross_section)
        
        all_coords = np.array([site.frac_coords for site in structure])
        chain_coords = all_coords[:, p_axis]
        
        features.append(np.mean(chain_coords))
        features.append(np.std(chain_coords))
        
        cross_coords = all_coords[:, np_axes]
        if len(cross_coords) > 1:
            cart_coords = structure.cart_coords[:, np_axes]
            max_spread = np.max(pdist(cart_coords, metric='euclidean'))
            features.append(max_spread)
        else:
            features.append(0.0)
        
        return features
    
    def _extract_1d_soap_features(self, structure) -> np.ndarray:
        """Extract SOAP features with r_cut constraint for 1D."""
        lattice = structure.lattice
        abc = [lattice.a, lattice.b, lattice.c]
        np_axes = self._get_non_periodic_axes()
        
        r_cut = self.feature_settings.get('soap_r_cut', 4.0)
        min_vacuum = min(abc[np_axes[0]], abc[np_axes[1]])
        if r_cut > min_vacuum / 2:
            r_cut = min_vacuum / 2 - 0.5
            if r_cut < 2.0:
                r_cut = 2.0
        
        return self._extract_soap_features(structure, periodic=True, r_cut=r_cut)
