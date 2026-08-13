"""
GEWUM 0D Structure Selector
Feature extraction optimized for 0D structures (molecules, clusters, nanoparticles)
"""
import numpy as np
from scipy.spatial.distance import pdist
try:
    from .base_selector import BaseStructureSelector
except ImportError:
    from base_selector import BaseStructureSelector


class StructureSelector0D(BaseStructureSelector):
    """
    Structure selector for 0D structures (molecules, clusters, nanoparticles).
    
    Features optimized for 0D non-periodic structures:
    - Geometric features: atom count, radius, diameter, shape parameters
    - SOAP descriptor with periodic=False
    - Coulomb Matrix for molecular fingerprinting
    
    Supported descriptors:
    - simple: Handcrafted geometric features (fast, no dependencies)
    - soap: SOAP descriptor (non-periodic mode)
    - coulomb: Coulomb Matrix (fast, designed for molecules)
    """
    
    @property
    def dimension_name(self) -> str:
        return "0D"
    
    def get_descriptor_choices(self) -> list:
        return ['simple', 'soap', 'coulomb']
    
    def extract_structure_features(self, structure) -> np.ndarray:
        """Extract features based on descriptor type."""
        if self.descriptor_type == 'soap':
            r_cut = self.feature_settings.get('soap_r_cut', 5.0)
            return self._extract_soap_features(structure, periodic=False, r_cut=r_cut)
        elif self.descriptor_type == 'coulomb':
            return self._extract_coulomb_features(structure)
        else:
            return self._extract_simple_features(structure)
    
    def _extract_simple_features(self, structure) -> np.ndarray:
        """Extract handcrafted geometric features for 0D structures."""
        features = []
        cart_coords = structure.cart_coords
        num_atoms = len(structure)
        
        features.append(num_atoms)
        
        centroid = np.mean(cart_coords, axis=0)
        distances_to_center = np.linalg.norm(cart_coords - centroid, axis=1)
        
        features.append(np.max(distances_to_center))    
        features.append(np.mean(distances_to_center))   
        features.append(np.std(distances_to_center))    
        
        if num_atoms > 1:
            pairwise_distances = pdist(cart_coords)
            features.append(np.max(pairwise_distances))    
            features.append(np.mean(pairwise_distances))   
            features.append(np.min(pairwise_distances))    
            features.append(np.std(pairwise_distances))   
            
            features.append(np.percentile(pairwise_distances, 25))
            features.append(np.percentile(pairwise_distances, 75))
        else:
            features.extend([0.0] * 6)
        
        if num_atoms > 1:
            centered_coords = cart_coords - centroid
            gyration_radius = np.sqrt(np.mean(np.sum(centered_coords**2, axis=1)))
            features.append(gyration_radius)
            
            inertia_tensor = np.dot(centered_coords.T, centered_coords) / num_atoms
            eigenvalues = np.linalg.eigvalsh(inertia_tensor)
            eigenvalues = np.sort(eigenvalues)[::-1]  
            
            if eigenvalues[0] > 0:
                features.append(eigenvalues[1] / eigenvalues[0])  
                features.append(eigenvalues[2] / eigenvalues[0])  
            else:
                features.extend([1.0, 1.0])
        else:
            features.extend([0.0, 1.0, 1.0])
        
        return np.array(features)
