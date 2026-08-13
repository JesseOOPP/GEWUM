"""
GEWUM Base Structure Selector Module
Provides common functionality for diversity-based structure selection across all dimensions (0D-3D)
"""
import os
import numpy as np
from abc import ABC, abstractmethod
from pymatgen.core import Structure
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from tqdm import tqdm
import shutil
import random
from scipy.spatial.distance import pdist

try:
    from dscribe.descriptors import SOAP, CoulombMatrix
    DSCRIBE_AVAILABLE = True
except ImportError:
    DSCRIBE_AVAILABLE = False

try:
    import hdbscan
    HDBSCAN_AVAILABLE = True
except ImportError:
    HDBSCAN_AVAILABLE = False


class BaseStructureSelector(ABC):
    """
    Abstract base class for structure selectors.
    
    Provides common functionality:
    - File loading and featurization pipeline
    - Clustering algorithms (random, kmeans, medoid, maxmin)
    - Feature normalization
    - File operations (copy, move)
    
    Subclasses must implement:
    - extract_structure_features(): Dimension-specific feature extraction
    - get_descriptor_choices(): Available descriptor types for this dimension
    """
    
    def __init__(self, cif_directory, feature_settings=None):
        """
        Initialize the selector.
        
        Args:
            cif_directory: Directory containing CIF files (or formula dir with structures.db)
            feature_settings: Dict with feature extraction options
                - descriptor_type: Feature descriptor type
                - use_pca: Whether to apply PCA for high-dim features (default: True for soap/coulomb)
                - pca_variance: Variance ratio to preserve (default: 0.95)
                - pca_n_components: Fixed number of components (overrides pca_variance)
                - Additional dimension-specific options
        """
        self.cif_directory = cif_directory
        self.structures = {}
        self.features = {}
        self.file_list = []
        self.feature_settings = feature_settings or {}
        self.descriptor_type = self.feature_settings.get('descriptor_type', 'simple')
        
        self.use_pca = self.feature_settings.get('use_pca', None)  # Auto-detect
        self.pca_variance = self.feature_settings.get('pca_variance', 0.95)
        self.pca_n_components = self.feature_settings.get('pca_n_components', None)
        
        parent_dir = os.path.dirname(os.path.abspath(cif_directory))
        self._db_path = os.path.join(parent_dir, 'structures.db')
        self._use_db = os.path.isfile(self._db_path)
        if self._use_db:
            try:
                self._sg_number = int(os.path.basename(cif_directory))
            except ValueError:
                self._sg_number = None
                self._use_db = False
        else:
            self._sg_number = None
        
        self._validate_descriptor()
    
    def _validate_descriptor(self):
        """Validate that the requested descriptor is available."""
        if self.descriptor_type in ['soap', 'coulomb'] and not DSCRIBE_AVAILABLE:
            raise ImportError(
                f"{self.descriptor_type.upper()} descriptor requires dscribe package. "
                "Install with: pip install dscribe"
            )
    
    @property
    @abstractmethod
    def dimension_name(self) -> str:
        """Return the dimension name (e.g., '3D', '2D', '1D', '0D')."""
        pass
    
    @abstractmethod
    def get_descriptor_choices(self) -> list:
        """Return list of available descriptor types for this dimension."""
        pass
    
    @abstractmethod
    def extract_structure_features(self, structure) -> np.ndarray:
        """
        Extract numerical features from a structure.
        
        Args:
            structure: pymatgen Structure object
        
        Returns:
            numpy array of features
        """
        pass
    
    def load_and_featurize(self):
        """
        Load all CIF files and extract features.
        Applies PCA for high-dimensional features (SOAP, Coulomb) automatically.
        
        Supports reading from structures.db when present in the parent directory.
        
        Returns:
            Tuple of (file_list, normalized_features)
        """
        if self._use_db and self._sg_number is not None:
            return self._load_and_featurize_from_db()
        return self._load_and_featurize_from_fs()
    
    def _load_and_featurize_from_db(self):
        """Load CIF files from structures.db in the parent directory."""
        import sqlite3
        
        conn = sqlite3.connect(self._db_path)
        conn.row_factory = sqlite3.Row
        try:
            rows = conn.execute(
                "SELECT sg_number, cif_name, cif_content FROM structures WHERE sg_number = ? AND stage = 'initial' ORDER BY cif_name",
                (self._sg_number,)
            ).fetchall()
        finally:
            conn.close()
        
        if not rows:
            raise ValueError(f"No CIF files found in DB for SG={self._sg_number}: {self.cif_directory}")
        
        all_features = []
        valid_files = []
        
        desc_name = self.descriptor_type.upper()
        for row in tqdm(rows, desc=f"Featurizing {self.dimension_name} structures ({desc_name})"):
            try:
                cif_name = row['cif_name']
                structure = Structure.from_str(row['cif_content'], fmt='cif')
                
                features = self.extract_structure_features(structure)
                self.features[cif_name] = features
                all_features.append(features)
                valid_files.append(cif_name)
            except Exception as e:
                print(f"Error processing SG={row['sg_number']}/{row['cif_name']}: {e}")
        
        self.feature_matrix = np.array(all_features)
        self.scaler = StandardScaler()
        self.normalized_features = self.scaler.fit_transform(self.feature_matrix)
        self.file_list = valid_files
        
        original_dim = self.feature_matrix.shape[1]
        print(f"Successfully processed {len(self.file_list)} {self.dimension_name} structures (from DB)")
        print(f"Original feature dimension: {original_dim}")
        
        self._apply_pca_if_needed(original_dim)
        
        return self.file_list, self.normalized_features
    
    def _load_and_featurize_from_fs(self):
        """Load CIF files from filesystem (original behavior)."""
        cif_files = [f for f in os.listdir(self.cif_directory) if f.endswith('.cif')]
        
        if not cif_files:
            raise ValueError(f"No CIF files found in directory: {self.cif_directory}")
        
        all_features = []
        valid_files = []
        
        desc_name = self.descriptor_type.upper()
        for cif_file in tqdm(cif_files, desc=f"Featurizing {self.dimension_name} structures ({desc_name})"):
            try:
                file_path = os.path.join(self.cif_directory, cif_file)
                structure = Structure.from_file(file_path)
                
                features = self.extract_structure_features(structure)
                self.features[cif_file] = features
                all_features.append(features)
                valid_files.append(cif_file)
            except Exception as e:
                print(f"Error processing {cif_file}: {e}")
        
        self.feature_matrix = np.array(all_features)
        self.scaler = StandardScaler()
        self.normalized_features = self.scaler.fit_transform(self.feature_matrix)
        self.file_list = valid_files
        
        original_dim = self.feature_matrix.shape[1]
        print(f"Successfully processed {len(self.file_list)} {self.dimension_name} structures")
        print(f"Original feature dimension: {original_dim}")
        
        self._apply_pca_if_needed(original_dim)
        
        return self.file_list, self.normalized_features
    
    def _apply_pca_if_needed(self, original_dim):
        """
        Apply PCA dimensionality reduction for high-dimensional features.
        Auto-enabled for SOAP/Coulomb descriptors (>50 dimensions).
        """
        use_pca = self.use_pca
        if use_pca is None:
            use_pca = (self.descriptor_type in ['soap', 'coulomb'] and original_dim > 50)
        
        if not use_pca:
            return
        
        n_samples = len(self.file_list)
        if self.pca_n_components is not None:
            n_components = min(self.pca_n_components, n_samples, original_dim)
        else:
            n_components = self.pca_variance
        
        self.pca = PCA(n_components=n_components, random_state=42)
        self.normalized_features = self.pca.fit_transform(self.normalized_features)
        
        reduced_dim = self.normalized_features.shape[1]
        variance_explained = sum(self.pca.explained_variance_ratio_) * 100
        
        print(f"PCA: {original_dim}D -> {reduced_dim}D (preserved {variance_explained:.1f}% variance)")
    
    def random_selection(self, target_count=50):
        """Random sampling selection."""
        if len(self.file_list) <= target_count:
            return self.file_list

        random.seed(42 + os.getpid())
        return random.sample(self.file_list, target_count)
    
    def cluster_based_selection(self, target_count=50):
        """K-means clustering with random selection from each cluster."""
        if len(self.file_list) <= target_count:
            return self.file_list
        
        kmeans = KMeans(n_clusters=target_count, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(self.normalized_features)
        
        selected_files = []
        for cluster_id in range(target_count):
            cluster_indices = np.where(cluster_labels == cluster_id)[0]
            if len(cluster_indices) > 0:
                selected_idx = np.random.choice(cluster_indices)
                selected_files.append(self.file_list[selected_idx])
        
        return selected_files
    
    def medoid_based_selection(self, target_count=50):
        """K-means clustering with medoid (closest to centroid) selection."""
        if len(self.file_list) <= target_count:
            return self.file_list
        
        kmeans = KMeans(n_clusters=target_count, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(self.normalized_features)
        
        selected_files = []
        for cluster_id in range(target_count):
            cluster_indices = np.where(cluster_labels == cluster_id)[0]
            if len(cluster_indices) > 0:
                cluster_features = self.normalized_features[cluster_indices]
                centroid = kmeans.cluster_centers_[cluster_id]
                distances = np.linalg.norm(cluster_features - centroid, axis=1)
                medoid_idx = cluster_indices[np.argmin(distances)]
                selected_files.append(self.file_list[medoid_idx])
        
        return selected_files
    
    def maxmin_diversity_selection(self, target_count=50, batch_size=2000):
        """Maximum-minimum distance selection for maximum diversity."""
        if len(self.file_list) <= target_count:
            return self.file_list
        
        selected_indices = []
        remaining_indices = list(range(len(self.file_list)))
        
        first_idx = np.random.choice(remaining_indices)
        selected_indices.append(first_idx)
        remaining_indices.remove(first_idx)
        
        for i in range(min(target_count - 1, len(remaining_indices))):
            if i % 10 == 0:
                print(f"Selection progress: {i}/{target_count-1}")
                
            max_min_distance = -1
            best_candidate = None
            
            for j in range(0, len(remaining_indices), batch_size):
                batch_indices = remaining_indices[j:j+batch_size]
                batch_features = self.normalized_features[batch_indices]
                
                selected_features = self.normalized_features[selected_indices]
                distances = self._batch_distances(batch_features, selected_features)
                min_distances = np.min(distances, axis=1)
                
                batch_best_idx = np.argmax(min_distances)
                batch_best_distance = min_distances[batch_best_idx]
                
                if batch_best_distance > max_min_distance:
                    max_min_distance = batch_best_distance
                    best_candidate = batch_indices[batch_best_idx]
            
            if best_candidate is not None:
                selected_indices.append(best_candidate)
                remaining_indices.remove(best_candidate)
        
        return [self.file_list[i] for i in selected_indices]
    
    def hdbscan_selection(self, min_cluster_size=6, min_samples=3,
                         cluster_selection_epsilon=0.0, alpha=1.0,
                         keep_noise=True):
        """
        HDBSCAN clustering with automatic cluster count determination.
        Selects medoid (closest to cluster center) from each cluster.
        
        Args:
            min_cluster_size: Minimum size of clusters (default: 6)
            min_samples: Min samples for core points (default: 3)
            cluster_selection_epsilon: Distance threshold for merging clusters (default: 0.0)
                                      Larger values -> fewer, larger clusters
            alpha: Density decay parameter (default: 1.0)
                  Larger values -> stricter distance sensitivity
            keep_noise: Whether to retain noise points in the final selection (default: True)
                        True  -> keep all noise points (more structures preserved)
                        False -> discard noise points (only cluster medoids kept)
        
        Returns:
            Tuple of (selected_files, cluster_count, noise_count)
        """
        if not HDBSCAN_AVAILABLE:
            raise ImportError(
                "HDBSCAN method requires hdbscan package. "
                "Install with: pip install hdbscan"
            )
        
        if min_samples is None:
            min_samples = 3
        
        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            metric='euclidean',
            cluster_selection_method='eom',
            cluster_selection_epsilon=cluster_selection_epsilon,
            alpha=alpha
        )
        cluster_labels = clusterer.fit_predict(self.normalized_features)
        
        unique_labels = set(cluster_labels)
        noise_label = -1
        valid_clusters = [l for l in unique_labels if l != noise_label]
        
        noise_count = int(np.sum(cluster_labels == noise_label))
        cluster_count = len(valid_clusters)
        
        print(f"HDBSCAN found {cluster_count} clusters, {noise_count} noise points")
        
        if cluster_count == 0:
            print("Warning: No clusters found. All points are noise - keeping all.")
            return self.file_list, 0, noise_count
        
        selected_files = []
        for cluster_id in valid_clusters:
            cluster_indices = np.where(cluster_labels == cluster_id)[0]
            cluster_features = self.normalized_features[cluster_indices]
            
            centroid = np.mean(cluster_features, axis=0)
            distances = np.linalg.norm(cluster_features - centroid, axis=1)
            medoid_idx = cluster_indices[np.argmin(distances)]
            selected_files.append(self.file_list[medoid_idx])
        
        noise_indices = np.where(cluster_labels == noise_label)[0]
        noise_files = [self.file_list[i] for i in noise_indices]
        
        if keep_noise:
            selected_files.extend(noise_files)
            print(f"Selected: {len(valid_clusters)} cluster medoids + {noise_count} noise points = {len(selected_files)} total")
        else:
            print(f"Selected: {len(valid_clusters)} cluster medoids (noise discarded: {noise_count} points)")
        
        return selected_files, cluster_count, noise_count
    
    def _batch_distances(self, features1, features2):
        """Calculate pairwise distances between two feature matrices."""
        return np.sqrt(((features1[:, np.newaxis, :] - features2[np.newaxis, :, :]) ** 2).sum(axis=2))
    
    def select(self, method='medoid', target_count=50, **kwargs):
        """
        Select structures using the specified method.
        
        Args:
            method: Selection method ('random', 'kmeans', 'medoid', 'maxmin', 'hdbscan')
            target_count: Number of structures to select (ignored for hdbscan)
            **kwargs: Additional method-specific parameters
                - hdbscan: min_cluster_size (default: 5), min_samples
        
        Returns:
            List of selected file names (for hdbscan, returns tuple with cluster info)
        """
        if method == 'random':
            return self.random_selection(target_count=target_count)
        elif method == 'kmeans':
            return self.cluster_based_selection(target_count=target_count)
        elif method == 'maxmin':
            return self.maxmin_diversity_selection(target_count=target_count)
        elif method == 'hdbscan':
            min_cluster_size = kwargs.get('min_cluster_size', 6)
            min_samples = kwargs.get('min_samples', 3)
            cluster_selection_epsilon = kwargs.get('cluster_selection_epsilon', 0.0)
            alpha = kwargs.get('alpha', 1.0)
            keep_noise = kwargs.get('keep_noise', True)
            return self.hdbscan_selection(
                min_cluster_size=min_cluster_size,
                min_samples=min_samples,
                cluster_selection_epsilon=cluster_selection_epsilon,
                alpha=alpha,
                keep_noise=keep_noise
            )
        else:  
            return self.medoid_based_selection(target_count=target_count)
    
    def copy_selected_structures(self, selected_files, output_dir):
        """Copy selected structures to output directory."""
        os.makedirs(output_dir, exist_ok=True)
        for file in selected_files:
            src_path = os.path.join(self.cif_directory, file)
            dst_path = os.path.join(output_dir, file)
            shutil.copy2(src_path, dst_path)
    
    def _extract_soap_features(self, structure, periodic=True, r_cut=6.0):
        """
        Extract SOAP descriptor features.
        
        Args:
            structure: pymatgen Structure object
            periodic: Whether to use periodic boundary conditions
            r_cut: SOAP cutoff radius
        
        Returns:
            numpy array of SOAP features
        """
        atoms = structure.to_ase_atoms()
        if not periodic:
            atoms.set_pbc([False, False, False])
        
        # sorted() so the SOAP feature-column order is deterministic across
        # processes/frames (set() iteration order is not stable).
        species = sorted(set(atoms.get_chemical_symbols()))
        
        n_max = self.feature_settings.get('soap_n_max', 4)
        l_max = self.feature_settings.get('soap_l_max', 4)
        r_cut = self.feature_settings.get('soap_r_cut', r_cut)
        
        soap = SOAP(
            species=species,
            r_cut=r_cut,
            n_max=n_max,
            l_max=l_max,
            periodic=periodic,
            average="inner",
            sparse=False
        )
        
        descriptor = soap.create(atoms)
        return descriptor.flatten()
    
    def _extract_coulomb_features(self, structure):
        """
        Extract Coulomb Matrix features (for 0D structures).
        
        Args:
            structure: pymatgen Structure object
        
        Returns:
            numpy array of Coulomb Matrix features
        """
        atoms = structure.to_ase_atoms()
        atoms.set_pbc([False, False, False])
        
        n_atoms_max = self.feature_settings.get('coulomb_n_atoms_max', 100)
        if len(atoms) > n_atoms_max:
            n_atoms_max = len(atoms)
        
        cm = CoulombMatrix(
            n_atoms_max=n_atoms_max,
            permutation="sorted_l2"
        )
        
        descriptor = cm.create(atoms)
        return descriptor.flatten()
