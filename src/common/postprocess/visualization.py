"""
GEWUM Structure Visualization Module
Extracts features from CIF files using existing selector infrastructure,
performs UMAP/t-SNE dimensionality reduction, and generates comparison plots.

Usage:
    python visualization.py --dim 3d --descriptor simple --total-dir ./total_cifs [options]
"""
import os
import sys
import re
import ast
import glob
import random
import argparse
import datetime
from pathlib import Path
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

from gewum.src.common.cif_archive import (
    find_cifs_db_aware, detect_stages_db_aware,
    load_structure, entry_to_path, entry_basename,
)

try:
    import umap
    UMAP_AVAILABLE = True
except ImportError:
    UMAP_AVAILABLE = False

try:
    from dscribe.descriptors import SOAP, CoulombMatrix
    DSCRIBE_AVAILABLE = True
except ImportError:
    DSCRIBE_AVAILABLE = False

def formula_to_mathtext(formula):
    """
    Convert a chemical formula string to matplotlib mathtext with subscripts.

    Example: 'Fe10Cr10Al1' -> '$\\mathrm{Fe_{10}Cr_{10}Al_{1}}$'
    """
    tokens = re.findall(r'([A-Z][a-z]?)(\d+)', formula)
    if not tokens:
        return formula
    inner = "".join(f"{elem}_{{{num}}}" for elem, num in tokens)
    return rf"$\mathrm{{{inner}}}$"


_ANGLES = [0, 45, -45, 90, 30, -30]
_ANGLE_CUM_WEIGHTS = []
_s = 0
for _w in [60, 15, 15, 5, 2.5, 2.5]:
    _s += _w
    _ANGLE_CUM_WEIGHTS.append(_s)


def _pick_angle(rng):
    """Return a random rotation angle using weighted selection."""
    r = rng.random() * _ANGLE_CUM_WEIGHTS[-1]
    for i, c in enumerate(_ANGLE_CUM_WEIGHTS):
        if r <= c:
            return _ANGLES[i]
    return _ANGLES[-1]


def parse_cifgen(filepath):
    """
    Parse cifgen.inp file and return a formula frequency dictionary.
    Each line format: ['Fe', 'Cr', 'Al'],(10, 10, 1),200
    Also returns a list of (formula, total_atoms) for sizing.
    """
    formula_counter = Counter()
    formula_atoms = {}
    total_lines = 0

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            total_lines += 1

            try:
                match = re.match(r"(\[.*?\])\s*,\s*\((.*?)\)\s*,\s*(\d+)", line)
                if not match:
                    continue

                elements_str = match.group(1)
                ratios_str = match.group(2)

                elements = ast.literal_eval(elements_str)
                ratios = [int(x.strip()) for x in ratios_str.split(",")]

                formula = "".join(f"{el}{r}" for el, r in zip(elements, ratios))
                formula_counter[formula] += 1
                formula_atoms[formula] = sum(ratios)

            except Exception as e:
                print(f"Warning: skipping unparseable line: {line}  ({e})")
                continue

    return formula_counter, formula_atoms, total_lines


def _generate_srss_mask(fig_w=16, fig_h=10, res_w=320, res_h=200):
    """
    Render the text 'SRSS' into a binary mask using matplotlib.

    Returns
    -------
    mask : np.ndarray of bool, shape (res_h, res_w)
        True where the letters are drawn.
    """
    dpi = 600
    fig_tmp = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
    ax_tmp = fig_tmp.add_axes([0, 0, 1, 1])
    ax_tmp.set_xlim(0, 1)
    ax_tmp.set_ylim(0, 1)
    ax_tmp.axis("off")
    fig_tmp.patch.set_facecolor("white")

    ax_tmp.text(
        0.5, 0.5, "SRSS",
        fontsize=430, fontweight="bold", fontfamily="Liberation Serif",
        ha="center", va="center", color="black",
    )

    fig_tmp.canvas.draw()
    buf = np.frombuffer(fig_tmp.canvas.buffer_rgba(), dtype=np.uint8)
    w_px, h_px = fig_tmp.canvas.get_width_height()
    buf = buf.reshape(h_px, w_px, 4)
    plt.close(fig_tmp)

    gray = buf[:, :, 0].astype(np.float32) * 0.299 \
         + buf[:, :, 1].astype(np.float32) * 0.587 \
         + buf[:, :, 2].astype(np.float32) * 0.114
    mask_full = gray < 128  

    bh = h_px // res_h
    bw = w_px // res_w
    cropped = mask_full[:bh * res_h, :bw * res_w]
    blocks = cropped.reshape(res_h, bh, res_w, bw)
    mask = blocks.any(axis=(1, 3))

    return mask


def _filter_min_distance(points, min_dist):
    """
    Greedily filter a list of (x, y) points so that no two accepted
    points are closer than *min_dist* (in normalised coordinates).
    """
    accepted = []
    for pt in points:
        too_close = False
        for ap in accepted:
            if (pt[0] - ap[0]) ** 2 + (pt[1] - ap[1]) ** 2 < min_dist ** 2:
                too_close = True
                break
        if not too_close:
            accepted.append(pt)
    return accepted


def generate_wordcloud(formula_counter, formula_atoms, output_path,
                       n_display=300, seed=42, cmap_name=None, dpi=150):
    """
    Generate a word cloud shaped as 'SRSS' using pure matplotlib.

    Parameters
    ----------
    formula_counter : Counter
        Mapping of formula string to frequency.
    formula_atoms : dict
        Mapping of formula string to total atom count (used for sizing).
    output_path : str
        File path for the saved PNG.
    n_display : int
        Target number of words to show (may be reduced by mask capacity).
    seed : int
        Random seed for reproducibility.
    cmap_name : str or None
        Colormap name (e.g. 'viridis'). If None, uses default pink-blue gradient.
    dpi : int
        Output DPI for the saved figure.
    """
    rng = random.Random(seed)
    np_rng = np.random.RandomState(seed)

    # --- Build the SRSS mask ------------------------------------------------
    fig_w, fig_h = 16, 10
    mask = _generate_srss_mask(fig_w=fig_w, fig_h=fig_h)
    ys, xs = np.where(mask)
    if len(ys) == 0:
        raise RuntimeError("SRSS mask is empty; check font availability.")

    # Normalise pixel coords to [0, 1]
    norm_x = xs / mask.shape[1]
    norm_y = 1.0 - ys / mask.shape[0]  

    order = np_rng.permutation(len(norm_x))
    candidates = [(float(norm_x[i]), float(norm_y[i])) for i in order]

    min_dist = 0.018
    candidates = _filter_min_distance(candidates, min_dist)
    print(f"Mask pixels: {len(ys)}, candidates after spacing filter: {len(candidates)}")

    formulas = list(formula_counter.keys())
    if len(formulas) > n_display:
        formulas = rng.sample(formulas, n_display)
    else:
        n_display = len(formulas)

    if len(candidates) < n_display:
        n_display = len(candidates)
        formulas = formulas[:n_display]
    else:
        candidates = candidates[:n_display]

    while len(formulas) < len(candidates):
        formulas.append(formulas[len(formulas) % len(formula_counter)])
    formulas = formulas[:len(candidates)]

    atom_counts = [formula_atoms.get(f, 10) for f in formulas]
    min_atoms = min(atom_counts) if atom_counts else 1
    max_atoms = max(atom_counts) if atom_counts else 1
    span = max_atoms - min_atoms if max_atoms != min_atoms else 1.0

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    if cmap_name is not None:
        cmap = plt.get_cmap(cmap_name)
    else:
        cmap = LinearSegmentedColormap.from_list("pink_blue", ["#E8879B", "#5BA8C8"])

    for idx, formula in enumerate(formulas):
        cx, cy = candidates[idx]

        normed = (formula_atoms.get(formula, 10) - min_atoms) / span
        fontsize = 6 + normed * 8 + rng.uniform(-1, 1)
        fontsize = max(6, min(fontsize, 14))

        colour = cmap(rng.random())
        rotation = _pick_angle(rng)

        ax.text(
            cx, cy, formula_to_mathtext(formula),
            fontsize=fontsize,
            color=colour,
            ha="center", va="center",
            rotation=rotation,
            fontweight="bold",
        )

    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight", pad_inches=0.1, facecolor="white")
    plt.close(fig)
    print(f"[Visualize] Word cloud saved to: {output_path}")

def _get_selector_class(dim):
    """Get the appropriate selector class for the given dimension."""
    try:
        from gewum.src.common.selection.selector_3d import StructureSelector3D
        from gewum.src.common.selection.selector_2d import StructureSelector2D
        from gewum.src.common.selection.selector_1d import StructureSelector1D
        from gewum.src.common.selection.selector_0d import StructureSelector0D
    except ImportError:
        sel_dir = os.path.join(os.path.dirname(__file__), '..', 'selection')
        sel_dir = os.path.abspath(sel_dir)
        if sel_dir not in sys.path:
            sys.path.insert(0, sel_dir)
        from selector_3d import StructureSelector3D
        from selector_2d import StructureSelector2D
        from selector_1d import StructureSelector1D
        from selector_0d import StructureSelector0D

    selector_map = {
        '3d': StructureSelector3D,
        '2d': StructureSelector2D,
        '1d': StructureSelector1D,
        '0d': StructureSelector0D,
    }
    if dim not in selector_map:
        raise ValueError(f"Unknown dimension: {dim}. Choose from {list(selector_map.keys())}")
    return selector_map[dim]


def _find_cif_files(directory, mode=None):
    """
    Find CIF files in directory with smart filtering based on RD workflow conventions.
    Transparently reads from structures.zip if present (no unpack needed).

    Args:
        directory: Root directory to search
        mode: Collection strategy
            - None or 'all': Find all CIF files recursively (legacy behavior)
            - 'total': Find all CIFs EXCLUDING those in relaxed/ subdirectories
            - 'selected': CIFs of the structures that were chosen by selection
                          and sent to relaxation. On DBs with a pre-relax
                          snapshot column, returns the snapshot (4-tuples
                          tagged ``'cif_content_initial'``). Otherwise falls
                          back to ``stage != 'removed'`` on ``cif_content``.
            - 'relaxed': Find CIFs ONLY inside relaxed/ subdirectories
            - 'initial': DB-only. Pre-relax CIF snapshot of stage='relaxed'
                          rows. Returns 4-tuples; mainly for DB-admin /
                          export workflows (``visualize`` uses 'selected'
                          which already prefers the snapshot when present).
    """
    mode_info = {
        'total': 'excluding relaxed/ subdirectories',
        'selected': 'pre-relax snapshot of selected rows (or stage != removed fallback)',
        'relaxed': 'only from relaxed/ subdirectories',
        'initial': 'pre-relax snapshot of relaxed rows (DB column cif_content_initial)',
    }
    effective_mode = mode or 'all'
    print(f"[Visualize] Searching for CIF files in: {directory} (mode: {effective_mode})")
    if effective_mode in mode_info:
        print(f"[Visualize] Filter: {mode_info[effective_mode]}")

    entries = find_cifs_db_aware(directory, mode=mode)
    print(f"[Visualize] Found {len(entries)} CIF files")
    return entries


def _featurize_worker(args):
    """Worker function for parallel feature extraction."""
    path, dim, descriptor, feature_settings = args
    try:
        structure = load_structure(path)

        SelectorClass = _get_selector_class(dim)
        selector = SelectorClass.__new__(SelectorClass)
        selector.feature_settings = feature_settings
        selector.descriptor_type = descriptor

        if dim == '1d':
            selector.periodic_direction = feature_settings.get('periodic_direction', 'z')

        if descriptor == 'soap':
            default_rcut = StructureVisualizer.DEFAULT_RCUT.get(dim, 6.0)
            r_cut = feature_settings.get('soap_r_cut', default_rcut)
            if dim == '1d':
                feat = selector._extract_1d_soap_features(structure)
            else:
                periodic = dim in ('3d', '2d', '1d')
                feat = selector._extract_soap_features(structure, periodic=periodic, r_cut=r_cut)
        elif descriptor == 'coulomb':
            feat = selector._extract_coulomb_features(structure)
        else:
            feat = selector._extract_simple_features(structure)

        return path, feat, None
    except Exception as e:
        return path, None, str(e)


class StructureVisualizer:
    """
    Feature extraction + dimensionality reduction + plotting for CIF datasets.

    Reuses the feature-extraction methods from the existing selector classes
    (base_selector / selector_0d-3d) without performing any selection.
    """

    DEFAULT_RCUT = {'3d': 6.0, '2d': 5.0, '1d': 4.0, '0d': 5.0}

    def __init__(self, dim, descriptor='simple', feature_settings=None):
        """
        Args:
            dim: Dimension string ('3d', '2d', '1d', '0d').
            descriptor: Feature descriptor type ('simple', 'soap', 'coulomb').
            feature_settings: Optional dict passed to the selector (soap_r_cut,
                soap_n_max, soap_l_max, pca_variance, ...).
        """
        self.dim = dim.lower()
        self.descriptor = descriptor.lower()
        self.feature_settings = dict(feature_settings or {})
        self.feature_settings.setdefault('descriptor_type', self.descriptor)

        self._validate()

        SelectorClass = _get_selector_class(self.dim)
        self._selector = SelectorClass.__new__(SelectorClass)
        self._selector.feature_settings = self.feature_settings
        self._selector.descriptor_type = self.descriptor

        if self.dim == '1d':
            self._selector.periodic_direction = self.feature_settings.get(
                'periodic_direction', 'z')

    def _validate(self):
        valid_descriptors = {
            '3d': ['simple', 'soap'],
            '2d': ['simple', 'soap'],
            '1d': ['soap'],
            '0d': ['simple', 'soap', 'coulomb'],
        }
        choices = valid_descriptors.get(self.dim, [])
        if self.descriptor not in choices:
            raise ValueError(
                f"Descriptor '{self.descriptor}' is not available for {self.dim.upper()}. "
                f"Choose from: {choices}"
            )
        if self.descriptor in ('soap', 'coulomb') and not DSCRIBE_AVAILABLE:
            raise ImportError(
                f"{self.descriptor.upper()} descriptor requires dscribe package. "
                "Install with: pip install dscribe"
            )

    def extract_features(self, cif_dir, mode=None, workers=1):
        """
        Extract features from all CIF files found (recursively) in *cif_dir*.

        Args:
            cif_dir: Directory containing CIF files.
            mode: Collection strategy passed to _find_cif_files.
                  See _find_cif_files for details.
            workers: Number of parallel workers (1 = serial).

        Returns:
            (features: np.ndarray, file_list: list[str])
        """
        cif_paths = _find_cif_files(cif_dir, mode=mode)
        if not cif_paths:
            raise FileNotFoundError(
                f"[Visualize] No CIF files found in directory: {cif_dir}"
            )
        return self._extract_features_from_paths(cif_paths, workers=workers)

    def _extract_features_from_paths(self, cif_paths, workers=1):
        """
        Extract features from a pre-collected list of CIF file paths.

        Handles featurization, normalization (StandardScaler), and optional
        PCA in one pass.  Used internally by both :meth:`extract_features`
        and :meth:`_auto_collect_features`.

        Args:
            cif_paths: list of absolute CIF file paths.
            workers: Number of parallel workers (1 = serial).

        Returns:
            (features: np.ndarray, valid_files: list[str])
        """
        all_features = []
        valid_files = []

        desc_label = self.descriptor.upper()
        if workers > 1:
            tasks = [(p, self.dim, self.descriptor, self.feature_settings) for p in cif_paths]
            with ProcessPoolExecutor(max_workers=workers) as pool:
                for path, feat, err in tqdm(pool.map(_featurize_worker, tasks, chunksize=max(1, len(tasks)//workers//4)),
                                             total=len(tasks), desc=f"[Visualize] Featurizing ({desc_label})"):
                    if feat is not None:
                        all_features.append(feat)
                        valid_files.append(path)
                    else:
                        print(f"[Visualize] Error processing {entry_basename(path)}: {err}")
        else:
            for path in tqdm(cif_paths, desc=f"[Visualize] Featurizing ({desc_label})"):
                try:
                    structure = load_structure(path)
                    feat = self._extract_single(structure)
                    all_features.append(feat)
                    valid_files.append(path)
                except Exception as e:
                    print(f"[Visualize] Error processing {entry_basename(path)}: {e}")

        if not all_features:
            raise RuntimeError(
                "[Visualize] No structures could be featurized from the given paths"
            )

        feature_matrix = np.array(all_features)
        print(f"[Visualize] Featurized {len(valid_files)} structures, "
              f"feature dim = {feature_matrix.shape[1]}")

        scaler = StandardScaler()
        normalised = scaler.fit_transform(feature_matrix)

        pca_variance = self.feature_settings.get('pca_variance', 0.95)
        if self.descriptor in ('soap', 'coulomb') and normalised.shape[1] > 50:
            n_components = min(pca_variance, normalised.shape[0], normalised.shape[1])
            pca = PCA(n_components=n_components, random_state=42)
            normalised = pca.fit_transform(normalised)
            var_explained = sum(pca.explained_variance_ratio_) * 100
            print(f"[Visualize] PCA: {feature_matrix.shape[1]}D -> "
                  f"{normalised.shape[1]}D ({var_explained:.1f}% variance)")

        return normalised, valid_files

    def _auto_collect_features(self, cif_dirs, workers=1):
        """
        Automatically collect and extract features for RD workflow stages
        from one or more CIF directories.

        For each directory the method checks for ``remove/`` and ``relaxed/``
        sub-directories (at any nesting level) and infers which workflow
        stages are present:

        * **total** - all CIF files excluding ``relaxed/``
        * **selected** - CIF files chosen by selection and sent to relaxation.
          When the DB carries a pre-relax CIF snapshot
          (``cif_content_initial`` populated by recent RD runs), this stage
          returns the *pre-relax* geometry of the rows currently at
          ``stage='relaxed'`` - i.e. exactly the structures that were
          selected. Together with the ``relaxed`` stage this gives a true
          per-row before/after comparison in the UMAP plot. Falls back to the
          legacy ``stage != 'removed'`` rule (reading ``cif_content``) when
          no snapshot is available.
        * **relaxed** - CIF files only inside ``relaxed/`` (the post-relax
          live CIF). Only when a ``relaxed/`` directory / DB rows exist.

        All CIF paths for a given stage are pooled across directories so
        that normalization (StandardScaler + optional PCA) is applied once
        on the merged set.

        Args:
            cif_dirs: list of directory paths.

        Returns:
            features_dict:    ``{'total': ndarray, ...}``
            file_lists_dict:  ``{'total': [path, ...], ...}``
            (only stages that were detected are included)
        """
        dir_info = {}  
        has_remove_any = False
        has_relaxed_any = False

        for cif_dir in cif_dirs:
            info = detect_stages_db_aware(cif_dir)
            dir_info[cif_dir] = info
            has_remove_any = has_remove_any or info['remove']
            has_relaxed_any = has_relaxed_any or info['relaxed']

            print(f"[Visualize] Auto-detecting workflow stages in: {cif_dir}")
            if info['remove']:
                print("[Visualize]   Found 'remove/' subdirectories "
                      "\u2192 will extract selected structures")
            else:
                print("[Visualize]   No 'remove/' subdirectories found")
            if info['relaxed']:
                print("[Visualize]   Found 'relaxed/' subdirectories "
                      "\u2192 will extract relaxed structures")
            else:
                print("[Visualize]   No 'relaxed/' subdirectories found")
            if info.get('initial', False):
                print("[Visualize]   Pre-relax CIF snapshot column detected "
                      "\u2192 'selected' will use the snapshot (true before/after)")
            stages_here = ['total']
            if info['remove']:
                stages_here.append('selected')
            if info['relaxed']:
                stages_here.append('relaxed')
            print(f"[Visualize]   Stages detected: {', '.join(stages_here)}")

        stages = ['total']
        if has_remove_any:
            stages.append('selected')
        if has_relaxed_any:
            stages.append('relaxed')

        stage_paths = {s: [] for s in stages}
        stage_dir_counts = {s: 0 for s in stages}

        for cif_dir in cif_dirs:
            info = dir_info[cif_dir]
            for stage in stages:
                if stage == 'selected' and not info['remove']:
                    continue
                if stage == 'relaxed' and not info['relaxed']:
                    continue
                if not has_remove_any and not has_relaxed_any:
                    mode = None  
                else:
                    mode = stage
                paths = _find_cif_files(cif_dir, mode=mode)
                if paths:
                    stage_paths[stage].extend(paths)
                    stage_dir_counts[stage] += 1

        stage_dir_counts['total'] = len(cif_dirs)

        print("[Visualize] --- Summary ---")
        for stage in stages:
            n_files = len(stage_paths[stage])
            n_dirs = stage_dir_counts[stage]
            print(f"[Visualize] {stage.capitalize()} structures: "
                  f"{n_files} files from {n_dirs} director{'y' if n_dirs == 1 else 'ies'}")

        # ------------------------------------------------------------------
        # Shared featurization pipeline.
        #
        # Previously each stage called _extract_features_from_paths
        # independently, which fit a **separate** StandardScaler (+ optional
        # PCA) per stage.  This destroyed the relative scale between stages:
        # features of 'selected' (a proper subset of 'total') ended up in a
        # different normalised space, so UMAP could not align them, and the
        # comparison plot showed selected and total as disjoint clouds even
        # though one is a subset of the other.
        #
        # Fix: collect all **unique** (by entry_to_path) structures across
        # all stages, featurize them **once** with a shared scaler + PCA,
        # then slice the resulting feature vectors back by stage.
        # ------------------------------------------------------------------
        unique_paths = []
        feat_key_by_path = {}  # entry_to_path(key) -> (stage, original_path)
        seen = set()
        for stage in stages:
            for p in stage_paths[stage]:
                key = entry_to_path(p) if not isinstance(p, str) else p
                if key not in seen:
                    seen.add(key)
                    unique_paths.append(p)

        if not unique_paths:
            return {}, {}

        print(f"\n[Visualize] Featurizing {len(unique_paths)} unique structures "
              "in one pass (shared StandardScaler + PCA)\n")
        pooled_features, pooled_valid = self._extract_features_from_paths(
            unique_paths, workers=workers
        )

        feat_by_key = {}
        for i, p in enumerate(pooled_valid):
            key = entry_to_path(p) if not isinstance(p, str) else p
            feat_by_key[key] = pooled_features[i]

        features_dict = {}
        file_lists_dict = {}

        for stage in stages:
            paths = stage_paths[stage]
            if not paths:
                print(f"[Visualize] Skipping '{stage}' \u2013 no CIF files found")
                continue
            flist = []
            feats = []
            for p in paths:
                key = entry_to_path(p) if not isinstance(p, str) else p
                vec = feat_by_key.get(key)
                if vec is None:
                    continue
                feats.append(vec)
                flist.append(p)
            if not feats:
                print(f"[Visualize] Skipping '{stage}' \u2013 no CIF files survived featurization")
                continue
            features_dict[stage] = np.array(feats)
            file_lists_dict[stage] = flist
            print(f"[Visualize] '{stage}' final: {len(flist)} structures "
                  f"({features_dict[stage].shape[1]}-D shared feature space)")

        return features_dict, file_lists_dict

    def _extract_single(self, structure):
        """Dispatch to the correct feature-extraction method."""
        if self.descriptor == 'soap':
            default_rcut = self.DEFAULT_RCUT.get(self.dim, 6.0)
            r_cut = self.feature_settings.get('soap_r_cut', default_rcut)
            if self.dim == '1d':
                return self._selector._extract_1d_soap_features(structure)
            periodic = self.dim in ('3d', '2d', '1d')
            return self._selector._extract_soap_features(
                structure, periodic=periodic, r_cut=r_cut)
        elif self.descriptor == 'coulomb':
            return self._selector._extract_coulomb_features(structure)
        else:
            return self._selector._extract_simple_features(structure)

    def reduce_dimensions(self, features_dict, method='umap', **kwargs):
        """
        Reduce feature matrices to 2D for plotting.

        Args:
            features_dict: dict mapping label -> np.ndarray of features.
                Convention: 'total' is the reference dataset;
                others ('selected', 'relaxed', ...) are projected into the
                same space.
            method: 'umap' or 'tsne'.
            **kwargs: Extra parameters forwarded to the reducer.

        Returns:
            dict[str, np.ndarray]  - same keys, each value is (N, 2).
        """
        method = method.lower()

        if method == 'umap':
            result = self._reduce_umap(features_dict, **kwargs)
        elif method == 'tsne':
            result = self._reduce_tsne(features_dict, **kwargs)
        else:
            raise ValueError(f"Unknown reduction method: {method}. Use 'umap' or 'tsne'.")

        self._last_reduction_method = method
        return result

    def _reduce_umap(self, features_dict, **kwargs):
        if not UMAP_AVAILABLE:
            raise ImportError(
                "umap-learn is required for UMAP visualization. "
                "Install with: pip install umap-learn"
            )
    
        params = dict(
            n_neighbors=30, min_dist=0.3, n_components=2,
            random_state=42, metric='euclidean', low_memory=True,
        )
        params.update(kwargs)
    
        keys = list(features_dict.keys())
        arrays = [features_dict[k] for k in keys]
        sizes = [a.shape[0] for a in arrays]
        combined = np.vstack(arrays)

        # Deduplicate: identical feature vectors (e.g. 'selected' entries
        # that also appear in 'total') MUST receive the exact same UMAP
        # coordinate.  Without deduplication UMAP's SGD treats each row
        # index independently and can drift duplicates apart.
        unique_feats, inverse = np.unique(combined, axis=0, return_inverse=True)
        print(f"[Visualize] Running UMAP on {unique_feats.shape[0]} unique "
              f"feature vectors (from {combined.shape[0]} total incl. duplicates) \u2026")

        reducer = umap.UMAP(**params)
        reduced_unique = reducer.fit_transform(unique_feats)
        all_reduced = reduced_unique[inverse]
    
        reduced = {}
        offset = 0
        for key, sz in zip(keys, sizes):
            reduced[key] = all_reduced[offset:offset + sz]
            offset += sz
            print(f"[Visualize]   '{key}': {sz} points")
    
        return reduced

    def _reduce_tsne(self, features_dict, **kwargs):
        from sklearn.manifold import TSNE

        params = dict(n_components=2, perplexity=30, random_state=42)
        params.update(kwargs)

        keys = list(features_dict.keys())
        if 'total' in keys:
            keys.remove('total')
            keys.insert(0, 'total')

        arrays = [features_dict[k] for k in keys]
        sizes = [a.shape[0] for a in arrays]
        combined = np.vstack(arrays)

        print(f"[Visualize] Running t-SNE on {combined.shape[0]} points ...")
        reducer = TSNE(**params)
        all_reduced = reducer.fit_transform(combined)

        reduced = {}
        offset = 0
        for key, sz in zip(keys, sizes):
            reduced[key] = all_reduced[offset:offset + sz]
            offset += sz

        return reduced

    DATASET_STYLES = {
        'total':    {'color': '#CCCCCC', 'alpha': 0.4, 's': 10, 'zorder': 1,
                     'label_prefix': 'Total'},
        'selected': {'color': '#F7B801', 'alpha': 0.7, 's': 14, 'zorder': 2,
                     'label_prefix': 'Selected'},
        'relaxed':  {'color': '#4ECDC4', 'alpha': 0.7, 's': 15, 'zorder': 3,
                     'label_prefix': 'Relaxed'},
    }

    def plot_comparison(self, reduced_dict, output_path, title=None, dpi=300):
        """
        Multi-dataset scatter plot.

        Args:
            reduced_dict: dict[str, np.ndarray] - 2-col arrays.
            output_path: File path for the saved figure.
            title: Optional plot title.
            dpi: Output resolution.
        """
        fig, ax = plt.subplots(figsize=(10, 8))

        method_label = self._guess_method_label()

        draw_order = sorted(
            reduced_dict.keys(),
            key=lambda k: self.DATASET_STYLES.get(k, {}).get('zorder', 10)
        )

        for key in draw_order:
            coords = reduced_dict[key]
            style = self.DATASET_STYLES.get(key, {
                'color': '#999999', 'alpha': 0.6, 's': 12, 'zorder': 5,
                'label_prefix': key.capitalize()
            })
            label = f"{style['label_prefix']} (n={len(coords)})"
            ax.scatter(
                coords[:, 0], coords[:, 1],
                s=style['s'], alpha=style['alpha'], c=style['color'],
                edgecolors='none', marker='o', rasterized=True,
                label=label, zorder=style['zorder'],
            )

        ax.set_xlabel(f"{method_label} Dimension 1", fontsize=12)
        ax.set_ylabel(f"{method_label} Dimension 2", fontsize=12)
        if title:
            ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
        ax.legend(loc='upper right', fontsize=10, markerscale=2,
                  framealpha=0.9)

        for spine in ax.spines.values():
            spine.set_linewidth(0.5)

        ax.set_facecolor('white')
        fig.set_facecolor('white')
        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        plt.close(fig)
        print(f"[Visualize] Saved comparison plot -> {output_path}")

    def plot_density(self, reduced_data, output_path, title=None, dpi=300):
        """
        2D density heatmap of a single reduced dataset using hexbin.

        Args:
            reduced_data: np.ndarray of shape (N, 2).
            output_path: File path for saved figure.
            title: Optional plot title.
            dpi: Output resolution.
        """
        fig, ax = plt.subplots(figsize=(10, 8))

        method_label = self._guess_method_label()

        hb = ax.hexbin(
            reduced_data[:, 0], reduced_data[:, 1],
            gridsize=60, cmap='YlOrRd', mincnt=1,
        )
        cb = fig.colorbar(hb, ax=ax, label='Count')

        ax.set_xlabel(f"{method_label} Dimension 1", fontsize=12)
        ax.set_ylabel(f"{method_label} Dimension 2", fontsize=12)
        if title:
            ax.set_title(title, fontsize=14, fontweight='bold', pad=15)

        for spine in ax.spines.values():
            spine.set_linewidth(0.5)

        ax.set_facecolor('white')
        fig.set_facecolor('white')
        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        plt.close(fig)
        print(f"[Visualize] Saved density plot  -> {output_path}")

    def _guess_method_label(self):
        """Heuristic label for axis (UMAP / t-SNE)."""
        if getattr(self, '_last_reduction_method', None) == 'tsne':
            return "t-SNE"
        return "UMAP"

    def save_features(self, features_dict, file_lists_dict, output_path):
        """
        Save feature matrices and file lists to a .npz archive.

        Args:
            features_dict:    {'total': ndarray, 'selected': ndarray, ...}
            file_lists_dict:  {'total': [path, ...], 'selected': [path, ...], ...}
            output_path:      Path for the .npz file.
        """
        save_kwargs = {
            'meta_dim': np.array([self.dim]),
            'meta_descriptor': np.array([self.descriptor]),
        }
        for key in features_dict:
            save_kwargs[f'features_{key}'] = features_dict[key]
        for key in file_lists_dict:
            # file_lists_dict entries can be 3-tuples or 4-tuples;
            # convert to path string first to avoid inhomogeneous-shape
            # errors from np.array(tuple, dtype=str).
            paths = [entry_to_path(f) if isinstance(f, tuple) else f
                     for f in file_lists_dict[key]]
            save_kwargs[f'files_{key}'] = np.array(paths, dtype=str)

        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        np.savez_compressed(output_path, **save_kwargs)
        print(f"[Visualize] Saved features -> {output_path}")

    @staticmethod
    def load_features(input_path):
        """
        Load a .npz archive previously saved by :meth:`save_features`.

        Returns:
            (features_dict, file_lists_dict, metadata)
        """
        data = np.load(input_path, allow_pickle=True)

        features_dict = {}
        file_lists_dict = {}
        metadata = {}

        for key in data.files:
            if key.startswith('features_'):
                label = key[len('features_'):]
                features_dict[label] = data[key]
            elif key.startswith('files_'):
                label = key[len('files_'):]
                file_lists_dict[label] = list(data[key])
            elif key.startswith('meta_'):
                label = key[len('meta_'):]
                metadata[label] = str(data[key][0]) if data[key].ndim > 0 else str(data[key])

        print(f"[Visualize] Loaded features from {input_path}")
        for k, v in features_dict.items():
            print(f"  {k}: {v.shape[0]} samples, {v.shape[1]} dims")

        return features_dict, file_lists_dict, metadata

    def generate_report(self, features_dict, file_lists_dict, output_path,
                        method, reduction_params):
        """
        Write a human-readable report summarising the visualisation run.

        Args:
            features_dict / file_lists_dict: Same dicts used elsewhere.
            output_path: Directory where output files were written.
            method: Reduction method ('umap' / 'tsne').
            reduction_params: Dict of reducer parameters.
        """
        report_path = os.path.join(output_path, 'visualization_report.txt')
        os.makedirs(output_path, exist_ok=True)

        lines = [
            "=" * 60,
            "GEWUM Visualization Report",
            "=" * 60,
            f"Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"Dimension: {self.dim.upper()}",
            f"Descriptor: {self.descriptor}",
            f"Reduction method: {method.upper()}",
            "",
            "--- Dataset Statistics ---",
        ]

        for key in features_dict:
            n_samples = features_dict[key].shape[0]
            n_dims = features_dict[key].shape[1]
            lines.append(f"  {key}: {n_samples} structures, {n_dims} feature dims")

        lines.append("")
        lines.append("--- Reduction Parameters ---")
        for k, v in (reduction_params or {}).items():
            lines.append(f"  {k}: {v}")

        lines.append("")
        lines.append("--- Output Files ---")
        for fname in sorted(os.listdir(output_path)):
            fpath = os.path.join(output_path, fname)
            if os.path.isfile(fpath):
                size_kb = os.path.getsize(fpath) / 1024
                lines.append(f"  {fname}  ({size_kb:.1f} KB)")

        lines.append("=" * 60)

        with open(report_path, 'w', encoding='utf-8') as fh:
            fh.write('\n'.join(lines) + '\n')
        print(f"[Visualize] Report saved -> {report_path}")

def _build_parser():
    parser = argparse.ArgumentParser(
        description="GEWUM Structure Visualization - feature extraction, "
                    "dimensionality reduction and comparison plotting",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--dim', choices=['0d', '1d', '2d', '3d'], required=True,
                        help='Structure dimension')
    parser.add_argument('--descriptor', choices=['simple', 'soap', 'coulomb'],
                        default='simple', help='Feature descriptor type (default: simple)')
    parser.add_argument('--cif-dir', action='append', nargs='+', default=None,
                        help='One or more CIF directories. The tool auto-detects '
                             'RD workflow stages (total/selected/relaxed) by checking '
                             'for remove/ and relaxed/ subdirectories. '
                             'When the structures.db carries pre-relax snapshots '
                             '(cif_content_initial column populated by recent RD '
                             'runs), the "selected" stage returns the pre-relax '
                             'geometry of those very rows, so the comparison '
                             'plot becomes a true per-row before/after view '
                             '(selected = pre-relax, relaxed = post-relax). '
                             'Accepts both "--cif-dir A B C" and repeated '
                             '"--cif-dir A --cif-dir B". '
                             'Mutually exclusive with --total-dir/--select-dir/--relax-dir.')
    parser.add_argument('--total-dir', required=False, default=None,
                        help='Directory with the total / reference CIF files '
                             '(manual mode, for fine-grained control)')
    parser.add_argument('--select-dir', default=None,
                        help='Directory with selected CIF files '
                             '(manual mode, optional)')
    parser.add_argument('--relax-dir', default=None,
                        help='Directory with relaxed CIF files '
                             '(manual mode, optional)')
    parser.add_argument('--reduction', choices=['umap', 'tsne'], default='umap',
                        help='Reduction method (default: umap)')
    parser.add_argument('--load-features', default=None,
                        help='Load previously saved .npz feature file instead of extracting')
    parser.add_argument('-o', '--output', default='./visualization',
                        help='Output directory (default: ./visualization)')
    parser.add_argument('--soap-r-cut', type=float, default=None,
                        help='SOAP cutoff radius (dimension-dependent default)')
    parser.add_argument('--soap-n-max', type=int, default=4,
                        help='SOAP radial basis n_max (default: 4)')
    parser.add_argument('--soap-l-max', type=int, default=4,
                        help='SOAP angular basis l_max (default: 4)')
    parser.add_argument('--pca-variance', type=float, default=0.95,
                        help='PCA variance ratio to preserve (default: 0.95)')
    parser.add_argument('--no-density', action='store_true',
                        help='Skip density heatmap generation')
    parser.add_argument('--dpi', type=int, default=300,
                        help='Figure DPI (default: 300)')
    parser.add_argument('--title', default=None,
                        help='Custom plot title')
    parser.add_argument('--font-family', type=str, default='Liberation Sans',
                        help='Font family for plots (default: Liberation Sans)')
    parser.add_argument('--font-size', type=int, default=None,
                        help='Global font size for plots (default: matplotlib default)')
    parser.add_argument('--workers', type=int, default=min(os.cpu_count() or 4, 4),
                        help='Number of parallel workers for feature extraction (default: min(cpu_count, 4))')
    parser.add_argument('--cifgen-inp', default=None,
                        help='Path to cifgen.inp file for composition word cloud generation')
    parser.add_argument('--wc-cmap', default=None,
                        help='Colormap name for word cloud (default: pink-blue gradient)')
    parser.add_argument('--no-wordcloud', action='store_true',
                        help='Skip word cloud generation even if cifgen.inp is provided')
    return parser


def main():
    parser = _build_parser()
    args = parser.parse_args()

    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['font.sans-serif'] = [args.font_family] + matplotlib.rcParams['font.sans-serif']
    if args.font_size is not None:
        matplotlib.rcParams['font.size'] = args.font_size

    feature_settings = {
        'descriptor_type': args.descriptor,
        'soap_n_max': args.soap_n_max,
        'soap_l_max': args.soap_l_max,
        'pca_variance': args.pca_variance,
        'include_lattice_params': True,
        'include_position_stats': True,
        'include_basic_features': True,
    }
    if args.soap_r_cut is not None:
        feature_settings['soap_r_cut'] = args.soap_r_cut

    viz = StructureVisualizer(dim=args.dim, descriptor=args.descriptor,
                              feature_settings=feature_settings)

    os.makedirs(args.output, exist_ok=True)

    features_dict = {}
    file_lists_dict = {}

    if args.load_features:
        features_dict, file_lists_dict, meta = StructureVisualizer.load_features(
            args.load_features)
    elif args.cif_dir:
        cif_dirs = [d for sub in args.cif_dir for d in sub]
        manual_set = any([args.total_dir, args.select_dir, args.relax_dir])
        if manual_set:
            print("[Visualize] Warning: --cif-dir is set together with "
                  "--total-dir/--select-dir/--relax-dir. "
                  "Using --cif-dir (auto-detect mode); manual dirs are ignored.")
        features_dict, file_lists_dict = viz._auto_collect_features(cif_dirs, workers=args.workers)
    else:
        if args.total_dir is None and args.select_dir is None and args.relax_dir is None:
            parser.error("At least one of --cif-dir, --total-dir, --select-dir, or "
                         "--relax-dir is required (or use --load-features).")

        dir_mode_map = [
            ('total',    args.total_dir,  'total'),
            ('selected', args.select_dir, 'selected'),
            ('relaxed',  args.relax_dir,  'relaxed'),
        ]
        for label, directory, mode in dir_mode_map:
            if directory is not None:
                print(f"\n[Visualize] Extracting features for '{label}' from {directory}")
                feat, flist = viz.extract_features(directory, mode=mode, workers=args.workers)
                features_dict[label] = feat
                file_lists_dict[label] = flist

    if not features_dict:
        print("[Visualize] Error: no feature data available.")
        sys.exit(1)

    npz_path = os.path.join(args.output, 'features.npz')
    viz.save_features(features_dict, file_lists_dict, npz_path)

    reduction_params = {}
    if args.reduction == 'umap':
        reduction_params = dict(n_neighbors=30, min_dist=0.3, n_components=2,
                                random_state=42, metric='euclidean', low_memory=True)
    else:
        reduction_params = dict(n_components=2, perplexity=30, random_state=42)

    reduced = viz.reduce_dimensions(features_dict, method=args.reduction,
                                    **reduction_params)

    import pandas as pd
    for key, coords in reduced.items():
        csv_path = os.path.join(args.output, f'reduced_{key}_{args.reduction}.csv')
        file_list = file_lists_dict.get(key, [])
        df = pd.DataFrame({
            'filename': [entry_basename(f) for f in file_list],
            'filepath': [entry_to_path(f) for f in file_list],
            f'{args.reduction}_dim1': coords[:, 0],
            f'{args.reduction}_dim2': coords[:, 1],
        })
        df.to_csv(csv_path, index=False)
        print(f"[Visualize] Saved {key} reduced coordinates: {csv_path} ({len(df)} points)")

    reduced_npz_path = os.path.join(args.output, f'reduced_{args.reduction}.npz')
    save_kwargs = {}
    for key, coords in reduced.items():
        save_kwargs[f'{key}_coords'] = coords
        if key in file_lists_dict:
            save_kwargs[f'{key}_files'] = np.array(
                [entry_to_path(f) for f in file_lists_dict[key]], dtype=object
            )
    save_kwargs['method'] = np.array(args.reduction)
    np.savez_compressed(reduced_npz_path, **save_kwargs)
    print(f"[Visualize] Saved all reduced coordinates: {reduced_npz_path}")

    method_upper = args.reduction.upper()
    scatter_path = os.path.join(args.output, f'comparison_{args.reduction}.png')
    title = args.title or f"{args.dim.upper()} {args.descriptor.upper()} - {method_upper}"
    viz.plot_comparison(reduced, scatter_path, title=title, dpi=args.dpi)

    if not args.no_density:
        for key, coords in reduced.items():
            if coords.shape[0] < 3:
                print(f"[Visualize] Skipping density for '{key}' (too few points)")
                continue
            density_path = os.path.join(args.output, f'density_{key}_{args.reduction}.png')
            density_title = (f"{args.dim.upper()} {args.descriptor.upper()} - "
                             f"{method_upper} Density ({key})")
            viz.plot_density(coords, density_path, title=density_title, dpi=args.dpi)

    if args.cifgen_inp and not args.no_wordcloud:
        if os.path.isfile(args.cifgen_inp):
            print(f"\n[Visualize] Generating composition word cloud from {args.cifgen_inp}")
            formula_freq, formula_atoms_wc, total_lines = parse_cifgen(args.cifgen_inp)
            print(f"[Visualize] Parsed {total_lines} lines, {len(formula_freq)} unique compositions")
            wc_path = os.path.join(args.output, 'composition_wordcloud.png')
            generate_wordcloud(formula_freq, formula_atoms_wc, wc_path,
                              cmap_name=args.wc_cmap, dpi=args.dpi)
        else:
            print(f"[Visualize] Warning: cifgen.inp not found: {args.cifgen_inp}")

    viz.generate_report(features_dict, file_lists_dict, args.output,
                        args.reduction, reduction_params)

    print(f"\n[Visualize] All outputs written to {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
