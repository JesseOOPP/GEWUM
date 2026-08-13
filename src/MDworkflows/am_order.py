"""
GEWUM SQ AM Crystallinity Gate

Screens a directory of relaxed configurations and separates crystalline
(recrystallized) frames from genuinely amorphous ones, so the downstream
amorphous ensemble is not contaminated by ordered structures.

Two complementary order criteria (a frame is flagged crystalline if EITHER
fires):
  1. Averaged Steinhardt bond-order parameter Q6 (Lechner-Dellago variant).
     Crystalline environments give high, sharp Q6; amorphous ones give low,
     broad Q6. A frame is crystalline when Q6 >= q6_threshold, where the
     threshold defaults to a fraction of the crystalline template's own Q6.
  2. SOAP distance to the crystalline template. A frame whose SOAP vector is
     within soap_threshold (cosine distance) of the template has relaxed back
     onto the original crystal and is flagged crystalline.

Flagged frames are moved into a sub-directory (default 'crystalline/'); the
amorphous survivors remain in place. An am_order.csv table records the per-
frame Q6, SOAP distance and verdict.

This module performs no relaxation; it is a pure structural classifier and is
used twice in the SQ pipeline: as a pre-selection filter (mode 'prefilter')
and as a post-refinement QC check (mode 'qc').
"""
import os
import csv
import sys
import shutil
import logging
import argparse
import numpy as np
# scipy renamed sph_harm -> sph_harm_y in 1.15 with a different argument order
# and angle convention; wrap both behind _ylm(m, l, azimuth, polar).
try:
    from scipy.special import sph_harm_y as _sph_harm_y

    def _ylm(m, l, azimuth, polar):
        return _sph_harm_y(l, m, polar, azimuth)
except ImportError:
    from scipy.special import sph_harm as _sph_harm

    def _ylm(m, l, azimuth, polar):
        return _sph_harm(m, l, azimuth, polar)
from ase.io import read
from ase.build import make_supercell
from ase.neighborlist import neighbor_list

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Fallback neighbour cutoff (Angstrom) when the RDF first-minimum auto-detect
# fails (e.g. too few atoms to build a histogram).
FALLBACK_CUTOFF = 3.5

# Below this template Q6 the relative threshold is unreliable (template may be
# mis-read or the cutoff too small); warn rather than silently emptying the pool.
MIN_TEMPLATE_Q6 = 0.2

CSV_FIELDS = ['frame', 'q6', 'q6_frac', 'soap_dist', 'verdict']


def _auto_cutoff(atoms, rmax=6.0, nbins=120):
    """Estimate a first-shell neighbour cutoff from the first RDF minimum.

    Builds a coarse pairwise-distance histogram (up to rmax), takes the tallest
    bin as the first coordination peak, then returns the first local minimum
    after it. Falls back to FALLBACK_CUTOFF when the estimate is unreliable.
    """
    try:
        d = neighbor_list('d', atoms, rmax)
    except Exception:
        return FALLBACK_CUTOFF
    if len(d) == 0:
        return FALLBACK_CUTOFF

    hist, edges = np.histogram(d, bins=nbins, range=(0.0, rmax))
    centers = 0.5 * (edges[:-1] + edges[1:])
    if hist.sum() == 0:
        return FALLBACK_CUTOFF

    peak_idx = int(np.argmax(hist))
    for k in range(peak_idx + 1, len(hist) - 1):
        if hist[k] <= hist[k - 1] and hist[k] <= hist[k + 1]:
            return float(centers[k])
    return FALLBACK_CUTOFF


def compute_q6(atoms, cutoff=None, l=6, return_atoms=False):
    """Global averaged Steinhardt order parameter Q_l (default l=6).

    Uses the Lechner-Dellago averaging: each atom's q_lm is first computed as
    the mean spherical harmonic over its bonds, then averaged again over the
    atom and its neighbours before taking the rotational invariant. This makes
    crystalline vs amorphous environments more separable than the raw q6.

    Args:
        atoms: ASE Atoms (periodic)
        cutoff: neighbour cutoff in Angstrom; auto-detected when None
        l: harmonic degree (6 for Q6)
        return_atoms: also return the per-atom q6 array (0.0 for atoms with no
            neighbours) so callers can compute a crystalline-atom fraction

    Returns:
        float global Q6 (mean over atoms that have neighbours; 0.0 when no bonds
        are found), or (global_q6, per_atom_q6) when return_atoms is True.
    """
    if cutoff is None:
        cutoff = _auto_cutoff(atoms)

    i_idx, j_idx, D = neighbor_list('ijD', atoms, cutoff)
    n_atoms = len(atoms)
    if len(i_idx) == 0:
        return (0.0, np.zeros(n_atoms)) if return_atoms else 0.0

    r = np.linalg.norm(D, axis=1)
    good = r > 1e-8
    r_safe = np.where(good, r, 1.0)
    # _ylm expects (m, l, azimuth, polar); build both angles here.
    azimuth = np.arctan2(D[:, 1], D[:, 0])
    polar = np.arccos(np.clip(D[:, 2] / r_safe, -1.0, 1.0))

    ms = np.arange(-l, l + 1)
    qlm = np.zeros((n_atoms, len(ms)), dtype=complex)
    counts = np.zeros(n_atoms)

    for mi, m in enumerate(ms):
        Y = _ylm(m, l, azimuth, polar)
        np.add.at(qlm[:, mi], i_idx, Y)
    np.add.at(counts, i_idx, 1.0)

    nz = counts > 0
    qlm[nz] /= counts[nz, None]

    # Lechner-Dellago averaging over the atom plus its neighbours.
    neigh_sum = np.zeros_like(qlm)
    for mi in range(len(ms)):
        np.add.at(neigh_sum[:, mi], i_idx, qlm[j_idx, mi])
    denom = counts + 1.0
    qlm_bar = (qlm + neigh_sum) / denom[:, None]

    q6_atom = np.sqrt(
        4.0 * np.pi / (2 * l + 1) * np.sum(np.abs(qlm_bar) ** 2, axis=1))

    if not np.any(nz):
        return (0.0, np.zeros(n_atoms)) if return_atoms else 0.0
    global_q6 = float(np.mean(q6_atom[nz]))
    if return_atoms:
        return global_q6, np.where(nz, q6_atom, 0.0)
    return global_q6


def _soap_vector(atoms, species, r_cut, n_max, l_max):
    """Averaged periodic SOAP descriptor for one structure (dscribe)."""
    from dscribe.descriptors import SOAP
    soap = SOAP(
        species=species,
        r_cut=r_cut,
        n_max=n_max,
        l_max=l_max,
        periodic=True,
        average="inner",
        sparse=False,
    )
    return soap.create(atoms).flatten()


def soap_distance(vec_a, vec_b):
    """Cosine distance (1 - cosine similarity) between two SOAP vectors."""
    na = float(np.linalg.norm(vec_a))
    nb = float(np.linalg.norm(vec_b))
    if na < 1e-12 or nb < 1e-12:
        return 1.0
    cos = float(np.dot(vec_a, vec_b) / (na * nb))
    cos = max(-1.0, min(1.0, cos))
    return 1.0 - cos


def _order_csv_path(frame_dir, mode):
    """Where to write the per-frame order table (parent of frame_dir)."""
    parent = os.path.dirname(os.path.abspath(frame_dir.rstrip(os.sep)))
    name = 'am_order.csv' if mode == 'qc' else 'am_order_prefilter.csv'
    return os.path.join(parent, name)


def screen_directory(frame_dir, template, supercell=(1, 1, 1),
                     q6_ref_frac=0.5, q6_abs=None, q6_cutoff=None,
                     q6_frac_thr=0.5, soap_thr=0.05, soap_params=None,
                     move_subdir='crystalline', mode='prefilter'):
    """Classify every top-level CIF in frame_dir as amorphous or crystalline.

    Crystalline frames are moved to frame_dir/<move_subdir>/; a CSV table with
    per-frame Q6 and SOAP distance is written to the parent of frame_dir.

    Returns:
        (n_amorphous, n_crystalline) or None on fatal error.
    """
    if soap_params is None:
        soap_params = {'r_cut': 6.0, 'n_max': 4, 'l_max': 4}

    if not os.path.isdir(frame_dir):
        logging.error(f"Frame directory not found: {frame_dir}")
        return None

    ref = template.copy()
    if any(s > 1 for s in supercell):
        ref = make_supercell(ref, np.diag(supercell))

    # Resolve ONE neighbour cutoff from the template and reuse it for every
    # frame, so the relative Q6 comparison uses a consistent bond definition.
    cutoff = q6_cutoff if q6_cutoff is not None else _auto_cutoff(ref)

    species = sorted(set(ref.get_chemical_symbols()))
    q6_template, _ = compute_q6(ref, cutoff=cutoff, return_atoms=True)
    try:
        tmpl_soap = _soap_vector(ref, species, **soap_params)
    except Exception as e:
        logging.error(f"Failed to build template SOAP vector: {e}")
        return None

    q6_thr = q6_abs if q6_abs is not None else q6_ref_frac * q6_template
    # An atom is locally crystalline when its own averaged q6 reaches the same
    # threshold; a frame whose crystalline-atom fraction is too high is flagged
    # even if its global-mean Q6 stays low (catches partial crystallization).
    q6_atom_thr = q6_thr
    if q6_abs is None and q6_template < MIN_TEMPLATE_Q6:
        logging.warning(
            f"[{mode}] Template Q6={q6_template:.4f} < {MIN_TEMPLATE_Q6}; "
            f"relative threshold may be unreliable. Check the template CIF or "
            f"pass an absolute --q6-threshold / --q6-cutoff.")
    logging.info(
        f"[{mode}] cutoff={cutoff:.3f} A, Q6_template={q6_template:.4f}, "
        f"q6_threshold={q6_thr:.4f}, q6_frac_threshold={q6_frac_thr:.2f}, "
        f"soap_threshold={soap_thr:.4f}")

    cif_files = sorted(f for f in os.listdir(frame_dir)
                       if f.endswith('.cif') and
                       os.path.isfile(os.path.join(frame_dir, f)))
    if not cif_files:
        logging.warning(f"No CIF files to screen in {frame_dir}")
        return 0, 0

    move_dir = os.path.join(frame_dir, move_subdir)
    records = []
    n_amorphous = 0
    n_crystalline = 0

    for name in cif_files:
        path = os.path.join(frame_dir, name)
        try:
            atoms = read(path)
        except Exception as e:
            logging.warning(f"Skipping unreadable frame {name}: {e}")
            continue

        try:
            q6, q6_atoms = compute_q6(atoms, cutoff=cutoff, return_atoms=True)
        except Exception as e:
            logging.warning(f"Q6 failed for {name}: {e}; treating as crystalline")
            q6 = float('inf')
            q6_atoms = None

        if q6_atoms is not None and q6_atoms.size:
            frac_cryst = float(np.mean(q6_atoms >= q6_atom_thr))
        else:
            frac_cryst = 1.0

        try:
            sdist = soap_distance(
                _soap_vector(atoms, species, **soap_params), tmpl_soap)
        except Exception as e:
            logging.warning(f"SOAP failed for {name}: {e}; skipping SOAP screen")
            sdist = 1.0

        is_crystalline = (q6 >= q6_thr) or (sdist <= soap_thr) or \
            (frac_cryst >= q6_frac_thr)
        verdict = 'crystalline' if is_crystalline else 'amorphous'
        records.append({
            'frame': name,
            'q6': f"{q6:.6f}" if np.isfinite(q6) else 'nan',
            'q6_frac': f"{frac_cryst:.4f}",
            'soap_dist': f"{sdist:.6f}",
            'verdict': verdict,
        })

        if is_crystalline:
            os.makedirs(move_dir, exist_ok=True)
            shutil.move(path, os.path.join(move_dir, name))
            n_crystalline += 1
        else:
            n_amorphous += 1

    csv_path = _order_csv_path(frame_dir, mode)
    try:
        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
            writer.writeheader()
            writer.writerows(records)
    except OSError as e:
        logging.warning(f"Failed to write {csv_path}: {e}")

    n_total = n_amorphous + n_crystalline
    purity = (n_amorphous / n_total * 100.0) if n_total else 0.0
    logging.info(
        f"[{mode}] Amorphicity gate: {n_amorphous}/{n_total} amorphous "
        f"({purity:.1f}%), {n_crystalline} moved to {move_subdir}/")
    logging.info(f"[{mode}] Order table: {csv_path}")

    return n_amorphous, n_crystalline


def _self_test():
    """Sanity check Q6 against analytic references (fcc ~0.575, bcc ~0.511)."""
    from ase.build import bulk
    for name, sym, ref in [('fcc', 'Cu', 0.575), ('bcc', 'Fe', 0.511)]:
        atoms = bulk(sym, name, a=3.6) * (3, 3, 3)
        q6 = compute_q6(atoms)
        print(f"{name} {sym}: Q6={q6:.4f} (reference ~{ref})")


def main():
    parser = argparse.ArgumentParser(
        description='GEWUM SQ crystallinity gate (Q6 + SOAP-to-template)',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('frame_dir', nargs='?',
                        help='Directory of relaxed CIF frames to screen')
    parser.add_argument('--template', type=str, default=None,
                        help='Crystalline template CIF (reference for Q6/SOAP)')
    parser.add_argument('--supercell', type=int, nargs=3, default=[1, 1, 1],
                        metavar=('NX', 'NY', 'NZ'),
                        help='Supercell expansion applied to the template '
                             'reference (match generation; default: 1 1 1)')
    parser.add_argument('--q6-ref-frac', type=float, default=0.5,
                        help='Crystalline if Q6 >= frac*Q6_template (default: 0.5)')
    parser.add_argument('--q6-threshold', type=float, default=None,
                        help='Absolute Q6 threshold (overrides --q6-ref-frac)')
    parser.add_argument('--q6-cutoff', type=float, default=None,
                        help='Neighbour cutoff in Angstrom (default: auto RDF min)')
    parser.add_argument('--q6-frac-threshold', type=float, default=0.5,
                        help='Crystalline if the crystalline-atom fraction '
                             '>= this value (default: 0.5)')
    parser.add_argument('--soap-threshold', type=float, default=0.05,
                        help='Crystalline if SOAP cosine distance to template '
                             '<= threshold (default: 0.05)')
    parser.add_argument('--soap-r-cut', type=float, default=6.0,
                        help='SOAP cutoff radius (default: 6.0)')
    parser.add_argument('--soap-n-max', type=int, default=4,
                        help='SOAP radial basis (default: 4)')
    parser.add_argument('--soap-l-max', type=int, default=4,
                        help='SOAP angular basis (default: 4)')
    parser.add_argument('--mode', type=str, default='prefilter',
                        choices=['prefilter', 'qc'],
                        help='prefilter (before selection) or qc (after refine)')
    parser.add_argument('--move-subdir', type=str, default='crystalline',
                        help='Sub-directory for flagged frames (default: crystalline)')
    parser.add_argument('--self-test', action='store_true',
                        help='Run the Q6 analytic self-check and exit')

    args = parser.parse_args()

    if args.self_test:
        _self_test()
        sys.exit(0)

    if not args.frame_dir or not args.template:
        parser.error('frame_dir and --template are required (unless --self-test)')

    try:
        template = read(args.template)
    except Exception as e:
        logging.error(f"Failed to read template {args.template}: {e}")
        sys.exit(1)

    soap_params = {
        'r_cut': args.soap_r_cut,
        'n_max': args.soap_n_max,
        'l_max': args.soap_l_max,
    }

    result = screen_directory(
        frame_dir=args.frame_dir,
        template=template,
        supercell=tuple(args.supercell),
        q6_ref_frac=args.q6_ref_frac,
        q6_abs=args.q6_threshold,
        q6_cutoff=args.q6_cutoff,
        q6_frac_thr=args.q6_frac_threshold,
        soap_thr=args.soap_threshold,
        soap_params=soap_params,
        move_subdir=args.move_subdir,
        mode=args.mode,
    )
    sys.exit(0 if result is not None else 1)


if __name__ == "__main__":
    main()
