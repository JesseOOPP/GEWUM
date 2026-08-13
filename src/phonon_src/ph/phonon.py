# -*- coding: utf-8 -*-
import datetime
import os
from typing import List, Tuple, Union

import numpy as np
from ase import Atoms
from phonopy import Phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from tqdm import tqdm

from gewum.src.phonon_src.utils.phonon_helpers import (
    get_primitive_cell,
    to_ase_atoms,
    to_phonopy_atoms,
)


class PhononWorkflow:
    """
    Phonon dispersion calculation workflow using phonopy.
    """

    def __init__(
        self,
        atoms: Atoms,
        find_prim: bool = False,
        work_dir: str = None,
        amplitude: float = 0.01,
        supercell_matrix: np.ndarray = None,
        qpoints_mesh: np.ndarray = None,
        band_path: List[Tuple[str, np.ndarray]] = None,
        band_npoints: int = 101,
        band_color: str = 'deepskyblue'
    ):
        """
        Args:
            atoms: ASE atoms object with attached calculator.
            find_prim: Use primitive cell for phonon calculation.
            work_dir: Directory for phonon results.
            amplitude: Displacement magnitude in Angstrom.
            supercell_matrix: Supercell matrix (3,) or (3,3). Required.
            qpoints_mesh: Q-point mesh for IBZ integration.
            band_path: Custom band path [(label, qpoint), ...].
            band_npoints: Points per band segment.
            band_color: Band line color.
        """
        assert atoms.calc is not None, "Atoms must have an attached calculator"
        
        if find_prim:
            self.atoms = get_primitive_cell(atoms)
            self.atoms.calc = atoms.calc
        else:
            self.atoms = atoms
        
        if work_dir:
            self.work_dir = work_dir
        else:
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
            self.work_dir = f"{timestamp}-{atoms.get_chemical_formula()}-phonon"
        
        self.amplitude = amplitude
        
        if supercell_matrix is None:
            raise ValueError("supercell_matrix is required")
        
        if supercell_matrix.shape == (3, 3):
            self.supercell_matrix = supercell_matrix
        elif supercell_matrix.shape == (3,):
            self.supercell_matrix = np.diag(supercell_matrix)
        else:
            raise ValueError("Supercell matrix must be (3,) or (3,3)")
        
        if qpoints_mesh is not None:
            assert qpoints_mesh.shape == (3,), "Qpoints mesh must be (3,)"
            self.qpoints_mesh = qpoints_mesh
        else:
            nrep = np.diag(self.supercell_matrix)
            self.qpoints_mesh = (6 * nrep if np.allclose(nrep, nrep[0]) else 3 * nrep).astype(int)
        
        self.band_path = band_path
        self.band_npoints = band_npoints
        self.band_color = band_color  

    def compute_force_constants(self, atoms: Atoms, supercell_matrix: np.ndarray):
        """Calculate 2nd-order force constants."""
        print(f"Supercell matrix: \n{supercell_matrix}")
        
        phonon = Phonopy(
            to_phonopy_atoms(atoms),
            supercell_matrix=supercell_matrix,
            primitive_matrix="auto",
            log_level=2,
        )
        phonon.generate_displacements(distance=self.amplitude)
        
        supercells = phonon.supercells_with_displacements
        forces = []
        print("\nCalculating forces for displaced structures...")
        for sc in tqdm(supercells):
            ase_sc = to_ase_atoms(sc)
            ase_sc.calc = self.atoms.calc
            forces.append(ase_sc.get_forces())
        
        phonon.forces = np.array(forces)
        phonon.produce_force_constants()
        phonon.symmetrize_force_constants()
        
        return phonon

    def compute_phonon_spectrum(self, atoms: Atoms, phonon: Phonopy, k_point_mesh):
        """Calculate phonon band structure."""
        print(f"Q-points mesh: {k_point_mesh}")
        phonon.run_mesh(k_point_mesh)
        print(f"Computing phonon dispersion for {atoms.symbols}...\n")
        
        output_name = f"{atoms.symbols}_phonon_band.png"
        
        if self.band_path:
            labels = [p[0] for p in self.band_path]
            qpoints = [p[1] for p in self.band_path]
            path_segments = [[qpoints[i], qpoints[i+1]] for i in range(len(qpoints)-1)]
            
            print("Using custom band path:")
            for i, (start, end) in enumerate(path_segments):
                print(f"  {labels[i]} -> {labels[i+1]}")
            
            band_qpoints, connections = get_band_qpoints_and_path_connections(
                path_segments, npoints=self.band_npoints
            )
            phonon.run_band_structure(
                band_qpoints,
                path_connections=connections,
                labels=labels,
                is_band_connection=True,
                with_eigenvectors=False,
                is_legacy_plot=True
            )
            fig = phonon.plot_band_structure()
            phonon.write_yaml_band_structure()
        else:
            print("Using auto band structure")
            fig = phonon.auto_band_structure(plot=True, write_yaml=True)
        
        ax = fig.gca()
        for line in ax.get_lines():
            if len(set(line.get_xdata())) > 1:
                line.set_color(self.band_color)
        
        fig.savefig(output_name, dpi=600, bbox_inches='tight')
        phonon.save(settings={"force_constants": True})

    @staticmethod
    def check_imaginary_freq(phonon: Phonopy):
        """Check for imaginary frequencies (threshold: -0.1 THz)."""
        band_dict = phonon.get_band_structure_dict()
        frequencies = np.concatenate(
            [np.array(f).flatten() for f in band_dict["frequencies"]]
        )
        has_imaginary = np.any(frequencies < -0.1)
        if has_imaginary:
            print("Warning! Imaginary frequencies found!")
        return has_imaginary

    def run(self):
        """Execute the phonon workflow."""
        current_path = os.path.abspath(".")
        
        try:
            os.makedirs(self.work_dir, exist_ok=True)
            os.chdir(self.work_dir)
            
            phonon = self.compute_force_constants(self.atoms, self.supercell_matrix)
            
            self.compute_phonon_spectrum(self.atoms, phonon, self.qpoints_mesh)
            
            has_imaginary = self.check_imaginary_freq(phonon)
            
            return has_imaginary, phonon
            
        finally:
            os.chdir(current_path)
