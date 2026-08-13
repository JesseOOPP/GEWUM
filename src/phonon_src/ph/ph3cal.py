import os
import datetime
import numpy as np
from ase import Atoms
from phono3py import Phono3py
from phono3py.interface.phono3py_yaml import Phono3pyYaml
from gewum.src.phonon_src.utils.phonon_helpers import (
    get_primitive_cell,
    to_ase_atoms,
    to_phonopy_atoms,
)
from typing import Iterable

class FC3CalculationWorkflow:
    """
    This class generates displacements and computes third-order force constants (FC3)
    using phono3py and mattersim.
    """

    def __init__(
        self,
        atoms: Atoms,
        find_prim: bool = False,
        work_dir: str = None,
        supercell_matrix: np.ndarray = None,
        max_atoms: int = None,
    ):
        """
        Initialize the FC3CalculationWorkflow.

        Args:
            atoms (Atoms): ASE Atoms object with structure and calculator.
            find_prim (bool, optional): Whether to find the primitive cell.
            work_dir (str, optional): Directory to store results.
            supercell_matrix (np.ndarray, optional): Supercell matrix for phono3py.
            max_atoms (int, optional): Maximum atoms in supercell.
        """
        assert atoms.calc is not None, "Atoms must have an attached calculator."
        
        if find_prim:
            self.atoms = get_primitive_cell(atoms)
            self.atoms.calc = atoms.calc
        else:
            self.atoms = atoms

        self.work_dir = work_dir or self._create_default_workdir()
        self.supercell_matrix = supercell_matrix
        self.max_atoms = max_atoms

    def _create_default_workdir(self):
        current_datetime = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
        return f"{current_datetime}-{self.atoms.get_chemical_formula()}-fc3"

    def run(self):
        """
        Execute the FC3 calculation workflow:
        1. Generate displacements with phono3py
        2. Compute forces for FC3
        3. Compute FC3
        """
        current_path = os.path.abspath(".")
        try:
            os.makedirs(self.work_dir, exist_ok=True)
            os.chdir(self.work_dir)

            unitcell_ph3 = to_phonopy_atoms(self.atoms)
            ph3 = Phono3py(
                unitcell_ph3,
                supercell_matrix=self.supercell_matrix,
                primitive_matrix='auto',
                log_level=2
            )
            ph3.generate_displacements()
            ph3.save("phono3py_disp.yaml")

            scs_with_disp = ph3.supercells_with_displacements
            forces_fc3 = []
            print(f"Calculating forces for {len(scs_with_disp)} supercells...")
            for sc in scs_with_disp:
                atoms_sc = to_ase_atoms(sc)
                atoms_sc.calc = self.atoms.calc
                forces_fc3.append(atoms_sc.get_forces())

            num_atoms = len(ph3.supercell)
            forces_fc3 = np.array(forces_fc3).reshape(-1, num_atoms, 3)
            np.savetxt("FORCES_FC3", forces_fc3.reshape(-1, 3))

            ph3yml = Phono3pyYaml()
            ph3yml.read("phono3py_disp.yaml")
            ph3 = Phono3py(
                ph3yml.unitcell,
                supercell_matrix=ph3yml.supercell_matrix,
                primitive_matrix=ph3yml.primitive_matrix
            )
            ph3.dataset = ph3yml.dataset
            ph3.forces = forces_fc3

            print("Computing third-order force constants (FC3)...")
            ph3.produce_fc3()
            print("FC3 calculation completed.")

            return ph3

        except Exception as e:
            print(f"Error in workflow: {str(e)}")
            raise
        finally:
            os.chdir(current_path)
