"""
GEWUM Third-Order Force Constants (FC3) Calculation Workflow
For thermal conductivity calculations using phono3py
"""
import os
import datetime
import numpy as np
from ase import Atoms
from phono3py import Phono3py
from phono3py.interface.phono3py_yaml import Phono3pyYaml
from typing import Optional
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms


# ============================================================
# Phonon utility functions (from phonon_helpers.py)
# ============================================================

def get_primitive_cell(atoms: Atoms):
    """Get primitive cell from ASE atoms object"""
    phonopy_atoms = Phonopy(
        to_phonopy_atoms(atoms), primitive_matrix="auto", log_level=2
    )
    primitive = phonopy_atoms.primitive
    return to_ase_atoms(primitive)


def to_phonopy_atoms(atoms: Atoms):
    """Transform ASE atoms object to Phonopy object"""
    phonopy_atoms = PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        cell=atoms.get_cell(),
        masses=atoms.get_masses(),
        positions=atoms.get_positions(),
    )
    return phonopy_atoms


def to_ase_atoms(phonopy_atoms):
    """Transform Phonopy object to ASE atoms object"""
    atoms = Atoms(
        symbols=phonopy_atoms.symbols,
        cell=phonopy_atoms.cell,
        masses=phonopy_atoms.masses,
        positions=phonopy_atoms.positions,
        pbc=True,
    )
    return atoms


# ============================================================


class FC3CalculationWorkflow:
    """
    Workflow for generating displacements and computing third-order force constants (FC3)
    using phono3py and MatterSim.
    """

    def __init__(
        self,
        atoms: Atoms,
        find_prim: bool = False,
        work_dir: Optional[str] = None,
        supercell_matrix: Optional[np.ndarray] = None,
        max_atoms: Optional[int] = None,
    ):
        """
        Initialize the FC3CalculationWorkflow.

        Args:
            atoms: ASE Atoms object with structure and calculator attached
            find_prim: Whether to find the primitive cell
            work_dir: Directory to store results (default: auto-generated)
            supercell_matrix: Supercell matrix for phono3py (default: 2x2x2)
            max_atoms: Maximum atoms in supercell (optional)
        """
        assert atoms.calc is not None, "Atoms must have an attached calculator."
        
        if find_prim:
            self.atoms = get_primitive_cell(atoms)
            self.atoms.calc = atoms.calc
        else:
            self.atoms = atoms

        self.work_dir = work_dir or self._create_default_workdir()
        self.supercell_matrix = supercell_matrix if supercell_matrix is not None else np.diag([2, 2, 2])
        self.max_atoms = max_atoms

    def _create_default_workdir(self) -> str:
        """Create default work directory name with timestamp"""
        current_datetime = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
        return f"{current_datetime}-{self.atoms.get_chemical_formula()}-fc3"

    def run(self) -> Phono3py:
        """
        Execute the FC3 calculation workflow:
        1. Generate displacements with phono3py
        2. Compute forces for FC3
        3. Compute FC3
        
        Returns:
            Phono3py object with computed FC3
        """
        current_path = os.path.abspath(".")
        try:
            os.makedirs(self.work_dir, exist_ok=True)
            os.chdir(self.work_dir)

            print(f"[FC3] Initializing phono3py with supercell matrix: {self.supercell_matrix.diagonal()}")
            unitcell_ph3 = to_phonopy_atoms(self.atoms)
            ph3 = Phono3py(
                unitcell_ph3,
                supercell_matrix=self.supercell_matrix,
                primitive_matrix='auto',
                log_level=2
            )
            ph3.generate_displacements()
            ph3.save("phono3py_disp.yaml")
            print(f"[FC3] Generated displacements: phono3py_disp.yaml")

            scs_with_disp = ph3.supercells_with_displacements
            forces_fc3 = []
            print(f"[FC3] Calculating forces for {len(scs_with_disp)} supercells...")
            
            for i, sc in enumerate(scs_with_disp):
                if (i + 1) % 10 == 0 or i == 0:
                    print(f"  Processing supercell {i + 1}/{len(scs_with_disp)}")
                atoms_sc = to_ase_atoms(sc)
                atoms_sc.calc = self.atoms.calc
                forces_fc3.append(atoms_sc.get_forces())

            num_atoms = len(ph3.supercell)
            forces_fc3 = np.array(forces_fc3).reshape(-1, num_atoms, 3)
            np.savetxt("FORCES_FC3", forces_fc3.reshape(-1, 3))
            print(f"[FC3] Forces saved: FORCES_FC3")

            ph3yml = Phono3pyYaml()
            ph3yml.read("phono3py_disp.yaml")
            ph3 = Phono3py(
                ph3yml.unitcell,
                supercell_matrix=ph3yml.supercell_matrix,
                primitive_matrix=ph3yml.primitive_matrix
            )
            ph3.dataset = ph3yml.dataset
            ph3.forces = forces_fc3

            print("[FC3] Computing third-order force constants...")
            ph3.produce_fc3()
            print("[FC3] FC3 calculation completed!")

            return ph3

        except Exception as e:
            print(f"[FC3] Error in workflow: {str(e)}")
            raise
        finally:
            os.chdir(current_path)
