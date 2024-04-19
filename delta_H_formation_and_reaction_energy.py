"""
calculates the heat of formation of aniline, water, and nitrobenzene
base on smiles as starting point
"""

import os
import sys

import psi4
from rdkit import Chem
from rdkit.Chem import AllChem

from mymodules import conformers

# Set the computational options
# psi4.set_options({"basis": "6-31G", "scf_type": "pk", "reference": "rhf"})
# uhf is "easier" to converge than rhf?
psi4.set_options({"basis": "6-31G", "scf_type": "pk", "reference": "uhf"})
psi4.core.set_num_threads(8)

HARTREE2KJMOL = psi4.constants.hartree2kJmol


def printLog(*args, **kwargs):
    """
    Prints the given arguments to the console and appends them to a log file in same path; filename hardcoded.

    Args:
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Returns:
        None
    """

    print(*args, **kwargs)
    with open(
        os.path.join(
            os.path.dirname(os.path.abspath(sys.argv[0])), "output_dHfrxn_struct.txt"
        ),
        "a",
    ) as file:
        print(*args, **kwargs, file=file)


def xyz_coordinates(smiles: str) -> str:
    """
    Generate XYZ coordinates for a molecule given its SMILES representation.

    Args:
        smiles (str): The SMILES representation of the molecule.

    Returns:
        str: The XYZ coordinates of the molecule.

    """
    rdkit_mol = Chem.MolFromSmiles(smiles)
    rdkit_mol = Chem.AddHs(rdkit_mol)
    AllChem.EmbedMolecule(rdkit_mol)
    AllChem.MMFFOptimizeMolecule(rdkit_mol)
    generator = conformers.ConformerGenerator(max_conformers=5, force_field="uff")
    # simplified for now; got to check case if multiple conformers are generated
    # actually perhaps not even necessary to generate conformers
    mol = generator.generate_conformers(rdkit_mol)
    return Chem.MolToXYZBlock(mol)


# Define the molecule using its SMILES string
mol1 = xyz_coordinates("Nc1ccccc1")
mol2 = xyz_coordinates("O")
mol3 = xyz_coordinates("O=N(=O)c1ccccc1")


# Create the Psi4 molecule object using the XYZ string
mol1_psi4 = psi4.geometry(mol1)
mol2_psi4 = psi4.geometry(mol2)
mol3_psi4 = psi4.geometry(mol3)

# calculate the energy of the molecules
energy1, _ = psi4.energy("B3LYP", molecule=mol1_psi4, return_wfn=True)
energy2, _ = psi4.energy("B3LYP", molecule=mol2_psi4, return_wfn=True)
energy3, _ = psi4.energy("B3LYP", molecule=mol3_psi4, return_wfn=True)
# unit returned is Hartree
heat_of_formation1 = energy1 * HARTREE2KJMOL
heat_of_formation2 = energy2 * HARTREE2KJMOL
heat_of_formation3 = energy3 * HARTREE2KJMOL

print("\n\n--------")

printLog(f"Heat of formation1: {heat_of_formation1:.2f} kJ/mol")
printLog(f"Heat of formation2: {heat_of_formation2:.2f} kJ/mol")
printLog(f"Heat of formation3: {heat_of_formation3:.2f} kJ/mol")
printLog(
    "\nknown reference values: 1: aniline +31, 2. Water: -258, 3. nitrobenzene +12 kJ/mol"
)
printLog("-> reaction energy (aniline+2H2O)-(nitrobenzene) = -590 kJ/mol")
printLog(
    f"Psi4 calculates: {(heat_of_formation1+2*heat_of_formation2)-(heat_of_formation3 + 3* 0.0):.2f} kJ/mol"
)
