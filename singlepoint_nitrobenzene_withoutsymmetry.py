""" some early testing, similar to water, result far off from reported """

import os
import sys

import psi4
from rdkit import Chem
from rdkit.Chem import AllChem

psi4.core.set_num_threads(8)


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
            os.path.dirname(os.path.abspath(sys.argv[0])), "output_nitrobenzene.txt"
        ),
        "a",
    ) as file:
        print(*args, **kwargs, file=file)


# os.environ["PSI_DISABLE_GPU"] = "True"
smiles = "O=N(=O)c1ccccc1"
rdkit_mol = Chem.MolFromSmiles(smiles)
rdkit_mol = Chem.AddHs(rdkit_mol)
AllChem.EmbedMolecule(rdkit_mol)
AllChem.MMFFOptimizeMolecule(rdkit_mol)
xyz = Chem.MolToXYZBlock(rdkit_mol)
# remove first two lines and insert symmetry information.
# lines = xyz.split("\n")
# xyz_coordinates = "\n".join(lines[2:])
# xyz_coordinates = "symmetry c1\n" + xyz_coordinates

# Create the Psi4 molecule object using the XYZ string
molecule = psi4.geometry(xyz)

# Set the computational options basis: 6-31G; dz
psi4.set_options({"basis": "6-31G", "scf_type": "pk", "reference": "uhf"})
# psi4.set_options({"basis": "dz", "scf_type": "pk", "reference": "uhf"})
energy = psi4.optimize("B3LYP", molecule=molecule)
# energy = psi4.energy("B3LYP", return_wfn=False, molecule=molecule)

# Calculate the total energy of reference species (H2 and O2)
H2 = psi4.geometry(
    """
    H
    H 1 0.74
    """
)
O2 = psi4.geometry(
    """
    O
    O 1 1.21
    """
)

N2 = psi4.geometry(
    """
    N
    N 1 1.09
    """
)

carbon = psi4.geometry(
    """
    C 0 0 0
    """
)


energy_carbon = psi4.energy("B3LYP", molecule=carbon, return_wfn=False)
energy_H2 = psi4.optimize("B3LYP", molecule=H2, return_wfn=False)
energy_O2 = psi4.optimize("B3LYP", molecule=O2, return_wfn=False)
energy_N2 = psi4.optimize("B3LYP", molecule=N2, return_wfn=False)

# Calculate the heat of formation
print("\n------------------------------------------ \n\n\n")
delta_energy = (
    energy - 6 * energy_carbon - 0.5 * energy_N2 - 5 * 0.5 * energy_H2 - energy_O2
)
heat_of_formation = delta_energy * psi4.constants.hartree2kJmol
printLog(f"Energy nitrobenzene: {energy*psi4.constants.hartree2kJmol:.2f} kJ/mol")
printLog(f"Heat of formation: {heat_of_formation:.2f} kJ/mol")
printLog("\n------------------------------------------ \n\n\n")
