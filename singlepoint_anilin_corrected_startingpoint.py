import os

import psi4
from rdkit import Chem
from rdkit.Chem import AllChem


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
            os.path.dirname(os.path.abspath(sys.argv[0])),
            "output_aniline_optimized.txt",
        ),
        "a",
    ) as file:
        print(*args, **kwargs, file=file)


# anilin
smiles = "Nc1ccccc1"
rdkit_mol = Chem.MolFromSmiles(smiles)
rdkit_mol = Chem.AddHs(rdkit_mol)
AllChem.EmbedMolecule(rdkit_mol)
AllChem.MMFFOptimizeMolecule(rdkit_mol)
xyz = Chem.MolToXYZBlock(rdkit_mol)

# Create the Psi4 molecule object using the XYZ string
molecule = psi4.geometry(xyz)
molecule = psi4.geometry(
    """
    N            0.677727093732    -2.236609343852     0.000000000000
    C            0.275765086528    -0.909928354252     0.000000000000
    C            0.066591160326    -0.219642308620    -1.213158709285
    C           -0.338200261606     1.116202119865    -1.205878072113
    C           -0.544811841825     1.797944007576     0.000000000000
    C           -0.338200261606     1.116202119865     1.205878072113
    C            0.066591160326    -0.219642308620     1.213158709285
    H            0.827002952168    -2.732132719414    -0.861176504387
    H            0.827002952168    -2.732132719414     0.861176504387
    H            0.223293040369    -0.736881840058    -2.156130375615
    H           -0.493124100728     1.627253046397    -2.151652894501
    H           -0.859421711109     2.835938255494     0.000000000000
    H           -0.493124100728     1.627253046397     2.151652894501
    H            0.223293040369    -0.736881840058     2.156130375615
"""
)
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
printLog("\n------------------------------------------ \n\n\n")
delta_energy = energy - 6 * energy_carbon - 0.5 * energy_N2 - 3 * energy_H2
heat_of_formation = delta_energy * psi4.constants.hartree2kJmol
printLog(f"Energy anilin: {energy*psi4.constants.hartree2kJmol:.2f} kJ/mol")
printLog(f"Heat of formation: {heat_of_formation:.2f} kJ/mol")
printLog("\n------------------------------------------ \n\n\n")
