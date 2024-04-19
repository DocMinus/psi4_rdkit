"""
calculates the heat of formation of aniline, water, and nitrobenzene
based on "sumformula", not smiles of complete molecule
"""

import os
import sys

import psi4
from rdkit import Chem
from rdkit.Chem import AllChem

# from mymodules import conformers

# Set the computational options
# psi4.set_options({"basis": "6-31G", "scf_type": "pk", "reference": "rhf"})
# uhf is "easier" to converge than rhf?
psi4.set_options({"basis": "6-31G", "scf_type": "pk", "reference": "uhf"})
psi4.core.set_num_threads(8)

HARTREE2KJMOL = psi4.constants.hartree2kJmol


def printLog(*args, **kwargs):
    """
    Prints the given arguments to the console and appends them to a log file in same path.

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
            "output_dHfrxn_sumformula.txt",
        ),
        "a",
    ) as file:
        print(*args, **kwargs, file=file)


def xyz_coordinates(smiles: str) -> str:
    rdkit_mol = Chem.MolFromSmiles(smiles)
    rdkit_mol = Chem.AddHs(rdkit_mol)
    AllChem.EmbedMolecule(rdkit_mol)
    AllChem.MMFFOptimizeMolecule(rdkit_mol)
    # skipping conformation generation in this case
    return Chem.MolToXYZBlock(rdkit_mol)


# define the molecules
# below works best for larger molecules
# smiles = {"molO": "[O]", "molC": "[C]", "molN": "[N]", "molH2": "[H][H]", "molH": "[H]"}
# mols = {k: xyz_coordinates(v) for k, v in smiles.items()}
# # Create the Psi4 molecule object using the XYZ string
# mols_psi4 = {k: psi4.geometry(v) for k, v in mols.items()}

H = psi4.geometry(
    """
H 0.0 0.0 0.0
"""
)

N = psi4.geometry(
    """
N 0.0 0.0 0.0
"""
)

O = psi4.geometry(
    """
O 0.0 0.0 0.0
"""
)

C = psi4.geometry(
    """
C 0.0 0.0 0.0
"""
)

H2 = psi4.geometry(
    """
H 0.0 0.0 0.0
H 0.0 0.0 0.74
"""
)

energy_name = "scf"  # "B3LYP"

energyO, _ = psi4.energy(name=energy_name, molecule=O, return_wfn=True)
energyC, _ = psi4.energy(name=energy_name, molecule=C, return_wfn=True)
energyN, _ = psi4.energy(name=energy_name, molecule=N, return_wfn=True)
energyH, _ = psi4.energy(name=energy_name, molecule=H, return_wfn=True)
energyH2, _ = psi4.energy(name=energy_name, molecule=H2, return_wfn=True)


energy_O_kjmol = energyO * HARTREE2KJMOL
energy_C_kjmol = energyC * HARTREE2KJMOL
energy_N_kjmol = energyN * HARTREE2KJMOL
energy_H_kjmol = energyH * HARTREE2KJMOL
energy_H2_kjmol = energyH2 * HARTREE2KJMOL

print("\n\n--------")

printLog(f"energy Oxygen: {energy_O_kjmol:.2f} kJ/mol")
printLog(f"energy Oxygen: {energyO:.2f} Hartree")
printLog(f"energy Carbon: {energy_C_kjmol:.2f} kJ/mol")
printLog(f"energy Carbon: {energyC:.2f} Hartree")
printLog(f"energy Nitrogen: {energy_N_kjmol:.2f} kJ/mol")
printLog(f"energy Nitrogen: {energyN:.2f} Hartree")
printLog(f"energy Hydrogen: {energy_H_kjmol:.2f} kJ/mol")
printLog(f"energy Hydrogen: {energyH:.2f} Hartree")
printLog(f"energy Hydrogen2: {energy_H2_kjmol:.2f} kJ/mol")
printLog(f"energy Hydrogen2: {energyH2:.2f} Hartree")

printLog("\nreferences: 1: aniline +31, 2. Water: -258, 3. nitrobenzene +12 kJ/mol")
printLog("-> reaction energy (aniline+3H2O)-(nitrobenzene) = -590 kJ/mol\n")

# aniline
aniline_Hartree = 6 * energyC + 6 * energyH + energyN
water_Hartree = 2 * energyH + energyO
nitrobenzene_Hartree = 6 * energyC + 5 * energyH + energyN + 2 * energyO
delta_hartree = aniline_Hartree + 2 * water_Hartree - nitrobenzene_Hartree + 6 * energyH
delta_kjmol = delta_hartree * 2625.50
printLog(f"aniline {aniline_Hartree * 2625.50:.2f} kJ/mol")
printLog(f"nitrobenzene {nitrobenzene_Hartree * 2625.50:.2f} kJ/mol")
printLog(f"water {water_Hartree * 2625.50:.2f} kJ/mol")
printLog(f"Psi4 calculates: {delta_kjmol:.2f} kJ/mol\n")
