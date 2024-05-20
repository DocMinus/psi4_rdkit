""" some early testing, heat of formation for water. result rather close to reported """

import psi4

from mymodules.psi_logger import Psi4Logger

Psi4Logger("testfile_output")


molecule = psi4.geometry(
    """
O
H 1 0.96
H 1 0.96 2 104.5
"""
)
# Set the computational options basis: 6-31G; dz
psi4.set_options({"basis": "6-31G", "scf_type": "pk", "reference": "uhf"})
# psi4.set_options({"basis": "dz", "scf_type": "pk", "reference": "uhf"})
optimize = psi4.optimize("B3LYP", molecule=molecule)
energy = psi4.energy("B3LYP", return_wfn=False, molecule=molecule)

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
# B3LYP vs hf
energy_H2 = psi4.energy("B3LYP", molecule=H2, return_wfn=False)
energy_O2 = psi4.energy("B3LYP", molecule=O2, return_wfn=False)

# Calculate the heat of formation
delta_energy = energy - energy_H2 - 0.5 * energy_O2
heat_of_formation = delta_energy * psi4.constants.hartree2kJmol
print(f"Energy water: {energy*psi4.constants.hartree2kJmol:.2f} kJ/mol")
print(f"Heat of formation: {heat_of_formation:.2f} kJ/mol")
print("reported value: -286 kJ/mol")
print("\n\n\n ------------------------------------------ \n\n\n")
