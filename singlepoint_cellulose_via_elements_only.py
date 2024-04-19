""" early testing; single sugar unit to 'simulate' cellulose by breaking down to elements instead of whole molecule"""

import psi4

# from rdkit import Chem
# from rdkit.Chem import AllChem

psi4.core.set_num_threads(8)


# smiles = "CO[C@@H]1O[C@H](CO)[C@@H](OC)[C@H](O)[C@H]1O"
# rdkit_mol = Chem.MolFromSmiles(smiles)
# rdkit_mol = Chem.AddHs(rdkit_mol)
# AllChem.EmbedMolecule(rdkit_mol)
# AllChem.MMFFOptimizeMolecule(rdkit_mol)
# xyz = Chem.MolToXYZBlock(rdkit_mol)

# molecule = psi4.geometry(xyz)
# # Set the computational options basis: 6-31G; dz
# psi4.set_options({"basis": "6-31G", "scf_type": "pk", "reference": "uhf"})
# energy = psi4.optimize("B3LYP", molecule=molecule)

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
    energy
    - 8 * energy_carbon
    - 0 * energy_N2
    - 16 * 0.5 * energy_H2
    - 6 * 0.5 * energy_O2
)
heat_of_formation = delta_energy * psi4.constants.hartree2kJmol
print(f"Heat of formation: {heat_of_formation:.2f} kJ/mol")
print("\n------------------------------------------ \n\n\n")
# outcome: -8650 kJ/mol; reality: -669 kJ/mol
