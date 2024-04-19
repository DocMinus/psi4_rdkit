from rdkit import Chem
from rdkit.Chem import AllChem

from mymodules import conformers


def print_bond_lengths(mol):
    conf = mol.GetConformer()
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtomIdx()
        end_atom = bond.GetEndAtomIdx()

        # Get the coordinates of the two atoms
        begin_pos = conf.GetAtomPosition(begin_atom)
        end_pos = conf.GetAtomPosition(end_atom)

        # Calculate the distance between the two atoms
        length = begin_pos.Distance(end_pos)

        # Print the bond length
        print(
            f"Bond length between atom {begin_atom} and atom {end_atom}: {length:.2f} Å"
        )


def main():
    smiles = "O=N(=O)c1ccccc1"
    rdkit_mol = Chem.MolFromSmiles(smiles)
    rdkit_mol = Chem.AddHs(rdkit_mol)
    AllChem.EmbedMolecule(rdkit_mol)
    AllChem.MMFFOptimizeMolecule(rdkit_mol)
    print_bond_lengths(rdkit_mol)

    ff = "uff"  # "mmff94s"
    generator = conformers.ConformerGenerator(max_conformers=5, force_field=ff)

    mol = generator.generate_conformers(rdkit_mol)
    print(generator.get_conformer_energies(rdkit_mol))
    print(generator.get_conformer_energies(mol))
    print(generator.get_conformer_rmsd(mol))
    print_bond_lengths(mol)


if __name__ == "__main__":
    main()
