#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import AllChem
import random
import pandas as pd

### DEFINITION OF USEFUL FONCTIONS FOR COMBINING SMILES

# Check the validity of the smiles
def sanitize_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        try:
            Chem.SanitizeMol(mol)
            return mol
        except:
            return None
    return None

# Find potential connection sites
def find_connection_sites(mol):
    free_valence = []
    for atom in mol.GetAtoms():
        if atom.GetImplicitValence() > 0:
            free_valence.append(atom.GetIdx())
    return free_valence

# Remove all the fictive atoms and try to "kekulizy" the molecule
def remove_dummy_atoms(mol):
    rw_mol = Chem.RWMol(mol)
    dummy_atoms = [atom for atom in rw_mol.GetAtoms() if atom.GetSymbol() == '*']
    for atom in reversed(sorted(dummy_atoms, key=lambda x: x.GetIdx())):
        rw_mol.RemoveAtom(atom.GetIdx())
    try:
        Chem.SanitizeMol(rw_mol, Chem.SANITIZE_KEKULIZE)
        return rw_mol
    except Chem.KekulizeException:
        return None

# Combining  two smiles together
def combine_fragments(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if not mol1 or not mol2:
        return None

    sites1 = find_connection_sites(mol1)
    sites2 = find_connection_sites(mol2)
    if not sites1 or not sites2:
        return None

    new_mol = Chem.RWMol(mol1)
    new_mol.InsertMol(mol2)
    # Here, we ensure that the new bond added does not violate the valence rule !
    if can_form_bond(mol1, sites1[0], mol2, sites2[0]):
        new_mol.AddBond(sites1[0], len(mol1.GetAtoms()) + sites2[0], Chem.BondType.SINGLE)
    else:
        return None

    final_mol = remove_dummy_atoms(new_mol)
    if final_mol is None:
        return None
    return Chem.MolToSmiles(final_mol)

# Check if we can create a new covalent bond or not
def can_form_bond(mol1, idx1, mol2, idx2):
    atom1 = mol1.GetAtomWithIdx(idx1)
    atom2 = mol2.GetAtomWithIdx(idx2)
    # Check the remaining valences !
    return atom1.GetImplicitValence() > 0 and atom2.GetImplicitValence() > 0

# Check if a SMILES is chemically valid
def validate_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    try:
        Chem.SanitizeMol(mol)
        return True
    except:
        return False

### LOAD THE FRAGMENTS AND GENERATE THE SMILES
with open('unique_fragments.txt', 'r') as file:
    fragments = [line.strip() for line in file if line.strip()]

# Random generation of new molecules
new_molecules = set()
num_mol = 1000
for _ in range(num_mol):
    frag1 = random.choice(fragments)
    frag2 = random.choice(fragments)
    new_smiles = combine_fragments(frag1, frag2)
    if new_smiles:
        new_molecules.add(new_smiles)

# Print and save the new generated molecules
print(f"Number of new unique and valid generated molecules: {len(new_molecules)}")
df = pd.DataFrame(new_molecules, columns=['SMILES'])
df.to_csv('smiles_generated.csv', index=False)

