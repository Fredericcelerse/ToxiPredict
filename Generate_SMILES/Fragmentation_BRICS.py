#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem import rdmolops
import pandas as pd

### DEFINE THE FRAGMENTATION STRATEGY (HERE WE USE BRICS FROM RDKIT)

def fragment_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    # Use BRICS to break the covalent bonds and obtain fragments
    fragmented_mol = BRICS.BreakBRICSBonds(mol)

    # Using GetMolFrags to extract the fragments as separated molecules
    fragments = rdmolops.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
    fragment_smiles = [Chem.MolToSmiles(frag) for frag in fragments]
    return fragment_smiles

### LOAD THE DATA

data_path = '../SMILES_and_toxicity.csv'
data_df = pd.read_csv(data_path)
data_df['Fragments'] = data_df['SMILES'].apply(fragment_smiles)

# Print some data
print(data_df['Fragments'].head())

# Extract all the fragments only once
all_fragments = [frag for sublist in data_df['Fragments'].tolist() for frag in sublist]

# Vanish duplicates by converting the list as an ensemble, then read again the list
unique_fragments = list(set(all_fragments))

# Save the unique fragments in an external file
with open('unique_fragments.txt', 'w') as file:
    for fragment in unique_fragments:
        file.write(f"{fragment}\n")

print(f"Number of unique fragments : {len(unique_fragments)}")

