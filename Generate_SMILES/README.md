# Generation of New SMILES

Through these two scripts, we aim to create a new small combinatorial database to test our AI model for toxicity.

1. **Load Data**: We first load the data from the [SMILES_and_toxicity.csv](../SMILES_and_toxicity.csv) file and store only the SMILES.
2. **Fragmentation**: We use the BRICS algorithm from RDKit to fragment the SMILES in order to create new fragments.
3. **Store Fragments**: The new fragments are stored in an external file named [unique_fragments.txt](unique_fragments.txt)
4. **Combine Fragments**: We then use this text file to randomly combine two fragments together and check the chemical validity of each result.
5. **Store Valid Molecules**: Finally, the valid molecules are stored in the external file [smiles_generated.csv](smiles_generated.csv)

```
python Fragmentation_BRICS.py
```
```
python Generate_new_smiles.py
```
