# Generation of new SMILES

Through these two scripts, we aimed at creating a new small comobinatorial database in  order to test our AI model for the toxicity.  

1. We first load the data from the [SMILES_and_toxicity.csv](../SMILES_and_toxicity.csv) file and stored only the SMILES
2. We then use the BRICS algorithm from RDKit to fragmentate the SMILES in order tto create new fragments
3. The new fragments are stored in an external file called [unique_fragments.txt](unique_fragments.txt)
4. We then use this txt file to combine randomly two fragments together, and check the  chemical validity of each result
5. Finally, the good molecules are stored in the external file [smiles_generated.csv](smiles_generated.csv)

```
python Fragmentation_BRICS.py
```
```
python Generate_new_smiles.py
```
