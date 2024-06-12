#!/usr/bin/env python

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from sklearn.ensemble import RandomForestClassifier
from joblib import load

### LOAD AND POST-TREAT THE NEW DATA

# We hidde the warning messages from rdkit
logger = RDLogger.logger()
logger.setLevel(RDLogger.CRITICAL)

# We load the data
data_path = "Generation_smiles/smiles_generated.csv"
data = pd.read_csv(data_path)

# Function to convert SMILES into fingerprints
def smiles_to_fingerprints(smiles, n_bits=2048):
	mol = Chem.MolFromSmiles(smiles)
	if mol is not None:
		fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
		return np.array(fingerprint)
	else:
		return None

# Creation of the fingerprints
data['Fingerprints'] = data['SMILES'].apply(lambda x: smiles_to_fingerprints(x))

# Cleaning of the new data
data = data.dropna(subset=['Fingerprints'])
print(data.head())

X = np.array(list(data['Fingerprints'].apply(lambda x: list(x))))

### LOAD THE MODEL FOR PREDICTION

rf_model = load('my_random_forst_model.joblib')
predictions = rf_model.predict(X)
data['Predicted_Toxicity'] = predictions
print(data[['SMILES', 'Predicted_Toxicity']])

