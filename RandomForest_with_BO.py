#!/usr/bin/env python

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score

### LOAD AND POST-TREAT THE DATABASE

# We hidde the warning messages from rdkit
logger = RDLogger.logger()
logger.setLevel(RDLogger.CRITICAL)

# We load the data
data_path = "./SMILES_and_toxicity.csv"
data = pd.read_csv(data_path)

# We then clean the data
data.dropna(inplace=True)
data.drop_duplicates(subset=['SMILES'], inplace=True)
print(data.info())

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

### SPLIT THE DATABASE

X = np.array(list(data['Fingerprints']))
y = data['Toxicity'].values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

### DEFINITION OF THE MODEL AND TRAINING

rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
rf_model.fit(X_train, y_train)
y_pred_rf = rf_model.predict(X_test)
y_prob_rf = rf_model.predict_proba(X_test)[:, 1]
accuracy_rf = accuracy_score(y_test, y_pred_rf)
auc_score_rf = roc_auc_score(y_test, y_prob_rf)

print("Random Forest - Accuracy:", accuracy_rf)
print("Random Forest - AUC Score:", auc_score_rf)

### BAYESIAN OPTIMIZATION

from skopt import BayesSearchCV
from skopt.space import Real, Categorical, Integer

# We define the hyperparameters space for optimization
param_grid = {
	'n_estimators': Integer(100, 500),
	'max_depth': Integer(3, 30),
	'min_samples_split': Integer(2, 20), 
	'min_samples_leaf': Integer(1, 10)
}

# We configure the Bayesian Optimization
opt = BayesSearchCV(
	rf_model, 
	param_grid, 
	n_iter = 32, 
	scoring = 'roc_auc', 
	cv = 3,
	n_jobs = -1, 
	random_state = 42
)

opt.fit(X_train, y_train)

print("Best parameters: ", opt.best_params_)
print("Best AUC score: ", opt.best_score_)

# We save the model
import joblib
joblib.dump(opt.best_estimator_, 'my_random_forst_model.joblib')

