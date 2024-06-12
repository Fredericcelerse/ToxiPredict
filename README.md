# ToxiPredict
This is a small project including AI-based approaches and extended chemical space exploration to predict the toxicity of any organic molecules

In this example, we developed a small Random Forest (optimized using Bayesian Optimization) which learns how to predict thte toxicity of any new organic molecules generated thanks to combinatorial approach. 

## Prerequisites

### Anaconda

To execute the code, we will set up a specific environment using Anaconda. To install it, visit [Anaconda Installation](https://docs.anaconda.com/free/anaconda/install/).

### Setup conda environment

We first create the conda environment: 
```
conda create -n ToxiPredict python=3.8
```

Then we activate the conda environment:
```
conda activate ToxiPredict
```

Once the environment is properly created, we will install the python libraries we need to execute the code:
```
conda install -c conda-forge rdkit pandas numpy scikit-learn scikit-optimize
```
```
pip install joblib
```

## Project architecture

The project is made in three parts:  
1. **Build and train an AI-based model to predict toxicity of organic molecules**  
2. **Generate new random molecules based on already existed organic molecules**  
3. **Apply the AI model on these new molecules**  

### Build and train an AI-based model to predict toxicity of organic molecules

The project starts by building and training an AI-based model that allows us to predict the toxicity of any organic molecule, which is indeed the kernel of this project. To do that, we employed the following pipeline and implemented it within the [RandomForest_with_BO.py](RandomForest_with_BO.py) file:  
**1.1**  We took an already existed database, coming from the work of Setyia et al. ([https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00566-4](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10826801/)). We called it here [SMILES_and_toxicity.csv](SMILES_and_toxicity.csv)  
**1.2**  We then load the data and convert each smiles by a fingerprint using the Morgan Fingerprints available in RDKit (https://www.rdkit.org/docs/GettingStartedInPython.html)  
**1.3**  The data have to be of the same size and are then used to train a simple Random Forest model with initial hyperparameters  
**1.4**  The hyperparameters are then optimized employing a Bayesian Optimization (BO, https://scikit-optimize.github.io/stable/auto_examples/bayesian-optimization.html)  
**1.5**  The model is finally saved into an  external file named [my_random_forest_model.joblib](my_random_forest_model.joblib)  

```
python RandomForest_with_BO.py
```

### Generate new random molecules based on already existed organic molecules

Once the model is saved, we used the already existed SMILES stored in [SMILES_and_toxicity.csv](SMILES_and_toxicity.csv) to create a new set of organic molecules to test our model. To do that, we used the BRICS approach implemented in RDKit to fragmentate the SMILES and stored them in a new external file. We then recombined them randomly to create a new set of random molecules. For more details, see the directory [Generate_SMILES](Generate_SMILES)

### Apply the AI model on these new molecules

The new molecules where we can test the code is stored in the [Generate_SMILES/smiles_generated.csv](Generate_SMILES/smiles_generated.csv) file. We can finally use the script [Predict_toxicity.py](Predict_toxicity.py) to see if the new generated molecules are toxic are not

```
python Predict_toxicity.py
```

## Conclusions

This small project was designed to show how AI and combinatorial approaches for database curation can be used to tackle biochemical issues. This project will be enhanced in the future, for instance by adding more restraints in the SMILES generation or by adding more properties to predict (such as stability or solubility) by using other models that can be combined together at the end during the prediction process.


