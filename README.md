# ToxiPredict
This is a small project including AI-based approach with Random Forest and Bayesian optimization to predict the toxicity of any organic molecules

In this example, we developed a small Message Passing Neural Network that considered each molecule of the dataset as a graph, and learns using a vector made of atomic and bond properties. 

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
1.1  
1.2  
1.3  

