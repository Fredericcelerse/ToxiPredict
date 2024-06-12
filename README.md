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
conda install -c conda-forge rdkit pandas
```
```
conda install -c pytorch pytorch
```
```
conda install -c pyg pyg pytorch-scatter
```
