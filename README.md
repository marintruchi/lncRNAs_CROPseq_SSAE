# lncRNAs_CROPseq_SSAE

This is the code from : 

Detecting subtle transcriptomic perturbations induced by lncRNAs Knock-Down in single-cell CRISPRi screening using a new sparse supervised autoencoder neural network

Marin Truchi, Caroline Lacoux, Cyprien Gille, Julien Fassy, Virginie Magnone, Rafael Lopez-Goncalvez, Cédric Girard-Riboulleau, Iris Manosalva-Pena, Marine Gautier-Isola, Kevin Lebrigand, Pascal Barbry, Salvatore Spicuglia, Georges Vassaux, Roger Rezzonico, Michel Barlaud, Bernard Mari
bioRxiv 2023.07.11.548494; doi: https://doi.org/10.1101/2023.07.11.548494


## Table of Contents
***
1. [Repository Contents](repository-contents)
2. [Installation](#installation)
3. [How to use](#how-to-use)
  
### **Repository Contents**
|File/Folder | Description |
|:---|:---:|
|`script_autoencoder.py`|Main script to train and evaluate the SAE|

|`script_PLSDA_RF_SVM.py`|Script to fit and evaluate classical methods|

|`datas`|Contains the  databases used in the paper|

|`functions`|Contains dedicated functions for the three main scripts|
    
### **Installation** 
---

To run this code, you will need :
- A version of python, 3.8 or newer. If you are new to using python, we recommend downloading anaconda ([here](https://www.anaconda.com/products/individual)) and using Spyder (available by default from the anaconda navigator) to run the code.
- [Pytorch](https://pytorch.org/get-started/locally/).
- The following packages, all of which except captum and shap are **usually included in the anaconda distribution** : [numpy](https://numpy.org/install/), [matplotlib](https://matplotlib.org/stable/users/installing/index.html), [scikit-learn](https://scikit-learn.org/stable/install.html), [pandas](https://pandas.pydata.org/getting_started.html), [shap](https://pypi.org/project/shap/), [captum](https://captum.ai/#quickstart). To install any package, you can use anaconda navigator's built-in environment manager.

See `requirements.txt` for the exact versions on which this code was developed.

### **How to use**

Everything is ready, you can just run the script you want using, for example, the run code button of your Spyder IDE. Alternatively, you can run the command `python [script_name].py` in the Anaconda Prompt from the root of this folder (i.e. where you downloaded and unzipped this repository).

Each script will produce results (statistical metrics, top features...) in a results folder.

You can change the database used, and other parameters, near the start of each script.
