# Detecting subtle transcriptomic perturbations in single-cell CRISPRi screening using a sparse supervised autoencoder neural network (SSAE)

This is the code from : 

Detecting subtle transcriptomic perturbations induced by lncRNAs Knock-Down in single-cell CRISPRi screening using a new sparse supervised autoencoder neural network

*Marin Truchi, Caroline Lacoux, Cyprien Gille, Julien Fassy, Virginie Magnone, Rafael Lopez-Goncalvez, Cédric Girard-Riboulleau, Iris Manosalva-Pena, Marine Gautier-Isola, Kevin Lebrigand, Pascal Barbry, Salvatore Spicuglia, Georges Vassaux, Roger Rezzonico, Michel Barlaud, Bernard Mari*

bioRxiv 2023.07.11.548494; doi: https://doi.org/10.1101/2023.07.11.548494


The original GitHub repository of the SSAE used in the manuscript is the property of [Michel Barlaud](https://github.com/MichelBarlaud/SAE-Supervised-Autoencoder-Omics/tree/main). The goal of this repository is to gather materials used for i) the preprocessing of CROP-seq data, ii) the detection of transcriptomic perturbation with the SSAE, iii) the figures of the manuscript. 

---

<figure>
  <img src="https://github.com/marintruchi/lncRNAs_CROPseq_SSAE/blob/main/SSAE_overview.jpg" alt="SSAE_overview"/>
  <figcaption>Two-step SSAE classification of perturbed cells among gRNA-targeted
cells
</figcaption>
</figure>

---

## **Repository Contents**
|Folder | Description |
|:----------|:----------|
|`data`|Contains the input matrices of the SSAE and the link to access the raw and preprocessed sequencing files of the CROP-seq library|
|`scripts`|Contains the scripts for the preprocessing of CROP-seq data, the preparation of SSAE input matrices, the execution of SSAE and the production of the manuscript's figures|
|`functions`|Contains dedicated functions for the execution of SSAE|

 ### **"scripts" repository Contents**   
|Script| Description |
|:----------|:----------|
|`Preprocess_CROPseq_lib.R`|R script to load, manipulate and prepare the count matrices for the SSAE |
|`Run_SSAE_script.py`|Main python script to run the SSAE|
|`Produce_figures.R`|R script to produce the figures of the manuscript|



## **How to run the SSAE** 

We recommend installing and running the SSAE in an anaconda environment.

### Requirements
- python >= 3.8 .
- [Pytorch](https://pytorch.org/get-started/locally/).
- The following packages : [numpy](https://numpy.org/install/), [matplotlib](https://matplotlib.org/stable/users/installing/index.html), [scikit-learn](https://scikit-learn.org/stable/install.html), [pandas](https://pandas.pydata.org/getting_started.html), [shap](https://pypi.org/project/shap/), [captum](https://captum.ai/#quickstart). 



### How to use it

For now, 
Everything is ready, you can just run the script you want using, for example, the run code button of your Spyder IDE. Alternatively, you can run the command `python [script_name].py` in the Anaconda Prompt from the root of this folder (i.e. where you downloaded and unzipped this repository).

Each script will produce results (statistical metrics, top features...) in a results folder.

You can change the database used, and other parameters, near the start of each script.
