## predicting_fixation

### Summary
The objective of this pipeline is to predict fixation of OTUs in a given community type or types from its **compositional data**. We establish successful fixation is reached when **only one OTU is present** in **each phylogenetic core group** (PCG) in **50%** of the instances of that community.

To do this, this pipeline trains a regression model with dilution-growth experiment datasets. Once trained, the model is capable of predicting **how many transfers are needed** for fixation to occur in a new, dataset given information about its **composition** and **PCGs**. 


### Datasets
In our case, our abundance data comes from [neutral mechanistic simulations of growth-dilution experiments](https://github.com/silvtal/dilgrowth/), so each iteration is an "instance" of a community type.

These neutral simulations are generated from two datasets:

1. Rhizosphere samples. We took 16S sequencing data from 8 samples and simulated multiple dilution-growth experiments for each sample. The [16S data was processed](https://github.com/silvtal/16S) with dada2 and abundance tables were obtained using QIIME 1. We [grouped the 3000+ 99% similarity OTUs](https://github.com/silvtal/BacterialCore) present in these samples into 12 PCGs. In order to reduce the total number of PCGs, we established that a phylogenetic tree branch needs to have a total abundance of 1% in each sample for it to be considered a PCG. We tried out different dilution factors and noted in which transfer fixation was reached for all PCGs.

These samples were submitted to an actual dilution-growth experiment, but it's important to note that our training data here does not take the actual fixation values but those generated with our mechanistic simulation model.

2. From artificial communities ("simcomms"). The initial abundance tables for this dataset don't come from 16S sequencing, but have instead been artificially generated. There's also a set of simulated "PCG tables" that allow for artificial functional groups to be simulated. 

The wrappers we used for running these simulations in batch for both datasets are in the `simulation_wrappers` folder.


### Main processing script[s]
For predicting fixation for new datasets, we need matrices containing the features or independent variables X (e.g. information about the initial community composition, before any dilution-growth experiment takes place, plus the dilution factor used in the simulated experiment) and the target information or dependent variable y (e.g. the number of transfers needed to reach successful fixation). This target information comes from running [simulations](https://github.com/silvtal/dilgrowth) (see `2_generating_training_data`).

The `create_data` scripts generate tables from the raw neutral simulations results, one per sample.

`create_data_simcomms.R` is designed for the multi-PCGs rhizosphere dataset and depends on the functions defined on `simul_fixation_functions.R`.

`create_data_simcomms.R` is a modification of this script tailored for the computationally generated communities. It depends on some of the functions defined on `simul_fixation_functions.R` but it also has some modifications to adapt for the specific input format of this dataset plus the fact that it doesn't have multiple PCGs. **This script is suitable for many community types**, with an option to include multiple functional groups separately or not.


### Analysis scripts (for one-group communities)

- `0__` scripts allow for preliminary visualization of the variables and their relationship with OTU fixation.

- `1__RF_success.R` goes a step further in order to quantify the effects of our variables. It defines the target variable, "success", as the dilution-growth cycle number in which one OTU is fixated, meaning it has reached a relative abundance of 50% or 90% from the total community size. Then it creates a series of random forest models to quantify the importance of each variable.

- In `2__GLM_failure.R` we define the target variable, "failure", as a boolean value which indicates whether fixation has failed (1, fixation did not occur) in a community before the Nth cycle or not (0). Since it's a binary variable, we can use a GLM which returns linear coefficients for each variable. We make multiple models for multiple N values ranging from 10 to 1000.

All figures are saved to the `figures` folder.
