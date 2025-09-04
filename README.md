## predicting_fixation

Data and _scripts_ in this repository were used to elaborate original article "Computational simulation of ecological drift for generating minimal microbiomes identifies key experimental and biotic factors influencing success" ([preprint](https://doi.org/10.1101/2025.07.10.664178)).

### Summary
This pipeline generates models that are capable of predicting fixation of OTUs in a given community type or types from its **compositional data**. We establish successful fixation is reached when **only one OTU is present** in **each phylogenetic core group** (PCG) in **50%** or **90%** of the instances of that community.

Models are trained with dilution-growth experiment data. Once trained, these models are capable of predicting **how many transfers are needed** for fixation to occur in a new, dataset given information about its **composition** and **PCGs**.

Specifically, our objective is using these models to **identify useful parameters** for prediction of fixation. Being able to predict how many cycles we need for a given starting community will help us design dilution-growth experiments.


### Dataset availability

The parsed datasets we have used are available [here](https://github.com/silvtal/predicting_fixation/tree/main/1_datasets/simulation_results).

There are two tables - one defines "success" (successful fixation) at a 50% threshold and the other defines it at a 90% threshold.


### Dataset generation

Our communities don't come from 16S sequencing, but have instead been **artificially generated with [this script](https://github.com/silvtal/predicting_fixation/blob/main/1_datasets/generate_simcomms.R)**. There's also a set of [artificially generated](https://github.com/silvtal/predicting_fixation/blob/main/1_datasets/generate_PCGtable.py) "PCG tables" that assign functional groups to the species in these simcomms.

In these datasets, each row corresponds to a microbial community that has been subjected to a dilution-growth experiment with a given dilution factor. The columns contain diversity metrics and other characteristics of our simcomms. There's also an additional column containing our target variable. The target variable is "successful fixation" - the number of dilution-growth cycles each community is expected to take until it reaches a state of successul fixation as defined above.

This target variable has been obtained from neutral simulations of growth-dilution experiments. The wrappers we used for running these simulations are in the `simulation_wrappers` folder and use our [`dilgrowth` R package](https://github.com/silvtal/dilgrowth/). The results from these simulations (which output many folders containing multiple abundance data tables) were parsed into the final available .csv files with custom scripts available at the [`2_generate_training_data`](https://github.com/silvtal/predicting_fixation/tree/main/2_generate_training_data) folder of this repository.



### Analysis

Analysis scripts for one-group communities can be found in the [`3_analysis`](https://github.com/silvtal/predicting_fixation/tree/main/3_analysis) folder of this repository.

- `0__` scripts allow for preliminary visualization of the variables and their relationship with OTU fixation.

- `1__RF_success.R` goes a step further in order to quantify the effects of our variables. It defines the target variable, "success", as the dilution-growth cycle number in which one OTU is fixated, meaning it has reached a relative abundance of 50% or 90% from the total community size. Then it creates a series of random forest models to quantify the importance of each variable.

- In `2__GLM_failure.R` we define the target variable, "failure", as a boolean value which indicates whether fixation has failed (1, fixation did not occur) in a community before the Nth cycle or not (0). Since it's a binary variable, we can use a GLM which returns linear coefficients for each variable. We make multiple models for multiple N values ranging from 10 to 1000.

All figures are saved to the `figures` folder.
