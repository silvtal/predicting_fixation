#!/bin/bash


nd=1000
ns=500
fa=0.5
cores=16
  
# submitdir==location of silvtal/simuls scripts/functions (for neutral simulation)
export submitdir="/home/silviatm/micro/varios/2021_06_28__null_models/null_model_v3"
export SCRIPT="v3.R" ## <- main simulation script
export resultsparent="/home/silviatm/micro/varios/2021_06_28__null_models/results_null_model"
# R==Rscript's location, can be empty, must add the final "/"
export R="/home/silviatm/R-4.0.5/bin/Rscript/"

## (1) input data
### (abundance table from 16S sequencing of tomato rhizosphere, obtained with qiime)
export abuntable="tomate_bc_0.01/Tree/0.99/table.from_biom_0.99.txt"
### PCG data, BacterialCore.py output
export pcgtable="tomate_bc_0.01/Tree/results.txt"

## create folder for slurm logs
mkdir $resultsdir"/"$(date -I)

## start loop

## each sample corresponds to a soil type
## -------------------------------------------------------------

for fdil in 0.5 0.1 0.05 0.01 0.008 0.005 0.001 0.25 0.05 0.04 0.025 0.005 0.004 0.0025 0.0005 0.0004 0.00025
do

# for each pcg
tail -n +2 $pcgtable | tr -d '"' | tr "\t" "|" | while IFS="|" read -r Core Prevalence Abundance Relative_abundances Min Max Average SD Leaves Taxonomy Leaves_number
do

  # for each sample
  for sa in eB1 eF7 eE1 eG3 eA1 eH3 eC2
  do

  echo -e \#\!/bin/bash > temp$sa$Core
  echo -e \#SBATCH \-o "$resultsdir/$(date -I)"/slurm.$sa.$Core.\%N.\%j.out \# STDOUT >> temp$sa$Core
  echo -e \#SBATCH \-e "$resultsdir/$(date -I)"/slurm.$sa.$Core.\%N.\%j.err \# STDERR >> temp$sa$Core
  echo '

  # 0. Variables
  # ============
  ## prepare

  # create workdir if it doesnt exist ; only do that for the first script that runs for each sample !!
  export workdir=/temporal/silviatm/'$sa'.'$Core'
  rm -rf $workdir #clean first

  if [ ! -d $workdir ]
  then
  mkdir -p -v $workdir
  # copy all files/folders in "submitdir" to "workdir"
  cp -r '$submitdir'/* $workdir
  fi

  # change directory to the temporary directory on the compute-node
  cd $workdir
  
  
  
  
  ## run simuls
  '$R'Rscript $workdir/'$SCRIPT' \
  -a $workdir/'$abuntable' -s '$sa' \
  --dilution '$fdil' --no_of_dil '$nd' --no_of_simulations '$ns' --subset "'$Leaves'" \
  --fixation_at '$fa' --fix_percentage TRUE --perc '$Average' \
  --outputname '$Core'_X'$sa' --outdir $workdir/'$fdil'_simuls_tomato_001trX'$sa' --cores '$cores' \
  --save_all TRUE
  
  ## copy MY OUTPUT DIR from the temporary directory on the compute-node
  ## I only want my outputdir log to be copied, not everything
  cp -rf --copy-contents $workdir/'$fdil'_simuls_tomato_001trX'$sa' '$resultsdir >> temp$sa$Core

  ## launch
  sbatch -A microbioma_serv -p biobis temp$i$Core

 done
done
done

