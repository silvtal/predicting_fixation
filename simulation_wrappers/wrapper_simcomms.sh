#!/bin/bash

########################################
##### 1000 simuls, fixation at 0.5 #####
##### dilution factors: 12 values  #####
########################################################################
##### Initial communities: SIMULATED with naive_simuls.R           #####
#####                      size: 10⁶, 10⁹                          #####
#####                      species: 10, 100, 1000                  #####
#####                      abundance distribs: uniform, log-normal #####
########################################################################
## --save_all TRUE :: necessary to detect when success happens

fa=0.5
cores=16
nd=100
ns=1000

# submitdir==location of silvtal/simuls scripts/functions (for neutral simulation)
export submitdir="/home/silviatm/micro/varios/2021_06_28__null_models/null_model_v3"
export SCRIPT="v3.R" ## <- main simulation script
export resultsparent="/home/silviatm/micro/varios/2021_06_28__null_models/results_null_model"
# R==Rscript's location, can be empty, must add the final "/"
export R="/home/silviatm/R-4.0.5/bin/Rscript/"

## (1) input data (initial communities simulated with naive_simuls.R)
export abuntableloc=/home/silviatm/micro/varios/2021_06_28__null_models/data/simcomms

## create folder for slurm logs
resultsdir=$resultsparent"/simcomms_"$(date -I)
mkdir $resultsdir

## start loop
for distrib in "uniform" "lognorm"
do
for size in 10000 1e+06
do
for sp in 10 100 1000
do

  export abuntable="$distrib"_"$sp"sp_size"$size".tsv

  for fdil in 0.8 0.5 0.1 0.05 0.01 0.008 0.005 0.001 0.25 0.05 0.04 0.025 0.005 0.004 0.0025 0.0005 0.0004 0.00025
  do

  # for each sample
  for sa in {1..30}
  do

  echo -e \#\!/bin/bash > temp$sa
  echo -e \#SBATCH \-o $resultsdir/slurm.$sa.$fdil.$distrib.$size.\%N.\%j.out \# STDOUT >> temp$sa
  echo -e \#SBATCH \-e $resultsdir/slurm.$sa.$fdil.$distrib.$size.\%N.\%j.err \# STDERR >> temp$sa
  echo '

  # 0. Variables
  # ============
  ## prepare

  # create workdir if it doesnt exist ; only do that for the first script that runs for each sample !!
  export workdir=/temporal/silviatm/'$sa'
  rm -rf $workdir #clean first

  if [ ! -d $workdir ]
  then
  mkdir -p -v $workdir
  # copy all files/folders in submitdir to workdir
  cp -r '$submitdir'/* $workdir
  # + the abundance table
  cp '$abuntableloc'/'$abuntable' $workdir
  fi

  # change directory to the temporary directory on the compute-node
  cd $workdir
  
  ## run simuls
  '$R'Rscript $workdir/'$SCRIPT' \
  -a $workdir/'$abuntable' -s '$sa' \
  --dilution '$fdil' --no_of_dil '$nd' --no_of_simulations '$ns' \
  --fixation_at '$fa' --fix_percentage TRUE \
  --outputname '$sa' --outdir $workdir/SIMCOMM_SIMS_'$fdil'_'$distrib'_'$size'_sp_'$sp'_'$sa' --cores '$cores' \
  --save_all TRUE
  
  ## copy MY OUTPUT DIR from the temporary directory on the compute-node
  ## I only want my outputdir+slurm log to be copied, not everything
  cp -r --copy-contents $workdir/SIMCOMM_SIMS_'$fdil'_'$distrib'_'$size'_sp_'$sp'_'$sa' '$resultsdir >> temp$sa$Core

  ## launch
  if compgen -G "${resultsparent}"/simcomms_*/SIMCOMM_SIMS_"${fdil}"_"${distrib}"_"${size}"_sp_"${sp}"_"${sa}" ]]; then
  echo 'skipping SIMCOMM_SIMS_'$fdil'_'$distrib'_'$size'_sp_'$sp'_'$sa', already exists'
  else
  echo 'running neutral simulations for SIMCOMM_SIMS_'$fdil'_'$distrib'_'$size'_sp_'$sp'_'$sa
  sbatch -A microbioma_serv -p bio temp$sa$Core
  echo ""
  fi

 done
done
done
done
