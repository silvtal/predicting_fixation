#!/bin/bash

########################################
##### 100 simuls, fixation at 50%  #####
##### dilution factors:  9 values  #####
##### only up to 200 cycles        #####
#########################################################################
##### Initial communities: SIMULATED with generate_simcomms.R       #####
##### (60 total)           size:    10⁶                             #####
#####                      species: 100, 1000                      #####
#####                      abundance distrib: log-normal            #####
#####                            ((30 communities, reps 1-30))      #####
#########################################################################
## --save_all TRUE :: necessary to detect when success happens

fa=1
cores=16
nd=200
ns=100
grow_step=0.01
is_grow_step_a_perc="TRUE"

# submitdir==location of silvtal/simuls scripts/functions (for neutral simulation)
export submitdir="/home/silviatm/micro/varios/2021_06_28__null_models/null_model_v3"
export SCRIPT="dilgrowth/scripts/dilgrowth.R" # _old_1g.R" # "v3.R" ## <- main simulation script
export resultsparent="/home/silviatm/micro/varios/2021_06_28__null_models/results_null_model"
# R==Rscript's location, can be empty, if not empty must add the final "/"
export R="/home/silviatm/R-4.0.5/bin/"

## (1) input data (initial communities simulated with naive_simuls.R)
export abuntableloc="/home/silviatm/micro/varios/2021_06_28__null_models/data/simcomms"
export pcgtableloc="/home/silviatm/micro/varios/2021_06_28__null_models/predicting_fixation/1_datasets/PCGtables/"

## create folder for slurm logs
resultsdir=$resultsparent"/simcomms_WITH_GROUPS_"$(date -I)
mkdir $resultsdir

fdils="0.00025 0.0005 0.001 0.0025 0.005 0.01 0.025 0.05 0.1"
## start loop
for distrib in "lognorm" #"uniform" # only one distrib now
do

for size in "1e+06"
do

for sp in 100 1000
do

  export abuntable="$distrib"_"$sp"sp_size"$size".tsv

for fdil in $fdils
do

# for each sample
for sa in {1..30}
do

for nichedist in EvenNiches SkewedNiches
do
for nicheN in 3 10
do
  pcgtable="$nichedist"_N"$nicheN"_"$sp"_sp.csv
  
  echo -e \#\!/bin/bash > temp$sa
  echo -e \#SBATCH \-o $resultsdir/slurm.$sa.$fdil.$distrib.$size.\%N.\%j.out \# STDOUT >> temp$sa
  echo -e \#SBATCH \-e $resultsdir/slurm.$sa.$fdil.$distrib.$size.\%N.\%j.err \# STDERR >> temp$sa
  echo -e \#SBATCH \-c $cores >> temp$sa
  echo '

  # 0. Variables
  # ============
  ## prepare

  # create workdir if it doesnt exist ; only do that for the first script that runs for each sample !!
  export workdir=/temporal/silviatm/'$sa'_WG
  rm -rf $workdir #clean first

  if [ ! -d $workdir ]
  then
  mkdir -p -v $workdir
  # copy all files/folders in submitdir to workdir
  cp -r '$submitdir'/* $workdir
  # + the abundance table
  cp '$abuntableloc'/'$abuntable' $workdir
  cp '$pcgtableloc'/'$pcgtable' $workdir
  fi

  # change directory to the temporary directory on the compute-node
  cd $workdir

  ## run simuls
  '$R'Rscript $workdir/'$SCRIPT' \
  -a $workdir/'$abuntable' -s '$sa' \
  -p $workdir/'$pcgtable' \
  --dilution '$fdil' --no_of_dil '$nd' --no_of_simulations '$ns' \
  --fixation_at '$fa' \
  --outputname '$sa' --outdir $workdir/SIMCOMM_SIMS_'$nicheN$nichedist'_'$fdil'_'$distrib'_sp_'$sp'_'$sa' --cores '$cores' \
  --save_all TRUE --grow_step '$grow_step' --is_grow_step_a_perc '$is_grow_step_a_perc'
  ## copy MY OUTPUT DIR from the temporary directory on the compute-node
  ## I only want my outputdir+slurm log to be copied, not everything
  cp -r --copy-contents $workdir/SIMCOMM_SIMS_'$nicheN$nichedist'_'$fdil'_'$distrib'_sp_'$sp'_'$sa' '$resultsdir >> temp$sa

  ## launch
  if compgen -G "${resultsparent}"/simcomms_*WITH*/SIMCOMM_SIMS_"$nicheN$nichedist"_"${fdil}"_"${distrib}"_sp_"${sp}"_"${sa}" ]]
  then
  echo 'skipping SIMCOMM_SIMS_'$nicheN$nichedist'_'$fdil'_'$distrib'_sp_'$sp'_'$sa', already exists'
  else
  echo 'running neutral simulations for SIMCOMM_SIMS_'$nicheN$nichedist'_'$fdil'_'$distrib'_sp_'$sp'_'$sa', PCG table: '$pcgtable
  sbatch -A microbioma_serv -p biobis temp$sa
  rm temp$sa
  echo ""
  fi
done
done
done
done
done
done
done
