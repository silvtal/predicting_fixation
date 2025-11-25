#!/bin/bash

########################################
##### Scenario 3: Interactions      #####
##### 100 simuls, fixation at 50%   #####
##### dilution factors: 9 values    #####
##### only up to 200 cycles          #####
#########################################################################
##### Initial communities: SIMULATED with generate_simcomms.R       #####
##### (1 community)        size:    10⁶                             #####
#####                      species: 100                            #####
#####                      abundance distrib: log-normal            #####
#####                      groups: 3 (even) or 10 (skewed)          #####
#####                      WITH INTER-GROUP INTERACTIONS            #####
#########################################################################
## --save_all TRUE :: necessary to detect when success happens
##
## NOTE: This wrapper uses qsubmit.pl (PBS system) instead of sbatch (SLURM)
##       It launches simuls_3.R for each dilution factor

# How This Works
#> PBS Reservation: By setting -n 16, you request 16 CPU cores on the compute node.
#  This will reserve those cores, preventing other jobs from using them and ensuring they are available for your job.
#> R Configuration: Inside your R script, mclapply with mc.cores=16 will attempt to use 16 cores
#> Memory Allocation: --mem 20 specifies the total RAM in gigabytes, not per core

# Others
# qstat -u <usuario> -n --> to check how my script is going
# qdel <jobID>          --> guess

CORES=16

# script location (in this repository)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/01_scripts"
SCRIPT="simuls_3.R"

for DILUTION in 0.00025 0.0005 0.001 0.0025 0.005 0.01 0.025 0.05 0.1
do
echo '#!/bin/bash' > launch_simuls_3.r
echo 'Rscript '$SCRIPT_DIR'/'$SCRIPT' --dilution '$DILUTION' --cores '$CORES >> launch_simuls_3.r

qsubmit.pl -q x86_64 -n $CORES --mem 20 -s launch_simuls_3
done

