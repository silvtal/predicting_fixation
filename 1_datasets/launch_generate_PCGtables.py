"""
This script generates PCG tables in batch. It consists of a loop that calls the `generate_table` function multiple
times, using different combinations of input arguments. The output tables are saved as CSV files with tab-separated
values.

Generated files are as follows:
    For each LN in LN_values (number of species; 10, 100 or 1000):
          |-- For each N in N_values[LN] (number of groups; different for each value of LN):

So, for every combination of number of total species LN and number of groups N, we will have files like :
|-- N<N>_LN.csv

Simulations (the next step):
    Each community (abundance table; see "simcomms" folder) will be paired with those PCG tables with the same number of
    species. This means multiple PCG tables per community. Each community will then be simulated with a varying number of
    groups + with a varying species distribution between groups (Skewed of Even)
    [We could also drop lognorm and use uniform communities only, or remake the lognorm ones; see Limitations]
    Each pairing will be simulated with many different dilution factors.
    Also, each abundance table actually contains 30 replicates, so everything will be simulated x30.

Limitations:
 - Leaves are randomly and evenly distributed between groups without checking their previous total abundance. If we have
   a lognormal community of 10 members and we choose to have 3 skewed groups (for instance, with final relative abundances
   of 70%, 20% and 10%), it might happen that this initial random distribution is very different from the final ones. In
   this example, group 1's 3 species might happen to sum only 11% of the total abundance, even when they have to grow to
   constitute the 70%.
      > This should not pose any serious problems.
      > However, maybe we should just use uniform communities for simplicity.
      > Which would mean crear una función nueva en el repo (predicting_fixation/1_datasets/generate_simcomms.R) similar a
      la de lognorm que ya hay (= un solo PCG) pero capaz de leer una pcg table existente.
 - No interactions yet.
      > Adding interactions means using a different script to create interaction tables. After that,
        we would need to make simulations for all combinations of:
   	  - Communities (distinct LN)
   	  - PCG tables  (with the same LN as the communities)
   	  - Interaction tables (with the same LN as the communities)
      > Furthermore, we would need to think about which kind of interactions we want (competitive, cooperative...) and
        between whom (inter-group or intra-group). In this case there would need to be more restrictions about which
        interaction table can be used for what community.

---------------------------------------------------------------------------------------------------------------------------
Input variables:
- LN_values: a list with values for the "Number of Leaves" argument to be used in the `generate_table` function.
- N_values: a list with values for the "Number of Groups" argument to be used in the `generate_table` function.
- input_folder: a string with the path to the folder containing the community files. Community files should be named
    with the format "<GroupDistr>_N<N>_<community_name>_<LN>sp.csv", where <GroupDistr> indicates the way leaves are
    distributed between groups, <N> is the number of groups, <LN> is the total number of leaves and <community_name>
    is a unique identifier for the community ("lognorm" or "uniform").
- output_folder: a string with the path to the folder where the output tables will be saved.

"""

# Import modules
import os
import random
import csv

# Set working directory
# os.chdir('/home/silvia/repos/predicting_fixation')

# Import main function
from generate_PCGtable import generate_PCGtable

# Set input variables
LN_values = [10, 100, 1000]
N_values = {
    10: [3],
    100: [3, 10],
    1000: [3, 10, 30]
}
Av_values = { 
    3: [0.404, 0.221, 0.375],
    10: [0.194, 0.048, 0.071, 0.056, 0.095, 0.052, 0.168, 0.135, 0.080, 0.102],
    30: [0.036, 0.031, 0.016, 0.009, 0.107, 0.028, 0.027, 0.027, 0.006, 0.036, 0.014, 0.025, 0.045, 0.022, 0.037, 0.034, 0.098, 0.052, 0.016, 0.019, 0.027, 0.062, 0.015, 0.007, 0.117, 0.015, 0.009, 0.043, 0.009, 0.010]
} # because I want niche sizes to be consistent; these are the "Skewed" type

av_type = ''

output_folder = 'PCGtables'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Start loop
for LN in LN_values:
    for N in N_values[LN]:
        Av = Av_values[N]
        # In case we wanted to generate random Av values, both uniformly and log-normally:
        # for distr in ['uniform', 'lognorm']:        
        #   if distr == "uniform":
        #       Av = [random.uniform(0, 1) for _ in range(N)]
        #       total = sum(Av)
        #       Av = [v/total for v in Av]
        #       av_type = 'EvenGroups'
        #   else:
        #       mu = random.uniform(-2, 2)
        #       sigma = random.uniform(0.1, 1)
        #       Av = list(random.lognormvariate(mu, sigma) for _ in range(N))
        #       total = sum(Av)
        #       Av = [v/total for v in Av]
        #       av_type = 'SkewedGroups'
                
        # Loop over community files
        outfile = f'N{N}_{LN}sp.csv'
        generate_PCGtable(N, Av, LN=LN, outfile=output_folder + '/' + outfile)
