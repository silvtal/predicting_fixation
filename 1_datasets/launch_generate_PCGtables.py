"""
This script generates PCG tables in batch. It consists of a loop that calls the `generate_table` function multiple
times, using different combinations of input arguments. The output tables are saved as CSV files with tab-separated
values.

Generated files are as follows:
    For each LN in LN_values (number of species; 100 or 1000):
          |
          |-- For each N in N_values[LN] (number of groups; 3 or 10):
          |        |
          |        |-- For each niche size distribution (Skewed or Even):
          |        |         |
          |        |         |-- <niche_size_distr>_N<N>_LN.csv

So, for every combination of number of total species LN and number of groups N, we will have files like :
|-- EvenGroups_N<N>_LN.csv
|-- SkewedGroups_N<N>_LN.csv


Simulations (the next step):
    Each community (abundance table; see "simcomms" folder) will be paired with those PCG tables with the same number of
    species. This means multiple PCG tables per community. Each community will then be simulated with a varying number of
    groups + with a varying species distribution between groups (Skewed of Even)
    Each pairing will be simulated with many different dilution factors.
    Also, each abundance table actually contains 30 replicates, so everything will be simulated x30.

Limitations:
 - Leaves are randomly and evenly distributed between groups without checking their previous total abundance. This is
   different to the distribution of the niche sizes -- here we refer to the richness of each niche. If we have
   a lognormal community of 10 members and we choose to have 3 skewed groups (for instance, with final relative abundances
   of 70%, 20% and 10%), it might happen that this initial random distribution is very different from the final ones. In
   this example, group 1's 3 species might happen to sum only 11% of the total abundance, even when they have to grow to
   constitute the 70%.
      > This should not pose any serious problems.
      > In any case, we use initial communities with 100 and 1000 species to avoid extreme cases like this.
 distribution communities. We are working with lognormal communities.
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
    is a unique identifier for the community ("Skewed" or "Even").
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
LN_values = [100, 1000]
N_values = { # we could make them different depending on size if we wanted
    100: [3, 10],
    1000: [3, 10]
}
Av_values = {# because I want niche sizes to be consistent
  "Skewed" : { # these are the "Skewed" type
    # 3: [0.404, 0.221, 0.375],
    3: [0.60, 0.30, 0.10],
    # 10: [0.194, 0.168, 0.135, 0.102, 0.095, 0.080, 0.071, 0.056, 0.052, 0.048], # Una log-normal
    # 10: [0.30, 0.20, 0.15, 0.12, 0.09, 0.07, 0.04, 0.02, 0.009, 0.001]          # Daniel
    10: [0.3, 0.2, 0.15, 0.1, 0.085, 0.06, 0.05, 0.035, 0.015, 0.005]
  },
  "Even" : { # these are the "Even" type
    3: [0.3333, 0.3333, 0.3334],
    10: [0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10] 
  } 
}

output_folder = 'PCGtables'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Start loop
for LN in LN_values:
    for N in N_values[LN]:
        for distr in ['Even', 'Skewed']:
          Av = Av_values[distr][N]
          if distr == "Even":
            av_type = 'EvenGroups'
          else:
            av_type = 'SkewedGroups'
                
          # Loop over community files
          outfile = f'{distr}Niches_N{N}_{LN}sp.csv'
          generate_PCGtable(N, Av, LN=LN, outfile=output_folder + '/' + outfile)
