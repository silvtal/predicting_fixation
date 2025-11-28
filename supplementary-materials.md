# Supplementary Materials: Detailed implementation of the dilution-growth simulation pipeline

## Overview

This document provides a **detailed technical description of the computational implementation** of the dilution-growth simulation pipeline used in "Computational simulation of ecological drift for generating functional minimal microbiomes identifies key experimental and biotic factors influencing success" (Talavera-Marcos & Aguirre de Cárcer).

All the scripts used for simulation and analysis are available at the `predicting_fixation` Github repository (!["github.com/silvtal/predicting_fixation"](github.com/silvtal/predicting_fixation)). The simulations are implemented using the `dilgrowth` R package (!["github.com/silvtal/dilgrowth"](github.com/silvtal/dilgrowth)), which uses C++ functions (via Rcpp) for computationally intensive operations.

There are three general steps to the pipeline: generation of initial communities, simulation of dilution-growth cycles and results analysis.

## Generation of initial communities

The first step of the pipeline involves generating initial community compositions that serve as starting points for the dilution-growth simulations. These communities are encoded as abundance tables where each column represents a distinct community with a specific combinations of parameters (richness, community size, abundance distribution).

The abundance tables are generated using the `generate_simcomms.R` script (available in the `predicting_fixation` repository under `1_datasets/`). This script runs two functions, one for communities with uniform distribution and another one for log-normal distributions. For uniform distributions, abundances are sampled with equal probability for all species. For log-normal distributions, abundances are sampled from a log-normal distribution (`meanlog = 0`, `sdlog = 1`) to create realistic skewed abundance patterns. The rest of the parameters are given as function arguments.

```
simulate_uniform_community <- function(num_species=10,
                                       comm_size=10**4,
                                       name="otu")
                                       {
                                        ...
                                       }
simulate_log_normal_community <- function(num_species=10,
                                          comm_size=10**4,
                                          name="otu",
                                          meanlog = 0,
                                          sdlog = 1)
                                          {
                                            ...
                                          }
```
Each parameter combination (richness: 10, 100, or 1000 species; community size: 10^4^ or 10^6^ individuals; distribution: uniform or log-normal) generates 30 replicate communities, resulting in 12 unique community types. The resulting abundance tables are available under `1_datasets/simcomms/`.

### Generation of functional group tables

Scenarios 2 and 3 involve exploring communities with multiple functional groups. These groups are defined by tables that dictate their relative niche sizes and their assigned populations.

The functional group tables are generated using the `generate_groups_table.py` function and the `launch_generate_groups_tables.py` script (available in the `predicting_fixation` repository under `1_datasets/`). This script creates tables where each row represents a functional group with a specific niche size (carrying capacity) and a list of populations assigned to that group.

The main function `generate_groups_table()` takes the following parameters:

```python
def generate_groups_table(N, Av, L=None, LN=None, outfile=None, seed=None):
    """
    Generate a table with N rows, where each row represents a group with an
    average value and a number of leaves.
    
    Args:
        N (int): Number of groups.
        Av (list of floats): List of niche sizes (size for each group). The length
            of the list must be equal to N, and the values must sum to 1.
        L (list of ints, optional): List of number of leaves for each group.
        LN (int, optional): Total number of leaves. If provided, the number of
            leaves for each group will be randomly generated.
        outfile (str, optional): Output file path.
        seed (int, optional): Random seed.
    """
```

The script `launch_generate_groups_tables.py` generates group tables in batch for different parameter combinations. For each combination of total number of species (LN: 100 or 1000) and number of groups (N: 3 or 10), it creates tables with either "Even" or "Skewed" niche size distributions. For Even distributions, all groups have equal niche sizes (e.g., 0.3333, 0.3333, 0.3334 for 3 groups). For Skewed distributions, groups have unequal niche sizes (0.60, 0.30, 0.10 for 3 groups and 0.3, 0.2, 0.15, 0.1, 0.085, 0.06, 0.05, 0.035, 0.015, 0.005 for 10 groups). Populations are randomly assigned to groups, with each group receiving approximately equal numbers of populations. The output tables are saved as .csv files. These tables are available in the `predicting_fixation` repository under `1_datasets/functional_groups_tables/`. 

## Simulation of dilution-growth processes

The second step of the pipeline involve performing the simulations themselves. All simulations were run using two scripts, both available in the `predicting_fixation` repository under `2_simulations/01_scripts`: `dilgrowth_1_2.R` for scenarios 1 and 2 and `simuls_3.R` for scenario 3. In this section, we will explain in detail the implementation of the simulations within these scripts.

![](process_fig.svg)

### Main function: `simulate_timeseries()`

The core of the simulation pipeline is the `simulate_timeseries()` function, which orchestrates the entire dilution-growth process. This function is implemented in the `dilgrowth` R package (github.com/silvtal/dilgrowth) and is called for each independent trajectory (replicate simulation) to manage the iterative cycle of dilution and growth.

In order to run multiple independent trajectories simultaneously, this function is called within a parallelized loop using `mclapply()`:

```r
abund_temp <- mclapply(X = 1:no_of_simulations,
                       FUN = function(iter) {
                         trajectory <- simulate_timeseries(
                           counts,
                           carrying_capacities = carrying_capacities,
                           interactions = interactions,
                           logistic = logistic,
                           dilution = dilution,
                           no_of_dil = no_of_dil,
                           fixation_at = fixation_at,
                           abun_total = abun_total,
                           growth_step = growth_step,
                           is_growth_step_a_perc = is_growth_step_a_perc,
                           keep_all_timesteps = save_all,
                           allow_group_extinctions = allow_group_extinctions,
                           force_continue = FALSE
                         )
                         return(trajectory)
                       }, mc.cores = cores)
```

#### Arguments

Arguments are as follow:

- **`counts`**: Initial abundance vector (named vector where names are species/population identifiers). This represents the starting community composition before any dilution-growth cycles, and corresponds with one column of the abundance tables created in step (1).
  - Internally, this vector is separated into a population identity vector (names) and an abundance vector. The latter will be transformed into a growth probability vector, which can be modified by interactions if an interaction matrix is present.
- **`carrying_capacities`**: Named vector indicating the carrying capacity (niche size) for each population. The names correspond to functional group assignments. For example:

    ```r
    carrying_capacities <- c(7000, 7000, 7000, 2000, 2000, 1000)
    names(carrying_capacities) <- c("Group1", "Group1", "Group1", 
                                    "Group2", "Group2", "Group3")
    ```

    This defines:
    - Group1: 3 populations, each with carrying capacity 7000 (total: 21000)
    - Group2: 2 populations, each with carrying capacity 2000 (total: 4000)
    - Group3: 1 population, with carrying capacity 1000 (total: 1000)

    This vector is extracted from the functional group tables generated in step (1). If `NULL`, simulations proceed without separate functional groups (scenario 1).

- **`interactions`**: Optional interaction matrix (A) where `A[i,j]` represents the effect of population j on population i. In Scenario 3, interaction matrices are generated on the go in the script `simuls_3.R` using the function `create_interact_table()`. Inside a nested loop, this function is called 100 times for each possible parameter combination (sign: positive, negative, or both; absolute value: 1, 0.1, 0.01, 0.001, or 0.0001; frequency: 1%, 10%, or 100% of possible inter-group pairs), once for each replicate. In other words, each replicate uses a different interactions matrix, but with the same parameters.
- **`dilution`**: Dilution factor (D).
- **`no_of_dil`**: Number of dilution-growth cycles to simulate.
- **`fixation_at`**: Relative abundance threshold (0-1) above which a population is considered "fixed".
- **`abun_total`**: Target total community abundance to reach after each growth phase.
- **`growth_step`**: Fixed percentage of the total community size that increases in each growth iteration / number of individuals that grow per iteration (if `is_growth_step_a_perc = FALSE`).
- **`is_growth_step_a_perc`**: Boolean indicating whether `growth_step` is a fixed value or a percentage of current community size. FALSE by default. In both cases, there will be as many iterations as needed until `abun_total` is reached by the community. 
- **`logistic`**: Boolean indicating whether to use logistic growth instead of fixed-step growth.
- **`keep_all_timesteps`**: Boolean indicating whether to save intermediate states (final population abundances) at each cycle.
- **`allow_group_extinctions`**: Boolean indicating whether simulations should continue if functional groups go extinct.
  - **`allow_group_extinctions = FALSE`**: Simulation stops with an error if any group goes extinct.
  - **`allow_group_extinctions = TRUE`**: Simulation continues, but a warning message is printed and the extinct group is recorded.
- **`force_continue`**: Boolean indicating whether to continue simulation if dilution would leave fewer than 1 individual. If `FALSE`, the simulation stops with an error. If `TRUE`, one random individual is selected to continue the simulation.

#### Structure

The function implements a `while` loop that continues until either (1) the maximum number of dilution cycles (`no_of_dil`) is reached, or (2) fixation is achieved in all functional groups (or the entire community if no groups are defined). Each iteration of the loop is an entire dilution-growtg cycle and consists of:

1. **Dilution step**: The current community is sampled proportionally to create a diluted community.
2. **Growth phase**: The diluted community grows back to `abun_total` through multiple iterations within a nested loop. A C++ function is called once per growth iteration (more details in the _Core growth functions_ section below).
3. **Fixation check**: The function checks whether fixation has occurred.
4. **State saving**: If `keep_all_timesteps = TRUE`, the current state is saved.

**!!** Each call to `simulate_timeseries()` is a whole dilution-growth process with multiple cycles. Each iteration of the `while` loop in `simulate_timeseries()` is a single cycle, and each call to the C++ growth functions corresponds to a single growth iteration within a cycle.

##### 1. Dilution

The dilution step reduces the community size by the dilution factor, D (encoded as `dilution`).

```r
this_timestep <- table(names(this_timestep)[.Internal(sample(
  x = length(this_timestep),
  size = sum(this_timestep)*dilution,
  replace = TRUE,
  prob = this_timestep))])
```
where `this_timestep` is the abundance vector of the whole community at each iteration (originally the `counts` argument).

This samples `sum(this_timestep)*dilution` individuals **with replacement**, where the probability of selecting each individual is proportional to its current abundance.

$$N_{\text{total}}(t,0) = \frac{N_{\text{total}}(t-1,\text{end})}{D}$$

where $N_{\text{total}}(t,0)$ is the total community abundance at the start of cycle $t$ (after dilution) and $N_{\text{total}}(t-1,\text{end})$ is the total community abundance at the end of the previous cycle (and at the start of the simulation). The result is a new abundance vector.


##### 2. Growth

Growth happens via C++ (Rcpp) functions. The growth phase operates differently depending on whether functional groups are defined. When `carrying_capacities = NULL` and therefore there are no separate groups, `simulate_timeseries()` uses the `growth_one_group()` helper function. Else, it uses the more general function `growth()` (or `growth_log()` if `logistic == TRUE`).

- Before entering the growth loop, `simulate_timeseries()` checks for extinct groups and handles them according to the `allow_group_extinctions` parameter.


##### 3. Fixation check

It's done via the `check_for_fixation()` R function, which determines whether fixation has occurred in the community.

- Without functional groups, fixation occurs when the relative abundance of any single population exceeds `fixation_at`:
  ```r
  fixated <- max(this_timestep)/sum(this_timestep) >= fixation_at
  ```

- With functional groups, fixation is checked separately for each group:
  ```r
  total_abundances <- aggregate(abundances ~ groups, df, sum)
  max_abundances <- aggregate(abundances ~ groups, df, max)
  fixated <- (max_abundances$abundances / total_abundances$abundances) >= fixation_at
  ```

`check_for_fixation()` returns a named logical vector indicating which groups have reached fixation. `simulate_timeseries()` will continue the simulations until all groups reach fixation, or until `no_of_dil` cycles have been completed.


### Core Growth Functions (C++ Implementation)

The computationally intensive growth operations are implemented in C++ for efficiency. All growth and supplementary functions are located in the `dilgrowth` repository as `src/growth.cpp`.

#### Setting the magnitude of the abundance increase: `check_step()`

First, the magnitude of growth ($I$) is defined via the helper function `check_step()`. This function is called before each growth iteration ($k$) within `simulate_timeseries()` to dynamically adjust the growth step.

If `growth_step` was defined as a percentage, then:

$$I_k = \text{growth_step} \times N_{TOTAL}(k)$$

This $I_k$ will be used by the next call to the growth functions.

   ```cpp
   step = trunc(growth_step * sum(this_timestep));
   step = max(1, step); // Never less than 1
   ```

- `check_step()` ensures that the total community size (or the maximum carrying capacity for each group) won't be exceeded.

   ```cpp
   if (trunc((sum(this_timestep)) + growth_step) > abun_total) {
     step = (abun_total - sum(this_timestep));
   }
   ```

- It also ensures that there's a growth of at least 1 unit.
   ```cpp
   else if (sum(this_timestep) < growth_step) {
     int half = trunc(sum(this_timestep)/2);
     step = std::max(half, 1);
   }
   ```

### Main growth function: `growth()` 

**!!** Each call to this function executes one single iteration inside a cycle.

When there are no functional groups or interactions, the probability for each population $i$ to be randomly selected for growth is equivalent to its relative abundance. If $N_i$ is the number of individuals of $i$, the probability $p_i(k)$ that an additional individual is assigned to population $i$ in iteration $k$ is:

$$P_i(k)={{N_i(k)}\over{N_{total}(k)}}$$

If functional groups are present, factors such as the number of functional groups, the assignment of populations to these groups and the relative abundance of each group (i.e. its niche size) need to be taken into account. As such, in scenarios with defined functional groups (scenarios 2 and 3), random selection occurs within each functional group, and all groups grow separately until reaching their respective niche sizes or carrying capacity ($K$). Let $G$ be the functional group to which population $i$ belongs. The conditional selection probability for population $i$, given the $K_G$ assigned to its group, is proportional to its relative abundance **within that group**:

$$P_{i\mid{G}}(k)={{p_i(k)}\over{\sum_{j \in G} p_j(k)}}$$

Where $N_i(k)$ is the abundance of population i, and $\sum_{j\in{G}}N_j(k)$ is the total abundance of group G. In order to keep the niche sizes stable over iterations, the magnitude of growth of each group, I_{k|G}, is weighted by its carrying capacity:

$$I_{ k\mid{G} } = I_k \times K_G$$

The increment for each selected individual is weighted by the group's relative carrying capacity. This ensures that groups with larger niche sizes grow proportionally.

Biotic interactions (present when `interactions` argument in `simulate_timeseries()` is not `NULL`) modify these probabilities as well. The provided interaction matrix ($A$) is multiplied by the relative abundances vector (which includes the $p_i$ values), and the result is added back to the original probability vector. The modified selection probability is:

$$p'(k)=p(k)+ A · p(k)$$

For each population $i$, the growth probability is calculated as the sum of multiple components: its own relative abundance and the relative abundances of the populations it interacts with, weighted by the strength (positive or negative) of each interaction:

$$p'_i(k)=p_i(k)+\sum_j{A_{ij} p_j(k)}$$

Where $A_{ij}$ is the effect of population j on population $I$ (i.e., the corresponding value in the interaction matrix. The process continues by sampling within the functional groups. Thus, the conditional selection probability for population $i$ to be picked for growth is:

$$P'_{i\mid{G}}(k)={{p'_i(k)}\over{\sum_{j \in{G}} p'_j(k)}}$$

The `growth()` function operates as follows:

1. Initialize probabilities for all populations using their relative abundances `prob = x/sum((x))`.
2. If interactions are provided, modify probabilities by multiplying by the interaction matrix: `prob = wrap(as<arma::vec>(prob) + (as<arma::mat>(interactions.get()) * as<arma::vec>(prob)))`. This implements the equation: **p'(k) = p(k) + A · p(k)**.
  - If all probabilities become zero after applying interactions, the function throws an error, as growth would be impossible.
3. Process each functional group separately. For each group $G$:
   - Extract probabilities for populations in group $G$. This is implemented by setting `group_prob[i] = prob[i]` for populations in group G, and `group_prob[i] = 0` for populations not in group G, then normalizing.
   - Calculate the number of individuals that can grow in this group (limited by current group abundance).
   - Sample `step` individuals from group $G$ using `pick_new_bugs()` (see section below) with group-specific probabilities.
   - For each selected individual $i$ in group $G$, increment abundance by: $1 × (K_G / \sum{K})$
     where $K_G$ is the carrying capacity of group $G$ and $\sum K$ is the sum of all carrying capacities.


#### Sampling Function (`pick_new_bugs()`)

This is a helper function used by all `growth` functions to randomly select which populations will grow in each iteration.

```cpp
NumericVector pick_new_bugs(NumericVector arr,
                            int size,
                            bool replace,
                            NumericVector prob)
```

Arguments are as follow:

- **`arr`**: Vector of indices [0, 1, 2, ..., n-1] where n is the number of populations. Populations can have various names when defined as the named vector `counts` (`simulate_timeseries()` argument), but `pick_new_bugs` will work with a simplified numeric vector.
- **`size`**: Number of individuals to sample (the aforementioned modified `growth_step`)
- **`replace`**: Always `true` (sampling with replacement)
- **`prob`**: Probability vector where `prob[i]` is the probability of selecting population i

This function uses `RcppArmadillo::sample()` to perform weighted random sampling and returns a vector of indices indicating which populations were selected for growth.


### Growth without functional groups: `growth_one_group()`

This function simulates growth when no functional groups are defined. It's a simplified version of `growth()`. It performs the following steps:

1. Calculate initial growth probabilities as relative abundance, according to the aforementioned formula: `prob = this_timestep/sum((this_timestep))`.
2. Apply interactions if provided.
3. Ensure no negative probabilities (set to 0 if negative).
4. Sample `growth_step` individuals using `pick_new_bugs()`. This creates an index vector, `new_bugs`.
5. For each index (selected population) in `new_bugs`, increment the corresponding abundance by 1. This in-place indexed update is computationally more efficient than constructing a dense increment vector and adding it to `this_timestep`.


### Supporting Functions

#### `roundVectorPreservingSum()`: Abundance Rounding

This function rounds abundance values to integers while preserving the total sum, which is essential for biological realism (populations cannot have fractional individuals). When `carrying_capacities` is provided, rounding is performed separately for each functional group to preserve group-specific sums.

First, it rounds all values (`rounded_vector <- round(vector)`). Then, it computes the difference (`diff_sum <- sum(original) - sum(rounded)`).

- If `diff_sum > 0`: Randomly select `diff_sum` individuals and increment by 1
- If `diff_sum < 0`: Randomly select `|diff_sum|` individuals (with abundance > 0) and decrement by 1

#### Logistic growth: `growth_log()`

This function implements logistic growth dynamics, where growth rates change as populations approach their carrying capacities. It was not used in the present work, but it works similarly to `growth()`, implementing a single iteration (growth step within a cycle) per call.

For each group G:
1. Calculates current group abundance.
2. For each population $i$ in group $G$, the function determines growth with `rbinom(1, 1, prob[i])`, computes the logistic growth step (`2 * (1 - sum_group / carrying_capacities[i])`) if selected, and updates the abundance accordingly.

**Key differences from `growth()`:**
- No fixed `growth_step` parameter.
- Growth magnitude is determined by a logistic function: `growth_step = 2 × (1 - sum_group / K_G)`
- Uses binomial sampling (`rbinom()`) instead of weighted sampling.
- Growth rate decreases as the group approaches its carrying capacity.
- If a group exceeds its carrying capacity, populations in that group can decrease in abundance.

---

There's also a `full_growth()` C++ function defined on the `dilgrowth` packages which triggers a loop in a similar but simpler way than `simulate_timeseries()`, finishing a whole cycle by calling `growth_one_group()`, `growth()` or `growth_log()`. This function, along with `growth_log()`, was not used for this work.

<!--TODO include?
## Computational Efficiency Considerations

### Parallelization

Multiple independent trajectories are run in parallel using `mclapply()` from the `parallel` package. Each trajectory is completely independent.

### C++ Implementation

The growth functions are implemented in C++ because the growth loops (which can run hundreds of times per cycle) benefit from compiled code

### Memory Management

- When `keep_all_timesteps = FALSE`, only the final state is stored, minimizing memory usage
- When `keep_all_timesteps = TRUE`, a data frame is pre-allocated with `(no_of_dil + 1)` rows
- Intermediate states during growth are not stored (only states after each complete dilution-growth cycle)
-->

## Analysis

The third step of the pipeline involves analyzing the simulation results to identify key factors that influence fixation success. All analysis scripts are available in the `3_analysis` folder of the `predicting_fixation` repository. Scripts are prefixed with `S1_`, `S2_`, or `S3_` to indicate which scenario they analyze (scenario 1: single communities without functional groups; scenario 2: communities with functional groups; scenario 3: communities with functional groups and explicit interactions).

### Scenario 1: Single-community analysis

#### Exploratory visualization (`S1_0__*` scripts)

- `S1_0__plot_target_variable_space.py`: Visualizes the distribution of the target variable (number of cycles until fixation) across different parameter combinations.
- `S1_0__plot_feature_effects.R`: Examines the relationship between individual features (diversity metrics, community size, dilution factor) and fixation.
- `S1_0__plot_dilfactor_effect.py`: Specifically explores how dilution factor affects fixation.

#### Predictive modeling

**Random Forest models (`S1_1__RF_success.R`)**: This script defines the target variable "success" as the dilution-growth cycle number in which one OTU reaches fixation (meaning it has achieved a relative abundance of 50% or 90% of the total community size).A series of random forest models are used to quantify the importance of each variable in predicting the number of cycles required for fixation.

**Generalized Linear Models (`S1_2__GLM_failure.R`)**: This script defines the target variable "failure" as a binary indicator (1 = fixation did not occur before cycle N, 0 = fixation occurred before cycle N). Multiple models are generated for different cycle thresholds (N values ranging from 10 to 1000).

All figures from scenario 1 analyses are saved to the `figures/` folder.

### Scenario 2: Multi-group analysis

#### Data visualization

The main scripts are:

- The `S2_g_plot_results_facet_by_dil*.r` parse simulation outputs (performing diversity calculations and recording group-level metrics) and generate preliminary visualizations.

- `S2_0g__plot_fixation_fixationPERgroups_gr_PREPARSED.R`: Calculates and visualizes fixation success rates per group.

- `S2_0g__plot_extinction_extinctionPERgroups_gr_PREPARSED.R`: Calculates and visualizes extinction rates per group.

All figures from scenario 2 analyses are saved to the `figures_groups/` folder.

### Scenario 3: Interactions analysis

#### Interaction effects (`S3_1g-i__*` scripts)

- `S3_1g-i__plot_interactions_successPERgroups_PREPARSED.R`: Examines how interactions affect success rates (fixation) within each functional group.
- `S3_1g-i__plot_interactions_extinctionPERgroups_PREPARSED.R`: Examines how interactions affect extinction rates within each functional group.

All figures from scenario 3 analyses are saved to the `figures_groups_INTERACTIONS/` folder.
