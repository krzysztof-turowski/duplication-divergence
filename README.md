# duplication-divergence

This repository provides software tools for the following studies of duplication-divergence random graph model on real data.
* Parameter estimation for duplication-divergence model\
  Associated paper: [1] [Revisiting Parameter Estimation in Biological Networks: Influence of Symmetries](https://www.cs.purdue.edu/homes/jithinks/files/publications/parameter_estimation_2019.pdf)
* Given a single snapshot of a dynamic network, provide algorithms to infer the arrival order of the nodes.
It implements the optimal and approximate solutions of the problem. The optimal solution is a result of an integer programming formulation with coefficients found by importance sampling techniques.\
  Associated paper: [2] [Temporal Ordered Clustering in Dynamic Networks](https://arxiv.org/abs/1905.00672)

##### Table of Contents
1. [Getting started](#gettingstarted)
2.

<a name="gettingstarted" />
## Getting started
Download the repository or clone it with the following command
```bash
git clone https://github.com/krzysztof-turowski/duplication-divergence.git DIR
cd DIR
git submodule update --init --recursive
```

**Configuration of libraries**:\
Compile some of the libraries and set environment variables using the following command
```bash
bash ./configure_libs.sh LIB1 LIB2 ...
```
where LIB1, LIB2 etc correspond to
- graph libraries: *koala* (default), *snap*, *networkit*
- optimization libraries: *glpk*, *gurobi*
- automorphisms library: *nauty*

One of the graph and optimization libraries should be configured before compiling the programs. The choice of libraries to use is done via variables **GRAPH_LIB** (default: *snap*) and **LP_SOLVER** (default: *glpk*) in the Makefile.

**Note:**
- Remember to invoke `source ~/.bashrc` or update enviroment *CPLUS_INCLUDE_PATH*, *LD_LIBRARY_PATH* and *LIBRARY_PATH* manually.
- The optimization libraries *glpk* and *gurobi* require separate download and installation. The *gurobi* need a license too (free license available for academic use).

## Overview of the library

**Available main programs**

| Program        | Function           |
| ------------- |-------------|
| `dd_automorphisms`    | Find the number of automorphisms and its p-value |
| `dd_ml_estimation`      | Maximum likelihood estimation (MLE) of the parameters      |
| `dd_recurrence_estimation` | Parameter estimation using our Recurrence-relation method      |
|`dd_temporal_bound` | Algorithms to approximate the upper bounds on optimal ordering to recover using LP relaxation |
| `dd_temporal_algorithms` | Temporal ordering algorithms |


To compile any of the above `program` issue
```bash
make program
```

The syntax and examples for running the compiled programs are provided at the beginning of the associated source file (`.cpp`) in `src` folder.

In addition to the above, the following files are available to reproduce the figures of the associated papers. See the source Python files for the syntax and examples.

**Available plotting programs**

| Program        | Function           |
| ------------- |-------------|
| `dd_ml_estimation_plot`      | Plot the log-likelihood as a function of the parameters `p` and `r` of the duplication-divergence model.
| `dd_recurrence_estimation_plot` | Plot estimated parameters via degree, open triangles and traingles recurrence-relations with our methods and find the confidence interval and convergence point |
|`dd_temporal_order_plot` | Plot the temporal order algorithms for recovering the node orders (ordered cluster labels) in the precision vs density curve|

## Examples for parameter estimation by automorphisms and recurrence-relations fitting

#### `dd_automorphisms`: finding the number of automorphisms and its p-value

**Compilation**

```bash
make dd_automorphisms
```

**Find logarithm of the number of automorphisms of a given real-world network**

```bash
./dd_automorphisms -action:real_graph -graph:G-a-thaliana.txt
```
A network is required to be in edge-list format.

**Find p-value for the log of automorphisms a given real-world network, its seed and a set of parameters**

```bash
./dd_automorphisms -action:real_seed -graph:G-a-thaliana.txt -st:100 -mode:pastor_satorras -p:0.98 -r:0.49
```
This would generate from the (required) seed file `files/G0-a-thaliana.txt` 100 sample graphs with `p = 0.98` and `r = 0.49` and then compute the empirical p-value.

**Reproduce cited results**

The script that was used to obtain number of automorphisms and test p-values in [1] can be run
```bash
bash scripts/run_automorphisms.sh
```

#### `dd_recurrence_estimation`: parameter estimation using recurrence-relation method

**Compilation**

```bash
make dd_recurrence_estimation
```

**Find estimated parameters for a given real-world network and a given model**

To find estimates of the parameters of the Pastor-Satorras model against a real-world network run:
```bash
./dd_ml_estimation -action:real_data -graph:G-100-20-PS-0.1-0.3.txt -mode:pastor_satorras -st:100
```
This generates a file `temp/G-100-20-PS-0.1-0.3-PS-RE.txt` containing results in text format.
Empirical tolerance intervals are computed from 100 sample graphs, generated from a given seed graph with found parameters.

**Plot estimated parameters as parameter curves**

To plot the results from this temporary file and save it as a PDF file run:
```bash
python -B ./src/dd_recurrence_estimation_plot.py G-100-20-PS-0.1-0.3-PS-RE.txt --export pdf
```

**Reproduce cited results**

The script that was used to obtain estimated parameters and plots in [1] can be run
```bash
bash scripts/run_recurrence_estimation.sh
```

#### `dd_ml_estimation`: maximum likelihood estimation (MLE) of the parameters

**Compilation**

```bash
make dd_ml_estimation
```

**Find MLE scores for a given real-world network and a given model**

To check MLE of the parameters of the Pastor-Satorras model against a real-world network run:
```bash
./dd_ml_estimation -action:real_data -graph:G-100-20-PS-0.1-0.3.txt -mode:pastor_satorras -st:100
```
This generates a file `temp/G-100-20-PS-0.1-0.3-PS-ML.txt` containing results in text format.

**Plot MLE scores on a heat map**

To plot the MLE result from this temporary file and save it as a PDF file run:
```bash
python -B ./src/dd_ml_estimation_plot.py G-100-20-PS-0.1-0.3-PS-ML.txt --export pdf
```

**Reproduce cited results**

The script that was used to obtain MLE scores and plots in [1] can be run
```bash
bash scripts/run_ml_estimation.sh
```

## Examples for emporal ordering algorithms

#### `dd_temporal_bound`: upper bound on temporal ordering

**Compilation**

```bash
make dd_temporal_bound
```

**Example runs**

To plot the exact upper bound from LP relaxation, averaged over 100 random graphs on 13 vertices, generated from Pastor-Satorras model run for example
```bash
./dd_temporal_bound -algorithm:exact -n:13 -n0:6 -mode:pastor_satorras -p:0.3 -r:0.6 -p0:0.6 -gt:100
```

To plot one of the example approximation algorithms (*local-unif-sampling*) based on 1000 permutation samples run
```bash
./dd_temporal_bound -algorithm:local-unif-sampling -n:13 -n0:6 -mode:pastor_satorras -p:0.3 -r:0.6 -p0:0.6 -gt:100 -st:1000
```

The results are **appended** to file `temp/synthetic-13-6-PS-0.300-0.60-TC.txt` containing data in text format.

#### `dd_temporal_algorithms`: temporal ordering algorithms

**Compilation**

```bash
make dd_temporal_algorithms
```

**Greedy algorithms**

```bash
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -algorithm:sort_by_degree -gt:100
```
The command run on 100 synthetic graphs on 50 vertices generated from Pastor-Satorras model with `p = 0.6`, `r = 1.0` and a random Erdos-Renyi seed graph on 10 vertices with edge probability 0.6. It returns the order given by *sort-by-degree* algorithm (see [2]).

The results are **appended** to file `temp/synthetic-50-10-PS-0.300-1.00-TA.txt` containing data in text format.

**Unsupervised learning**

```bash
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:p_uv_threshold -gt:100 -st:100000 -threshold:0.7
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:sort_by_p_uv_sum -gt:100 -st:100000 -binsize:1
```
First command executes the *$p_{uv}$-threshold* algorithm with $\tau = 0.7$, second executes *sort-by-$p_{uv}$-sum* algorithm with bin size 1 (for details, see [2]).

**Supervised learning**

```bash
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:p_uv_threshold -gt:100 -st:100000 -threshold:0.7 -perfect:0.001
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:sort_by_p_uv_sum -gt:100 -st:100000 -binsize:1 -perfect:0.001
```
These are exactly the same algorithms as above, but this time each uses 0.1% of the total number of recognizable pairs from the original order.

#### `dd_temporal_order_plot`: plot temporal ordering results

Suppose we have both temporary files `G-test-TC.txt` (containing data related to the temporal bounds) and `G-test-TA.txt` (containing results of various algorithms).
To plot results from files and save it as a PDF file run:
```bash
python -B ./src/dd_ml_estimation_plot.py G-test --export pdf
```
If you want to plot not only the mean result, but also the distribution of the results for each algorithm add `--detailed` option.
In case any of the files is missing, the script will report it in a command line, but plot the existing data nevertheless.

**Reproduce cited results**

The scripts that were used to find all the data collected in [2] are as following:
```bash
bash scripts/run_temporal_curve.sh
bash scripts/run_temporal_reinforced.sh
```

## Data

All the data files we have used in our experiments are provided in `files` folder: arXiv, Simple English Wikipedia and CollegeMsg networks, and the protein-protein interaction (PPI) networks.
We collect the PPI data of seven species from [BioGRID](https://thebiogrid.org/).
The age of the proteins are based on phylogentic tree data gathered from [ProteinHistorian](http://lighthouse.ucsf.edu/ProteinHistorian/) and are also provided.
