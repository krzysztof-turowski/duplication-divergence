# duplication-divergence

This repository provides software tools for the following studies of duplication-divergence random graph model on real data.
* Parameter estimation for duplication-divergence model\
  Associated paper: [1] [Revisiting Parameter Estimation in Biological Networks: Influence of Symmetries](https://www.biorxiv.org/content/10.1101/674739v1)
* Given a single snapshot of a dynamic network, provide algorithms to infer the arrival order of the nodes.
It implements the optimal and approximate solutions of the problem. The optimal solution is a result of an integer programming formulation with coefficients found by importance sampling techniques.\
  Associated paper: [2] [Recovery of Vertex Orderings in Dynamic Graphs: Algorithms for a General Solution](https://arxiv.org/abs/1905.00672)

#### Table of Contents
1. [Getting started](#getting-started)
2. [Overview of the library](#overview)
3. [Examples for parameter estimation by automorphisms and recurrence-relations fitting](#examples-estimation)
4. [Examples for temporal ordering algorithms](#examples-temporal-order)
5. [Data](#data)

## <a name="getting-started"></a>Getting started
Download the repository or clone it with the following command
```bash
git clone https://github.com/krzysztof-turowski/duplication-divergence.git DIR
cd DIR
git submodule update --init --recursive
```

**Configuration of libraries**:\
Set environment variables and compile libraries that reside inside the `lib` folder using the following command
```bash
bash ./configure_libs.sh LIB1 LIB2 ...
```
where *LIB1*, *LIB2* etc correspond to.
- graph libraries: *koala*, *snap*, *networkit*,
- LP optimization libraries: *glpk*, *gurobi*,
- automorphisms library: *nauty*.

One of the graph and optimization libraries should be configured before compiling the programs in this repo. The choice of libraries to use can be changed via variables **GRAPH_LIB** (default: *snap*) and **LP_SOLVER** (default: *glpk*) in the Makefile.

**Note:**
- Remember to invoke `source ~/.bashrc` or update enviroment variables *CPLUS_INCLUDE_PATH*, *LD_LIBRARY_PATH* and *LIBRARY_PATH* manually after running the bash script.
- The optimization libraries *glpk* and *gurobi* require separate download and installation. The *gurobi* need a license too (free license available for academic use).

## <a name="overview"></a>Overview of the library

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

The syntax and usage for running the compiled programs are provided at the beginning of the associated source files (`.cpp`) in `src` folder.

In addition to the above, the following files are available to reproduce the figures of the associated papers. See the Python source files for the syntax and usage.

**Available plotting programs**

| Program        | Function           |
| ------------- |-------------|
| `dd_ml_estimation_plot`      | Plot the log-likelihood as a function of the parameters `p` and `r` of the duplication-divergence model.
| `dd_recurrence_estimation_plot` | Plot estimated parameters via degree, open triangles and traingles recurrence-relations with our methods and find the confidence interval and convergence point |
|`dd_temporal_order_plot` | Plot the temporal order algorithms for recovering the node orders (ordered cluster labels) in the precision vs density curve|

## <a name="examples-estimation"></a>Examples for parameter estimation by automorphisms and recurrence-relations fitting

#### `dd_automorphisms`: finding the number of automorphisms and its p-value

Compilation

```bash
make dd_automorphisms
```

Find logarithm of the number of automorphisms of a given real-world network (edge-list format):

```bash
./dd_automorphisms -action:real_graph -graph:G-a-thaliana.txt
```
A network is required to be in edge-list format.

Find p-value of the log of automorphisms a given real-world network, its seed graph and the set of graph parameters:

```bash
./dd_automorphisms -action:real_seed -graph:G-a-thaliana.txt -st:100 -mode:pastor_satorras -p:0.98 -r:0.49
```
This will first generate 100 random graphs from the seed graph `files/G0-a-thaliana.txt` with `p = 0.98` and `r = 0.49`, and then will compute the empirical p-value.

Reproduce the results in [1] to obtain the number of automorphisms and p-values of the graph examples:

```bash
bash scripts/run_automorphisms.sh
```

### `dd_recurrence_estimation`: parameter estimation using recurrence-relation method

Compilation

```bash
make dd_recurrence_estimation
```

Find estimated parameters for a given real-world network and a given model: to find estimates of the parameters of the fitted Pastor-Satorras model for a given real-world network.
```bash
./dd_ml_estimation -action:real_data -graph:G-100-20-PS-0.1-0.3.txt -mode:pastor_satorras -st:100
```
This generates a file `temp/G-100-20-PS-0.1-0.3-PS-RE.txt` containing results in text format. Empirical tolerance intervals are computed from 100 sample graphs, generated from a given seed graph with estimated parameters.

Plot estimated parameters from the outputted above temporary file and save it as a PDF file:
```bash
python -B ./src/dd_recurrence_estimation_plot.py G-100-20-PS-0.1-0.3-PS --export pdf
```

Reproduce cited results in [1] to obtain estimated parameters for various graphs:
```bash
bash scripts/run_recurrence_estimation.sh
```

#### `dd_ml_estimation`: maximum likelihood estimation (MLE) of the parameters

Compilation

```bash
make dd_ml_estimation
```

Find log-likelihood scores of various parameters for a given real-world network with Pastor-Satorras model
```bash
./dd_ml_estimation -action:real_data -graph:G-100-20-PS-0.1-0.3.txt -mode:pastor_satorras -st:100
```
This generates results in a text file `temp/G-100-20-PS-0.1-0.3-PS-ML.txt`.

Plot MLE scores on a heat map from the above temporary file and save it as a PDF file:
```bash
python -B ./src/dd_ml_estimation_plot.py G-100-20-PS-0.1-0.3-PS-ML.txt --export pdf
```

Reproduce cited results in [1] of the MLE scores and plots.
```bash
bash scripts/run_ml_estimation.sh
```

## <a name="examples-temporal-order"></a>Examples for temporal ordering algorithms

#### `dd_temporal_bound`: upper bound on temporal ordering

Compilation
```bash
make dd_temporal_bound
```

Plot the exact upper bound of the optimization using linear programming (LP) relaxation, averaged over 100 random graphs with 13 vertices generated from Pastor-Satorras model:
```bash
./dd_temporal_bound -algorithm:exact -n:13 -n0:6 -mode:pastor_satorras -p:0.3 -r:0.6 -p0:0.6 -gt:100
```

Plot one of the approximation algorithms (*local-unif-sampling*) based on 1000 permutation samples:
```bash
./dd_temporal_bound -algorithm:local-unif-sampling -n:13 -n0:6 -mode:pastor_satorras -p:0.3 -r:0.6 -p0:0.6 -gt:100 -st:1000
```
The results are *appended* to a text file `temp/synthetic-13-6-PS-0.300-0.60-TC.txt`.

#### `dd_temporal_algorithms`: temporal ordering algorithms

Compilation
```bash
make dd_temporal_algorithms
```

Greedy algorithms:
```bash
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -algorithm:sort_by_degree -gt:100
```
The command run on 100 synthetic graphs with 50 vertices generated from the Pastor-Satorras model having paremeters `p = 0.6`, `r = 1.0` and a random Erdos-Renyi seed graph with 10 vertices and edge probability 0.6. It returns the order given by *sort-by-degree* algorithm (see [2]).

The results are appended to a text file `temp/synthetic-50-10-PS-0.300-1.00-TA.txt`.

**Unsupervised learning** solution:
```bash
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:p_uv_threshold -gt:100 -st:100000 -threshold:0.7
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:sort_by_p_uv_sum -gt:100 -st:100000 -binsize:1
```
The first command executes the *$p_{uv}$-threshold* algorithm with $\tau = 0.7$. The second command executes *sort-by-$p_{uv}$-sum* algorithm with bin size 1 (for details, see [2]).

**Supervised learning** solution:
```bash
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:p_uv_threshold -gt:100 -st:100000 -threshold:0.7 -perfect:0.001
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:sort_by_p_uv_sum -gt:100 -st:100000 -binsize:1 -perfect:0.001
```
These are exactly the same algorithms as above, but each uses 0.1% of the original order for training the model.

#### `dd_temporal_order_plot`: plot temporal ordering results

Suppose we have both temporary files `G-test-TC.txt` (containing data related to the temporal bounds) and `G-test-TA.txt` (containing results of various algorithms). To plot results from these files and save it as a PDF, run:
```bash
python -B ./src/dd_ml_estimation_plot.py G-test --export pdf
```
If you wish to plot not only the average result, but also the distribution of the results for each algorithm add `--detailed` option. In case any of the files is missing, the script will report it in the command prompt, but plots rest of the data.

Reproduce cited results in [2]:
```bash
bash scripts/run_temporal_curve.sh
bash scripts/run_temporal_reinforced.sh
```

## <a name="data"></a>Data

All the data files we have used in our experiments are provided in `files` folder: arXiv, Simple English Wikipedia and CollegeMsg networks, and the protein-protein interaction (PPI) networks.
We collect the PPI data of seven species from [BioGRID](https://thebiogrid.org/).
The age of the proteins are based on phylogentic tree data gathered from [ProteinHistorian](http://lighthouse.ucsf.edu/ProteinHistorian/) and are also provided.
