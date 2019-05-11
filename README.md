# duplication-divergence

This repository provides software tools for the following studies of duplication-divergence random graph model on real data:
* Parameter estimation for duplication-divergence model\
  Associated paper: [1] [Revisiting Parameter Estimation in Biological Networks: Influence of Symmetries](https://www.cs.purdue.edu/homes/jithinks/files/publications/parameter_estimation_2019.pdf)
* Given a single snapshot of a dynamic network, provide algorithms to infer the arrival order of the nodes.
It implements the optimal and approximate solutions of the problem. The optimal solution is a result of an integer programming formulation with coefficients found by importance sampling techniques.\
  Associated paper: [2] [Temporal Ordered Clustering in Dynamic Networks](https://arxiv.org/abs/1905.00672)

## Getting started
Download the repository or clone it with the following command
```bash
git clone https://github.com/krzysztof-turowski/duplication-divergence.git DIR
cd DIR
git submodule update --init --recursive
```

## Installation
Install libraries and set environment variables using the following command
```bash
bash ./configure_libs.sh LIB1 LIB2 ...
```
where
- available graph libraries: *koala*, *snap*, *networkit* ,
- available LP libraries: *glpk*, *gurobi* (note: need separate download and license),
- available automorphisms library: *nauty*.

The choice of libraries to use is done via variables **GRAPH_LIB** (default: *snap*) and **LP_SOLVER** (default: *glpk*) in Makefile.

## Overview of the library

### Available programs

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

| Program        | Function           |
| ------------- |-------------|
| `dd_ml_estimation_plot`      | Plot the log-likelihood as a function of the parameters `p` and `r` of the duplication-divergence model.
| `dd_recurrence_estimation_plot` | Plot estimated parameters via degree, open triangles and traingles recurrence-relations with our methods and find the confidence interval and convergence point |
|`dd_temporal_order_plot` | Plot the temporal order algorithms for recovering the node orders (ordered cluster labels) in the precision vs density curve|

## Parameter estimation by automorphisms and recurrence-relations estimation

### `dd_automorphisms`: finding the number of automorphisms and its p-value

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

### `dd_recurrence_estimation`: parameter estimation using recurrence-relation method

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

### `dd_ml_estimation`: maximum likelihood estimation (MLE) of the parameters

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

## Temporal ordering algorithms

### `dd_temporal_bound`: upper bound on temporal ordering

**Compilation**

```bash
make dd_temporal_bound
```

### `dd_temporal_algorithms`: temporal ordering algorithms

**Compilation**

```bash
make dd_temporal_algorithms
```

**Note:** the results from `dd_temporal_bound` and `dd_temporal_algorithms` are **appended** to the temporary files.

### `dd_temporal_order_plot`: plot temporal ordering results

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

All the data files we have used in our experiments are provided in `files` folder - Arxiv, Simple English Wikipedia and CollegeMsg networks, and the protein-protein interaction (PPI) networks.
We collect the PPI data of seven species from [BioGRID](https://thebiogrid.org/).
The phylogentic tree based age of the proteins are gathered from [ProteinHistorian](http://lighthouse.ucsf.edu/ProteinHistorian/) and are also provided.
