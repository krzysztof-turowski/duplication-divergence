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
- available graph libraries: *koala*, *snap*, *networkit*,
- available LP libraries: *glpk*, *gurobi* (note: need separate download and license.),
- available automorphisms library: *nauty*.

The choice of libraries to use is done via variables **GRAPH_LIB** and **LP_SOLVER** in Makefile.

## Compiling and Running the Algorithms

**Available programs**

| Program        | Function           |
| ------------- |-------------|
| `dd_automorphisms`    | Find the number of automorphisms and its p-value |
| `dd_ml_estimation`      | Maximum likelihood estimation (MLE) of the parameters      |
| `dd_recurrence_estimation` | Parameter estimation using our Recurrence-relation method      |
|`dd_temporal_bound` | Algorithms to approximate the integer programming optimization solution with the linear programming optimization(uses Gurobi)|
| `dd_temporal_algorithms` | Temporal ranking algorithms |


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

<!-- **Finding the number of automorphisms and its p-value**
```bash
make dd_automorphisms
```

**Maximum likelihood estimation (MLE) of the parameters:**
```bash
make dd_ml_estimation
```
An example run:
```bash
# Write MLE of the parameters of the pastor_satorras model that fits the real data into a file
./dd_ml_estimation -action:real_data -graph:G-test.txt -mode:pastor_satorras -st:1000
# Plot the MLE result and save it as a PDF file
python -B ./src/dd_ml_estimation_plot.py G-test <ML_file> --export pdf

```

**Parameter estimation using our Recurrence-relation method:**
```bash
make dd_recurrence_estimation
```
An example run:
```bash
# Write estimates of the parameters of the pastor_satorras model that fits the real data into a file
./dd_recurrence_estimation -action:real_data -graph:G-test.txt -mode:pastor_satorras -st:1000
# Plot the estimate result and save it as a PDF file
python -B ./src/dd_recurrence_estimation_plot.py G-test <estimate_file> --export pdf
```

#### Temporal ranking algorithms

**Temporal ranking bound**\
Algorithms to approximate the integer programming optimization solution with the linear programming optimization(uses Gurobi)

```bash
make dd_temporal_bound
```
**Temporal ranking algorithms**

```bash
make dd_temporal_algorithms
``` -->

## Data
All the data files we have used in our experiments are provided in `files` folder - Arxiv, Simple English Wikipedia and CollegeMsg networks, and the protein-protein interaction (PPI) networks. We collect the PPI data of seven species from [BioGRID](https://thebiogrid.org/). The phylogentic tree based age of the proteins is gathered from [ProteinHistorian](http://lighthouse.ucsf.edu/ProteinHistorian/) and is also provided.
