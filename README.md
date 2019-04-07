# duplication-divergence
Software for duplication-divergence model analysis

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
- available LP libraries: *glpk*, *gurobi* (note: need separate download),
- available automorphisms library: *nauty*.

The choice of libraries to use is done via variables **GRAPH_LIB** and **LP_SOLVER** in Makefile.
