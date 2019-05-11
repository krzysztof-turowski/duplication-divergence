#!/bin/bash

make dd_ml_estimation

./dd_ml_estimation -action:real_data -graph:G-100-20-PS-0.1-0.3.txt -mode:pastor_satorras -st:100
python -B ./dd_ml_estimation_plot.py G-100-20-PS-0.1-0.3-ML.txt --export pdf

./dd_ml_estimation -action:real_data -graph:G-100-20-PS-0.99-3.0.txt -mode:pastor_satorras -st:100
python -B ./dd_ml_estimation_plot.py G-100-20-PS-0.99-3.0-ML.txt --export pdf
