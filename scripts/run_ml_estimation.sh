#!/bin/bash

make dd_ml_estimation

./dd_ml_estimation -action:real_data -graph:G-100-20-PS-0.1-0.3.txt -mode:pastor_satorras -st:100
python -B ./src/dd_ml_estimation_plot.py G-100-20-PS-0.1-0.3-PS --export pdf

./dd_ml_estimation -action:real_data -graph:G-100-20-PS-0.99-3.txt -mode:pastor_satorras -st:100
python -B ./src/dd_ml_estimation_plot.py G-100-20-PS-0.99-3-PS --export pdf
