#!/bin/bash

make dd_recurrence_estimation

./dd_recurrence_estimation -action:real_data -graph:G-100-20-PS-0.1-0.3.txt -st:100 -mode:pastor_satorras
python -B ./src/dd_recurrence_estimation_plot.py G-100-20-PS-0.1-0.3-PS --export pdf

./dd_recurrence_estimation -action:real_data -graph:G-100-20-PS-0.99-3.txt -st:100 -mode:pastor_satorras
python -B ./src/dd_recurrence_estimation_plot.py G-100-20-PS-0.99-3-PS --export pdf

./dd_recurrence_estimation -action:real_data -graph:G-a-thaliana.txt -st:100 -mode:pastor_satorras
python -B ./src/dd_recurrence_estimation_plot.py G-a-thaliana-PS --export pdf

./dd_recurrence_estimation -action:real_data -graph:G-c-elegans.txt -st:100 -mode:pastor_satorras
python -B ./src/dd_recurrence_estimation_plot.py G-c-elegans-PS --export pdf

./dd_recurrence_estimation -action:real_data -graph:G-d-melanogaster.txt -st:100 -mode:pastor_satorras
python -B ./src/dd_recurrence_estimation_plot.py G-d-melanogaster-PS --export pdf

./dd_recurrence_estimation -action:real_data -graph:G-homo-sapiens.txt -st:100 -mode:pastor_satorras
python -B ./src/dd_recurrence_estimation_plot.py G-homo-sapiens-PS --export pdf

./dd_recurrence_estimation -action:real_data -graph:G-mus-musculus.txt -st:100 -mode:pastor_satorras
python -B ./src/dd_recurrence_estimation_plot.py G-mus-musculus-PS --export pdf

./dd_recurrence_estimation -action:real_data -graph:G-s-cerevisiae.txt -st:100 -mode:pastor_satorras
python -B ./src/dd_recurrence_estimation_plot.py G-s-cerevisiae-PS --export pdf

./dd_recurrence_estimation -action:real_data -graph:G-s-pombe.txt -st:100 -mode:pastor_satorras
python -B ./src/dd_recurrence_estimation_plot.py G-s-pombe-PS --export pdf
