#!/bin/bash

make dd_temporal_bound
make dd_temporal_algorithms

./dd_temporal_bound -algorithm:exact -n:13 -n0:6 -mode:pastor_satorras -p:0.3 -r:0.6 -p0:0.6 -gt:100
./dd_temporal_bound -algorithm:local-unif-sampling -n:13 -n0:6 -mode:pastor_satorras -p:0.3 -r:0.6 -p0:0.6 -gt:100 -st:10
./dd_temporal_bound -algorithm:local-unif-sampling -n:13 -n0:6 -mode:pastor_satorras -p:0.3 -r:0.6 -p0:0.6 -gt:100 -st:100
./dd_temporal_bound -algorithm:local-unif-sampling -n:13 -n0:6 -mode:pastor_satorras -p:0.3 -r:0.6 -p0:0.6 -gt:100 -st:1000
./dd_temporal_bound -algorithm:high-prob-sampling -n:13 -n0:6 -mode:pastor_satorras -p:0.3 -r:0.6 -p0:0.6 -gt:100 -st:10
./dd_temporal_bound -algorithm:high-prob-sampling -n:13 -n0:6 -mode:pastor_satorras -p:0.3 -r:0.6 -p0:0.6 -gt:100 -st:100
./dd_temporal_bound -algorithm:high-prob-sampling -n:13 -n0:6 -mode:pastor_satorras -p:0.3 -r:0.6 -p0:0.6 -gt:100 -st:1000
python -B ./src/dd_temporal_order_plot.py synthetic-13-6-PS-0.300-0.60-TA --export pdf

./dd_temporal_bound -algorithm:exact -n:13 -n0:6 -mode:pastor_satorras -p:0.6 -r:0.6 -p0:0.6 -gt:100
./dd_temporal_bound -algorithm:local-unif-sampling -n:13 -n0:6 -mode:pastor_satorras -p:0.6 -r:0.6 -p0:0.6 -gt:100 -st:10
./dd_temporal_bound -algorithm:local-unif-sampling -n:13 -n0:6 -mode:pastor_satorras -p:0.6 -r:0.6 -p0:0.6 -gt:100 -st:100
./dd_temporal_bound -algorithm:local-unif-sampling -n:13 -n0:6 -mode:pastor_satorras -p:0.6 -r:0.6 -p0:0.6 -gt:100 -st:1000
./dd_temporal_bound -algorithm:high-prob-sampling -n:13 -n0:6 -mode:pastor_satorras -p:0.6 -r:0.6 -p0:0.6 -gt:100 -st:10
./dd_temporal_bound -algorithm:high-prob-sampling -n:13 -n0:6 -mode:pastor_satorras -p:0.6 -r:0.6 -p0:0.6 -gt:100 -st:100
./dd_temporal_bound -algorithm:high-prob-sampling -n:13 -n0:6 -mode:pastor_satorras -p:0.6 -r:0.6 -p0:0.6 -gt:100 -st:1000
python -B ./src/dd_temporal_order_plot.py synthetic-13-6-PS-0.600-0.60-TA --export pdf

./dd_temporal_bound -algorithm:local-unif-sampling -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -gt:100 -st:10000
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -algorithm:sort_by_degree -gt:100
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -algorithm:peel_by_degree -gt:100
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -algorithm:peel_by_neighborhood -gt:100
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -algorithm:p_uv_threshold -gt:100 -st:100000 -threshold:0.5
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -algorithm:p_uv_threshold -gt:100 -st:100000 -threshold:0.6
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -algorithm:p_uv_threshold -gt:100 -st:100000 -threshold:0.7
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -algorithm:p_uv_threshold -gt:100 -st:100000 -threshold:0.8
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -algorithm:sort_by_p_uv_sum -gt:100 -st:100000 -binsize:1
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -algorithm:sort_by_p_uv_sum -gt:100 -st:100000 -binsize:5
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -algorithm:sort_by_p_uv_sum -gt:100 -st:100000 -binsize:10
python -B ./src/dd_temporal_order_plot.py synthetic-50-10-PS-0.300-1.00 --export pdf

./dd_temporal_bound -algorithm:local-unif-sampling -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -gt:100 -st:10000
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:sort_by_degree -gt:100
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:peel_by_degree -gt:100
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:sort_by_neighborhood -gt:100
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:peel_by_neighborhood -gt:100
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:p_uv_threshold -gt:100 -st:100000 -threshold:0.5
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:p_uv_threshold -gt:100 -st:100000 -threshold:0.6
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:p_uv_threshold -gt:100 -st:100000 -threshold:0.7
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:p_uv_threshold -gt:100 -st:100000 -threshold:0.8
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:sort_by_p_uv_sum -gt:100 -st:100000 -binsize:1
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:sort_by_p_uv_sum -gt:100 -st:100000 -binsize:5
./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:sort_by_p_uv_sum -gt:100 -st:100000 -binsize:10
python -B ./src/dd_temporal_order_plot.py synthetic-50-10-PS-0.600-1.00 --export pdf
