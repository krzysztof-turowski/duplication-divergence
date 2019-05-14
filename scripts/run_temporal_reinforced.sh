#!/bin/bash

make dd_temporal_algorithms

PERFECT=(0.001 0.01 0.1)
BINSIZE=(1 5 10)
THRESHOLD=(0.5 0.7 0.9)
for perfect in ${PERFECT[@]}; do
  for binsize in ${BINSIZE[@]}; do
    ./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -algorithm:sort_by_p_uv_sum -gt:100 -st:100000 -binsize:$binsize -perfect:$perfect
  done
  for threshold in ${THRESHOLD[@]}; do
    ./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.3 -r:1.0 -p0:0.6 -algorithm:p_uv_threshold -gt:100 -st:100000 -threshold:$threshold -perfect:$perfect
  done
done
python -B ./src/dd_temporal_order_plot.py synthetic-50-10-PS-0.300-1.00 --export pdf

for perfect in ${PERFECT[@]}; do
  for binsize in ${BINSIZE[@]}; do
    ./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:sort_by_p_uv_sum -gt:100 -st:100000 -binsize:$binsize -perfect:$perfect
  done
  for threshold in ${THRESHOLD[@]}; do
    ./dd_temporal_algorithms -action:synthetic -n:50 -n0:10 -mode:pastor_satorras -p:0.6 -r:1.0 -p0:0.6 -algorithm:p_uv_threshold -gt:100 -st:100000 -threshold:$threshold -perfect:$perfect
  done
done
python -B ./src/dd_temporal_order_plot.py synthetic-50-10-PS-0.600-1.00 --export pdf

PERFECT=(0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1)
./dd_temporal_algorithms -action:real_data -graph:G-hep-th-citations-scc.txt -mode:pastor_satorras -p:0.72 -r:1.0 -algorithm:sort_by_degree
./dd_temporal_algorithms -action:real_data -graph:G-hep-th-citations-scc.txt -mode:pastor_satorras -p:0.72 -r:1.0 -algorithm:peel_by_degree
./dd_temporal_algorithms -action:real_data -graph:G-hep-th-citations-scc.txt -mode:pastor_satorras -p:0.72 -r:1.0 -algorithm:sort_by_neighborhood
./dd_temporal_algorithms -action:real_data -graph:G-hep-th-citations-scc.txt -mode:pastor_satorras -p:0.72 -r:1.0 -algorithm:peel_by_neighborhood
for perfect in ${PERFECT[@]}; do
  ./dd_temporal_algorithms -action:real_data -graph:G-hep-th-citations-scc.txt -mode:pastor_satorras -p:0.72 -r:1.0 -algorithm:sort_by_p_uv_sum -st:100000 -binsize:1 -perfect:$perfect
done
for perfect in ${PERFECT[@]}; do
  ./dd_temporal_algorithms -action:real_data -graph:G-hep-th-citations-scc.txt -mode:pastor_satorras -p:0.72 -r:1.0 -algorithm:p_uv_threshold -st:100000 -threshold:0.5 -perfect:$perfect
done

./dd_temporal_algorithms -action:real_data -graph:G-dynamic-simplewiki-10k.txt -mode:pastor_satorras -p:0.66 -r:0.5 -algorithm:sort_by_degree
./dd_temporal_algorithms -action:real_data -graph:G-dynamic-simplewiki-10k.txt -mode:pastor_satorras -p:0.66 -r:0.5 -algorithm:peel_by_degree
./dd_temporal_algorithms -action:real_data -graph:G-dynamic-simplewiki-10k.txt -mode:pastor_satorras -p:0.66 -r:0.5 -algorithm:sort_by_neighborhood
./dd_temporal_algorithms -action:real_data -graph:G-dynamic-simplewiki-10k.txt -mode:pastor_satorras -p:0.66 -r:0.5 -algorithm:peel_by_neighborhood
for perfect in ${PERFECT[@]}; do
  ./dd_temporal_algorithms -action:real_data -graph:G-dynamic-simplewiki-10k.txt -mode:pastor_satorras -p:0.66 -r:0.5 -algorithm:sort_by_p_uv_sum -st:100000 -binsize:1 -perfect:$perfect
done
for perfect in ${PERFECT[@]}; do
  ./dd_temporal_algorithms -action:real_data -graph:G-dynamic-simplewiki-10k.txt -mode:pastor_satorras -p:0.66 -r:0.5 -algorithm:p_uv_threshold -st:100000 -threshold:0.5 -perfect:$perfect
done

./dd_temporal_algorithms -action:real_data -graph:G-college-msg.txt -mode:pastor_satorras -p:0.65 -r:0.45 -algorithm:sort_by_degree
./dd_temporal_algorithms -action:real_data -graph:G-college-msg.txt -mode:pastor_satorras -p:0.65 -r:0.45 -algorithm:peel_by_degree
./dd_temporal_algorithms -action:real_data -graph:G-college-msg.txt -mode:pastor_satorras -p:0.65 -r:0.45 -algorithm:sort_by_neighborhood
./dd_temporal_algorithms -action:real_data -graph:G-college-msg.txt -mode:pastor_satorras -p:0.65 -r:0.45 -algorithm:peel_by_neighborhood
for perfect in ${PERFECT[@]}; do
  ./dd_temporal_algorithms -action:real_data -graph:G-college-msg.txt -mode:pastor_satorras -p:0.65 -r:0.45 -algorithm:sort_by_p_uv_sum -st:100000 -binsize:1 -perfect:$perfect
done
for perfect in ${PERFECT[@]}; do
  ./dd_temporal_algorithms -action:real_data -graph:G-college-msg.txt -mode:pastor_satorras -p:0.65 -r:0.45 -algorithm:p_uv_threshold -st:100000 -threshold:0.5 -perfect:$perfect
done
