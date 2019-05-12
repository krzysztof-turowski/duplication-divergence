#!/bin/bash

make dd_automorphisms

./dd_automorphisms -action:real_graph -graph:G-a-thaliana.txt
./dd_automorphisms -action:real_seed -graph:G-a-thaliana.txt -st:100 -mode:pastor_satorras -p:0.98 -r:0.49
./dd_automorphisms -action:real_seed -graph:G-a-thaliana.txt -st:100 -mode:pastor_satorras -p:0.44 -r:0.43

./dd_automorphisms -action:real_graph -graph:G-c-elegans.txt
./dd_automorphisms -action:real_seed -graph:G-c-elegans.txt -st:100 -mode:pastor_satorras -p:0.85 -r:0.35
./dd_automorphisms -action:real_seed -graph:G-c-elegans.txt -st:100 -mode:pastor_satorras -p:0.47 -r:0.14

./dd_automorphisms -action:real_graph -graph:G-d-melanogaster.txt
./dd_automorphisms -action:real_seed -graph:G-d-melanogaster.txt -st:100 -mode:pastor_satorras -p:0.53 -r:0.92
./dd_automorphisms -action:real_seed -graph:G-d-melanogaster.txt -st:100 -mode:pastor_satorras -p:0.44 -r:0.75

./dd_automorphisms -action:real_graph -graph:G-homo-sapiens.txt
./dd_automorphisms -action:real_seed -graph:G-homo-sapiens.txt -st:100 -mode:pastor_satorras -p:0.64 -r:0.49
./dd_automorphisms -action:real_seed -graph:G-homo-sapiens.txt -st:100 -mode:pastor_satorras -p:0.43 -r:2.39

./dd_automorphisms -action:real_graph -graph:G-mus-musculus.txt
./dd_automorphisms -action:real_seed -graph:G-mus-musculus.txt -st:100 -mode:pastor_satorras -p:0.96 -r:0.32
./dd_automorphisms -action:real_seed -graph:G-mus-musculus.txt -st:100 -mode:pastor_satorras -p:0.48 -r:0.12

./dd_automorphisms -action:real_graph -graph:G-s-cerevisiae.txt
./dd_automorphisms -action:real_seed -graph:G-s-cerevisiae.txt -st:100 -mode:pastor_satorras -p:0.98 -r:0.35
./dd_automorphisms -action:real_seed -graph:G-s-cerevisiae.txt -st:100 -mode:pastor_satorras -p:0.28 -r:38.25

./dd_automorphisms -action:real_graph -graph:G-s-pombe.txt
./dd_automorphisms -action:real_seed -graph:G-s-pombe.txt -st:100 -mode:pastor_satorras -p:0.983 -r:0.85
./dd_automorphisms -action:real_seed -graph:G-s-pombe.txt -st:100 -mode:pastor_satorras -p:0.46 -r:1.02

for p in $(seq 0.0 1.0 0.1);
do for r in $(seq 0.0 0.4 5.0);
  do ./dd_automorphisms -action:synthetic -n:200 -n0:20 -p0:1.0 -st:100 -mode:pastor_satorras -p:$p -r:$r
  done
done

for p in $(seq 0.0 1.0 0.1);
do for r in $(seq 0.0 0.4 5.0);
  do ./dd_automorphisms -action:synthetic -n:1000 -n0:20 -p0:1.0 -st:100 -mode:pastor_satorras -p:$p -r:$r
  done
done

for p in $(seq 0.0 1.0 0.1);
do for r in $(seq 0.0 0.4 5.0);
  do ./dd_automorphisms -action:synthetic -n:5000 -n0:20 -p0:1.0 -st:100 -mode:pastor_satorras -p:$p -r:$r
  done
done
