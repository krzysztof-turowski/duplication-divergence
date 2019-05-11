#!/bin/bash

make dd_automorphisms

./dd_automorphisms -action:real_graph -graph:G-a-thaliana.txt
./dd_automorphisms -action:real_seed -graph:G-a-thaliana.txt -st:100 -mode:pastor_satorras -p:0.98 -r:0.49

./dd_automorphisms -action:real_graph -graph:G-c-elegans.txt
./dd_automorphisms -action:real_seed -graph:G-c-elegans.txt -st:100 -mode:pastor_satorras -p:0.85 -r:0.35

./dd_automorphisms -action:real_graph -graph:G-d-melanogaster.txt
./dd_automorphisms -action:real_seed -graph:G-d-melanogaster.txt -st:100 -mode:pastor_satorras -p:0.53 -r:0.92

./dd_automorphisms -action:real_graph -graph:G-homo-sapiens.txt
./dd_automorphisms -action:real_seed -graph:G-homo-sapiens.txt -st:100 -mode:pastor_satorras -p:0.64 -r:0.49  

./dd_automorphisms -action:real_graph -graph:G-mus-musculus.txt
./dd_automorphisms -action:real_seed -graph:G-mus-musculus.txt -st:100 -mode:pastor_satorras -p:0.96 -r:0.32

./dd_automorphisms -action:real_graph -graph:G-s-cerevisiae.txt
./dd_automorphisms -action:real_seed -graph:G-s-cerevisiae.txt -st:100 -mode:pastor_satorras -p:0.98 -r:0.35

./dd_automorphisms -action:real_graph -graph:G-s-pombe.txt
./dd_automorphisms -action:real_seed -graph:G-s-pombe.txt -st:100 -mode:pastor_satorras -p:0.983 -r:0.85
