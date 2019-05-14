#!/bin/bash

make dd_graphlets

./dd_graphlets -action:real_seed -graph:G-a-thaliana.txt -gt:100 -mode:pastor_satorras -p:0.98 -r:0.49
./dd_graphlets -action:real_seed -graph:G-a-thaliana.txt -gt:100 -mode:pastor_satorras -p:0.44 -r:0.43

./dd_graphlets -action:real_seed -graph:G-c-elegans.txt -gt:100 -mode:pastor_satorras -p:0.85 -r:0.35
./dd_graphlets -action:real_seed -graph:G-c-elegans.txt -gt:100 -mode:pastor_satorras -p:0.47 -r:0.14

./dd_graphlets -action:real_seed -graph:G-d-melanogaster.txt -gt:100 -mode:pastor_satorras -p:0.53 -r:0.92
./dd_graphlets -action:real_seed -graph:G-d-melanogaster.txt -gt:100 -mode:pastor_satorras -p:0.44 -r:0.75

./dd_graphlets -action:real_seed -graph:G-homo-sapiens.txt -gt:100 -mode:pastor_satorras -p:0.64 -r:0.49
./dd_graphlets -action:real_seed -graph:G-homo-sapiens.txt -gt:100 -mode:pastor_satorras -p:0.43 -r:2.39

./dd_graphlets -action:real_seed -graph:G-mus-musculus.txt -gt:100 -mode:pastor_satorras -p:0.96 -r:0.32
./dd_graphlets -action:real_seed -graph:G-mus-musculus.txt -gt:100 -mode:pastor_satorras -p:0.48 -r:0.12

./dd_graphlets -action:real_seed -graph:G-s-cerevisiae.txt -gt:100 -mode:pastor_satorras -p:0.98 -r:0.35
./dd_graphlets -action:real_seed -graph:G-s-cerevisiae.txt -gt:100 -mode:pastor_satorras -p:0.28 -r:38.25

./dd_graphlets -action:real_seed -graph:G-s-pombe.txt -gt:100 -mode:pastor_satorras -p:0.983 -r:0.85
./dd_graphlets -action:real_seed -graph:G-s-pombe.txt -gt:100 -mode:pastor_satorras -p:0.46 -r:1.02
