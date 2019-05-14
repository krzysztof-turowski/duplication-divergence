#!/bin/bash

python -B -W ignore ./src/dd_power_law.py G-a-thaliana -action real_graph
python -B -W ignore ./src/dd_power_law.py G-a-thaliana -action real_seed -gt 100 -mode pastor_satorras -p 0.98 -r 0.49
python -B -W ignore ./src/dd_power_law.py G-a-thaliana -action real_seed -gt 100 -mode pastor_satorras -p 0.44 -r 0.43

python -B -W ignore ./src/dd_power_law.py G-c-elegans -action real_graph
python -B -W ignore ./src/dd_power_law.py G-c-elegans -action real_seed -gt 100 -mode pastor_satorras -p 0.85 -r 0.35
python -B -W ignore ./src/dd_power_law.py G-c-elegans -action real_seed -gt 100 -mode pastor_satorras -p 0.47 -r 0.14

python -B -W ignore ./src/dd_power_law.py G-d-melanogaster -action real_graph
python -B -W ignore ./src/dd_power_law.py G-d-melanogaster -action real_seed -gt 100 -mode pastor_satorras -p 0.53 -r 0.92
python -B -W ignore ./src/dd_power_law.py G-d-melanogaster -action real_seed -gt 100 -mode pastor_satorras -p 0.44 -r 0.75

python -B -W ignore ./src/dd_power_law.py G-homo-sapiens -action real_graph
python -B -W ignore ./src/dd_power_law.py G-homo-sapiens -action real_seed -gt 100 -mode pastor_satorras -p 0.64 -r 0.49
python -B -W ignore ./src/dd_power_law.py G-homo-sapiens -action real_seed -gt 100 -mode pastor_satorras -p 0.43 -r 2.39

python -B -W ignore ./src/dd_power_law.py G-mus-musculus -action real_graph
python -B -W ignore ./src/dd_power_law.py G-mus-musculus -action real_seed -gt 100 -mode pastor_satorras -p 0.96 -r 0.32
python -B -W ignore ./src/dd_power_law.py G-mus-musculus -action real_seed -gt 100 -mode pastor_satorras -p 0.48 -r 0.12

python -B -W ignore ./src/dd_power_law.py G-s-cerevisiae -action real_graph
python -B -W ignore ./src/dd_power_law.py G-s-cerevisiae -action real_seed -gt 100 -mode pastor_satorras -p 0.98 -r 0.35
python -B -W ignore ./src/dd_power_law.py G-s-cerevisiae -action real_seed -gt 100 -mode pastor_satorras -p 0.28 -r 38.25

python -B -W ignore ./src/dd_power_law.py G-s-pombe -action real_graph
python -B -W ignore ./src/dd_power_law.py G-s-pombe -action real_seed -gt 100 -mode pastor_satorras -p 0.983 -r 0.85
python -B -W ignore ./src/dd_power_law.py G-s-pombe -action real_seed -gt 100 -mode pastor_satorras -p 0.46 -r 1.02
