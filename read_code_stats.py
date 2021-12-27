#!/usr/bin/env python3

import pstats
import os

file="rand_circ_stats.txt"
# output="stats.txt"
# os.system("python -m cProfile -o " + output + " " + file + " 5 1 2 2")
p = pstats.Stats(file)
p.sort_stats('cumulative').print_stats(10)