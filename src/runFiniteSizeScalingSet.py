import numpy as np
import sys
import os
sizes = [
    8,
    9,
    10,
    14,
    16,
    20,
    27,
    32,
    40,
    50,
    64
]
for n in sizes:
    print(n)
    os.system("julia runParallelSet.jl -N {0} -A 4 ".format(n))
for path in os.listdir("../outputs"):
    path = "../outputs/" + path
    print(path)
    os.system("julia analysisMetropolis.jl -D \"{0}\" ".format(path))