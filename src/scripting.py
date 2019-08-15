import sys
import os
sizes = [
    8,
    10,
    14,
    16,
    20,
    26,
    32,
    40,
    50,
    64
]
for n in sizes:
    print(n)
    os.system("julia runParallelSet.jl -N {0} -A 7 -E PenalizedEnergy -M LatticeGas".format(n))
for path in os.listdir("../outputs"):
    path = "../outputs/" + path
    print(path)
    os.system("julia analysisMetropolis.jl -D \"{0}\" ".format(path))
os.system("echo \" Tu simulaciones ya terminaron \" | mail -s \"Simulaciones terminadas\" sayeg@ciencias.unam.mx ")