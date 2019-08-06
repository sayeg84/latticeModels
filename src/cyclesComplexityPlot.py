import pandas as pd
import numpy as no
import matplotlib.pyplot as plt
plt.figure(figsize=(6,4))
df = pd.read_csv("../cyclesComplexity.csv")
df["ours"] *= 1000
df["LightGraphs"] *= 1000
#plt.scatter(df["n"].tolist(),df["ours"].tolist(),label="Nosotros")
plt.loglog(df["n"].tolist(),df["ours"].tolist(),ms=2)

#plt.scatter(df["n"].tolist(),df["LightGraphs"].tolist(),label="LightGraphs")
plt.loglog(df["n"].tolist(),df["LightGraphs"].tolist(),ms=2)
plt.grid()
plt.legend()
plt.xlabel("lado del lattice (cuadrada)")
plt.ylabel("ms de ejecucion")
plt.tight_layout()
plt.savefig("../cyclesComplexity.png",dpi=300)