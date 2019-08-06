import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os
from matplotlib import rc
plt.rc('text', usetex=True)
plt.rc('text.latex', unicode=True)
plt.rc('font', family='serif')
parser = argparse.ArgumentParser(description = "Make macroscopic variables plots")
parser.add_argument("--path", type = str, required = True)
args = parser.parse_args()
os.chdir(args.path)
df = pd.read_csv("Macroscopic.csv")
df = df.sort_values(by="kT")
var = [["M","E"],["Xi","CV"]]
fig, axs  = plt.subplots(nrows=2,ncols=2,sharex = True,figsize=(12,8))

for i in range(2):
    for j in range(2):
        if i==1:
            axs[i][j].set_xlabel("kT")
        axs[i][j].set_ylabel(var[i][j])
        try:
            axs[i][j].errorbar(df["kT"].tolist(),df[var[i][j]].tolist(),yerr = df["Msigma"].tolist())
        except KeyError:
            print("no "+ var[i][j] +" STD")
        axs[i][j].plot(df["kT"].tolist(),df[var[i][j]].tolist())
        axs[i][j].scatter(df["kT"].tolist(),df[var[i][j]].tolist())
        axs[i][j].grid()
plt.tight_layout()
fig.savefig("tempPlots.png",dpi=300)