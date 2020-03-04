import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os
import numpy as np
from matplotlib import rc
plt.rc('text', usetex=True)
plt.rc('text.latex', unicode=True)
plt.rc('font', family='serif')
parser = argparse.ArgumentParser(description = "Make macroscopic variables plots")
parser.add_argument("--path", type = str, required = True)
args = parser.parse_args()
os.chdir(args.path)
c1 = "royalblue"    
c2 = "darkred"

al = 0.6

def arrowPlot(axes,xs,ys,inverse = False,f=4,color="k"):
    if inverse:
         xs = [xs[-(i+1)] for i in range(len(xs))]
         ys = [ys[-(i+1)] for i in range(len(ys))]
    for i in range(len(xs)-1):
        if i % f ==1:
            axes.annotate("",xy=(xs[i+1],ys[i+1]),xytext=(xs[i],ys[i]),arrowprops=dict(arrowstyle="-|>,head_width=0.5,head_length=1", color = color,alpha=al))
        else:
            axes.plot([xs[i],xs[i+1]],[ys[i],ys[i+1]],color = color,alpha=al)

df = pd.read_csv("exothermic/Macroscopic.csv")
df = df.sort_values(by="kT")
var = [["M","E"],["Xi","CV"]]
varLabel = [[r"$\rho $",r"$U / N$"],[r"$\chi$",r"$C_v /N$"]]
fig, axs  = plt.subplots(nrows=2,ncols=2,sharex = True,figsize=(12,8))

fig.suptitle("Mu = {0}, J = {1}, C = {2}".format(df["Mu"].to_list()[0],df["J"].to_list()[0],df["C"].to_list()[0]),fontsize=24)

for i in range(2):
    for j in range(2):
        if i==1:
            axs[i][j].set_xlabel("kT")
        axs[i][j].set_ylabel(varLabel[i][j])
        try:
            axs[i][j].errorbar(df["kT"].tolist(),df[var[i][j]].tolist(),yerr = df["Msigma"].tolist(),color = c1,alpha=al)
        except KeyError:
            print("no "+ var[i][j] +" STD")
        arrowPlot(axs[i][j],df["kT"].tolist(),df[var[i][j]].tolist(),inverse=True,f=int(df["kT"].shape[0]/5),color = c1)
        #axs[i][j].scatter(df["kT"].tolist(),df[var[i][j]].tolist())

df = pd.read_csv("endothermic/Macroscopic.csv")
df = df.sort_values(by="kT")
var = [["M","E"],["Xi","CV"]]
varLabel = [[r"$\rho $",r"$U / N$"],[r"$\chi$",r"$C_v /N$"]]


for i in range(2):
    for j in range(2):
        if i==1:
            axs[i][j].set_xlabel("kT")
        axs[i][j].set_ylabel(varLabel[i][j])
        try:
            axs[i][j].errorbar(df["kT"].tolist(),df[var[i][j]].tolist(),yerr = df["Msigma"].tolist(),color = c2,alpha=al)
        except KeyError:
            print("no "+ var[i][j] +" STD")
        #axs[i][j].plot(df["kT"].tolist(),df[var[i][j]].tolist())
        #axs[i][j].scatter(df["kT"].tolist(),df[var[i][j]].tolist())
        arrowPlot(axs[i][j],df["kT"].tolist(),df[var[i][j]].tolist(),f=int(df["kT"].shape[0]/5),color = c2)
        axs[i][j].grid()
#plt.tight_layout()
fig.savefig("tempHystPlots.png",dpi=300)