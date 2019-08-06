import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import os
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib import rc
plt.rc('text', usetex=True)
plt.rc('text.latex', unicode=True)
plt.rc('font', family='serif')
parser = argparse.ArgumentParser(description = "Make macroscopic variables plots")
parser.add_argument("--path", type = str, required = True)
args = parser.parse_args()
os.chdir(args.path)
df = pd.read_csv("Macroscopic.csv")
mus = sorted(list(set(df["Mu"].tolist())))
cs = sorted(list(set(df["C"].tolist())))

indexes=[[(mus[i],cs[j]) for i in range(len(mus))] for j in range(len(cs))]
mag=[[df[df["Mu"]==mus[i]][df["C"]==cs[j]]["M"].values[0] for i in range(len(mus))] for j in range(len(cs))]
ener=[[df[df["Mu"]==mus[i]][df["C"]==cs[j]]["E"].values[0] for i in range(len(mus))] for j in range(len(cs))]
cv=[[df[df["Mu"]==mus[i]][df["C"]==cs[j]]["CV"].values[0] for i in range(len(mus))] for j in range(len(cs))]
chi=[[df[df["Mu"]==mus[i]][df["C"]==cs[j]]["Xi"].values[0] for i in range(len(mus))] for j in range(len(cs))]
mags=[[df[df["Mu"]==mus[i]][df["C"]==cs[j]]["Msigma"].values[0] for i in range(len(mus))] for j in range(len(cs))]
eners=[[df[df["Mu"]==mus[i]][df["C"]==cs[j]]["Esigma"].values[0] for i in range(len(mus))] for j in range(len(cs))]
cvs=[[df[df["Mu"]==mus[i]][df["C"]==cs[j]]["CVsigma"].values[0] for i in range(len(mus))] for j in range(len(cs))]
chis=[[df[df["Mu"]==mus[i]][df["C"]==cs[j]]["Xisigma"].values[0] for i in range(len(mus))] for j in range(len(cs))]

def groupedPlots(color2,colorbar="double"):

    f1 = 20
    f2 = 24
    f3 = 32
    size = 50
    
    fig, axs = plt.subplots(2,2,figsize=(16,12),sharex=True)
    
    normalize = mcolors.Normalize(vmin=min(cs), vmax=max(cs))
    colormap = color2
    
    fig.subplots_adjust(hspace=0.05)
    
    for c in range(len(cs)):

        
        axs[0][0].plot(mus,mag[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
        
        axs[0][0].errorbar(mus,mag[c],yerr=mags[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
        
        axs[0][0].scatter(mus,mag[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), s=size)
        
        axs[1][0].plot(mus,chi[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
        
        axs[1][0].errorbar(mus,chi[c],yerr=chis[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
        
        axs[1][0].scatter(mus,chi[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), s=size)
        
        axs[0][1].plot(mus,ener[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
        
        axs[0][1].errorbar(mus,ener[c],yerr=eners[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
        
        axs[0][1].scatter(mus,ener[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), s=size)
        
        axs[1][1].plot(mus,cv[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
        
        axs[1][1].errorbar(mus,cv[c],yerr=cvs[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
        
        axs[1][1].scatter(mus,cv[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), s=size)
    
    
    

    axs[1][0].set_xlabel(r"$\mu$",fontsize=f2)
    axs[0][0].set_ylabel(r"$\rho $",fontsize=f2)
    axs[1][0].set_ylabel(r"$\chi$",fontsize=f2)
    axs[1][1].set_xlabel(r"$\mu$",fontsize=f2)
    axs[0][1].set_ylabel(r"$U / N$",fontsize=f2)
    axs[1][1].set_ylabel(r"$C_v /N$",fontsize=f2)
    
    ys=np.linspace(0.0,1.0,21)
    xticksindex = range(0,len(mus),int(len(mus)/5.0))
    yticksindex = range(0,len(ys),int(len(ys)/6.0))
    xticks = [round(mus[i],2) for i in xticksindex]
    yticks = [ys[i] for i in yticksindex]
    axs[0][0].set_yticks(yticks)
    axs[0][0].set_xticks(xticks)
    axs[1][0].set_xlabel(r"$\mu$",fontsize=f2)
    axs[0][0].set_ylabel(r"$\rho$",fontsize=f2)
    axs[1][0].set_ylabel(r"$\chi$",fontsize=f2)

    ys=np.linspace(-3.0,0.5,21)
    xticksindex = range(0,len(mus),int(len(mus)/5.0))
    yticksindex = range(0,len(ys),int(len(ys)/6.0))
    xticks = [round(mus[i],2) for i in xticksindex]
    yticks = [ys[i] for i in yticksindex]
    axs[0][1].set_yticks(yticks)
    axs[0][1].set_xticks(xticks)
    axs[1][1].set_xlabel(r"$\mu$",fontsize=f2)
    axs[0][1].set_ylabel(r"$U / N$",fontsize=f2)
    axs[1][1].set_ylabel(r"$C_v /N$",fontsize=f2)

    
    for i in range(len(axs)):
        for j in range(len(axs[i])):
            axs[i][j].grid()
            axs[i][j].tick_params(axis="both",labelsize=f1)
    fig.subplots_adjust(wspace=0.28)
    
    if colorbar == "double":
    
        scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
        scalarmappaple.set_array(cs)
        for ax in axs:
            cbar = fig.colorbar(scalarmappaple,ax=ax)
            csaux = [cs[i] for i in range(0,len(cs),2)]
            cbar.set_ticks(csaux)
            cbar.set_ticks(csaux)
            cbar.set_ticklabels([str(x) for x in csaux])
            cbar.set_label(r"$C$",size=f2)
            cbar.ax.tick_params(labelsize=f1)
            plt.tight_layout()

    
    elif colorbar == "single":
        
        cbar_ax = fig.add_axes([0.93,0.15,0.03,0.7])
        scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=color2)
        scalarmappaple.set_array(cs)
        cbar = fig.colorbar(scalarmappaple,cax=cbar_ax)
        csaux = [cs[i] for i in range(0,len(cs),2)]
        cbar.set_ticks(csaux)
        cbar.set_ticklabels([str(x) for x in csaux])
        cbar.set_label(r"$C$",size=f2,labelpad=40,rotation=0)
        cbar.ax.tick_params(labelsize=22)

    else:
        print("Error: color bar not supported")
    plt.savefig("mu-C.png",dpi=300)

groupedPlots(cm.terrain,"single")

