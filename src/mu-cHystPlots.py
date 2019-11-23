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


df1 = pd.read_csv("mu-increasing/Macroscopic.csv")
mus = sorted(list(set(df1["Mu"].tolist())))
cs = sorted(list(set(df1["C"].tolist())))

indexes=[[(mus[i],cs[j]) for i in range(len(mus))] for j in range(len(cs))]
mag1=[[df1[df1["Mu"]==mus[i]][df1["C"]==cs[j]]["M"].values[0] for i in range(len(mus))] for j in range(len(cs))]
ener1=[[df1[df1["Mu"]==mus[i]][df1["C"]==cs[j]]["E"].values[0] for i in range(len(mus))] for j in range(len(cs))]
cv1=[[df1[df1["Mu"]==mus[i]][df1["C"]==cs[j]]["CV"].values[0] for i in range(len(mus))] for j in range(len(cs))]
chi1=[[df1[df1["Mu"]==mus[i]][df1["C"]==cs[j]]["Xi"].values[0] for i in range(len(mus))] for j in range(len(cs))]
mags1=[[df1[df1["Mu"]==mus[i]][df1["C"]==cs[j]]["Msigma"].values[0] for i in range(len(mus))] for j in range(len(cs))]
eners1=[[df1[df1["Mu"]==mus[i]][df1["C"]==cs[j]]["Esigma"].values[0] for i in range(len(mus))] for j in range(len(cs))]
cvs1=[[df1[df1["Mu"]==mus[i]][df1["C"]==cs[j]]["CVsigma"].values[0] for i in range(len(mus))] for j in range(len(cs))]
chis1=[[df1[df1["Mu"]==mus[i]][df1["C"]==cs[j]]["Xisigma"].values[0] for i in range(len(mus))] for j in range(len(cs))]

df2 = pd.read_csv("mu-decreasing/Macroscopic.csv")
mus = sorted(list(set(df2["Mu"].tolist())))
cs = sorted(list(set(df2["C"].tolist())))

indexes=[[(mus[i],cs[j]) for i in range(len(mus))] for j in range(len(cs))]
mag2=[[df2[df2["Mu"]==mus[i]][df2["C"]==cs[j]]["M"].values[0] for i in range(len(mus))] for j in range(len(cs))]
ener2=[[df2[df2["Mu"]==mus[i]][df2["C"]==cs[j]]["E"].values[0] for i in range(len(mus))] for j in range(len(cs))]
cv2=[[df2[df2["Mu"]==mus[i]][df2["C"]==cs[j]]["CV"].values[0] for i in range(len(mus))] for j in range(len(cs))]
chi2=[[df2[df2["Mu"]==mus[i]][df2["C"]==cs[j]]["Xi"].values[0] for i in range(len(mus))] for j in range(len(cs))]
mags2=[[df2[df2["Mu"]==mus[i]][df2["C"]==cs[j]]["Msigma"].values[0] for i in range(len(mus))] for j in range(len(cs))]
eners2=[[df2[df2["Mu"]==mus[i]][df2["C"]==cs[j]]["Esigma"].values[0] for i in range(len(mus))] for j in range(len(cs))]
cvs2=[[df2[df2["Mu"]==mus[i]][df2["C"]==cs[j]]["CVsigma"].values[0] for i in range(len(mus))] for j in range(len(cs))]
chis2=[[df2[df2["Mu"]==mus[i]][df2["C"]==cs[j]]["Xisigma"].values[0] for i in range(len(mus))] for j in range(len(cs))]


def arrowPlot(axes,xs,ys,inverse = False,f=4,color="k"):
    if inverse:
         xs = [xs[-(i+1)] for i in range(len(xs))]
         ys = [ys[-(i+1)] for i in range(len(ys))]
    for i in range(len(xs)-1):
        if i % f ==0:
            axes.arrow(xs[i],ys[i],xs[i+1]-xs[i],ys[i+1]-ys[i],head_width = 0.06,head_length = 0.13,zorder=1,color = color)
        else:
            axes.plot([xs[i],xs[i+1]],[ys[i],ys[i+1]],color = color)

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

        
        arrowPlot(axs[0][0],mus,mag1[c],color=colormap(normalize(cs[c])))
        
        #axs[0][0].errorbar(mus,mag1[c],yerr=mags1[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
        
        arrowPlot(axs[1][0],mus,chi1[c],color=colormap(normalize(cs[c])))
        
        #axs[1][0].errorbar(mus,chi1[c],yerr=chis1[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
        
        arrowPlot(axs[0][1],mus,ener1[c],color=colormap(normalize(cs[c])))
        
        #axs[0][1].errorbar(mus,ener1[c],yerr=eners1[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
        
        arrowPlot(axs[1][1],mus,cv1[c],color=colormap(normalize(cs[c])))
        
        #axs[1][1].errorbar(mus,cv1[c],yerr=cvs1[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)

        ##other process

        arrowPlot(axs[0][0],mus,mag2[c],inverse=True,color=colormap(normalize(cs[c])))
        
        #axs[0][0].errorbar(mus,mag2[c],yerr=mags2[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
        
        arrowPlot(axs[1][0],mus,chi2[c],inverse=True,color=colormap(normalize(cs[c])))
        
        #axs[1][0].errorbar(mus,chi2[c],yerr=chis2[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
        
        arrowPlot(axs[0][1],mus,ener2[c],inverse=True,color=colormap(normalize(cs[c])))
        
        #axs[0][1].errorbar(mus,ener2[c],yerr=eners2[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
        
        arrowPlot(axs[1][1],mus,cv2[c],inverse=True,color=colormap(normalize(cs[c])))
        
        #axs[1][1].errorbar(mus,cv2[c],yerr=cvs2[c],label = r"$C = {0}$".format(cs[c]),color=colormap(normalize(cs[c])), lw = 2)
    
    
    

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

