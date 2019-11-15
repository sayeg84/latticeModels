import numpy as np
import pandas as pd
import scipy.optimize as opt
import scipy.interpolate as inter
import matplotlib.pyplot as plt
import os
sizes = []
dataFrames = []

for name in os.listdir("../finiteSize"):
    sizes.append(int(name.split(".")[0]))
    df = pd.read_csv("../finiteSize/"+name,header=0)
    df = df.sort_values(by="kT")
    dataFrames.append(df)

iters = range(len(sizes))
tempArr = dataFrames[1]["kT"].tolist()
plt.figure(figsize=(12,8))
gammas = []
for i in iters:
#    if sizes[i]==64: 
        xi = dataFrames[i]["Xi"].tolist()
        m = dataFrames[i]["M"].tolist()
        ic = np.argmax(xi)
        ic = tempArr.index(2.26)
        #tc  = tempArr[ic]
        tc=2.26
        newTemp = [ np.abs((tempArr[j+ic+1] - tc) / tc) for j in range(len(tempArr) - ic-1)]
        newTemp = np.log(newTemp)
        newXi = np.log(xi[(ic+1):])
        poly = np.polyfit(newTemp,newXi,deg=1)
        gamma = poly[0]
        gammas.append(gamma)
        plt.scatter(newTemp,newXi, label = "L = {0}".format(sizes[i]))
        plt.plot(newTemp,[gamma*x + poly[1] for x in newTemp])
plt.legend()
plt.grid()
plt.show()
print np.mean(gammas)


betas = []
for i in iters:
#    if sizes[i]==64: 
        m = dataFrames[i]["M"].tolist()
        ic = np.argmax(xi)
        ic = tempArr.index(2.26)
        #tc  = tempArr[ic]
        tc = 2.26
        newTemp = [ np.abs((tempArr[j] - tc) / tc) for j in range(ic)]
        newTemp = np.log(newTemp)
        newXi = np.log(m[0:ic])
        poly = np.polyfit(newTemp,newXi,deg=1)
        beta = poly[0]
        betas.append(beta)
        plt.scatter(newTemp,newXi, label = "L = {0}".format(sizes[i]))
        plt.plot(newTemp,[beta*x + poly[1] for x in newTemp])
plt.legend()
plt.grid()
plt.show()

print np.mean(betas)


plt.figure(figsize=(12,8))
alphas = []
for i in iters:
#    if sizes[i] == 64:
        cv = dataFrames[i]["CV"].tolist()
        ic = np.argmax(cv)
        tc = tempArr[ic]
        newCv = np.log(cv[0:ic])
        newTemp = [ np.abs((tempArr[j] - tc) / tc) for j in range(ic)]
        newTemp = np.log(newTemp)
        print newTemp
        #poly = np.polyfit(newTemp,newCv,deg=1)
        plt.scatter(newTemp,newCv, label = "L = {0}".format(sizes[i]))
        #plt.plot(newTemp,[alpha*x + poly[1] for x in newTemp])
plt.legend()
plt.grid()
plt.show()


plt.figure(figsize=(12,8))
alphas = []
for i in iters:
#    if sizes[i] == 64:
        cv = dataFrames[i]["CV"].tolist()
        ic = np.argmax(xi)
        tc  = tempArr[ic]
        newTemp = [ np.abs((tempArr[j] - tc) / tc) for j in range(ic)]
        newTemp = np.log(newTemp)
        newcv = np.log(cv[0:ic])
        poly = np.polyfit(newTemp,newcv,deg=1)
        alpha = poly[0]
        alphas.append(alpha)
        plt.scatter(newTemp,newcv, label = "L = {0}".format(sizes[i]))
        plt.plot(newTemp,[alpha*x + poly[1] for x in newTemp])
plt.legend()
plt.grid()
plt.show()
print np.mean(alphas)
#Finite Size Scaling method

maxtemps = []
for i in iters:
    xi = dataFrames[i]["Xi"].tolist()
    ic = np.argmax(xi)
    tc = tempArr[ic]
    maxtemps.append(tc)
def distance(nu):
    newSizes =[float(np.float_power(sizes[j],-1/nu)) for j in iters]
    poly = np.polyfit(newSizes,maxtemps,deg=1)
    ys = [poly[0]*y + poly[1] for y in newSizes]
    R = np.corrcoef(maxtemps,ys)[0][1]
    return 1-R
A = opt.minimize(distance,2.0)
print A