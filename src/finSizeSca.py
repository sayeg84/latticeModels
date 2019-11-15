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
for i in iters:
        plt.plot(tempArr,dataFrames[i]["Xi"].tolist(),label="L = {0}".format(sizes[i]),marker='o',markersize=6)

        plt.legend()

plt.grid()

plt.show()

mint = min(tempArr)

maxt = max(tempArr)

newTempArr = np.linspace(mint,maxt,100)

plt.figure(figsize=(12,8))

for i in range(len(sizes)):
        tck = inter.splrep(dataFrames[i]["kT"].tolist(),dataFrames[i]["Xi"].tolist())
        
        ys = inter.splev(newTempArr,tck,der=0)
        
        plt.plot(newTempArr,ys,label="L = {0}".format(sizes[i]))
        
        plt.scatter(dataFrames[i]["kT"].tolist(),dataFrames[i]["Xi"].tolist())
        
        plt.legend()
        
plt.grid()

plt.show()

"""
def distance(arr):
        Xis = [np.array(df["Xi"].tolist()) for df in dataFrames]
        Xis = [Xis[i]*np.power(sizes[i],-arr[0]/arr[1]) for i in range(len(sizes))]
        return np.linalg.norm([np.linalg.norm(xi1 - xi2,2) for xi1 in Xis for xi2 in Xis],2)
def distance2(arr):
        iters = range(len(sizes))
        Xis = [np.array(df["Xi"].tolist()) for df in dataFrames]
        Xis = [Xis[i]*np.power(sizes[i],-arr[0]/arr[1]) for i in iters]
        #maxindxs  = [dataFrames[i]["kT"].tolist().index(2.26) for i in iters]
        maxindxs  = [np.argmax(Xis[i]) for i in iters]
        scales = [[(temp-arr[2])*np.power(sizes[i],1/arr[1]) for temp in df["kT"].tolist()] for i in iters]
        maxpoints = [np.array([scales[i][maxindxs[i]],Xis[i][maxindxs[i]]]) for i in iters]
        vect = [np.linalg.norm(m1-m2,2) for m1 in maxpoints for m2 in maxpoints]
        return np.linalg.norm(vect,2)
"""
def functionNorm(f1,f2,xmin,xmax,n=101):
        xs = np.linspace(xmin,xmax,n)
        
        f1s = np.array(f1(xs))
        
        f2s = np.array(f2(xs))
        
        return np.linalg.norm(f1s-f2s,np.Inf)
def distance(arr):
        iters = range(len(sizes))
        
        # finding spline interpolation functions fo all data
        fxis = [inter.splrep(dataFrames[i]["kT"].tolist(),dataFrames[i]["Xi"].tolist()) for i in iters]

        # building domain of all points
        scales = np.array([[(temp-arr[2])*np.power(sizes[i],1/arr[1]) for temp in newTempArr] for i in iters])
        
        # building universal domain
        tmin = np.min(scales)
        tmax = np.min(scales)
        newTemp = np.linspace(tmin,tmax,101)

        sizedTemps = [x*np.power(sizes[i],-1/arr[1])+arr[2] for x in newTemp]
        xis =  [np.array(inter.splev(sizedTemps,fxis[i],der=0)*np.power(sizes[i],-arr[0]/arr[1])) for i in iters]
        
        vect = [np.linalg.norm(xi1-xi2,2) for xi1 in xis for xi2 in xis]
        return max(vect)
#distance(np.array([1.0,1.0,1.0]))
A = opt.minimize(distance,np.array([1.5,1.0,2.6]))
print(A)
#gamma = A.x[0]
#nu = A.x[1]
#tc = A.x[2]
gamma = 1.75
nu = 1.0
tc=2.26

# building domain of all points
scales = np.array([[(temp-tc)*np.power(sizes[i],1/nu) for temp in newTempArr] for i in iters])

# building universal domain
tmin = np.min(scales)
tmax = np.max(scales)
newTemp = np.linspace(tmin,tmax,101)
print tmin
print tmax 

plt.figure(figsize=(14,8))
for i in iters:
        fxis = [inter.splrep(dataFrames[j]["kT"].tolist(),dataFrames[j]["Xi"].tolist()) for j in iters]
        sizedTemps = [x*np.power(sizes[i],-1/nu)+tc for x in newTemp]
        Ys = np.array(inter.splev(sizedTemps,fxis[i],der=0))*np.power(sizes[i],-gamma/nu)
        plt.plot(newTemp,Ys,label="L = {0}".format(sizes[i]),marker='o',ms=6)
        plt.ylim(0,0.05)
        plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(14,8))
for i in iters:
        Ys = np.array(dataFrames[i]["Xi"].tolist())
        Ys = Ys * np.power(sizes[i],-gamma/nu)
        Xs = [(temp-tc)*np.power(sizes[i],1/nu) for temp in tempArr]
        plt.plot(Xs,Ys,label="L = {0}".format(sizes[i]),marker='o',ms=6)
        plt.ylim()
        plt.legend()
plt.grid()
plt.show()
