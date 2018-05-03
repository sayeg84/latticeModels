include("inOut.jl")
include("algorithms.jl")
include("statEnsemble.jl")
using Algorithms,InOut,StatEnsemble

println("Begining")
println()

#### parameter asignation

#square lattice size
N=4

param=[N,
# coupling constant J
1,
# magnetic field B J
0,
# maximum steps
10^10,
# save frecuency (Metropolis) or final iterations (Wang Landau)
0,
# energy bins (WangLandau)
convert(Int64,ceil(N^2/2)-N),
#cycle constant
0.3]

#making test run for compilation
#initLatt=ones(param[1],param[1])
#neigLatt=Auxiliar.NeighborIndexLattice(initLatt,Auxiliar.SquareLatticeNeighborsIndex)
#Algorithms.WangLandauCycle(param,initLatt,neigLatt,printLog=false)


println("Simulation")
println()
#start counting time of execution
time=Dates.time()
N=10
param=[N,1,0,0,0,convert(Int64,ceil(N^2/2)-N),0.8]

initLatt=ones(param[1],param[1])
neigLatt=Auxiliar.NeighborIndexLattice(initLatt,Auxiliar.SquareLatticeNeighborsIndex)
println("Simulating")

X=Algorithms.WangLandauCycle(param,initLatt,neigLatt)
energyIntervals=X[1]
s=X[2]
tempArray=linspace(0.1,20,200)
ener=[]
cv=[]
mag=[]
for temp in tempArray
    push!(ener,StatEnsemble.DOSEnergy(s,energyIntervals,temp,param)/(param[1]^2))
    push!(cv,StatEnsemble.DOSCV(s,energyIntervals,temp,param)/(param[1]^2))
    push!(mag,StatEnsemble.DOSMag(s,energyIntervals,X[3],StatEnsemble.PenalizedEnergy,temp,param)/(param[1]^2))

end

param[4]=X[5]
println()
println("Writing Output")
InOut.MakeDirectories()
InOut.MakePlots(tempArray,mag,ener,cv,param)
time=Dates.time()-time
InOut.MakeTable(tempArray,mag,ener,cv,param,time)
InOut.MakeDOSTable(s,energyIntervals,param,time)