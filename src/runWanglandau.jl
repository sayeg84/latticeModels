time=Dates.time()
include("inOut.jl")
include("algorithms.jl")
include("statEnsemble.jl")
using Algorithms,InOut,StatEnsemble

println("Begining")
println()
println("Simulating")
param=[]
#asignamos los parámetros
N=16
push!(param,N)
push!(param,1)
push!(param,0)
push!(param,0)
push!(param,0)
push!(param,convert(Int64,ceil(N^2/2)-N))

#initLatt=rand([-1,1],param[1],param[1])
initLatt=ones(param[1],param[1])
neigLatt=Auxiliar.NeighborIndexLattice(initLatt,Auxiliar.SquareLatticeNeighborsIndex)
X=Algorithms.WangLandauSimple(param,initLatt,neigLatt,test=false)
energyIntervals=X[1]
s=X[2]
tempArray=linspace(0.1,5,50)
ener=[]
cv=[]
mag=[]
for temp in tempArray
    push!(ener,StatEnsemble.DOSEnergy(s,energyIntervals,temp,param)/(param[1]^2))
    push!(cv,StatEnsemble.DOSCV(s,energyIntervals,temp,param)/(param[1]^2))
    push!(mag,StatEnsemble.DOSMag(s,energyIntervals,temp,param)/(param[1]^2))

end
param[4]=1
param[4]=X[4]
println()
println("Writing Output")
InOut.MakeDirectories()
InOut.MakePlots(tempArray,mag,ener,cv,param)
time=Dates.time()-time
InOut.MakeTable(tempArray,mag,ener,cv,param,time)
InOut.MakeDOSTable(s,energyIntervals,param,time)