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
push!(param,10)
push!(param,1)
push!(param,0)
push!(param,10^6)
push!(param,10^4)

initLatt=rand([-1,1],param[1],param[1])
X=Algorithms.WangLandau(param,initLatt,10,test=true)
energyIntervals=X[1]
s=X[2]
tempArray=linspace(0.1,5,50)
ener=[]
cv=[]
for temp in tempArray
    push!(ener,StatEnsemble.DOSEnergy(s,energyIntervals,temp,param))
    push!(cv,StatEnsemble.DOSCV(s,energyIntervals,temp,param))
end

mag=zeros(50)
println()
println("Writing Output")
InOut.MakeDirectories()
InOut.MakePlots(tempArray,mag,ener,cv,param)
time=Dates.time()-time
InOut.MakeTable(tempArray,mag,ener,cv,param,time)