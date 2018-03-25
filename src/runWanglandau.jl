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
N=30
push!(param,N)
push!(param,1)
push!(param,0)
push!(param,10^6)
push!(param,10^4)
push!(param,N^2/2)

#initLatt=rand([-1,1],param[1],param[1])
initLatt=ones(param[1],param[1])
X=Algorithms.WangLandau(param,initLatt,test=false)
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

println()
println("Writing Output")
InOut.MakeDirectories()
InOut.MakePlots(tempArray,mag,ener,cv,param)
time=Dates.time()-time
InOut.MakeTable(tempArray,mag,ener,cv,param,time)
InOut.MakeDOSTables(s,energyIntervals,param,time)