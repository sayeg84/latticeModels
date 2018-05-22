include("inOut.jl")
include("algorithms.jl")
include("statEnsemble.jl")
using Algorithms,InOut,StatEnsemble

function ParseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-N", "--Nlatt" 
            help = "Lattice size"
            arg_type=Int64
            default = 10
         "-J", "--Jconst"
            help = "Coupling constant of Ising model"
            arg_type = Float64
            default = 1.0
         "-B", "--Bfield"
            help = "Magnetic field"
            arg_type = Float64
            default = 0.0
         "-C", "--Cconst"
            help = "Cycle constant"
            arg_type = Float64
            default = 1.0
    end
    return parse_args(s)
end

parsedArgs = ParseCommandline()

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
1]

#making test run for compilation
#initLatt=ones(param[1],param[1])
#neigLatt=Auxiliar.NeighborIndexLattice(initLatt,Auxiliar.SquareLatticeNeighborsIndex)
#Algorithms.WangLandauCycle(param,initLatt,neigLatt,printLog=false)


#start counting time of execution
time=Dates.time()
N=parsedArgs["Nlatt"]
param=[N,parsedArgs["Jconst"],parsedArgs["Bfield"],0,0,convert(Int64,ceil(N^2/2)-N),parsedArgs["Cconst"]]

initLatt=-1*ones(param[1],param[1])
neigLatt=Auxiliar.NeighborIndexLattice(initLatt,Auxiliar.SquareLatticeNeighborsIndex)
println("Simulating")
println()
X=Algorithms.WangLandauCycle(param,initLatt,neigLatt,printLog=false)
energyIntervals=X[1]
s=X[2]
tempArray=linspace(0.1,10,100)
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