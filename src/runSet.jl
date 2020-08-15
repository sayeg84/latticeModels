println("Ising model")
println()
println()
println("Importing libraries")
println()

using ArgParse, Statistics, Dates

include("lattices.jl")
include("algorithms.jl")
include("inOut.jl")



println()
println("Parsing Arguments")
println()
function ParseCommandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "-N", "--Nlatt" 
            help = "Lattice size"
            arg_type = Int64
            default = 10
        "-D", "--dim" 
            help = "Dimension"
            arg_type = Int64
            default = 2
        "-G", "--Geometry" 
            help = "Lattice geometry"
            arg_type = String
            default = "PeriodicSquare"
        "-E", "--EnerFunc" 
            help = "Energy function"
            arg_type = String
            default = "NormalEnergy"
        "-M", "--Model" 
            help = "Lattice geometry"
            arg_type = String
            default = "SpinLattice"
        "-S", "--sweeps"
            help = "logarithm base 10 of total sweeps"
            arg_type = Float64
            default = 2.0
        "-A", "--Systems"
            help = "Number of independent systems simulated for same simulation parameters"
            arg_type = Int64
            default = 2
    end
    return parse_args(s)
end

parsedArgs = ParseCommandline()

algoParam=Array{Int64,1}([
    floor(10^parsedArgs["sweeps"])*(parsedArgs["Nlatt"]^parsedArgs["dim"]),
    parsedArgs["Systems"]
])
metaParam=[
    parsedArgs["Nlatt"],
    parsedArgs["dim"],
    string(parsedArgs["Geometry"],"LatticeNeighbors"),
    parsedArgs["EnerFunc"],
    parsedArgs["Model"]
]



println()
println("Initializing Lattice")
println()
# Initializing first lattice
neigFunc = getfield(Lattices,Symbol(metaParam[3]))
sysFunc = getfield(Main,Symbol(metaParam[5])) 
initSys  = [ sysFunc(neigFunc,metaParam[1],metaParam[2]) for i in 1:algoParam[2]] 
enerFunc = getfield(StatEnsemble,Symbol( metaParam[4]))

#initializing parameters

bArray = [0,2,4,6,8]
jArray = [1]
cArray = [0]
tArray = [0,1,2,3]

#bArray = range(-3.3,stop = -1.5,length = 31)
#jArray = [2.0]
#cArray = [0.8]
#tArray = [0.5]



params1 = [Array{Float64,1}([bArray[i1],jArray[i2],cArray[i3],tArray[end-i4]])  for i3 in 1:length(cArray) for i1 in 1:length(bArray) for i4 in 0:(length(tArray)-1) for i2 in 1:length(jArray)]

params2 = [vcat(par,[iter]) for iter in 1:algoParam[2] for par in params1]

paramDict = Dict(params1[i] => i for i in 1:length(params1))


InOut.MakeAndEnterDirectories()
InOut.WriteAlgoParamTable(algoParam,"metropolis")
InOut.WriteMetaParamTable(metaParam)
InOut.WriteAdjMat(initSys[1])
InOut.WriteSimulParamDict(paramDict)

println()
println("Running simulations")
println()


for i in 1:length(params2)
    iter = Int(params2[i][end])
    simulParam = params2[i][1:4]
    name = string(paramDict[simulParam],"_",iter)
    print("Running ")
    println(name)
    println(simulParam)
    @time res = Algorithms.MetropolisFast(initSys[iter],enerFunc,simulParam,algoParam)
    InOut.MetropolisAllOut(initSys[iter],res,name)
    initSys[iter] = copy(res[2])
end

#=
@time for simulParam in params
    current="B, J, C, T = $(simulParam) "
    println()
    println(current)
    println()
    InOut.MakeAndEnterDirectories()
    InOut.WriteAlgoParamTable(algoParam,"metropolis")
    InOut.WriteMetaParamTable(metaParam)
    InOut.WriteAdjMat(initSys[1])
    println()
    println("Making simulation")
    println()
    # making simulation
    for i in 1:algoParam[2]
        println(i)
        global initSys
        res = Algorithms.MetropolisOptimal(initSys[i],enerFunc,simulParam,algoParam)
        name=string(simulParam,"_",i)
        mkdir(name)
        cd(name)
        InOut.MetropolisAllOut(initSys[i],res[1],algoParam)
        InOut.WriteSimulParamTable(simulParam)
        initSys[i] = copy(res[2])
        cd("..")
    end
    InOut.ExitDirectories()
end    
=#