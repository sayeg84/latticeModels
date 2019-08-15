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

    @add_arg_table s begin
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
        "-S", "--steps"
            help = "logarithm base 10 of total steps"
            arg_type = Float64
            default = 5.0
            "-A", "--Systems"
            help = "Number of independent systems simulated for same simulation parameters"
            arg_type = Int64
            default = 2
    end
    return parse_args(s)
end

parsedArgs = ParseCommandline()

algoParam=Array{Int64,1}([
    floor(10^parsedArgs["steps"]),
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
#bArray = [0]
#jArray = [1]
#cArray = [0]
#tArray = range(0.1 , stop = 5 , length = 21)

bArray = range(-3.3,stop = -1.5,length = 31)
jArray = [2.0]
cArray = [0.8]
tArray = [0.5]


println()
println("Running simulations")
println()

params=[Array{Float64,1}([bArray[i1],jArray[i2],cArray[i3],tArray[end-i4]]) for i1 in 1:length(bArray), i2 in 1:length(jArray), i3 in 1:length(cArray), i4 in 0:(length(tArray)-1)]


@time for simulParam in params
    current="B, J, C, T = $(simulParam) "
    println()
    println(current)
    println()
    InOut.MakeAndEnterDirectories()
    println()
    println("Making simulation")
    println()
    # making simulation
    for i in 1:algoParam[2]
        println(i)
        global initSys
        res = Algorithms.MetropolisFast(initSys[i],enerFunc,simulParam,algoParam)
        name=string(simulParam,"_",i)
        mkdir(name)
        cd(name)
        InOut.MetropolisAllOut(initSys[i],res[1],algoParam)
        InOut.WriteAlgoParamTable(algoParam,"metropolis")
        InOut.WriteSimulParamTable(simulParam)
        InOut.WriteMetaParamTable(metaParam)
        InOut.WriteAdjMat(initSys[i])
        initSys[i] = copy(res[2])
        cd("..")
    end
    InOut.ExitDirectories()
end    