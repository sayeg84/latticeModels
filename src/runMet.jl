println("Metropolis ising model")
println()
println()
println("Importing libraries")
println()

include("lattices.jl")
include("algorithms.jl")
include("inOut.jl")

using ArgParse, Statistics, Dates

println()
println("Parsing Arguments")
println()
function ParseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-N", "--Nlatt" 
            help = "Lattice side"
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
        "-B", "--Bfield"
            help = "Magnetic field"
            arg_type = Float64
            default = 0.0
         "-J", "--Jconst"
            help = "Coupling constant of Ising model"
            arg_type = Float64
            default = 1.0
        "-C", "--Cconst"
            help = "Cycle constant"
            arg_type = Float64
            default = 0.0
        "-T", "--temp"
            help = "Temperature"
            arg_type = Float64
            default = 5.0
        "-S", "--sweeps"
            help = "logarithm base 10 of total sweeps"
            arg_type = Float64
            default = 4.0
        "-A", "--Systems"
            help = "Number of independent systems simulated for same simulation parameters"
            arg_type = Int64
            default = 2
    end
    return parse_args(s)
end

parsedArgs = ParseCommandline()

algoParam = Array{Int64,1}([
    floor(10^parsedArgs["sweeps"])*parsedArgs["Nlatt"]^parsedArgs["dim"],
    parsedArgs["Systems"]
])

simulParam = Array{Float64,1}([
    parsedArgs["Bfield"],
    parsedArgs["Jconst"],
    parsedArgs["Cconst"],
    parsedArgs["temp"]
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
    res = Algorithms.MetropolisFast(initSys[i],enerFunc,simulParam,algoParam)
    name=string(simulParam,"_",i)
    mkdir(name)
    cd(name)
    InOut.MetropolisAllOut(initSys[i],res[1],algoParam)
    InOut.WriteSimulParamTable(simulParam)
    initSys[i] = copy(res[2])
    cd("..")
end
InOut.ExitDirectories()


