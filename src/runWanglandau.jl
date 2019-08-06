println("Wang-Landau ising model")
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
            default = "square"
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
            default = 1.0
        "-T", "--temp"
            help = "Temperature"
            arg_type = Float64
            default = 5.0
        "-P", "--flatness"
            help = "percentage of the average that the minimum must overcome for the histogram to be flat"
            arg_type = Float64
            default = 0.5
        "-F", "--Ffactor"
            help = "Modification factor for histogram"
            arg_type = Float64
            default = 1.0
        #=
        "-M", "--MaxSteps"
            help = "Logarithm base 10 of the maximum number of steps "
            arg_type = Float64
            default = 6.5
        =#
    end
    return parse_args(s)
end

parsedArgs = ParseCommandline()



simulParam=Array{Float64,1}([
    parsedArgs["Bfield"],
    parsedArgs["Jconst"],
    parsedArgs["Cconst"],
    parsedArgs["temp"]
])

metaParam=[
    parsedArgs["Nlatt"],
    parsedArgs["dim"],
    parsedArgs["Geometry"],
    parsedArgs["EnerFunc"],
    parsedArgs["Model"]
]

algoParam=Array{Float64,1}([
    ceil((metaParam[1])*(metaParam[1]/2 - 1)),
    parsedArgs["flatness"],
    parsedArgs["Ffactor"],
    10^parsedArgs["MaxSteps"]
])

println()
println("Initializing Lattice")
println()
# Initializing first lattice



initSys  = getfield(Main,Symbol(metaParam[5]))(Lattices.PeriodicSquareLatticeNeighbors,metaParam[1],metaParam[2])
println()
println("Making simulation")
println()

X = Algorithms.WangLandauOptimal(simulParam,algoParam,metaParam,initSys,printLog=false)
energyIntervals=X[1]
s=X[2]
mag=X[3]



println()
println("Writing Output")
println()


InOut.MakeAndEnterDirectories()
name="Single"
if ~(isdir(name))
    mkdir(name)
end
cd(name)
InOut.WriteDOSTable(s,mag,energyIntervals)
InOut.WriteAlgoParamTable(algoParam,"WL")
InOut.WriteSimulParamTable(simulParam)
InOut.WriteMetaParamTable(metaParam)