println("Metropolis ising model")
println()
println()
println("Importing libraries")
println()

include("inOut.jl")
include("algorithms.jl")
include("statEnsemble.jl")
include("auxiliar.jl")
include("geometry.jl")
#using Algorithms,InOut,StatEnsemble, ArgParse
using ArgParse, Statistics, Dates

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
            default = "square"
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
        "-S", "--steps"
            help = "logarithm base 10 of total steps"
            arg_type = Float64
            default = 6.0
        "-F", "--frequency"
            help = "logarithm base 10 of saving frecuency. Must be less than steps"
            arg_type = Float64
            default = 4.0
        "-A", "--averages"
            help = "Number of averages performed"
            arg_type = Int64
            default = 5
    end
    return parse_args(s)
end

parsedArgs = ParseCommandline()

algoParam=Array{Int64,1}([
    floor(10^parsedArgs["steps"]),
    floor(10^parsedArgs["frequency"]),
    parsedArgs["averages"]
])

simulParam=Array{Float64,1}([
    parsedArgs["Bfield"],
    parsedArgs["Jconst"],
    parsedArgs["Cconst"],
    parsedArgs["temp"]
])

geoParam=[
    parsedArgs["Nlatt"],
    parsedArgs["dim"],
    parsedArgs["Geometry"]
]

println()
println("Initializing Lattice")
println()
# Initializing first lattice

latt , neigLatt =  Geometry.BuildLattices(geoParam,"cycle")

println()
println("Making simulation")
println()
# making simulation
resul=[]
for i in 1:algoParam[3]
    X=Algorithms.Metropolis(simulParam,algoParam,latt,neigLatt,"cycle")
    push!(resul,X)
    latt=copy(X[end])
end

println()
println("Writing output")
println()

InOut.MakeAndEnterDirectories()
name="Single"
mkdir(name)
cd(name)
InOut.MetropolisOut(mean(resul),algoParam)
InOut.WriteAlgoParamTable(algoParam,"metropolis")
InOut.WriteSimulParamTable(simulParam)
InOut.WriteGeoParamTable(geoParam)
cd("..")
InOut.ExitDirectories()



