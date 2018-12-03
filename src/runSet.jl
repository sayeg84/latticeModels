println("Ising model")
println()
println()
println("Importing libraries")
println()
using ArgParse, Statistics, Dates, Distributed, SharedArrays

include("inOut.jl")
include("algorithms.jl")
include("statEnsemble.jl")
include("auxiliar.jl")
include("geometry.jl")
#using Algorithms,InOut,StatEnsemble, ArgParse


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
geoParam=[
    parsedArgs["Nlatt"],
    parsedArgs["dim"],
    parsedArgs["Geometry"]
]
model="cycle"
latt , neigLatt =  Geometry.BuildLattices(geoParam,model)


println()
println("Initializing Lattice")
println()
# Initializing first lattice

#initializing parameters
bArray=[-1.0]
jArray=[2.0]
cArray=[0.8]
tArray=[0.5]


println()
println("Running simulations")
println()

params=[Array{Float64,1}([bArray[i1],jArray[i2],cArray[i3],tArray[end-i4]]) for i1 in 1:length(bArray), i2 in 1:length(jArray), i3 in 1:length(cArray), i4 in 0:(length(tArray)-1)]

for simulParam in params
    current="B, J, C, T = $(simulParam) "
    println()
    println(current)
    println()
    InOut.MakeAndEnterDirectories()
    for i in 1:algoParam[3]
        println(i)
        X=Algorithms.Metropolis(simulParam,algoParam,latt,neigLatt,"cycle")
        latt=copy(X[end])
        name=string(simulParam,"_",i)
        mkdir(name)
        cd(name)
        InOut.MetropolisOut(X,algoParam)
        InOut.WriteAlgoParamTable(algoParam,"metropolis")
        InOut.WriteSimulParamTable(simulParam)
        InOut.WriteGeoParamTable(geoParam)
        cd("..")
    end
    InOut.ExitDirectories()    
    
end
