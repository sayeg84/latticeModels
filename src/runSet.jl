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
    10000,
    100,
    10
])

geoParam=[
    10,
    parsedArgs["dim"],
    parsedArgs["Geometry"]
]

println()
println("Initializing Lattice")
println()
# Initializing first lattice
model="Cycle"
latt , neigLatt =  Geometry.BuildLattices(geoParam,model)
#initializing parameters
bArray=range(-3,stop=-2,length=21)
#bArray=[-3.0]
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

    resul=[]
    for i in 1:algoParam[3]
        println(i)
        X=Algorithms.Metropolis(simulParam,algoParam,latt,neigLatt,model)
        push!(resul,X)
        latt=deepcopy(X[end])
    end

    InOut.MakeAndEnterDirectories()
    mkdir(string(simulParam))
    cd(string(simulParam))
    InOut.MetropolisOut(mean(resul),algoParam)
    InOut.WriteAlgoParamTable(algoParam,"metropolis")
    InOut.WriteSimulParamTable(simulParam)
    InOut.WriteGeoParamTable(geoParam)
    cd("..")
    InOut.ExitDirectories()
end
    