println("Ising model")
println()
println()
println("Importing libraries")
println()

using ArgParse, Statistics, Dates, Distributed

Distributed.addprocs(4)

@sync @everywhere include("inOut.jl")
@everywhere include("algorithms.jl")
include("statEnsemble.jl")
include("auxiliar.jl")
include("geometry.jl")

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

#initializing parameters
#bArray=range(-3.0,stop=-1.0,length=21)
jArray=[2.0]
cArray=[0.8]
cArray=range(0.5,stop=1.5,length=21)
tArray=[0.5]

println()
println("Initializing Lattice")
println()

model="cycle"
latt , neigLatt =  Geometry.BuildLattices(geoParam,model)

params=[Array{Float64,1}([bArray[i1],jArray[i2],cArray[i3],tArray[end-i4]]) for i1 in 1:length(bArray), i2 in 1:length(jArray), i3 in 1:length(cArray), i4 in 0:(length(tArray)-1)]

println()
println("Running simulations")
println()

@time for s in params
    simulParam = s
    current="B, J, C, T = $(simulParam) "
    println()
    println(current)
    println()
    
    @sync @everywhere InOut.MakeAndEnterDirectories()
    @sync @distributed for i in 1:algoParam[3]
        X=Algorithms.Metropolis(simulParam,algoParam,latt,neigLatt,model)
        name=string(simulParam,"_",i)
        mkdir(name)
        cd(name)
        InOut.MetropolisOut(X,algoParam)
        InOut.WriteAlgoParamTable(algoParam,"metropolis")
        InOut.WriteSimulParamTable(simulParam)
        InOut.WriteGeoParamTable(geoParam)
        cd("..")
    end
    @everywhere InOut.ExitDirectories()
    
end
