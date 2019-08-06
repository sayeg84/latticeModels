println("Ising model")
println()
println()
println("Importing libraries")
println()

using  Distributed
@everywhere using ArgParse, Statistics, Dates

Distributed.addprocs(4)


@sync @everywhere include("inOut.jl")
@everywhere include("algorithms.jl")
@everywhere include("lattices.jl")


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
        "-S", "--steps"
            help = "logarithm base 10 of total steps"
            arg_type = Float64
            default = 5.0
        "-F", "--frequency"
            help = "logarithm base 10 of saving frecuency. Must be less than steps"
            arg_type = Float64
            default = 3.5
        "-A", "--averages"
            help = "Number of averages performed"
            arg_type = Int64
            default = 8
    end
    return parse_args(s)
end

parsedArgs = ParseCommandline()



algoParam=Array{Int64,1}([
    floor(10^parsedArgs["steps"]),
    floor(10^parsedArgs["frequency"]),
    parsedArgs["averages"]
])

metaParam=[
    parsedArgs["Nlatt"],
    parsedArgs["dim"],
    parsedArgs["Geometry"],
    parsedArgs["EnerFunc"],
    parsedArgs["Model"]
]

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
println("Initializing Lattice")
println()



let y  = getfield(Main,Symbol(metaParam[5]))(Lattices.PeriodicSquareLatticeNeighbors,metaParam[1],metaParam[2])
    @sync @everywhere initSys = $y
end

let z  = getfield(StatEnsemble,Symbol( metaParam[4]))
    @sync @everywhere enerFunc = $z
end


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
        global initSys
        res = Algorithms.MetropolisFast(initSys,enerFunc,simulParam,algoParam)
        name=string(simulParam,"_",i)
        mkdir(name)
        cd(name)
        InOut.MetropolisAllOut(initSys,res[1],algoParam)
        InOut.WriteAlgoParamTable(algoParam,"metropolis")
        InOut.WriteSimulParamTable(simulParam)
        InOut.WriteMetaParamTable(metaParam)
        InOut.WriteAdjMat(initSys)
        initSys = copy(res[2])
        cd("..")
    end
    @everywhere InOut.ExitDirectories()
    
end
