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
        "-B", "--Barray"
            help = "String describing the begining, end and length of the Barray. Reverse order (begining > end) is supported"
            arg_type = String
            default = "-3.5,0.0,21"
        "-J", "--Jarray"
            help = "String describing the begining, end and length of the Jarray. Reverse order (begining > end) is supported"
            arg_type = String
            default = "2.0,2.0,1"
        "-C", "--Carray"
            help = "String describing the begining, end and length of the Carray. Reverse order (begining > end) is supported"
            arg_type = String
            default = "0.2,1.8,21"
        "-T", "--kTarray"
            help = "String describing the begining, end and length of the kTarray. Reverse order (begining > end) is supported"
            arg_type = String
            default = "0.5,0.5,1"
        "-O", "--order"
            help = "String describing the order of parameters to retrive. String must be a permutation of B,J,C,kT separated by commas"
            arg_type = String
            default = "B,J,C,kT"
    end
    return parse_args(s)
end

function TupleToRange(str)
    numbers = split(str,',')
    if length(numbers) != 3
        error("String $(str) non parsable")
    else
        tup = (parse(Float64,numbers[1]),parse(Float64,numbers[2]),parse(Int64,numbers[3]))
        if tup[3] == 1
            return tup[1]
        else
            return range(tup[1],stop = tup[2], length = tup[3])
        end
    end
end

parsedArgs = ParseCommandline()

algoParam=NTuple{2,Int64}((
    floor(10^parsedArgs["sweeps"])*(parsedArgs["Nlatt"]^parsedArgs["dim"]),
    parsedArgs["Systems"]
))

metaParam=(
    parsedArgs["Nlatt"],
    parsedArgs["dim"],
    string(parsedArgs["Geometry"],"LatticeNeighbors"),
    parsedArgs["EnerFunc"],
    parsedArgs["Model"]
)



println()
println("Initializing Lattice")
println()
# Initializing first lattice
neigFunc = getfield(Lattices,Symbol(metaParam[3]))
sysFunc = getfield(Main,Symbol(metaParam[5])) 
initSys  = [ sysFunc(neigFunc,metaParam[1],metaParam[2]) for i in 1:algoParam[2]] 
enerFunc = getfield(StatEnsemble,Symbol( metaParam[4]))

#initializing parameters

arrays = []
varsOrder = split(parsedArgs["order"],',')
for var in varsOrder
    println(var)
    array = TupleToRange(parsedArgs[string(var,"array")])
    println(array)
    push!(arrays,array)
end
orderDict = Dict(varsOrder[i] => i for i in 1:length(varsOrder))
params1 = [simulParam for simulParam in Iterators.product(arrays...)]
params1 = reshape(params1,length(params1))

#params1 = [collect(simulParam[[orderDict["B"],orderDict["J"],orderDict["C"],orderDict["kT"]]]) for simulParam in params1]
#params2 = [vcat(par,[iter]) for iter in 1:algoParam[2] for par in params1]
#display(params2)

# combination of parameters in desired order
params1 = [simulParam[[orderDict["B"],orderDict["J"],orderDict["C"],orderDict["kT"]]] for simulParam in params1]
# adding number of iterations
params2 = [(par...,iter) for iter in 1:algoParam[2] for par in params1]
simulParamDict = Dict(params1[i] => i for i in 1:length(params1))


InOut.MakeAndEnterDirectories()
InOut.WriteAlgoParamTable(algoParam,"metropolis")
InOut.WriteMetaParamTable(metaParam)
InOut.WriteAdjMat(initSys[1])
InOut.WriteSimulParamDict(simulParamDict)

println()
println("Running simulations")
println()


for i in 1:length(params2)
    iter = Int(params2[i][end])
    simulParam = params2[i][1:4]
    name = string(simulParamDict[simulParam],"_",iter)
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
