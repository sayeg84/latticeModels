println("Ising model")
println()

using  Distributed
using DelimitedFiles

cores = 3
Distributed.addprocs(cores)

@everywhere using ArgParse, Statistics, Dates
@sync @everywhere include("inOut.jl")
@everywhere include("algorithms.jl")
@everywhere include("lattices.jl")



println()
println("Parsing Arguments")
println()
function ParseCommandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "-N", "--Nlatt" 
            help = "Lattice side"
            arg_type = Int64
            default = 10
        "-D", "--Dim" 
            help = "Dimension for square lattice"
            arg_type = Int64
            default = 2
        "-L", "--Network" 
            help = "Network for running the model. Must be direction of a CSV file containing an adjacency matrix or a lattice-creator function of Geometry class"
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
        "-S", "--Sweeps"
            help = "logarithm base 10 of total sweeps"
            arg_type = Float64
            default = 3.0
        "-A", "--Systems"
            help = "Number of independent systems simulated for same simulation parameters"
            arg_type = Int64
            default = cores
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
            return range(tup[1], stop = tup[2], length = tup[3])
        end
    end
end

const parsedArgs = ParseCommandline()

const algoParam=NTuple{2,Int64}((
    floor(10^parsedArgs["Sweeps"])*(parsedArgs["Nlatt"]^parsedArgs["Dim"]),
    parsedArgs["Systems"]
))

const metaParam=(
    parsedArgs["Nlatt"],
    parsedArgs["Dim"],
    parsedArgs["Network"],
    parsedArgs["EnerFunc"],
    parsedArgs["Model"]
)


println()
println("Importing libraries")
println()



println()
println("Initializing Lattice")
println()

# get the function of the model
const sysFunc = getfield(Main,Symbol(parsedArgs["Model"])) 

if isfile(parsedArgs["Network"])
    adjMat = DelimitedFiles.readdlm(parsedArgs["Network"],',',Int64)
    initSys = [sysFunc(adjMat) for i in 1:cores]
else
    neigFunc = getfield(Lattices,Symbol(string(parsedArgs["Network"],"LatticeNeighbors")))
    initSys = [sysFunc(neigFunc,parsedArgs["Nlatt"],parsedArgs["Dim"]) for i in 1:cores]
end


let z  = getfield(StatEnsemble,Symbol( parsedArgs["EnerFunc"]))
    @sync @everywhere const enerFunc = $z
end

println()
println("Initializing parameters")
println()
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

@everywhere InOut.MakeAndEnterDirectories()
InOut.WriteAlgoParamTable(algoParam,"metropolis")
InOut.WriteMetaParamTable(metaParam)
InOut.WriteAdjMat(initSys[1])
InOut.WriteSimulParamDict(simulParamDict)

println()
println("Running simulations")
println()


@time @sync @distributed for i in 1:length(params2)
    iter = Int(params2[i][end])
    simulParam = params2[i][1:4]
    nwork = myid()-1
    if nwork==1
        per = round((i-1)*100*cores/length(params2); digits= 2)
        prog = "Progress: $(per) % "
        #current="Current: B, J, C, T = $(simulParam) "
        println()
        println(prog)
        #println(current)
        println()
    end
    name = string(simulParamDict[simulParam],"_",iter)
    res = Algorithms.MetropolisFast(initSys[iter],enerFunc,simulParam,algoParam)
    InOut.MetropolisAllOut(initSys[iter],res,name)
    initSys[iter] = copy(res[2])
    if i==length(params2)
        println()
        println("Progress: 100 %")
        println()
    end
end

# Legacy code

#=
@time @sync @distributed  for i in 1:length(params)
    simulParam = params[i][1:4]
    nwork = myid()-1
    if nwork==1
        per = round((i-1)*100*cores/length(params); digits= 2)
        prog = "Progress: $(per) % "
        #current="Current: B, J, C, T = $(simulParam) "
        println()
        println(prog)
        #println(current)
        println()
    end
    InOut.MakeAndEnterDirectories()
    InOut.WriteAlgoParamTable(algoParam,"metropolis")
    InOut.WriteMetaParamTable(metaParam)
    InOut.WriteAdjMat(initSys[1])
    res = Algorithms.MetropolisFast(initSys[nwork],enerFunc,simulParam,algoParam)
    name = string(simulParam,"_",string(nwork))
    mkdir(name)
    cd(name)
    InOut.MetropolisAllOut(initSys[nwork],res[1],algoParam)
    InOut.WriteSimulParamTable(simulParam)
    #initSys[nwork] = sysFunc(neigFunc,metaParam[1],metaParam[2])  
    initSys[nwork] = copy(res[2])  
    cd("..")
    InOut.ExitDirectories()
end
=#


