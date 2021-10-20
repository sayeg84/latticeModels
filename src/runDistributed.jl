println("Ising model")
println()

using Distributed
using DelimitedFiles
include("init.jl")

println()
println("Parsing Arguments")
println()


const parsedArgs = setParseCommandline()

const cores = Sys.CPU_THREADS
Distributed.addprocs(parsedArgs["processes"])

@everywhere using Statistics, Dates
@sync @everywhere include("inOut.jl")
@everywhere include("algorithms.jl")
@everywhere include("lattices.jl")
@everywhere include("init.jl")







println()
println("Initializing Lattice")
println()

const metaParam,algoParam = getParams(parsedArgs)

# get the function of the model
const sysFunc = getfield(Main,Symbol(parsedArgs["Model"])) 

if isfile(parsedArgs["Network"])
    adjMat = DelimitedFiles.readdlm(parsedArgs["Network"],',',Int64)
    initSys = [sysFunc(adjMat) for i in 1:cores]
else
    initSys = [sysFunc(parsedArgs["Network"],metaParam[1],periodic = parsedArgs["Periodic"]) for i in 1:cores]
end


let z  = getfield(StatEnsemble,Symbol( parsedArgs["EnerFunc"]))
    @sync @everywhere const enerFunc = $z
end



println()
println("Initializing parameters")
println()


const simulParams, simulParamDict = getSimulParams(parsedArgs)

@everywhere InOut.MakeAndEnterDirectories()
InOut.WriteAlgoParamTable(algoParam,"metropolis")
InOut.WriteMetaParamTable(metaParam)
InOut.WriteAdjMat(initSys[1])
InOut.WriteSimulParamDict(simulParamDict)

println()
println("Running simulations")
println()


@time @sync @distributed for i in 1:length(simulParams)
    iter = Int(simulParams[i][end])
    simulParam = simulParams[i][1:4]
    nwork = myid()-1
    if nwork==1
        per = round((i-1)*100*cores/length(simulParams); digits= 2)
        prog = "Progress: $(per) % "
        #current="Current: B, J, C, T = $(simulParam) "
        println()
        println(prog)
        #println(current)
        println()
    end
    name = string(simulParamDict[simulParam],"_",iter)
    res = Algorithms.NPTFast(initSys[iter],enerFunc,simulParam,algoParam)
    InOut.MetropolisAllOut(initSys[iter],res,name)
    initSys[iter] = copy(res[2])
    if i==length(simulParams)
        println()
        println("Progress: 100 %")
        println()
    end
end

