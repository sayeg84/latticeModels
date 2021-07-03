println("Ising model")
println()
println()
println("Importing libraries")
println()


include("initialization.jl")



println()
println("Parsing Arguments")
println()

const parsedArgs = setParseCommandline()
const metaParam,algoParam = getParams(parsedArgs)


println()
println("Importing libraries")
println()



println()
println("Initializing Lattice")
println()

# get the function of the model
const sysFunc = getfield(Main,Symbol(parsedArgs["Model"])) 
const sizes = StringToTuple(parsedArgs["Nlatt"])
const enerFunc = getfield(StatEnsemble,Symbol( metaParam[4]))


if isfile(parsedArgs["Network"])
    adjMat = DelimitedFiles.readdlm(parsedArgs["Network"],',',Int64)
    initSys = [sysFunc(adjMat) for i in 1:algoParam[2]]
else
    initSys = [sysFunc(parsedArgs["Network"],metaParam[1],periodic = parsedArgs["Periodic"]) for i in 1:algoParam[2]]
end

const simulParams, simulParamDict = getSimulParams(parsedArgs)

InOut.MakeAndEnterDirectories()
InOut.WriteAlgoParamTable(algoParam,"metropolis")
InOut.WriteMetaParamTable(metaParam)
InOut.WriteAdjMat(initSys[1])
InOut.WriteSimulParamDict(simulParamDict)

println()
println("Running simulations")
println()


for i in 1:length(simulParams)
    iter = Int(simulParams[i][end])
    simulParam = simulParams[i][1:4]
    name = string(simulParamDict[simulParam],"_",iter)
    print("Running ")
    println(name)
    println(simulParam)
    @time res = Algorithms.MetropolisFast(initSys[iter],enerFunc,simulParam,algoParam)
    InOut.MetropolisAllOut(initSys[iter],res,name)
    initSys[iter] = copy(res[2])
end

