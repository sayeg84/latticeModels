println("Ising model")
println()

using  Distributed
cores = Sys.CPU_THREADS-1
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
        "-W", "--BDirection" 
            help = "Whether to go increasing or decreasing on B param"
            arg_type = String
            default = "increasing"
        "-S", "--sweeps"
            help = "logarithm base 10 of total sweeps"
            arg_type = Float64
            default = 3.0
        "-A", "--Systems"
            help = "Number of independent systems simulated for same simulation parameters"
            arg_type = Int64
            default = cores
    end
    return parse_args(s)
end

parsedArgs = ParseCommandline()

algoParam=Array{Int64,1}([
    floor(10^parsedArgs["sweeps"])*parsedArgs["Nlatt"]^parsedArgs["dim"],
    parsedArgs["Systems"]
])

metaParam=[
    parsedArgs["Nlatt"],
    parsedArgs["dim"],
    string(parsedArgs["Geometry"],"LatticeNeighbors"),
    parsedArgs["EnerFunc"],
    parsedArgs["Model"]
]

println()
println("Importing libraries")
println()



bArray = range(-3.5,stop = 0.0,length = 41)
jArray = [4.0]
cArray = [0.9]
#cArray = [0.0]
tArray = range(0.2,stop=1.6,length=31)
#tArray = range(0.1,10.0,length = 41)

println()
println("Initializing Lattice")
println()


neigFunc = getfield(Lattices,Symbol(metaParam[3]))
sysFunc = getfield(Main,Symbol(metaParam[5])) 
initSys = [sysFunc(neigFunc,metaParam[1],metaParam[2]) for i in 1:cores]


let z  = getfield(StatEnsemble,Symbol( metaParam[4]))
    @sync @everywhere enerFunc = $z
end



if parsedArgs["BDirection"]=="increasing"
    params1 = [Array{Float64,1}([bArray[i1],jArray[i2],cArray[i3],tArray[end-i4]])  for i3 in 1:length(cArray) for i1 in 1:length(bArray) for i4 in 0:(length(tArray)-1) for i2 in 1:length(jArray)]
else 
    params1 = [Array{Float64,1}([bArray[end-i1],jArray[i2],cArray[i3],tArray[end-i4]])  for i3 in 1:length(cArray) for i1 in 1:length(bArray) for i4 in 0:(length(tArray)-1) for i2 in 1:length(jArray)]
end
params2 = [vcat(par,[iter]) for iter in 1:algoParam[2] for par in params1]
paramDict = Dict(params1[i] => i for i in 1:length(params1))


@everywhere InOut.MakeAndEnterDirectories()
InOut.WriteAlgoParamTable(algoParam,"metropolis")
InOut.WriteMetaParamTable(metaParam)
InOut.WriteAdjMat(initSys[1])
InOut.WriteSimulParamDict(paramDict)

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
    name = string(paramDict[simulParam],"_",iter)
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


