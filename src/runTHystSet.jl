println("Ising model")
println()
println()
println("Importing libraries")
println()

using  Distributed
@everywhere using ArgParse, Statistics, Dates

cores = Sys.CPU_THREADS-1
Distributed.addprocs(cores)


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
            default = 3.0
        "-A", "--Systems"
            help = "Number of independent systems simulated for same simulation parameters"
            arg_type = Int64
            default = cores
        "-B","--bconst"
            help = "value of B constant"
            arg_type = Float64
            default = -2.0
        "-J","--jconst"
            help = "value of J constant"
            arg_type = Float64
            default = 2.0
        "-C","--cconst"
            help = "value of C constant"
            arg_type = Float64
            default = 0.9
        "-T","--temp"
            help = "value of center temperature"
            arg_type = Float64
            default = 0.5
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

#initializing parameters

bArray = [parsedArgs["bconst"]]
jArray = [parsedArgs["jconst"]]
cArray = [parsedArgs["cconst"]]

tcenter = parsedArgs["temp"]
tArray = range(tcenter/2, stop = tcenter*3/2, length = 31)

println()
println("Initializing Lattice")
println()

neigFunc = getfield(Lattices,Symbol(metaParam[3]))
sysFunc = getfield(Main,Symbol(metaParam[5])) 
initSys = [sysFunc(neigFunc,metaParam[1],metaParam[2]) for i in 1:cores]


let z  = getfield(StatEnsemble,Symbol( metaParam[4]))
    @sync @everywhere enerFunc = $z
end


params=[Array{Float64,1}([bArray[i1],jArray[i2],cArray[i3],tArray[end-i4],iter]) for iter in 1:algoParam[2] for i3 in 1:length(cArray) for i1 in 1:length(bArray) for i4 in 0:(length(tArray)-1) for i2 in 1:length(jArray)]

println()
println("Running exotermic simulations")
println()

@time @sync @distributed  for i in 1:length(params)
    simulParam = params[i][1:4]
    nwork = myid()-1
    if nwork==1
        per = round((i-1)*100*algoParam[2]/length(params); digits= 2)
        prog = "Progress: $(per) % "
        #current="Current: B, J, C, T = $(simulParam) "
        println()
        println(prog)
        #println(current)
        println()
    end
    InOut.MakeAndEnterDirectories()
    if ~(isdir("exothermic"))
        try mkdir("exothermic") catch SystemError end
    end
    cd("exothermic")
    InOut.WriteAlgoParamTable(algoParam,"metropolis")
    InOut.WriteMetaParamTable(metaParam)
    InOut.WriteAdjMat(initSys[1])
    res = Algorithms.MetropolisFast(initSys[nwork],enerFunc,simulParam,algoParam)
    name = string(simulParam,"_",string(nwork))
    mkdir(name)
    cd(name)
    InOut.MetropolisAllOut(initSys[nwork],res[1],algoParam)
    InOut.WriteSimulParamTable(simulParam)
    initSys[nwork] = copy(res[2])
    cd("..")
    cd("..")
    InOut.ExitDirectories()
end


params=[Array{Float64,1}([bArray[i1],jArray[i2],cArray[i3],tArray[i4],iter]) for iter in 1:algoParam[2] for i3 in 1:length(cArray) for i1 in 1:length(bArray) for i4 in 1:(length(tArray)) for i2 in 1:length(jArray)]

println()
println("Running endothermic simulations")
println()

@time @sync @distributed  for i in 1:length(params)
    simulParam = params[i][1:4]
    nwork = myid()-1
    if nwork==1
        per = round((i-1)*100*algoParam[2]/length(params); digits= 2)
        prog = "Progress: $(per) % "
        #current="Current: B, J, C, T = $(simulParam) "
        println()
        println(prog)
        #println(current)
        println()
    end
    InOut.MakeAndEnterDirectories()
    if ~(isdir("endothermic"))
        try mkdir("endothermic") catch SystemError end
    end
    cd("endothermic")
    InOut.WriteAlgoParamTable(algoParam,"metropolis")
    InOut.WriteMetaParamTable(metaParam)
    InOut.WriteAdjMat(initSys[1])
    res = Algorithms.MetropolisFast(initSys[nwork],enerFunc,simulParam,algoParam)
    name = string(simulParam,"_",string(nwork))
    mkdir(name)
    cd(name)
    InOut.MetropolisAllOut(initSys[nwork],res[1],algoParam)
    InOut.WriteSimulParamTable(simulParam)
    initSys[nwork] = copy(res[2])
    cd("..")
    cd("..")
    InOut.ExitDirectories()
end