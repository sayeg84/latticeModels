println()
println("Metropolis Data analyzer")
println()

println()
println("Importing libraries")
println()

include("statEnsemble.jl")
include("inOut.jl")

using DelimitedFiles, ArgParse, Statistics

original=pwd()

println()
println("Parsing arguments")
println()

function ParseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-D", "--path" 
            help = "Directory path to analyse set"
            arg_type = String
            required=true
        "-N", "--ncut" 
            help = "Directory path to analyse set"
            arg_type = Float64
            default = 1/2
        "-F", "--frequency" 
            help = "Frecuency of elements loaded in chain"
            arg_type = Int64
            default = 100

    end
    return parse_args(s)
end
parsedArgs = ParseCommandline()

function NormalDistribution(x ; mu = 0 , sigma = 1)
    return 1/sqrt(2*pi*sigma^2) * exp(-(x-mu)^2/(2*sigma^2))
end

function MacroscopicVariables(name,model;cut=1/2)
    M=[]
    E=[]
    CV=[]
    Xi=[]
    #getting list of directories
    original=pwd()
    cd(name)
    dirs=InOut.Folders()
    for simulParam in dirs
        @time data = InOut.ReadSingleSimul(simulParam,parsedArgs["frequency"])
        ncut=Int64(floor(length(data[1])*cut))
        X=data[1][ncut:end]
        simulParam = InOut.ParseArray(string(split(simulParam,"_")[1]))
        push!(M,StatEnsemble.Magnetization(X))
        push!(E,StatEnsemble.InternalEnergy(X,StatEnsemble.PenalizedEnergy,data[3]))
        push!(CV,StatEnsemble.HeatCapacity(X,StatEnsemble.PenalizedEnergy,data[3]))
        push!(Xi,StatEnsemble.MagneticSucep(X,data[3])) 
    end
    cd(original)
    return ( dirs , M , E , CV , Xi )
end

function AverageMacroscopicVariables(vals)
    sims = [x for x in Set([string(split(d,"_")[1]) for d in vals[1]])]
    M = []
    E = []
    CV = []
    Xi = []
    Msigma = []
    Esigma = []
    CVsigma = []
    Xisigma = []
    for sim in sims
        indexes = [ i for i in 1:length(vals[1]) if string(split(vals[1][i],"_")[1]) == sim ]
        push!(M,Statistics.mean(vals[2][indexes]))
        push!(E,Statistics.mean(vals[3][indexes]))
        push!(CV,Statistics.mean(vals[4][indexes]))
        push!(Xi,Statistics.mean(vals[5][indexes]))
        push!(Msigma,Statistics.std(vals[2][indexes]))
        push!(Esigma,Statistics.std(vals[3][indexes]))
        push!(CVsigma,Statistics.std(vals[4][indexes]))
        push!(Xisigma,Statistics.std(vals[5][indexes]))
    end
    return ( sims , M , E , CV , Xi, Msigma , Esigma , CVsigma , Xisigma )
end

function MacroscopicTables(name,vals)
    original = pwd()
    cd(name)
    open("Macroscopic.csv","w") do f
        write(f,"Mu,J,C,kT,M,E,CV,Xi,Msigma,Esigma,CVsigma,Xisigma\n")
        for i in 1:length(vals[1])
            simulParam = InOut.ParseArray(vals[1][i])
            write(f,"$(simulParam[1]),$(simulParam[2]),$(simulParam[3]),$(simulParam[4]),$(vals[2][i]),$(vals[3][i]),$(vals[4][i]),$(vals[5][i]),$(vals[6][i]),$(vals[7][i]),$(vals[8][i]),$(vals[9][i])\n")
        end
    end
    cd(original)
end

println()
println("Making Analysis")
println("ncut = $(parsedArgs["ncut"])")
println()



model = "cycle"
@time vals = MacroscopicVariables(parsedArgs["path"],model,cut = parsedArgs["ncut"])
vals = AverageMacroscopicVariables(vals)
MacroscopicTables(parsedArgs["path"],vals)

