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



function MacroscopicVariables(initSys,changes,metaParam,simulParam,frequency,cut=1/2)
    sys = copy(initSys)
    M = Array{Float64,1}()
    E = Array{Float64,1}()
    enerFunc = enerFunc = getfield(StatEnsemble,Symbol( metaParam[4]))
    push!(M,StatEnsemble.Magnetization(sys))
    push!(E,enerFunc(sys,simulParam))
    for i in 1:length(changes)
        if changes[i] != 0
            ChangeSpin!(sys,changes[i])
        end
        if mod(i,frequency) == 0
            push!(M,StatEnsemble.Magnetization(sys))
            push!(E,enerFunc(sys,simulParam))
        end
    end
    ncut = Int64(floor(length(M)*cut))
    M = M[ncut:end]
    m = Statistics.mean(M)
    E = E[ncut:end]
    e = Statistics.mean(E)
    cv = Statistics.var(M,corrected=false) / ((simulParam[4]^2)*(N(initSys)))
    xi = Statistics.var(E.*N(initSys),corrected=false) / (simulParam[4]*(N(initSys)))
    return m,e,cv,xi
end

function GetMacroscopicVariables(name;cut=1/2)
    M = Array{Float64,1}()
    E = Array{Float64,1}()
    CV = Array{Float64,1}()
    Xi = Array{Float64,1}()
    #getting list of directories
    original = pwd()
    cd(name)
    dirs = InOut.Folders()
    for folder in dirs
        @time data = InOut.ReadSingleSimul(folder)
        res = MacroscopicVariables(data[1], data[2], data[3], data[4], parsedArgs["frequency"],cut)
        push!(M,res[1])
        push!(E,res[2])
        push!(CV,res[3])
        push!(Xi,res[4])
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



@time vals = GetMacroscopicVariables(parsedArgs["path"],cut = parsedArgs["ncut"])
vals = AverageMacroscopicVariables(vals)
MacroscopicTables(parsedArgs["path"],vals)

