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

    @add_arg_table! s begin
        "-D", "--path" 
            help = "Directory path to analyse set"
            arg_type = String
            required=true
        "-S", "--fromSystem" 
            help = "Whether to recalculate M and E from systems by rebuilding chain"
            arg_type = Bool
            default = false
        "-C", "--cut" 
            help = "Cut of chain for macroscopic variables"
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



function MacroscopicVariablesFromChain(initSys,changes,metaParam,simulParam,frequency,cut=1/2)
    sys = copy(initSys)
    M = Array{Float64,1}()
    E = Array{Float64,1}()
    enerFunc = getfield(StatEnsemble,Symbol(metaParam[4]))
    push!(M,StatEnsemble.Magnetization(sys))
    push!(E,enerFunc(sys,simulParam/N(initSys)))
    for i in 1:length(changes)
        if changes[i] != 0
            ChangeSpin!(sys,changes[i])
        end
        if mod(i,frequency) == 0
            push!(M,StatEnsemble.Magnetization(sys))
            push!(E,enerFunc(sys,simulParam)/N(initSys))
        end
    end
    ncut = Int64(floor(length(M)*cut))
    M = M[ncut:end]
    m = Statistics.mean(M)
    E = E[ncut:end]
    e = Statistics.mean(E)
    chi = Statistics.var(M.*N(initSys),corrected=false) / (simulParam[4]*(N(initSys)))
    cv = Statistics.var(E.*N(initSys),corrected=false) / ((simulParam[4]^2)*(N(initSys)))
    return m , e , cv , chi
end

function MacroscopicVariables(initSys,M,E,simulParam,frequency,cut=1/2)
    ncut = Int64(floor(length(M)*cut))
    M = M[ncut:frequency:end]./N(initSys)
    m = Statistics.mean(M)
    E = E[ncut:frequency:end]./N(initSys)
    e = Statistics.mean(E)
    chi = Statistics.var(M.*N(initSys),corrected=false) / (simulParam[4]*(N(initSys)))
    cv = Statistics.var(E.*N(initSys),corrected=false) / ((simulParam[4]^2)*(N(initSys)))
    return m , e , cv , chi
end



function GetMacroscopicVariables(; cut=1/2, fromSystem = false)
    M = Array{Float64,1}()
    E = Array{Float64,1}()
    CV = Array{Float64,1}()
    Chi = Array{Float64,1}()
    indexes = []
    #getting list of directories

    simulParamDict = InOut.ReadSimulParamDict()
    metaParam = InOut.ReadMetaParamTable()
    algoParam = InOut.ReadAlgoParamTable()
    adjMat = InOut.ReadAdjMat()     
    for i in 1:size(simulParamDict)[1]
        simulParamIndex = Int(simulParamDict[i,1])
        for j in 1:algoParam[2]
            if fromSystem
                # read calculating variables by rebuilding chain
                @time data = InOut.ReadSingleSimul(simulParamIndex,j,metaParam,adjMat) ; res = MacroscopicVariablesFromChain(data[1], data[2], metaParam, simulParamDict[i,2:end], parsedArgs["frequency"],cut)
            else
                # read using calculations already made
                @time data = InOut.ReadSingleSimul(simulParamIndex,j,metaParam,adjMat) ; res = MacroscopicVariables(data[1],data[3],data[4],simulParamDict[i,2:end],parsedArgs["frequency"],cut)
            end
            push!(M,res[1])
            push!(E,res[2])
            push!(CV,res[3])
            push!(Chi,res[4])
            push!(indexes,(simulParamIndex,j))
        end
    end
    return ( indexes , M , E , CV , Chi )
end

function AverageMacroscopicVariables(macros)
    simulParamIndexes = [x for x in Set([var[1] for var in macros[1]])]
    sort!(simulParamIndexes)
    M = []
    E = []
    CV = []
    Chi = []
    Msigma = []
    Esigma = []
    CVsigma = []
    Chisigma = []
    for i in 1:length(simulParamIndexes)
        indexes = [ j for j in 1:length(macros[1]) if macros[1][j][1] == simulParamIndexes[i]]
        push!(M,Statistics.mean(macros[2][indexes]))
        push!(E,Statistics.mean(macros[3][indexes]))
        push!(CV,Statistics.mean(macros[4][indexes]))
        push!(Chi,Statistics.mean(macros[5][indexes]))
        push!(Msigma,Statistics.std(macros[2][indexes]))
        push!(Esigma,Statistics.std(macros[3][indexes]))
        push!(CVsigma,Statistics.std(macros[4][indexes]))
        push!(Chisigma,Statistics.std(macros[5][indexes]))
    end
    return ( M , E , CV , Chi, Msigma , Esigma , CVsigma , Chisigma )
end

function MacroscopicTables(avgMacros)
    simulParamDict = InOut.ReadSimulParamDict()
    open("macroscopic.csv","w") do f
        write(f,"Index,B,J,C,kT,M,E,CV,Chi,Msigma,Esigma,CVsigma,Chisigma\n")
        DelimitedFiles.writedlm(f,hcat(sortslices(simulParamDict,dims=1),avgMacros...),',')
    end
    cd(original)
end

println()
println("Making Analysis")
println("cut = $(parsedArgs["cut"])")
println()


original = pwd()
cd(parsedArgs["path"])
@time macros = GetMacroscopicVariables(cut = parsedArgs["cut"],fromSystem=parsedArgs["fromSystem"])
avgMacros = AverageMacroscopicVariables(macros)
MacroscopicTables(avgMacros)
cd(original)
