println()
println("Metropolis Data analyzer")
println()

println()
println("Importing libraries")
println()

include("statEnsemble.jl")
include("auxiliar.jl")
include("inOut.jl")
include("geometry.jl")

using Plots, Rsvg, DelimitedFiles, ArgParse, Statistics, LaTeXStrings
pyplot()

original=pwd()

println()
println("Parsing arguments")
println()

function ParseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-D", "--Dirname" 
            help = "Directory path to analyse set"
            arg_type = String
            required=true
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
        @time data = InOut.ReadSingleSimul(simulParam)
        ncut=Int64(floor(length(data[1])*cut))
        X=data[1][ncut:end]
        simulParam = InOut.ParseArray(string(split(simulParam,"_")[1]))
        push!(M,StatEnsemble.Magnetization(X,data[2]))
        push!(E,StatEnsemble.InternalEnergy(X,model,data[2],data[3]))
        push!(CV,StatEnsemble.HeatCapacity(X,model,data[2],data[3]))
        push!(Xi,StatEnsemble.MagneticSucep(X,data[2],data[3])) 
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
    for sim in sims
        indexes = [ i for i in 1:length(vals[1]) if string(split(vals[1][i],"_")[1]) == sim ]
        push!(M,Statistics.mean(vals[2][indexes]))
        push!(E,Statistics.mean(vals[3][indexes]))
        push!(CV,Statistics.mean(vals[4][indexes]))
        push!(Xi,Statistics.mean(vals[5][indexes]))
    end
    return ( sims , M , E , CV , Xi )
end

function MacroscopicTables(name,vals)
    original = pwd()
    cd(name)
    open("Macroscopic.csv","w") do f
        write(f,"Mu,J,C,kT,M,E,CV,Xi\n")
        for i in 1:length(vals[1])
            simulParam = InOut.ParseArray(vals[1][i])
            write(f," $(simulParam[1]) , $(simulParam[2]), $(simulParam[3]), $(simulParam[4]) , $(vals[2][i]), $(vals[3][i]), $(vals[4][i]), $(vals[5][i]) \n")
        end
    end
    cd(original)
end



function ThermalizationPlots(X,model,geoParam,simulParam,algoParam,cut=1/2)
    if model=="normal"
        enerFunc=StatEnsemble.NormalEnergy
        n1 = L"Magnetizaci\acute{o}n"
        n2 = L"Energ\acute{i}a"
        s1 = "Magnetizacion"
        s2 = "Energia"
    else
        enerFunc=StatEnsemble.PenalizedEnergy2
        n1 = L"Densidad"
        n2 = L"Energ\acute{i}a"
        s1 = "Densidad"
        s2 = "Energia"
    end
    neigLatt=Geometry.BuildLattices(geoParam,model)[2]
    steps=[(i-1)*(algoParam[2]) for i in 1:length(X)]
    mag=[abs(sum(x)) for x in X]
    en=[enerFunc(x,neigLatt,simulParam) for x in X]
    #making first plot
    scatter(
        steps,mag,xlabel=L"Pasos \; de \; Monte \; Carlo",ylabel=n1,legend=:none
    )
    plot!(
        steps,mag,legend=:none
    )
    Plots.savefig(string("Termalization",s1,".png"))
    Plots.closeall()
    #making second plot
    scatter(
        steps,en,xlabel=L"Pasos \; de \; Monte \; Carlo",ylabel=n2,legend=:none
    )
    plot!(
        steps,en,legend=:none
    )
    Plots.savefig(string("Termalization",s2,".png"))
    Plots.closeall()

    ncut=Int64(floor(length(X)*cut))
    epsilon = 0.1

    Plots.histogram(mag,xlabel = L"Valor", ylabel = L"Probabilidad",title=LaTeXString(string(s1," antes del corte")),normed=true,label=L"Histograma",legend=:topleft)
    sig = Statistics.std(mag)
    mu = Statistics.mean(mag)
    xs = range(sign(minimum(mag)) * (minimum(mag)*(1-epsilon)),stop = sign(maximum(mag))*(maximum(mag)*(1+epsilon)),length=200)
    Plots.plot!(xs, [NormalDistribution(x,mu=mu,sigma=sig) for x in xs],label = LaTeXString("N ( $(round(mu,digits=2)) , $(round(sig,digits=2)) )"))
    Plots.savefig(string(s1,"Antes",".png"))
    Plots.closeall()

    
    mag2 = mag[ncut:end]
    epsilon = 0.05



    Plots.histogram(mag2,xlabel = L"Valor", ylabel = L"Probabilidad",title = LaTeXString(string(s1," después del corte")),label = L"Histograma",normed=true,legend=:topleft)
    sig = Statistics.std(mag2)
    mu = Statistics.mean(mag2)
    xs = range(sign(minimum(mag2)) * (minimum(mag2)*(1-epsilon)),stop = sign(maximum(mag2))*(maximum(mag2)*(1+epsilon)),length=200)
    Plots.plot!(xs, [NormalDistribution(x,mu=mu,sigma=sig) for x in xs],label = LaTeXString("N ( $(round(mu,digits=2)) , $(round(sig,digits=2)) )"))
    Plots.savefig(string(s1,"Despues",".png"))
    Plots.closeall()

end



function ThermalizationAnalysis(name)
    original = pwd()
    cd(name)
    dirs = InOut.Folders()
    for folder in dirs
        data = InOut.ReadSingleSimul(folder)
        cd(folder)
        ThermalizationPlots(data[1],"cycle",data[2],data[3],data[4])
        cd("..")
    end
    cd(original)
end


#=
function MuCycPlots(name,cut= 1/2)
    #getting list of directories
    original=pwd()
    cd(name)
    dirs=InOut.Folders()
    sims=Set([split(d,"_")[1] for d in dirs])
    mus = []
    cs = []
    rhos = []
    for simulParam in sims
        @time X = InOut.ReadSingleMeanSimul(simulParam)
        ex = InOut.ParseArray(string(simulParam))
        ncut=Int64(floor(length(X)*cut))
        push!(mus,ex[1])
        push!(cs,ex[3])
        push!(rhos,sum(X[ncut:end]))
    end
    #saving the steps simulation as the last one 
    return mus, cs, rhos
end
=#
println()
println("Making Analysis")
println()

#ThermalizationAnalysis(parsedArgs["Dirname"])
#=
model = "cycle"
vals = MacroscopicVariables(parsedArgs["Dirname"],model)
vals = AverageMacroscopicVariables(vals)
MacroscopicTables(parsedArgs["Dirname"],vals)

=#

#=
mus=[]
rhos=[]
for iter in a 
    simulParam=InOut.ParseArray(iter[1])
    b=iter[2]
    ncut=Int64(length(b*3/4))
    X=[b[i] for i in ncut:length(b)]
    push!(mus,simulParam[1])
    push!(rhos,Statistics.mean([sum(x) for x in X]))
end

Plots.scatter(mus,rhos,
title="Fase Intermedia",
xlabel="Mu",
ylabel="Densidad",
label="simulación con N=100)",
legend=:topleft
)
Plots.savefig(joinpath(parsedArgs["Dirname"],"intermediate.png"))
=#


