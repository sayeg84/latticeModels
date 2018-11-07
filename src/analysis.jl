println()
println("Data analyzer")
println()

println()
println("Importing libraries")
println()

include("statEnsemble.jl")
include("auxiliar.jl")
include("inOut.jl")
include("geometry.jl")

using Plots, Rsvg, DelimitedFiles, ArgParse, Statistics
gr()

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



function MacroscopicVariables(dic,model,geoParam;cut=1/2)
    indexes=[]
    M=[]
    E=[]
    CV=[]
    Xi=[]
    for iter in dic
        simulParam=iter[1]
        ncut=Int64(floor(length(iter[2])*cut))
        X=[iter[2][i] for i in ncut:length(iter[2])]
        push!(indexes,simulParam)
        push!(M,StatEnsemble.Magnetization(X,geoParam))
        push!(E,StatEnsemble.InternalEnergy(X,model,geoParam,simulParam))
        push!(CV,StatEnsemble.HeatCapacity(X,model,geoParam,simulParam))
        push!(Xi,StatEnsemble.MagneticSucep(X,geoParam,simulParam)) 
    end
    return (indexes,M,E,CV,Xi)
end

function ThermalizationAnalysis(X,model,geoParam,simulParam,algoParam)
    
    if model=="normal"
        enerFunc=StatEnsemble.NormalEnergy
        n1="Magnetizacion"
        n2="Energia"
    else
        enerFunc=StatEnsemble.PenalizedEnergy2
        n1="Densidad"
        n2="Energia"
    end
    neigLatt=Geometry.BuildLattices(geoParam,model)[2]

    steps=[(i-1)*(algoParam[2]) for i in 1:length(X)]
    mag=[abs(sum(x)) for x in X]
    en=[enerFunc(x,neigLatt,simulParam) for x in X]

    #making first plot
    scatter(
        steps,mag,xlabel="Pasos de Monte Carlo",ylabel=n1,title="Termalizacion"
    )
    plot!(
        steps,mag
    )

    Plots.savefig("Termalization1.png")
    #making second plot
    scatter(
        steps,en,xlabel="Pasos de Monte Carlo",ylabel=n2,title="Termalizacion"
    )
    plot!(
        steps,en
    )
    Plots.savefig("Termalization2.png")
end

function PlotMag(vals)
    M=[]
    for element in vals[1]
    end
end

a=InOut.MetropolisIn(parsedArgs["Dirname"])
#=
cd(parsedArgs["Dirname"])
for iter in a
    cd(iter[1])
    simulParam=InOut.ParseArray(iter[1])
    geoParam=InOut.ReadGeoParamTable()
    algoParam=InOut.ReadAlgoParamTable()
    ThermalizationAnalysis(iter[2],"cycle",geoParam,simulParam,algoParam)
    cd("..")
end
cd(original)
=#
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
