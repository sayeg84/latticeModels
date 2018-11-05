println()
println("Data analyzer")
println()

println()
println("Importing libraries")
println()

include("statEnsemble.jl")
include("auxiliar.jl")
include("inOut.jl")

using Plots, Rsvg, DelimitedFiles, ArgParse, Statistics
gr()

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

dic=InOut.MetropolisIn(parsedArgs["Dirname"])

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

function PlotMag(vals)
    
end

mus=[]
rhos=[]
for iter in a 
    simulParam=iter[1]
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
label="simulación con N=400)"
)
Plots.savefig("finally.png")

#=
for i in -3.0:0.05:-2.0
    println(i)
    indx=[i,2.0,0.8,0.5]
    b=a[indx]
    simulParam=indx[1]
    dens=[sum(x) for x in b]
    Plots.scatter(dens,title="termalizacion",ylim=(180,220),xlabel="pasos",ylabel="densidad")
    Plots.plot!(dens)
    Plots.savefig("termalizacion$(Int64(round(i*100))).png")
end
#=
using InteractiveUtils
println(varinfo(r"a"))
=#