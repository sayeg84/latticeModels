println()
println("Metropolis Dynamic Plots Maker")
println()

println()
println("Importing libraries")
println()

include("statEnsemble.jl")
include("auxiliar.jl")
include("inOut.jl")
include("geometry.jl")
include("cyclesAta.jl")
#=

Refactor all

using ArgParse, Statistics, PyPlot
PyPlot.rc("text", usetex=true)
PyPlot.rc("font", family="serif")
PyPlot.rc("text.latex", unicode=true)

function ParseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-D", "--Dirname" 
            help = "Directory path to make gif"
            arg_type = String
            required=true
        "-X"
            help = "Whether or not to show box containing system"
            action = :store_true
        "-L"
            help = "Whether or not to show box containing system"
            action = :store_true
    end
    return parse_args(s)
end
parsedArgs = ParseCommandline()
original = pwd()

function MakeFiguresSequence(X,metaParam,algoParam)
    neigLatt = Geometry.BuildLattices(metaParam,"cycle")[2]
    for i in 1:length(X)
        f1 = 16
        f2 = 12
        s = 20
        println("Making step $i plot")
        cyc = CyclesAta.ciclos2(X[i])
        if parsedArgs["L"]
            PyPlot.figure()
            PyPlot.scatter(size(X[1])[2]+3,size(X[1])[1]+3,c="k",label="Atomos")
            PyPlot.scatter(size(X[1])[2]+3,size(X[1])[1]+3,c="blueviolet",label="Componentes \n rigidas")
            PyPlot.legend(bbox_to_anchor=(1,1),fontsize=f1)
        else
            PyPlot.figure(figsize=(6,4))
        end
        PyPlot.title(string("Paso de Monte Carlo: ",(i-1)*algoParam[2]),fontsize=f1)
        PyPlot.xlim(0,size(X[1])[2]+1)
        PyPlot.ylim(0,size(X[1])[1]+1)
        for indx in CartesianIndices(size(X[1]))
            if X[i][indx] == 1
                PyPlot.scatter([indx[2]],[indx[1]],color="k")
                for vec in neigLatt[indx]
                    if X[i][vec] == 1 && maximum([abs((vec-indx)[i]) for i in 1:length(vec)]) <= 1
                        PyPlot.plot([indx[2],vec[2]],[indx[1],vec[1]],color="k")
                    end    
                end
            end
        end
        for indx in CartesianIndices(size(X[1]))
            if cyc[indx] == 1
                PyPlot.scatter([indx[2]],[indx[1]],color="blueviolet",zorder=5)
                for vec in neigLatt[indx]
                    if cyc[vec] == 1 && maximum([abs((vec-indx)[i]) for i in 1:length(vec)]) <= 1
                        PyPlot.plot([indx[2],vec[2]],[indx[1],vec[1]],color="blueviolet",zorder=5)
                    end    
                end
            end
        end
        if parsedArgs["X"]
            PyPlot.xticks([])
            PyPlot.yticks([])
        else
            PyPlot.axis("off")
        end
        step = lpad((i-1)*algoParam[2],Int(ceil(log10(algoParam[1])))+1,'0')
        PyPlot.tight_layout()
        PyPlot.savefig(string("plots/",step,".png"),dpi=300)
        PyPlot.close_figs()
    end
end


data = InOut.ReadSingleSimul(parsedArgs["Dirname"])
cd(parsedArgs["Dirname"])
if ~isdir("plots")
    mkdir("plots")
end
println()
println("Making figures")
println()
MakeFiguresSequence(data[1],data[2],data[4])
cd(original)
=#