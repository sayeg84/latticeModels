println()
println("Metropolis Dynamic Plots Maker")
println()

println()
println("Importing libraries")
println()

include("statEnsemble.jl")
include("inOut.jl")
include("lattices.jl")

# importing function from lattice module




using ArgParse, Statistics, PyPlot
PyPlot.rc("text", usetex=true)
PyPlot.rc("font", family="serif")
PyPlot.rc("text.latex", unicode=true)

function ParseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-D", "--path" 
            help = "Directory path to make gif"
            arg_type = String
            required=true
        "-F","--freq"
            help = "Frequency of plots"
            arg_type = Int64
            default = 100 
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

dirs = InOut.Folders()
splitedPath = splitpath(parsedArgs["Path"])
afterpath = splitedPath[end]
path = splitedPath[1]
for i in 2:(length(splitedPath)-1)
    global path
    path = joinpath(path,splitedPath[i])
end
cd(path)
simul = InOut.ReadSingleSimul(afterpath)
sys = simul[1]
changes = simul[2]
sysArray = Array{typeof(sys),1}([copy(sys)])
for i in 1:length(changes)
    if changes[i] != 0
        ChangeSpin!(sys,changes[i])
    end
    if mod(i,parsedArgs["Freq"]) == 0
        push!(sysArray,copy(sys))
    end
end
#cycArray = [Cycles.CycSubSys(x) for x in sysArray]


function plotSystem(sys,c;zpos=0)
    sizes = fill(sys.shape[1],sys.shape[2])
    sizes = tuple(sizes...)
    positions = Lattices.SquarePositionLattice(sizes)
    positions = reshape(positions,length(positions))
    ones = [i for i in 1:N(sys) if sys.sites[i]==1]
    xpositions = [positions[i][1] for i in ones ]
    ypositions = [positions[i][2] for i in ones ]
    PyPlot.scatter(xpositions,ypositions,s=140,color=c,zorder=zpos)
    visited = falses(N(sys))
    for i in 1:N(sys)
        if sys.sites[i] == 1
            for j in sys.edgList[i]
                if sys.sites[j]==1
                    PyPlot.plot([positions[i][1],positions[j][1]],[positions[i][2],positions[j][2]],color=c,zorder=zpos-1)
                end
            end
        end
    end
end
cd(afterpath)
if !(isdir("animPlots"))
    mkdir("animPlots")
end
cd("animPlots")
for i in 1:length(sysArray)
    PyPlot.figure(figsize=(8,8))
    plotSystem(sysArray[i],"k")
    plotSystem(Cycles.CycSubSys(sysArray[i]),(0.6,0.4,0.8))
    if parsedArgs["X"]
        PyPlot.xticks([])
        PyPlot.yticks([])
    else
        PyPlot.axis("off")
    end
    PyPlot.tight_layout()
    PyPlot.savefig("$i.pdf")    
end
cd("..")
cd("..")



#=
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


data = InOut.ReadSingleSimul(parsedArgs["path"])
cd(parsedArgs["path"])
if ~isdir("plots")
    mkdir("plots")
end
println()
println("Making figures")
println()
MakeFiguresSequence(data[1],data[2],data[4])
cd(original)
=#