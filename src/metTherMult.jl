println()
println("Metropolis Thermalization Multi Data analyzer")
println()

println()
println("Importing libraries")
println()

include("statEnsemble.jl")
include("inOut.jl")

using DelimitedFiles, ArgParse, Statistics, PyPlot
##remove references to algoParam[2]

function NormalDistribution(x ; mu = 0 , sigma = 1)
    return 1/sqrt(2*pi*sigma^2) * exp(-(x-mu)^2/(2*sigma^2))
end

function polyfit(xs::Vector, ys::Vector, deg::Int) 
    return collect([v ^ p for v in xs, p in 0:deg]) \ ys
end

#PyCall.PyDict(matplotlib["rcParams"])["text.latex.unicode"] = "True"
PyPlot.rc("text", usetex=true)
PyPlot.rc("font", family="serif")
PyPlot.rc("text.latex", unicode=true)
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
        "-N", "--Ncut" 
            help = "Fraction to make cut"
            arg_type = Float64
            default = 1/2
         "-P" 
            help = "Plot or not normal distribution adjusting histograms"
            action = :store_true
        "-C" 
            help = "Center histograms in mean or not"
            action = :store_true
        "-S" 
            help = "Put scatter in thermalization plots"
            action = :store_true
    end
    return parse_args(s)
end
parsedArgs = ParseCommandline()
original = pwd()


function MacroscopicVariables(initSys,changes,metaParam,simulParam,frequency)
    sys = copy(initSys)
    M = Array{Float64,1}()
    E = Array{Float64,1}()
    enerFunc = enerFunc = getfield(StatEnsemble,Symbol( metaParam[4]))
    push!(M,StatEnsemble.Magnetization(sys))
    push!(E,enerFunc(sys,simulParam)/N(initSys))
    for i in 1:length(changes)
        if changes[i] != 0
            ChangeSpin!(sys,changes[i])
        end
        if mod(i,frequency) == 0
            push!(M,StatEnsemble.Magnetization(sys))
            push!(E,enerFunc(sys,simulParam)/N(initSys))
        end
    end
    return M, E
end

dirs = InOut.Folders()
splitedPath = splitpath(parsedArgs["Dirname"])
afterpath = splitedPath[end]
path = splitedPath[1]
for i in 2:(length(splitedPath)-1)
    global path
    path = joinpath(path,splitedPath[i])
end
cd(path)
dirs = InOut.Folders()
simils = [dir for dir in dirs if split(dir,"_")[1] == split(afterpath,"_")[1]]
Ms = []
Es = []
for path in simils
    global data = InOut.ReadSingleSimul(path)
    M,E = MacroscopicVariables(data[1],data[2],data[3],data[4],1)
    push!(Ms,copy(M))
    push!(Es,copy(E))
end
cd(original)
open("$(split(afterpath,"_")[1])_chains.csv","w") do io
    for i in 1:(length(simils)-1)
        write(io,string("M",i,","))
        write(io,string("E",i,","))
    end
    write(io,string("M",length(simils),","))
    write(io,string("E",length(simils),"\n"))
    for i in 1:length(Ms[1])
        for j in 1:(length(simils)-1)
            write(io,"$(Ms[j][i]),")
            write(io,"$(Es[j][i]),")
        end
        write(io,"$(Ms[end][i]),")
        write(io,"$(Es[end][i])\n")
    end
end

#T, XM, XE = StatEnsemble.AutoCorrelations(M,E,1/3,1)
#=

if data[3][5] == "SpinLattice"
    n1 = "\$ \\mathcal{M} / N \$"
    n1c = "\$ Autocorr \\left(  \\mathcal{M} / N  \\right) \$"
    n2 = "\$ U / N \$"
    n2c = "\$ Autocorr \\left(  U / N  \\right) \$"
    s1 = "Magnetizacion"
    s2 = "Energia"
elseif data[3][5] == "LatticeGas"
    n1 = "\$ \\rho  \$"
    n1c = "\$ Autocorr \\left(  \\rho  \\right) \$"
    n2 = "\$ U / N \$"
    n2c = "\$ Autocorr \\left(  U / N  \\right) \$"
    s1 = "Densidad"
    s2 = "Energia"
end

f1=22
f2=26
f3=30
siz=20

steps = collect(1:length(Ms[1]))
PyPlot.figure(figsize=(10,6.66))
for i in 1:length(simils)
    PyPlot.plot(steps,Ms[i],lw= 3,alpha = 0.7)
    if parsedArgs["S"]
        PyPlot.scatter(steps,Ms[i],s=siz,alpha = 0.5)
    end
end
PyPlot.title("Termalización",fontsize=f3)
PyPlot.xlabel("Pasos de Monte Carlo",fontsize=f2)
PyPlot.ylabel(n1,fontsize=f2)
PyPlot.grid()
PyPlot.xticks(fontsize=f1)
PyPlot.yticks(fontsize=f1)
PyPlot.xlim(0.0,2000)
#PyPlot.ylim(-0.1,1.1)
PyPlot.tight_layout()
PyPlot.show()

corrTimeM = 1
PyPlot.figure(figsize=(6,4))

for i in 1:length(simils)
    if parsedArgs["C"]
        M1 = [m - Statistics.mean(Ms[i][1:corrTimeM:end]) for m in Ms[i][1:corrTimeM:end]]
        PyPlot.xlabel(string(n1,"\$ - \\langle \$",n1,"\$ - \\rangle \$"),fontsize=f2)
    else
        M1 = Ms[i][1:corrTimeM:end]
        PyPlot.xlabel(n1,fontsize=f2)
    end
    lim1 = minimum(M1)
    lim2 = maximum(M1)
    sigma = Statistics.std(M1)
    mu = Statistics.mean(M1)
    if parsedArgs["P"]
        xs = range(mu-3*sigma,stop=mu+3*sigma,length=200)
        PyPlot.hist(M1,density=true,bins=20,alpha = 0.5)
        PyPlot.plot(xs, [NormalDistribution(x,mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=4)) , $(round(sigma,digits=4)) )",color="C1")
        PyPlot.ylabel("Probabilidad",fontsize=f2)
    else
        PyPlot.hist(M1,label="\$ \\mu_{X} = $(round(mu,digits=3)) \$ \n \$ \\sigma_{X} = $(round(sigma,digits=3)) \$",bins=20,alpha = 0.5)
        PyPlot.ylabel("Frecuencia",fontsize=f2)
    end
end
PyPlot.title(string(s1," antes del corte"),fontsize=f3)
#PyPlot.xlim(lim1,lim2)
PyPlot.ylabel("Probabilidad",fontsize=f2)
PyPlot.xticks(fontsize=f1)
PyPlot.yticks(fontsize=f1)
PyPlot.grid()
PyPlot.legend(fontsize=f2,loc=0)
PyPlot.tight_layout()
#PyPlot.savefig(string(s1,"Antes",".png"),dpi=300)
PyPlot.show()
#PyPlot.close_figs()








ncut = Int64(floor(length(Ms[1])*parsedArgs["Ncut"]))
PyPlot.figure(figsize=(6,4))
PyPlot.title(string(s1," despues del corte"),fontsize=f3)
for i in 1:length(simils)
    M2 = Ms[i][ncut:corrTimeM:end]
    if parsedArgs["C"]
        M2 = [m - Statistics.mean(M2) for m in M2]
        PyPlot.xlabel(string(n1,"\$ - \\langle \$",n1,"\$ - \\rangle \$"),fontsize=f2)
    else
        PyPlot.xlabel(n1,fontsize=f2)
    end 
    sigma = Statistics.std(M2)
    mu = Statistics.mean(M2)
    if parsedArgs["P"]
        xs = range(mu-3*sigma,stop = mu+3*sigma,length=200)
        PyPlot.hist(M2,density=true,bins=20)
        PyPlot.plot(xs, [NormalDistribution(x,mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=4)) , $(round(sigma,digits=4)) )",color="C1")
        PyPlot.ylabel("Probabilidad",fontsize=f2)
    else
        PyPlot.hist(M2,label="\$ \\mu_{X} = $(round(mu,digits=3)) \$ \n \$ \\sigma_{X} = $(round(sigma,digits=3)) \$",bins=20)
        PyPlot.ylabel("Frecuencia",fontsize=f2)
    end
end
#PyPlot.xlim(lim1,lim2)
PyPlot.xticks(fontsize=f1)
PyPlot.yticks(fontsize=f1)
PyPlot.grid()
PyPlot.legend(fontsize=f2,loc=0)
PyPlot.tight_layout()
#PyPlot.savefig(string(s1,"Despues",".png"),dpi=300)
PyPlot.show()
#PyPlot.close_figs()

=#