println()
println("Metropolis Thermalization Single Data analyzer")
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
        "-D", "--path" 
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


splitedPath = splitpath(parsedArgs["path"])
afterpath = splitedPath[end]
path=splitedPath[1]
for i in 2:(length(splitedPath)-1)
    global path
    path = joinpath(path,splitedPath[i])
end
cd(path)
data = InOut.ReadSingleSimul(afterpath)
M,E = MacroscopicVariables(data[1],data[2],data[3],data[4],1)


T, XM, XE = StatEnsemble.AutoCorrelations(M,E,1/3,1)


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


steps = collect(1:length(M))
PyPlot.figure(figsize=(10,6.66))
PyPlot.plot(steps,M,lw= 3)
if parsedArgs["S"]
    PyPlot.scatter(steps,M,s=siz)
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

PyPlot.figure(figsize=(10,6.66))

if parsedArgs["S"]
    PyPlot.scatter(T,XM,s=siz)
end
PyPlot.plot(T,XM,lw=3)
b,m = polyfit(T,log.(XM),1)
corrTimeM = Int64(floor(-1/m))
#PyPlot.semilogy(T,[exp(m*t+b) for t in T],label="Ajuste lineal con T = $(-1/m)")
PyPlot.title("Termalización",fontsize=f3)
PyPlot.xlabel("Pasos de Monte Carlo",fontsize=f2)
PyPlot.ylabel(n1c,fontsize=f2)
PyPlot.grid()
PyPlot.xticks(fontsize=f1)
PyPlot.yticks(fontsize=f1)
PyPlot.legend()
PyPlot.tight_layout()
PyPlot.show()

PyPlot.figure(figsize=(10,6.66))
PyPlot.plot(steps,E,lw= 3,color="C2")
if parsedArgs["S"]
    PyPlot.scatter(steps,E,s=siz,color="C2")
end
PyPlot.title("Termalización",fontsize=f3)
PyPlot.xlabel("Pasos de Monte Carlo",fontsize=f2)
PyPlot.ylabel(n2,fontsize=f2)
PyPlot.grid()
PyPlot.xticks(fontsize=f1)
PyPlot.yticks(fontsize=f1)
PyPlot.xlim(0.0,2000)
#PyPlot.ylim(-4.1,1.1)
PyPlot.tight_layout()
PyPlot.savefig(string("Termalization",s2,".png"),dpi=300)
PyPlot.show()


PyPlot.figure(figsize=(10,6.66))
if parsedArgs["S"]
    PyPlot.scatter(T,XE,s=siz,color="C2")
end
PyPlot.plot(T,XE,color="C2",lw=3)
b,m = polyfit(T,log.(XE),1)
corrTimeE = Int64(floor(-1/m))
#PyPlot.semilogy(T,[exp(m*t+b) for t in T],label="Ajuste lineal con T = $(-1/m)")
PyPlot.title("Termalización",fontsize=f3)
PyPlot.xlabel("Pasos de Monte Carlo",fontsize=f2)
PyPlot.ylabel(n2c,fontsize=f2)
PyPlot.grid()
PyPlot.xticks(fontsize=f1)
PyPlot.yticks(fontsize=f1)
PyPlot.legend()
PyPlot.tight_layout()
PyPlot.show()

print("Correlation time for sum: ")
println(corrTimeM)
print("Correlation time for energy: ")
println(corrTimeE)
if corrTimeM> 1000
    corrTimeM = 10
end
if corrTimeE> 1000
    corrTimeE = 10
end

    
PyPlot.title(string(s1," antes del corte"),fontsize=f3)
if parsedArgs["C"]
    M1 = [m - Statistics.mean(M) for m in M[1:corrTimeM:end]]
    PyPlot.xlabel(string(n1,"\$ - \\langle \$",n1,"\$ - \\rangle \$"),fontsize=f2)
else
    M1 = M[1:corrTimeM:end]
    PyPlot.xlabel(n1,fontsize=f2)
end
lim1 = minimum(M1)
lim2 = maximum(M1)
PyPlot.xlim(lim1,lim2)
PyPlot.ylabel("Probabilidad",fontsize=f2)

PyPlot.xticks(fontsize=f1)
PyPlot.yticks(fontsize=f1)
sigma = Statistics.std(M1)
mu = Statistics.mean(M1)
if parsedArgs["P"]
    xs = range(mu-3*sigma,stop=mu+3*sigma,length=200)
    PyPlot.hist(M1,density=true,color="C0",bins=20)
    PyPlot.plot(xs, [NormalDistribution(x,mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=4)) , $(round(sigma,digits=4)) )",color="C1")
    PyPlot.ylabel("Probabilidad",fontsize=f2)
else
    PyPlot.hist(M1,label="\$ \\mu_{X} = $(round(mu,digits=3)) \$ \n \$ \\sigma_{X} = $(round(sigma,digits=3)) \$",color="C0",bins=20)
    PyPlot.ylabel("Frecuencia",fontsize=f2)
end
PyPlot.grid()
PyPlot.legend(fontsize=f2,loc=0)
PyPlot.tight_layout()
#PyPlot.savefig(string(s1,"Antes",".png"),dpi=300)
PyPlot.show()
#PyPlot.close_figs()


ncut = Int64(floor(length(M)*parsedArgs["Ncut"]))

M2 = M[ncut:corrTimeM:end]
PyPlot.figure(figsize=(6,4))
PyPlot.title(string(s1," despues del corte"),fontsize=f3)
if parsedArgs["C"]
    M2 = [m - Statistics.mean(M2) for m in M2]
    PyPlot.xlabel(string(n1,"\$ - \\langle \$",n1,"\$ - \\rangle \$"),fontsize=f2)
else
    PyPlot.xlabel(n1,fontsize=f2)
end 

PyPlot.xlim(lim1,lim2)

PyPlot.xticks(fontsize=f1)
PyPlot.yticks(fontsize=f1)
sigma = Statistics.std(M2)
mu = Statistics.mean(M2)
if parsedArgs["P"]
    xs = range(mu-3*sigma,stop = mu+3*sigma,length=200)
    PyPlot.hist(M2,density=true,color="C0",bins=20)
    PyPlot.plot(xs, [NormalDistribution(x,mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=4)) , $(round(sigma,digits=4)) )",color="C1")
    PyPlot.ylabel("Probabilidad",fontsize=f2)
else
    PyPlot.hist(M2,label="\$ \\mu_{X} = $(round(mu,digits=3)) \$ \n \$ \\sigma_{X} = $(round(sigma,digits=3)) \$",color="C0",bins=20)
    PyPlot.ylabel("Frecuencia",fontsize=f2)
end
PyPlot.grid()
PyPlot.legend(fontsize=f2,loc=0)
PyPlot.tight_layout()
#PyPlot.savefig(string(s1,"Despues",".png"),dpi=300)
PyPlot.show()
#PyPlot.close_figs()


    
PyPlot.title(string(s2," antes del corte"),fontsize=f3)
if parsedArgs["C"]
    E1 = [m - Statistics.mean(E) for m in E[1:corrTimeE:end]]
    PyPlot.xlabel(string(n2,"\$ - \\langle \$",n2,"\$ - \\rangle \$"),fontsize=f2)
else
    E1 = E[1:corrTimeE:end]
    PyPlot.xlabel(n2,fontsize=f2)
end
lim1 = minimum(E1)
lim2 = maximum(E1)
PyPlot.xlim(lim1,lim2)
PyPlot.ylabel("Probabilidad",fontsize=f2)

PyPlot.xticks(fontsize=f1)
PyPlot.yticks(fontsize=f1)
sigma = Statistics.std(E1)
mu = Statistics.mean(E1)
if parsedArgs["P"]
    xs = range(mu-3*sigma,stop=mu+3*sigma,length=200)
    PyPlot.hist(E1,density=true,color="C2",bins=20)
    PyPlot.plot(xs, [NormalDistribution(x,mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=4)) , $(round(sigma,digits=4)) )",color="C1")
    PyPlot.ylabel("Probabilidad",fontsize=f2)
else
    PyPlot.hist(E1,label="\$ \\mu_{X} = $(round(mu,digits=3)) \$ \n \$ \\sigma_{X} = $(round(sigma,digits=3)) \$",color="C2",bins=20)
    PyPlot.ylabel("Frecuencia",fontsize=f2)
end
PyPlot.grid()
PyPlot.legend(fontsize=f2,loc=0)
PyPlot.tight_layout()
#PyPlot.savefig(string(s2,"Antes",".png"),dpi=300)
PyPlot.show()
#PyPlot.close_figs()


ncut = Int64(floor(length(E)*parsedArgs["Ncut"]))

E2 = E[ncut:corrTimeE:end]
PyPlot.figure(figsize=(6,4))
PyPlot.title(string(s2," despues del corte"),fontsize=f3)
if parsedArgs["C"]
    E2 = [m - Statistics.mean(E2) for m in E2]
    PyPlot.xlabel(string(n2,"\$ - \\langle \$",n2,"\$ - \\rangle \$"),fontsize=f2)
else
    PyPlot.xlabel(n2,fontsize=f2)
end 

#PyPlot.xlim(lim1,lim2)

PyPlot.xticks(fontsize=f1)
PyPlot.yticks(fontsize=f1)
sigma = Statistics.std(E2)
mu = Statistics.mean(E2)
if parsedArgs["P"]
    xs = range(mu-3*sigma,stop = mu+3*sigma,length=200)
    PyPlot.hist(E2,density=true,color="C2",bins=20)
    PyPlot.plot(xs, [NormalDistribution(x,mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=4)) , $(round(sigma,digits=4)) )",color="C1")
    PyPlot.ylabel("Probabilidad",fontsize=f2)
else
    PyPlot.hist(E2,label="\$ \\mu_{X} = $(round(mu,digits=3)) \$ \n \$ \\sigma_{X} = $(round(sigma,digits=3)) \$",color="C2",bins=20)
    PyPlot.ylabel("Frecuencia",fontsize=f2)
end
PyPlot.legend(fontsize=f2,loc=1)
PyPlot.grid()
PyPlot.tight_layout()
#PyPlot.savefig(string(s2,"Despues",".png"),dpi=300)
PyPlot.show()
#PyPlot.close_figs()




#=
function ThermalizationPlots(X,model,metaParam,simulParam,algoParam,cut=1/2)
    if model=="normal"
        enerFunc=StatEnsemble.NormalEnergy
        n1 = "\$ \\mathcal{M} / N \$"
        n2 = "\$ U / N \$"
        s1 = "Magnetizacion"
        s2 = "Energia"
    else
        enerFunc=StatEnsemble.PenalizedEnergy2
        n1 = "\$ \\rho  \$"
        n2 = "\$ U / N \$"
        s1 = "Densidad"
        s2 = "Energia"
    end
    f1=16
    f2=20
    f3=22
    size=20

    neigLatt=Geometry.BuildLattices(metaParam,model)[2]
    steps=[(i-1)*(algoParam[2]) for i in 1:length(X)]
    den=[abs(sum(x))/(metaParam[1]^metaParam[2]) for x in X]
    print("Calculating energy ... ")
    en=[enerFunc(x,neigLatt,simulParam)/(metaParam[1]^metaParam[2]) for x in X]
    println("done")

    println("Making first plot")


    PyPlot.figure(figsize=(6,4))
    if parsedArgs["S"]
        PyPlot.scatter(steps,den,s=size)
    end
    PyPlot.plot(steps,den)
    PyPlot.title("Termalización",fontsize=f3)
    PyPlot.xlabel("Pasos de Monte Carlo",fontsize=f2)
    PyPlot.ylabel(n1,fontsize=f2)
    PyPlot.grid()
    PyPlot.xticks(fontsize=f1)
    PyPlot.yticks(fontsize=f1)
    PyPlot.ylim(0,1)
    PyPlot.tight_layout()
    PyPlot.savefig(string("Termalization",s1,".png"),dpi=300)
    PyPlot.close_figs()


    println("Making second plot")
    
    
    PyPlot.figure(figsize=(6,4))
    if parsedArgs["S"]
        PyPlot.scatter(steps,en,color="C2",s=size)
    end
    PyPlot.plot(steps,en,color="C2")
    PyPlot.title("Termalización",fontsize=f3)
    PyPlot.xlabel("Pasos de Monte Carlo",fontsize=f2)
    PyPlot.ylabel(n2,fontsize=f2)
    PyPlot.grid()
    PyPlot.xticks(fontsize=f1)
    PyPlot.yticks(fontsize=f1)
    PyPlot.ylim(-4,1)
    PyPlot.tight_layout()
    PyPlot.savefig(string("Termalization",s2,".png"),dpi=300)
    PyPlot.close_figs()

    ncut=Int64(floor(length(X)*cut))
    println("Making third plot")


    PyPlot.figure(figsize=(6,4))
    
    PyPlot.title(string(s1," antes del corte"),fontsize=f3)
    if parsedArgs["C"]
        den1 = [m - Statistics.mean(den) for m in den]
        PyPlot.xlabel(string(n1,"\$ - \\langle \$",n1,"\$ - \\rangle \$"),fontsize=f2)
    else
        den1 = den
        PyPlot.xlabel(n1,fontsize=f2)
    end
    lim1 = minimum(den1)
    lim2 = maximum(den1)
    PyPlot.xlim(lim1,lim2)
    PyPlot.ylabel("Probabilidad",fontsize=f2)
    PyPlot.grid()
    PyPlot.xticks(fontsize=f1)
    PyPlot.yticks(fontsize=f1)
    sigma = Statistics.std(den1)
    mu = Statistics.mean(den1)
    if parsedArgs["P"]
        xs = range(mu-3*sigma,stop=mu+3*sigma,length=200)
        PyPlot.hist(den1,density=true,color="C0",bins=20)
        PyPlot.plot(xs, [NormalDistribution(x,mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=4)) , $(round(sigma,digits=4)) )",color="C1")
        PyPlot.ylabel("Probabilidad",fontsize=f2)
    else
        PyPlot.hist(den1,label="\$ \\mu_{X} = $(round(mu,digits=3)) \$ \n \$ \\sigma_{X} = $(round(sigma,digits=3)) \$",color="C0",bins=20)
        PyPlot.ylabel("Frecuencia",fontsize=f2)
    end
    PyPlot.legend(fontsize=f2,loc=0)
    PyPlot.tight_layout()
    PyPlot.savefig(string(s1,"Antes",".png"),dpi=300)
    PyPlot.close_figs()


    
    println("Making fourth plot")
    
    den2 = den[ncut:end]
    PyPlot.figure(figsize=(6,4))
    PyPlot.title(string(s1," despues del corte"),fontsize=f3)
    if parsedArgs["C"]
        den2 = [m - Statistics.mean(den2) for m in den2]
        PyPlot.xlabel(string(n1,"\$ - \\langle \$",n1,"\$ - \\rangle \$"),fontsize=f2)
    else
        PyPlot.xlabel(n1,fontsize=f2)
    end 

    PyPlot.xlim(lim1,lim2)
    PyPlot.grid()
    PyPlot.xticks(fontsize=f1)
    PyPlot.yticks(fontsize=f1)
    sigma = Statistics.std(den2)
    mu = Statistics.mean(den2)
    if parsedArgs["P"]
        xs = range(mu-3*sigma,stop = mu+3*sigma,length=200)
        PyPlot.hist(den2,density=true,color="C0",bins=20)
        PyPlot.plot(xs, [NormalDistribution(x,mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=4)) , $(round(sigma,digits=4)) )",color="C1")
        PyPlot.ylabel("Probabilidad",fontsize=f2)
    else
        PyPlot.hist(den2,label="\$ \\mu_{X} = $(round(mu,digits=3)) \$ \n \$ \\sigma_{X} = $(round(sigma,digits=3)) \$",color="C0",bins=20)
        PyPlot.ylabel("Frecuencia",fontsize=f2)
    end
    PyPlot.legend(fontsize=f2,loc=0)
    PyPlot.tight_layout()
    PyPlot.savefig(string(s1,"Despues",".png"),dpi=300)
    PyPlot.close_figs()
    
    println("Making fifth plot")


    PyPlot.figure(figsize=(6,4))
    PyPlot.title(string(s2," antes del corte"),fontsize=f3)
    if parsedArgs["C"]
        en1 = [e - Statistics.mean(en) for e in en]
        PyPlot.xlabel(string(n2,"\$ - \\langle \$",n2,"\$ - \\rangle \$"),fontsize=f2)
    else
        en1 = en
        PyPlot.xlabel(string(n2),fontsize=f2)
    end

    lim1 = minimum(en1)
    lim2 = maximum(en1)
    PyPlot.xlim(lim1,lim2)
    PyPlot.grid()
    PyPlot.xticks(fontsize=f1)
    PyPlot.yticks(fontsize=f1)
    sigma = Statistics.std(en1)
    mu = Statistics.mean(en1)
    if parsedArgs["P"]
        xs = range(mu-3*sigma,stop=mu+3*sigma,length=200)
        PyPlot.hist(en1,density=true,color="C2",bins=20)
        PyPlot.plot(xs, [NormalDistribution(big(x),mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=4)) , $(round(sigma,digits=4)) )",color="C1")
        PyPlot.ylabel("Probabilidad",fontsize=f2)
    else
        PyPlot.hist(en1,label="\$ \\mu_{X} = $(round(mu,digits=3)) \$ \n \$ \\sigma_{X} = $(round(sigma,digits=3)) \$",color="C2",bins=20)
        PyPlot.ylabel("Frecuencia",fontsize=f2)
    end
    PyPlot.legend(fontsize=f2,loc=0)
    PyPlot.tight_layout()
    PyPlot.savefig(string(s2,"Antes",".png"),dpi=300)
    PyPlot.close_figs()


    
    println("Making sixth plot")
    en2 = en[ncut:end]
    PyPlot.figure(figsize=(6,4))
    PyPlot.title(string(s2," despues del corte"),fontsize=f3)
    if parsedArgs["C"]
        en2 = [e - Statistics.mean(en2) for e in en2]
        PyPlot.xlabel(string(n2,"\$ - \\langle \$",n2,"\$ - \\rangle \$"),fontsize=f2)
    else
        PyPlot.xlabel(string(n2),fontsize=f2)
    end

    PyPlot.xlim(lim1,lim2)
    PyPlot.ylabel("Probabilidad",fontsize=f2)
    PyPlot.grid()
    PyPlot.xticks(fontsize=f1)
    PyPlot.yticks(fontsize=f1)
    sigma = Statistics.std(en2)
    mu = Statistics.mean(en2)
    if parsedArgs["P"]
        xs = range(mu-3*sigma,stop=mu+3*sigma,length=200)
        PyPlot.hist(en2,density=true,color="C2",bins=20)
        PyPlot.plot(xs, [NormalDistribution(big(x),mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=4)) , $(round(sigma,digits=4)) )",color="C1")
        PyPlot.ylabel("Probabilidad",fontsize=f2)
    else
        PyPlot.hist(en2,label="\$ \\mu_{X} = $(round(mu,digits=3)) \$ \n \$ \\sigma_{X} = $(round(sigma,digits=3)) \$",color="C2",bins=20)
        PyPlot.ylabel("Frecuencia",fontsize=f2)
    end
    PyPlot.legend(fontsize=f2,loc=0)
    PyPlot.tight_layout()
    PyPlot.savefig(string(s2,"Despues",".png"),dpi=300)
    PyPlot.close_figs()

    #=

        en ceros: mu = -3.325, C = 0.7

        enmedio: mu = -1.575, C = 0.9

        enmediosiempre: mu = -0.525, C = 1.5

        arriba: mu = -0.525, C = 0.5

    =#

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


println()
println("Making Analysis")
println()

data = InOut.ReadSingleSimul(parsedArgs["path"])
cd(parsedArgs["path"])
ThermalizationPlots(data[1],"cycle",data[2],data[3],data[4])
cd(original)
=#