println()
println("Metropolis Thermalization Single Data analyzer")
println()

println()
println("Importing libraries")
println()

include("statEnsemble.jl")
include("auxiliar.jl")
include("inOut.jl")
include("geometry.jl")

using PyPlot, DelimitedFiles, ArgParse, Statistics

function NormalDistribution(x ; mu = 0 , sigma = 1)
    return 1/sqrt(2*pi*sigma^2) * exp(-(x-mu)^2/(2*sigma^2))
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
            help = "Directory path to analyse set"
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

function ThermalizationPlots(X,model,geoParam,simulParam,algoParam,cut=1/2)
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

    neigLatt=Geometry.BuildLattices(geoParam,model)[2]
    steps=[(i-1)*(algoParam[2]) for i in 1:length(X)]
    den=[abs(sum(x))/(geoParam[1]^geoParam[2]) for x in X]
    print("Calculating energy ... ")
    en=[enerFunc(x,neigLatt,simulParam)/(geoParam[1]^geoParam[2]) for x in X]
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
        PyPlot.hist(den1,density=true,label="Histograma",color="C0")
        PyPlot.plot(xs, [NormalDistribution(x,mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=4)) , $(round(sigma,digits=4)) )",color="C1")
        PyPlot.ylabel("Probabilidad",fontsize=f2)
    else
        PyPlot.hist(den1,label="\$ \\mu_{X} = $(round(mu,digits=3)) \$ \n \$ \\sigma_{X} = $(round(sigma,digits=3)) \$",color="C0")
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
        PyPlot.hist(den2,density=true,label="histograma",color="C0")
        PyPlot.plot(xs, [NormalDistribution(x,mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=4)) , $(round(sigma,digits=4)) )",color="C1")
        PyPlot.ylabel("Probabilidad",fontsize=f2)
    else
        PyPlot.hist(den2,label="\$ \\mu_{X} = $(round(mu,digits=3)) \$ \n \$ \\sigma_{X} = $(round(sigma,digits=3)) \$",color="C0")
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
        PyPlot.hist(en1,density=true,label="histograma",color="C2")
        PyPlot.plot(xs, [NormalDistribution(big(x),mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=4)) , $(round(sigma,digits=4)) )",color="C1")
        PyPlot.ylabel("Probabilidad",fontsize=f2)
    else
        PyPlot.hist(en1,label="\$ \\mu_{X} = $(round(mu,digits=3)) \$ \n \$ \\sigma_{X} = $(round(sigma,digits=3)) \$",color="C2")
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
        PyPlot.hist(en2,density=true,label="histograma",color="C2")
        PyPlot.plot(xs, [NormalDistribution(big(x),mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=4)) , $(round(sigma,digits=4)) )",color="C1")
        PyPlot.ylabel("Probabilidad",fontsize=f2)
    else
        PyPlot.hist(en2,label="\$ \\mu_{X} = $(round(mu,digits=3)) \$ \n \$ \\sigma_{X} = $(round(sigma,digits=3)) \$",color="C2")
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

data = InOut.ReadSingleSimul(parsedArgs["Dirname"])
cd(parsedArgs["Dirname"])
ThermalizationPlots(data[1],"cycle",data[2],data[3],data[4])
cd(original)