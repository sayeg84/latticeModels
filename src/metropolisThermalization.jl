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
    end
    return parse_args(s)
end
parsedArgs = ParseCommandline()

function ThermalizationPlots(X,model,geoParam,simulParam,algoParam,cut=1/2)
    if model=="normal"
        enerFunc=StatEnsemble.NormalEnergy
        n1 = "\$ \\mathcal{M} / N \$"
        n2 = "\$ E / N \$"
        s1 = "Magnetizacion"
        s2 = "Energia"
    else
        enerFunc=StatEnsemble.PenalizedEnergy2
        n1 = "\$ \\rho / N \$"
        n2 = "\$ E / N \$"
        s1 = "Densidad"
        s2 = "Energia"
    end
    f1=20
    f2=24
    f3=28
    size=20

    neigLatt=Geometry.BuildLattices(geoParam,model)[2]
    steps=[(i-1)*(algoParam[2]) for i in 1:length(X)]
    mag=[abs(sum(x))/(geoParam[1]^geoParam[2]) for x in X]
    en=[enerFunc(x,neigLatt,simulParam)/(geoParam[1]^geoParam[2]) for x in X]
    
    println("Making first plot")

    PyPlot.figure(figsize=(6,4))
    PyPlot.scatter(
        steps,mag,s=size
    )
    PyPlot.plot(
        steps,mag
    )
    PyPlot.title("Termalización")
    PyPlot.xlabel("Pasos de Monte Carlo")
    PyPlot.ylabel(n1)
    PyPlot.grid()

    PyPlot.savefig(string("Termalization",s1,".png"),dpi=300)
    PyPlot.close_figs()

    println("Making second plot")
    PyPlot.figure(figsize=(6,4))
    PyPlot.scatter(
        steps,en,color="C2",s=size
    )
    PyPlot.plot(
        steps,en,color="C2"
    )
    PyPlot.title("Termalización")
    PyPlot.xlabel("Pasos de Monte Carlo")
    PyPlot.ylabel(n2)
    PyPlot.grid()
    PyPlot.savefig(string("Termalization",s2,".png"),dpi=300)
    PyPlot.close_figs()

    ncut=Int64(floor(length(X)*cut))
    epsilon = 0.2

    
    println("Making third plot")
    PyPlot.figure(figsize=(6,4))
    PyPlot.hist(mag,density=true,label="Histograma",color="C0")


    PyPlot.title(string(s1," antes del corte"))
    PyPlot.xlabel("Valor")
    PyPlot.ylabel("Probabilidad")
    PyPlot.grid()

    sigma = Statistics.std(mag)
    mu = Statistics.mean(mag)
    xs = range(sign(minimum(mag)) * (minimum(mag)*(1-epsilon)),stop = sign(maximum(mag))*(maximum(mag)*(1+epsilon)),length=200)
    PyPlot.plot(xs, [NormalDistribution(x,mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=2)) , $(round(sigma,digits=2)) )",color="C1")
    PyPlot.legend(loc=2)
    PyPlot.savefig(string(s1,"Antes",".png"),dpi=300)
    PyPlot.close_figs()

    mag = mag[ncut:end]
    epsilon = 0.05

    
    println("Making fourth plot")
    PyPlot.figure(figsize=(6,4))
    PyPlot.hist(mag,density=true,label="Histograma",color="C0")


    PyPlot.title(string(s1," despues del corte"))
    PyPlot.xlabel("Valor")
    PyPlot.ylabel("Probabilidad")
    PyPlot.grid()

    sigma = Statistics.std(mag)
    mu = Statistics.mean(mag)
    xs = range(sign(minimum(mag)) * (minimum(mag)*(1-epsilon)),stop = sign(maximum(mag))*(maximum(mag)*(1+epsilon)),length=200)
    PyPlot.plot(xs, [NormalDistribution(x,mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=2)) , $(round(sigma,digits=2)) )",color="C1")
    PyPlot.legend(loc=2)
    PyPlot.savefig(string(s1,"Despues",".png"),dpi=300)
    PyPlot.close_figs()


    
    println("Making fifth plot")
    PyPlot.figure(figsize=(6,4))
    PyPlot.hist(en,density=true,label="Histograma",color="C2")


    PyPlot.title(string(s2," antes del corte"))
    PyPlot.xlabel("Valor")
    PyPlot.ylabel("Probabilidad")
    PyPlot.grid()

    sigma = Statistics.std(en)
    mu = Statistics.mean(en)
    xs = range(mu-3*sigma,stop=mu+3*sigma,length=200)
    PyPlot.plot(xs, [NormalDistribution(big(x),mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=2)) , $(round(sigma,digits=2)) )",color="C1")
    PyPlot.legend(loc=2)
    PyPlot.savefig(string(s2,"Antes",".png"),dpi=300)
    PyPlot.close_figs()

    en = en[ncut:end]
    epsilon = 0.05

    
    println("Making sixth plot")
    PyPlot.figure(figsize=(6,4))
    PyPlot.hist(en,density=true,label="Histograma",color="C2")


    PyPlot.title(string(s2," despues del corte"))
    PyPlot.xlabel("Valor")
    PyPlot.ylabel("Probabilidad")
    PyPlot.grid()

    sigma = Statistics.std(en)
    mu = Statistics.mean(en)
    xs = range(mu-3*sigma,stop=mu+3*sigma,length=200)
    PyPlot.plot(xs, [NormalDistribution(x,mu=mu,sigma=sigma) for x in xs],label ="N ( $(round(mu,digits=2)) , $(round(sigma,digits=2)) )",color="C1")
    PyPlot.legend(loc=2)
    PyPlot.savefig(string(s2,"Despues",".png"),dpi=300)
    PyPlot.close_figs()

    
    #=
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
    =#
end

println()
println("Making Analysis")
println()

data = InOut.ReadSingleSimul(parsedArgs["Dirname"])
cd(parsedArgs["Dirname"])
ThermalizationPlots(data[1],"cycle",data[2],data[3],data[4])
cd(original)