include("inOut.jl")
include("algorithms.jl")
include("statEnsemble.jl")
include("auxiliar.jl")
#using Algorithms,InOut,StatEnsemble,ArgParse
using ArgParse, Statistics, Dates
time=Dates.time()
function ParseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-N", "--Nlatt" 
            help = "Lattice size"
            arg_type=Int64
            default = 10
         "-J", "--Jconst"
            help = "Coupling constant of Ising model"
            arg_type = Float64
            default = 1.0
         "-B", "--Bfield"
            help = "Magnetic field"
            arg_type = Float64
            default = 0.0
         "-C", "--Cconst"
            help = "Cycle constant"
            arg_type = Float64
            default = 1.0
    end
    return parse_args(s)
end

parsedArgs = ParseCommandline()


println("Begining")
println()



N=parsedArgs["Nlatt"]
param=[N,
# coupling constant J
parsedArgs["Jconst"],
# magnetic field B J
parsedArgs["Bfield"],
# maximum steps
10^5,
# save frecuency (Metropolis) or final iterations (Wang Landau)
10^3,
# energy bins (WangLandau)
convert(Int64,ceil(N^2/2)-N),
#cycle constant
parsedArgs["Cconst"]]

t=[0.6]
mus=range(-3,stop=-1,length=50)
println("Simulating")
println()

    mag=zeros(length(mus))
    ener=zeros(length(mus))
    cv=zeros(length(mus))
    init=rand(Array{Int8}([1,0]),convert(Int64,param[1]),convert(Int64,param[1]))
    for j in 0:(length(mus)-1)
        temp=t[1]
        temp1=[]
        temp2=[]
        temp3=[]
        X=[]
        param[3]=mus[end-j]
        @time for i in 1:5
            println(param[1])
            neigLatt=Auxiliar.NeighborIndexLattice(init,Auxiliar.SquareLatticeNeighborsIndex)
            X=Algorithms.Metropolis(temp,param,StatEnsemble.PenalizedEnergy2,init,neigLatt)
            initial=convert(Int64,floor(length(X)*3/4))
            M=[abs(sum(X[k])) for k in initial:length(X)]
            push!(temp1,Statistics.mean(M))
            en=[StatEnsemble.PenalizedEnergy2(x,param,neigLatt) for x in X]
            temp4=[en[k] for k in initial:length(X)]
            push!(temp2,Statistics.mean(temp4))
            push!(temp3,Statistics.mean([x^2 for x in temp4]))
            print(mus[end-j])
            print("-")
            println(i)
        end
        init=copy(X[end])
        mag[end-j]=Statistics.mean(temp1)
        ener[end-j]=Statistics.mean(temp2)
        cv[end-j]=Statistics.mean(temp3)
    end
    cv=[(cv[i]-ener[i]^2)/(t[1]^2) for i in 1:length(mus)]
    #normalizating energy
    ener=ener*1/(param[1]^2)
    mag=mag*1/(param[1]^2)
    cv=cv*1/(param[1]^2)

println()
println("Writing Output")
InOut.MakeDirectories()
InOut.MakePlots(mus,mag,ener,cv,param)
time=Dates.time()-time
InOut.MakeTable(mus,mag,ener,cv,param,time)