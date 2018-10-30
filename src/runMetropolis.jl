include("inOut.jl")
include("algorithms.jl")
include("statEnsemble.jl")
include("auxiliar.jl")
#using Algorithms,InOut,StatEnsemble, ArgParse
using ArgParse, Statistics, Dates

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

#=
function metropolisMontecarlo(initLatt::Array{Int8,2},neigLatt::Array{Array{CartesianIndex{2},1},2},param::Array{Float64,1},tempArray::Array{Float64,1},nAveg::Int64)

    time=Dates.time()

    mag=zeros(length(tempArray))
    ener=zeros(length(tempArray))
    cv=zeros(length(tempArray))
    for j in 0:(length(tempArray)-1)
        temp=tempArray[end-j]

        tempMag=[]
        tempEner=[]
        tempCV=[]
        X=[]


        for i in 1:nAveg

            X=Algorithms.Metropolis(temp,param,StatEnsemble.Energy,initLatt,neigLatt)
            initial=convert(Int64,floor(length(X)*3/4))
            M=[abs(sum(X[k])) for k in initial:length(X)]
            push!(tempMag,Statistics.mean(M))
            en=[StatEnsemble.Energy(X[k],param,neigLatt) for k in initial:length(X)]
            push!(tempEner,Statistics.mean(en))
            push!(tempCV,Statistics.mean([x^2 for x in en]))

            print(tempArray[end-j])
            print("-")
            println(i)
        end
        initLatt=copy(X[end])
        mag[end-j]=Statistics.mean(tempMag)
        ener[end-j]=Statistics.mean(tempEner)
        cv[end-j]=Statistics.mean(tempCV)
    end
    cv=[(cv[i]-ener[i]^2)/(tempArray[i]^2) for i in 1:length(tempArray)]
    #normalizating energy
    ener=ener*1/(param[1]^2)
    mag=mag*1/(param[1]^2)
    cv=cv*1/(param[1]^2)

    println()
    println("Writing Output")
    InOut.MakeDirectories()
    InOut.MakePlots(tempArray,mag,ener,cv,param)
    time=Dates.time()-time
    InOut.MakeTable(tempArray,mag,ener,cv,param,time)


end
=#
println("Begining")
println()

tempArray=collect(range(0.1,stop=5,length=50))
N=parsedArgs["Nlatt"]

param=[N,
# coupling constant J
parsedArgs["Jconst"],
# magnetic field B J
parsedArgs["Bfield"],
# maximum steps
10^6,
# save frecuency (Metropolis) or final iterations (Wang Landau)
10^4,
# energy bins (WangLandau)
ceil(N^2/2)-N,
#cycle constant
parsedArgs["Cconst"]]

println("Simulating")
println()

init=rand(Array{Int8}([1,-1]),Int64(param[1]),Int64(param[1]))
neigLatt=Auxiliar.NeighborIndexLattice(init,Auxiliar.SquareLatticeNeighborsIndex)

t=[0.6]
mus=range(-3,stop=-1,length=50)
println("Simulating")
println()

    mag=zeros(length(mus))
    ener=zeros(length(mus))
    cv=zeros(length(mus))
    init=rand(Array{Int8}([1,-1]),convert(Int64,param[1]),convert(Int64,param[1]))
    for j in 0:(length(mus)-1)
        temp=t[1]
        temp1=[]
        temp2=[]
        temp3=[]
        X=[]
        param[3]=mus[j+1]
        @time for i in 1:5
            println(param[1])
            neigLatt=Auxiliar.NeighborIndexLattice(init,Auxiliar.SquareLatticeNeighborsIndex)
            X=Algorithms.Metropolis(temp,param,StatEnsemble.Energy,init,neigLatt)
            initial=convert(Int64,floor(length(X)*3/4))
            M=[abs(sum(X[k])) for k in initial:length(X)]
            push!(temp1,Statistics.mean(M))
            en=[StatEnsemble.Energy(x,param,neigLatt) for x in X]
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

#=
metropolisMontecarlo(init,neigLatt,param,tempArray,5)
=#

