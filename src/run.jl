time=Dates.time()
include("inOut.jl")
include("algorithms.jl")
include("statEnsemble.jl")
using Algorithms,InOut,StatEnsemble
println("Begining")
println()
println("Simulating")
t=linspace(0.1,5,50)
param=[]
#asignamos los parámetros
push!(param,5)
push!(param,1)
push!(param,0)
push!(param,10^6)
push!(param,10^4)
#=
mag=zeros(length(t))
ener=zeros(length(t))
cv=zeros(length(t))
init=rand([1,-1],param[1],param[1])
for j in 0:(length(t)-1)
    temp=t[end-j]
    temp1=[]
    temp2=[]
    temp3=[]
    X=[]
    for i in 1:5
        X=Algorithms.Metropolis(temp,param,init)
        initial=convert(Int64,floor(length(X)*3/4))
        M=[abs(sum(X[k])) for k in initial:length(X)]
        push!(temp1,mean(M))
        en=[StatEnsemble.Energy(x,param,StatEnsemble.SquareLatticeNeighbors) for x in X]
        temp4=[en[k] for k in initial:length(X)]
        push!(temp2,mean(temp4))
        push!(temp3,mean([x^2 for x in temp4]))
        print(t[end-j])
        print("-")
        println(i)
    end
    init=copy(X[end])
    mag[end-j]=mean(temp1)
    ener[end-j]=mean(temp2)
    cv[end-j]=mean(temp3)
end
cv=[(cv[i]-ener[i]^2)/(t[i]^2) for i in 1:length(t)]
#normalizating energy
ener=ener*1/(param[1]^2)
mag=mag*1/(param[1]^2)

println()
println("Writing Output")
InOut.MakeDirectories()
InOut.MakePlots(t,mag,ener,cv,param)
time=Dates.time()-time
InOut.MakeTable(t,mag,ener,cv,param,time)
=#
initLatt=rand([-1,1],param[1],param[1])
X=Algorithms.WangLandau(param,initLatt,20;test=true)
println(X[1])
println(X[2])
println(X[3])
