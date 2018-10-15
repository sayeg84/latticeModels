
include("statEnsemble.jl")
include("auxiliar.jl")
include("algorithms.jl")
using Auxiliar
using Algorithms
using StatEnsemble
param=[6,1,0,10^6,10^4,10]
#=
x=StatEnsemble.MagArray(collect(linspace(-2,2,11)),param,test=true)
println(x)
#push!(param,convert(Int64,ceil(N^2/2)-N))
#initLatt=rand([-1,1],param[1],param[1])
initLatt=ones(param[1],param[1])
X=Algorithms.WangLandauSimple(param,initLatt,test=false)
println(X[1])
println(X[2])
println()
println("partition")
println(StatEnsemble.Partition(X[2],X[1],5,param))
println()
println()
a=StatEnsemble.DOSMag(X[2],X[1],100,param,test=false)
=#
#m=ones(param[1],param[1])
m=[1.0 -1.0 1.0 -1.0 1.0 -1.0;
-1.0 1.0 -1.0 1.0 -1.0 1.0;
1.0 -1.0 1.0 -1.0 1.0 -1.0;
-1.0 1.0 -1.0 1.0 -1.0 1.0;
1.0 -1.0 1.0 -1.0 1.0 -1.0;
-1.0 1.0 -1.0 1.0 -1.0 1.0]
@time idx=Auxiliar.NeighborIndexLattice(m,Auxiliar.SquareLatticeNeighborsIndex)
pos=CartesianIndex((1,1))
x=0
@time for neigPos in idx[pos]
    x=x+m[neigPos]
end
println(x)
x=0
@time x=Auxiliar.SquareLatticeNeighbors(m,[1,1])
println(x)
@time e=StatEnsemble.Energy(m,param,idx)
println(e)
#m[6,4]=-1
@time e=StatEnsemble.Energy(m,param,idx)
println(e)
println(Auxiliar.NeighborSum(m,idx,CartesianIndex((6,4))))