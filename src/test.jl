include("statEnsemble.jl")
include("auxiliar.jl")
using Auxiliar
#=
using StatEnsemble
param=[10,1,0,10^6,10^4]
x=StatEnsemble.MagArray(collect(linspace(-2,2,11)),param,test=false)
println(x)
=#
L=collect(1:2:20)
Auxiliar.MirrorList!(L)
println(L)