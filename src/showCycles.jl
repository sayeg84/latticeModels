include("graphTheory.jl")
include("auxiliar.jl")
using GraphTheory,Auxiliar
M= [
    0 1 1 1 1 1 1 1 1 0 0 ;
    0 0 0 0 1 0 0 0 1 0 0 ;
    0 0 1 1 1 1 1 0 1 1 1 ;
    0 0 1 0 1 0 1 0 1 0 1 ;
    0 0 1 0 1 1 1 1 1 0 1 ;
    0 0 1 0 0 0 1 0 0 0 0 ;
    0 0 1 1 1 1 1 0 0 0 0 ;
    0 0 0 0 0 0 0 0 0 0 0 ;
    ]
neigLatt=Auxiliar.NeighborIndexLattice(M,Auxiliar.SquareLatticeNeighborsIndex)
@time C=GraphTheory.Cycles(M,neigLatt)
A=zeros(size(M))
for pos in C 
    A[pos]=1
end
@show M
@show A