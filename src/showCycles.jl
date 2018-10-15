include("graphTheory.jl")
include("auxiliar.jl")
function PrintByLines(M)
    println()
    println()
    for i in 1:(size(M)[1])
        @show M[i,:]
    end
    println()
    println()
end
using GraphTheory,Auxiliar
M=Array{Int8,2}([-1 1 1 1 -1 1 -1 -1; -1 1 1 1 -1 1 -1 1; 1 1 1 1 1 1 1 -1; 1 -1 -1 1 1 -1 -1 -1; 1 -1 1 -1 1 -1 1 -1; -1 -1 1 1 -1 -1 -1 1; -1 -1 1 1 -1 -1 -1 1; 1 -1 -1 1 1 -1 1 -1])

neigLatt=Auxiliar.NeighborIndexLattice(M,Auxiliar.SquareLatticeNeighborsIndex)
@time C=GraphTheory.Cycles(M,neigLatt)
A=zeros(Int8,size(M))
for pos in C 
    A[pos]=1
end
PrintByLines(M)
PrintByLines(A)