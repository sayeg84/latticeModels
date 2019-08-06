
using Distributed
addprocs(1)
@everywhere include("structs.jl")
@everywhere include("lattices.jl")
#=
let y  = SpinLattice(Lattices.PeriodicSquareLatticeNeighbors,2,2)
@everywhere x = copy($y)
end
=#
println()
println()
println("distribuido")
println()
println()
@everywhere x =SpinLattice(Lattices.PeriodicSquareLatticeNeighbors,2,2) 
@show x.linearLatt
@sync @distributed for i in 1:1
    global x
    X = []
    S = []
    push!(X,sum(copy(x).linearLatt))
    push!(S,sum(x.linearLatt))
    for j in 1:3
        ChangeSpin!(x,rand(1:4))
        push!(X,sum(copy(x).linearLatt))
        push!(S,sum(x.linearLatt))
    end
    @show X
    @show S
end

println()
println()
println("sin distribuido")
println()
println()
x =SpinLattice(Lattices.PeriodicSquareLatticeNeighbors,2,2) 
@show x.linearLatt
for i in 1:1
    global x
    X = []
    S = []
    push!(X,sum(copy(x).linearLatt))
    push!(S,sum(x.linearLatt))
    for j in 1:3
        ChangeSpin!(x,rand(1:4))
        push!(X,sum(copy(x).linearLatt))
        push!(S,sum(x.linearLatt))
    end
    @show X
    @show S
end
