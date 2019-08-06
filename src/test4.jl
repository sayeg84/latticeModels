
using Distributed
addprocs(1)
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
@everywhere x = rand([-1,1],2,2)
@show x
@sync @distributed for i in 1:1
    global x
    X = []
    S = []
    push!(X,sum(copy(x)))
    push!(S,sum(x))
    for j in 1:3
        x[rand(1:4)] *= -1
        push!(X,sum(copy(x)))
        push!(S,sum(x))
    end
    @show X
    @show S
end

println()
println()
println("sin distribuido")
println()
println()
x = rand([-1,1],2,2)
@show x
for i in 1:1
    global x
    X = []
    S = []
    push!(X,sum(copy(x)))
    push!(S,sum(x))
    for j in 1:3
        x[rand(1:4)] *= -1
        push!(X,sum(copy(x)))
        push!(S,sum(x))
    end
    @show X
    @show S
end
