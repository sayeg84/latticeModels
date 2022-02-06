
include("hamiltonians.jl")

abstract type AbstractModel end


struct Model{T1<:AbstractSystem,T2<:AbstractHamiltonian} <: AbstractModel
    system::T1
    hamiltonian::T2
    variables::Tuple{Vararg{Function}}
end

abstract type AbstractVariable end



if abspath(PROGRAM_FILE) == @__FILE__

    const system = LatticeGas("square",(4,4),random=true)
    energy1 = IsingHamiltonian(2,3)
    energy2 = IsingHamiltonian(2,4)
    m1 = Model(system,energy1,(SitesSum,))
    m2 = Model(system,energy1,(SitesSum,))
    @show m1.system.sites
    m1.system.sites[1] = 2
    @show m2.system.sites
end