include("systems.jl")

abstract type AbstractHamiltonian end

struct IsingHamiltonian <: AbstractHamiltonian
    B::Float64
    J::Float64
end

# benchmarks seems to indicate that this is the fastest way to calculate energy
# this must be related to not making memory allocations for any object
function Energy(system::AbstractSystem,hamil::IsingHamiltonian)::Float64
    ener = 0
    for i in 1:N(system)
        ener += -system.sites[i]*hamil.B 
        for j in system.edgList[i]
            ener -= hamil.J/2*system.sites[i]*system.sites[j]
        end
    end
    return ener
end


struct PenalizedHamiltonian <: AbstractHamiltonian
    B::Float64
    J::Float64
    C::Float64
end

function Energy(system::AbstractSystem,hamil::PenalizedHamiltonian)::Float64
    ener = - hamil.B*sum(system.sites) 
           - hamil.J*sum(system.neigSum.*system.sites)/2
           + hamil.C*EdgInCyc(sys)
    return ener
end

if abspath(PROGRAM_FILE) == @__FILE__
    #   Testing
    H1 = IsingHamiltonian(rand(),rand())
    H2 = PenalizedHamiltonian(rand(),rand(),rand())
    system = LatticeGas("square",(100,100),random=true)
end



