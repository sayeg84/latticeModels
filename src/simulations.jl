include("model.jl")

abstract type AbstractSimulation end

struct NVTMetropolis{F<:AbstractMutation} <: AbstractSimulation
    kT::Float64       #ignored for some models
    mcSweeps::Integer
    mutation::F
end



struct NPTMetropolis{F<:AbstractMutation} <: AbstractSimulation
    kT::Float64       #ignored for some models
    mcSweeps::Integer
    mutation::F
end

function advance!(model::AbstractModel,simul::NVTMetropolis,oldEnergy::Real)
    newSystem = mutate(model.system,simul.mutation)
    newEnergy = Energy(newSystem,model.hamiltonian)
    prob = (-newEnergy+oldEnergy)/simul.kT
    r = log(rand())
    if r < prob
        model.system = newSystem
        return true, newEnergy
    else
        return false, oldEnergy
    end
end

function runSimulation(initialModel::AbstractModel,simul::NVTMetropolis;saveChain::bool=true)
    n_steps = simul.mcSweeps*initialModel.system.n
    vars = zeros(n_steps,length(model.variables)+1)
    if saveChain
        chain = zeros(n_steps,length(simul.mutation))
        for t in 2:n_steps
            accepted, vars[t,1] = (model,simul,vars[t-1,1])
            if accepted

            else

            end
        end
        return model, vars, chain
    else
        return model, vars
    end
end