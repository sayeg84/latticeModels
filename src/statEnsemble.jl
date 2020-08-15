
include("cycles.jl")
module StatEnsemble
    
    import ..SpinLattice
    import ..LatticeGas
    import ..IsingModel
    import ..ChangeSpin
    import ..N
    import ..Cycles.lgEdgesInCycles
    import ..Cycles.lgNewCycles
    import ..Cycles.ciclos2
    import ..Cycles.EdgInCyc
    import ..Cycles.CycSubSys


    using Statistics,LinearAlgebra
    """
        NormalEnergy(latt,neigLatt,simulParam)

        Calculates normal energy of Ising model for a lattice `latt` with simulation parameters `simulParam` and neighbor lattice `neigLatt`
    """
    function NormalEnergy(x::IsingModel,simulParam;printLog=false)
        e=0
        B=simulParam[1]
        J=simulParam[2]
        e = -B*sum(x.linearLatt) - J * sum(x.linearNeigSumLatt .* x.linearLatt ) /2
        return e
    end

    """
        PenalizedEnergy(latt,neigLatt,simulParam)

        Calculates energy of the modified Ising model for a lattice `latt` with simulation parameters `simulParam` and neighbor lattice `neigLatt` searching for cycles using the algorithm contained in "graphTheory.jl"
    """
    function PenalizedEnergy(x::IsingModel,simulParam;printLog=false)
        e=0
        B=simulParam[1]
        J=simulParam[2]
        C=simulParam[3]
        e = -B*sum(x.linearLatt) - J * sum(x.linearNeigSumLatt .* x.linearLatt ) /2 + C * EdgInCyc(x)
        return e
    end

    function ProbAta(ener,x::LatticeGas,pos::Integer,simulParam;printLog=false)
        B=simulParam[1]
        J=simulParam[2]
        C=simulParam[3]
        x2 =ChangeSpin(x,pos)
        newcyc = lgNewCycles(x,x2)
        deltae = -B-J*x.linearNeigSumLatt[pos] + 2*((B+J*x.linearNeigSumLatt[pos])*x.linearLatt[pos]) + C * newcyc
        return exp(-big(deltae)/simulParam[4])
    end



    function Prob(x::IsingModel,pos,enerFunc,simulParam)
        e1=enerFunc(x,simulParam)
        e2=enerFunc(ChangeSpin(x,pos),simulParam)
        return exp(big((-e2+e1)/simulParam[4]))
    end

    function ProbOptimal(x::SpinLattice,pos,enerFunc,simulParam)
        if enerFunc == NormalEnergy
            deltaener = 2*x.linearLatt[pos]*(simulParam[1] + simulParam[2]* sum(x.linearLatt[x.linearNeigLatt[pos]]))
        end
        return exp(-deltaener/simulParam[4])
    end


    function ProbOptimal(x::LatticeGas,pos,enerFunc,simulParam)
        if enerFunc == NormalEnergy
            deltaener = (2*x.linearLatt[pos]-1)*(simulParam[1] + simulParam[2]* sum(x.linearLatt[x.linearNeigLatt[pos]]))
        end
        return exp(-deltaener/simulParam[4])
    end


    function Magnetization(sys::IsingModel;absolute=true)
        if absolute
            a = abs(sum(sys.linearLatt))
        else
            a = sum(sys.linearLatt)
        end
        return a/N(sys)
    end
    
    function Magnetization(X::Array{IsingModel,1};absolute=true)
        if absolute
            a = Statistics.mean([abs(sum(x.linearLatt)) for x in X])
        else
            a = Statistics.mean([sum(x.linearLatt) for x in X])
        end
        return a/N(X[1])
    end

    function MagnetizationUpdate(sys::SpinLattice,pos::Integer)
        return -2*sys.linearLatt[pos]
    end

    function MagnetizationUpdate(sys::LatticeGas,pos::Integer)
        return 1 - 2*sys.linearLatt[pos]
    end


    function MagneticSucep(X::Array{IsingModel,1},simulParam;absolute=true)
        if absolute
            a=Statistics.mean([abs(sum(x.linearLatt))^2 for x in X])
            b=(Statistics.mean([abs(sum(x.linearLatt)) for x in X]))^2
        else
            a=Statistics.mean([sum(x.linearLatt)^2 for x in X])
            b=(Statistics.mean([sum(x.linearLatt) for x in X]))^2
        end
        return 1/(simulParam[4]*(N(X[1])))*(a-b)
    end

    function InternalEnergy(sys::IsingModel,enerFunc,simulParam)
        a = enerFunc(sys,simulParam)
        return a/N(X[1])
    end

    function InternalEnergy(X::Array{IsingModel,1},enerFunc,simulParam)
        a = Statistics.mean([enerFunc(x,simulParam) for x in X])
        return a/N(X[1])
    end


    function HeatCapacity(X::Array{IsingModel,1},enerFunc,simulParam)
        a=Statistics.mean([enerFunc(x,simulParam)^2 for x in X])
        b=(Statistics.mean([enerFunc(x,simulParam) for x in X]))^2
        return 1/((simulParam[4]^2)*(N(X[1])))*(a-b)
    end


    function Neigs(r,pos,posLatt,delta = 1e-1)
        center = posLatt[pos]
        neigs = Array{CartesianIndex{length(size(posLatt))},1}()
        for indx in CartesianIndices(size(posLatt))
            if abs(norm(posLatt[indx] - center)-r) < delta
                push!(neigs,indx)
            end
        end
        return neigs
    end


    function PairCorrelation(x::IsingModel,r,posLatt,delta)
        m2 = (sum(x.linearLatt)/N(x))^2
        s = 0
        for indx in CartesianIndices(size(posLatt))
            neigs = Neigs(r,indx,posLatt,delta)
            neigs = LinearIndices(size(posLatt))[neigs]
            indx = LinearIndices(size(posLatt))[indx]
            aux = x.linearLatt[indx]*prod(x.linearLatt[neigs])/length(neigs)
            s += aux - m2
        end
        return s/N(x)
    end

    function AutoCorrelations(M,E,cut=1/2,dt=10)
        res = []
        n = length(M)
        ncut = Int64(floor(n*cut))
        T = 1:n
        push!(res,collect(1:dt:ncut))
        for var in [M,E]
            corr = [sum([var[j]*var[j+i] for j in 1:dt:(n-i)])/(n-i+1) for i in 1:dt:ncut]
            push!(res,corr)
        end
        return res
    end

    function Partition(s,energyIntervals,temp,metaParam)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*((metaParam[1]^metaParam[2]))/2
            x=x+exp(big(s[i])-ener/temp)
        end
        return x
    end

    function DOSMagnetizationNoSave(s,energyIntervals,enerFunc,temp,metaParam;printLog=false)
        m=MagArray(energyIntervals,enerFunc,metaParam)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*(metaParam[1]^metaParam[2])/2
           x=x+m[i]*exp(big(s[i])-ener/temp)
            if printLog
                println(energyIntervals[i])
                println(m[i])
                println(s[i])
                println(ener)
                println(ener/temp)
                println(big(s[i])-ener/temp)
                println(exp(big(s[i])-ener/temp))
                println()
            end
        end
        return x/(Partition(s,energyIntervals,temp,metaParam))
    end

    function DOSMagnetization(s,energyIntervals,mag,temp,metaParam;printLog=false)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*(metaParam[1]^metaParam[2])/2
           x=x+mag[i]*exp(big(s[i])-ener/temp)
            if printLog
                println(energyIntervals[i])
                println(m[i])
                println(s[i])
                println(ener)
                println(ener/temp)
                println(big(s[i])-ener/temp)
                println(exp(big(s[i])-ener/temp))
                println()
            end
        end
        return abs(x)/(Partition(s,energyIntervals,temp,metaParam)*(metaParam[1]^metaParam[2]))
    end

    function DOSMagneticSucep(s,energyIntervals,mag,temp,metaParam;printLog=false)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*(metaParam[1]^metaParam[2])/2
           x=x+(mag[i]^2)*exp(big(s[i])-ener/temp)
            if printLog
                println(energyIntervals[i])
                println(m[i])
                println(s[i])
                println(ener)
                println(ener/temp)
                println(big(s[i])-ener/temp)
                println(exp(big(s[i])-ener/temp))
                println()
            end
        end
        x=x/(Partition(s,energyIntervals,temp,metaParam)*(metaParam[1]^metaParam[2]))
        return (    (x - (metaParam[1]^metaParam[2])*DOSMagnetization(s,energyIntervals,mag,temp,metaParam)^2))/temp
    end

    function DOSInternalEnergy(s,energyIntervals,temp,metaParam)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*(metaParam[1]^metaParam[2])/2
            x=x+ener*exp(big(s[i])-ener/temp)
        end
        return x/(Partition(s,energyIntervals,temp,metaParam)*(metaParam[1]^metaParam[2]))
    end
        
    function DOSHeatCapacity(s,energyIntervals,temp,metaParam)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*(metaParam[1]^metaParam[2])/2
            x = x + ener^2 *exp(big(s[i])-ener/temp)
        end
        x=x/(Partition(s,energyIntervals,temp,metaParam)*(metaParam[1]^metaParam[2]))
        return (x-(metaParam[1]^metaParam[2])*DOSInternalEnergy(s,energyIntervals,temp,metaParam)^2)/temp^2
    end

    
end
