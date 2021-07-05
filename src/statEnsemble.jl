
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
    function NormalEnergy(sys::IsingModel,simulParam;printLog=false)
        e=0
        B=simulParam[1]
        J=simulParam[2]
        e = -B*sum(sys.sites) - J * sum(sys.neigSum .* sys.sites ) /2
        return e
    end

    """
        PenalizedEnergy(latt,neigLatt,simulParam)

        Calculates energy of the modified Ising model for a lattice `latt` with simulation parameters `simulParam` and neighbor lattice `neigLatt` searching for cycles using the algorithm contained in "graphTheory.jl"
    """
    function PenalizedEnergy(sys::IsingModel,simulParam;printLog=false)
        e=0
        B=simulParam[1]
        J=simulParam[2]
        C=simulParam[3]
        e = -B*sum(sys.sites) - J * sum(sys.neigSum .* sys.sites ) /2 + C * EdgInCyc(sys)
        return e
    end

    """
    LogProbAta(ener,sys::LatticeGas,pos::Integer,simulParam;printLog=false)

    Idea for a new possible probability acceptance ratio that approximates all new cycles as the new C
    """
    function LogProbAta(ener,sys::LatticeGas,pos::Integer,simulParam;printLog=false)
        B=simulParam[1]
        J=simulParam[2]
        C=simulParam[3]
        x2 =ChangeSpin(sys,pos)
        newcyc = lgNewCycles(sys,x2)
        deltae = -B-J*sys.neigSum[pos] + 2*((B+J*sys.neigSum[pos])*sys.sites[pos]) + C * newcyc
        return -deltae/simulParam[4]
    end


    """
    """
    function LogProb(sys::IsingModel,pos,enerFunc,simulParam)
        e1=enerFunc(sys,simulParam)
        e2=enerFunc(ChangeSpin(sys,pos),simulParam)
        return (-e2+e1)/simulParam[4]
    end

    function LogProbOptimal(sys::SpinLattice,pos,enerFunc,simulParam)
        if enerFunc == NormalEnergy
            deltaener = 2*sys.sites[pos]*(simulParam[1] + simulParam[2]* sum(sys.sites[sys.edgList[pos]]))
        end
        return -deltaener/simulParam[4]
    end


    function LogProbOptimal(sys::LatticeGas,pos,enerFunc,simulParam)
        if enerFunc == NormalEnergy
            deltaener = (2*sys.sites[pos]-1)*(simulParam[1] + simulParam[2]* sum(sys.sites[sys.edgList[pos]]))
        end
        return -deltaener/simulParam[4]
    end


    function Magnetization(sys::IsingModel;absolute=true)
        if absolute
            a = abs(sum(sys.sites))
        else
            a = sum(sys.sites)
        end
        return a/N(sys)
    end
    
    function Magnetization(sysArr::Array{IsingModel,1};absolute=true)
        if absolute
            a = Statistics.mean([abs(sum(sys.sites)) for sys in sysArr])
        else
            a = Statistics.mean([sum(sys.sites) for sys in sysArr])
        end
        return a/N(sysArr[1])
    end

    function MagnetizationUpdate(sys::SpinLattice,pos::Integer)
        return -2*sys.sites[pos]
    end

    function MagnetizationUpdate(sys::LatticeGas,pos::Integer)
        return 1 - 2*sys.sites[pos]
    end


    function MagneticSucep(sysArr::Array{IsingModel,1},simulParam;absolute=true)
        if absolute
            a=Statistics.mean([abs(sum(sys.sites))^2 for sys in sysArr])
            b=(Statistics.mean([abs(sum(sys.sites)) for sys in sysArr]))^2
        else
            a=Statistics.mean([sum(sys.sites)^2 for sys in sysArr])
            b=(Statistics.mean([sum(sys.sites) for sys in sysArr]))^2
        end
        return 1/(simulParam[4]*(N(sysArr[1])))*(a-b)
    end

    function InternalEnergy(sys::IsingModel,enerFunc,simulParam)
        a = enerFunc(sys,simulParam)
        return a/N(sysArr[1])
    end

    function InternalEnergy(sysArr::Array{IsingModel,1},enerFunc,simulParam)
        a = Statistics.mean([enerFunc(sys,simulParam) for sys in sysArr])
        return a/N(sysArr[1])
    end


    function HeatCapacity(sysArr::Array{IsingModel,1},enerFunc,simulParam)
        a=Statistics.mean([enerFunc(sys,simulParam)^2 for sys in sysArr])
        b=(Statistics.mean([enerFunc(sys,simulParam) for sys in sysArr]))^2
        return 1/((simulParam[4]^2)*(N(sysArr[1])))*(a-b)
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


    function PairCorrelation(sys::IsingModel,r,posLatt,delta)
        m2 = (sum(sys.sites)/N(sys))^2
        s = 0
        for indx in CartesianIndices(size(posLatt))
            neigs = Neigs(r,indx,posLatt,delta)
            neigs = LinearIndices(size(posLatt))[neigs]
            indx = LinearIndices(size(posLatt))[indx]
            aux = sys.sites[indx]*prod(sys.sites[neigs])/length(neigs)
            s += aux - m2
        end
        return s/N(sys)
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
        sum=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*((metaParam[1]^metaParam[2]))/2
            sum=sum+exp(big(s[i])-ener/temp)
        end
        return sum
    end

    function DOSMagnetizationNoSave(s,energyIntervals,enerFunc,temp,metaParam;printLog=false)
        m=MagArray(energyIntervals,enerFunc,metaParam)
        sum=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*(metaParam[1]^metaParam[2])/2
           sum=sum+m[i]*exp(big(s[i])-ener/temp)
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
        return sum/(Partition(s,energyIntervals,temp,metaParam))
    end

    function DOSMagnetization(s,energyIntervals,mag,temp,metaParam;printLog=false)
        sum=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*(metaParam[1]^metaParam[2])/2
            sum=sum+mag[i]*exp(big(s[i])-ener/temp)
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
        return abs(sum)/(Partition(s,energyIntervals,temp,metaParam)*(metaParam[1]^metaParam[2]))
    end

    function DOSMagneticSucep(s,energyIntervals,mag,temp,metaParam;printLog=false)
        sum=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*(metaParam[1]^metaParam[2])/2
           sum=sum+(mag[i]^2)*exp(big(s[i])-ener/temp)
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
        sum=sum/(Partition(s,energyIntervals,temp,metaParam)*(metaParam[1]^metaParam[2]))
        return (    (sum - (metaParam[1]^metaParam[2])*DOSMagnetization(s,energyIntervals,mag,temp,metaParam)^2))/temp
    end

    function DOSInternalEnergy(s,energyIntervals,temp,metaParam)
        sum=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*(metaParam[1]^metaParam[2])/2
            sum=sum+ener*exp(big(s[i])-ener/temp)
        end
        return sum/(Partition(s,energyIntervals,temp,metaParam)*(metaParam[1]^metaParam[2]))
    end
        
    function DOSHeatCapacity(s,energyIntervals,temp,metaParam)
        sum=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*(metaParam[1]^metaParam[2])/2
            sum = sum + ener^2 *exp(big(s[i])-ener/temp)
        end
        sum=sum/(Partition(s,energyIntervals,temp,metaParam)*(metaParam[1]^metaParam[2]))
        return (sum-(metaParam[1]^metaParam[2])*DOSInternalEnergy(s,energyIntervals,temp,metaParam)^2)/temp^2
    end

    
end
