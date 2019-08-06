
include("cycles.jl")
module StatEnsemble
    import ..SpinLattice
    import ..LatticeGas
    import ..IsingModel
    import ..ChangeSpin
    import ..N
    import ..Cycles.lgEdgesInCycles
    import ..Cycles.ciclos2


    using Statistics
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
        e = -B*sum(x.linearLatt) - J * sum(x.linearNeigSumLatt .* x.linearLatt ) /2 + C * lgEdgesInCycles(x)
        return e
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

    function Magnetization(X::Array{IsingModel,1},absolute=true)
        if absolute
            a = Statistics.mean([abs(sum(x.linearLatt)) for x in X])
        else
            a = Statistics.mean([sum(x.linearLatt) for x in X])
        end
        return a/N(X[1])
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


    function InternalEnergy(X::Array{IsingModel,1},enerFunc,simulParam)
        a = Statistics.mean([enerFunc(x,simulParam) for x in X])
        return a/N(X[1])
    end


    function HeatCapacity(X::Array{IsingModel,1},enerFunc,simulParam)
        a=Statistics.mean([enerFunc(x,simulParam)^2 for x in X])
        b=(Statistics.mean([enerFunc(x,simulParam) for x in X]))^2
        return 1/((simulParam[4]^2)*(N(X[1])))*(a-b)
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

    #=

    ########### Wang-Landau

    function MagArray(energyIntervals,enerFunc,param;printLog=false)
        m=ones(param[1],param[1])
        neigLatt=Geometry.IndexLattice(m,Geometry.SquareLatticeNeighbors)
        aux=[]
        for i in 1:(length(energyIntervals)-1)
            push!(aux,[])
        end

        for i in 1:param[1]
            k=Auxiliar.Modl(i,2)
            for j in k:2:param[1]
                mag=abs(sum(m))
                e=enerFunc(m,param,neigLatt)/(param[1]^2)
                pos=Auxiliar.SearchSortedMod(energyIntervals,e)
                push!(aux[pos],mag)
                if printLog
                    println(m)
                    println(mag)
                    println(e)
                    println(pos)
                end
                Auxiliar.ChangeSpin!(m,[i,j])
            end
        end
        fin=[Auxiliar.MeanMod(aux[i]) for i in 1:length(aux)]
        return fin
    end




    
=#