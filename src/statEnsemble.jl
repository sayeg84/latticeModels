module StatEnsemble
    #using Auxiliar,GraphTheory
    include("auxiliar.jl")
    include("geometry.jl")
    include("graphTheory.jl")
    include("cyclesAta.jl")
    using Statistics
    """
        NormalEnergy(latt,neigLatt,simulParam)

        Calculates normal energy of Ising model for a lattice `latt` with simulation parameters `simulParam` and neighbor lattice `neigLatt`
    """
    function NormalEnergy(latt,neigLatt,simulParam;printLog=false)
        e=0
        s=size(latt)
        B=simulParam[1]
        J=simulParam[2]
        for pos in CartesianIndices(s)
            e=e-B*latt[pos]-J*latt[pos]*Geometry.NeighborSum(latt,neigLatt,pos)*1/2
        end
        return e
    end

    """
        PenalizedEnergy(latt,neigLatt,simulParam)

        Calculates energy of the modified Ising model for a lattice `latt` with simulation parameters `simulParam` and neighbor lattice `neigLatt` searching for cycles using the algorithm contained in "graphTheory.jl"
    """
    function PenalizedEnergy(latt,neigLatt,simulParam;printLog=false)
        e=0
        dim=length(s)
        B=simulParam[1]
        J=simulParam[2]
        for pos in CartesianIndices(s)
            e=e-(B+J*Geometry.NeighborSum(latt,neigLatt,pos)*1/2)*latt[pos]
        end
        cycles=GraphTheory.Cycles(latt,neigLatt,printLog=printLog)
        for pos in cycles
            e=e+simulParam[3]*latt[pos]
        end
        return e
    end

    """
        PenalizedEnergy2(latt,neigLatt,simulParam)

        Calculates energy of the modified Ising model for a lattice `latt` with simulation parameters `simulParam` and neighbor lattice `neigLatt` searching for cycles using the algorithm contained in "cyclesAta.jl"
    """
    function PenalizedEnergy2(latt,neigLatt,simulParam;printLog=false)
        e=0
        s=size(latt)
        dim=length(s)
        B=simulParam[1]
        J=simulParam[2]
        for pos in CartesianIndices(size(latt))
            e=e-(B+J*Geometry.NeighborSum(latt,neigLatt,pos)*1/2)*latt[pos]
        end
        #println("hay energia penalizada")
        cycles=CyclesAta.ciclos2(latt)
        e=e+simulParam[3]*sum(cycles)
        return e
    end

    """
        Prob(latt,neigLatt,pos,simulParam,enerFunc)

        Calculates the probability of changing spin in position `pos` in an Ising model with energy function `enerFunc`for a lattice `latt` with params `param` and neighbor lattice `neigLatt` searching for cycles using the algorithm contained in "cyclesAta.jl"
    """
    function Prob(latt,neigLatt,pos,simulParam,enerFunc,nrml)
        e1=enerFunc(latt,neigLatt,simulParam)
        e2=enerFunc(Auxiliar.ChangeSpin(latt,pos,normal=nrml),neigLatt,simulParam)
        return exp((-e2+e1)/simulParam[4])
    end


    function Magnetization(array,geoParam;ab=true)
        if ab
            a=Statistics.mean([abs(sum(x)) for x in array])
        else 
            a=Statistics.mean([sum(x) for x in array])
        end
        return a/(geoParam[1]^geoParam[2])
    end

    function InternalEnergy(array,model,geoParam,simulParam)
        neigLatt=BuildLattices(geoParam,model)
        if model=="normal"
            enerFunc=NormalEnergy
        else
            enerFunc=PenalizedEnergy2
        end
        a=Statistics.mean([enerFunc(x,neigLatt,simulParam) for x in array])
        return a/(geoParam[1]^geoParam[2])
    end

    function HeatCapacity(array,model,geoParam,simulParam)
        neigLatt=BuildLattices(geoParam,model)
        if model=="normal"
            enerFunc=NormalEnergy
        else
            enerFunc=PenalizedEnergy2
        end
        a=Statistics.mean([enerFunc(x,neigLatt,simulParam)^2 for x in array])
        b=(Statistics.mean([enerFunc(x,neigLatt,simulParam) for x in array]))^2
        return 1/((simulParam[4]^2)*(geoParam[1]^geoParam[2]))*(a-b)
    end

    function MagneticSucep(array,geoParam,simulParam;ab=true)
        if ab
            a=Statistics.mean([sum(x)^2 for x in array])
            b=(Statistics.mean([sum(x) for x in array]))^2
        else
            a=Statistics.mean([abs(sum(x))^2 for x in array])
            b=(Statistics.mean([abs(sum(x)) for x in array]))^2
        end
        return 1/(simulParam[4]*(geoParam[1]^geoParam[2]))*(a-b)
    end




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

    function Partition(s,energyIntervals,temp,geoParam)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*((geoParam[1]^geoParam[2]))/2
            x=x+exp(big(s[i])-ener/temp)
        end
        return x
    end

    function DOSMag(s,energyIntervals,mag,enerFunc,temp,geoParam;printLog=false)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*(geoParam[1]^geoParam[2])/2
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
        return abs(x)/(Partition(s,energyIntervals,temp,geoParam))
    end

    function DOSMagNoSave(s,energyIntervals,enerFunc,temp,geoParam;printLog=false)
        m=MagArray(energyIntervals,enerFunc,geoParam)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*(geoParam[1]^geoParam[2])/2
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
        return x/(Partition(s,energyIntervals,temp,geoParam))
    end

    function DOSEnergy(s,energyIntervals,temp,geoParam)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*(geoParam[1]^geoParam[2])/2
            x=x+ener*exp(big(s[i])-ener/temp)
        end
        return x/(Partition(s,energyIntervals,temp,geoParam))
    end
        
    function DOSCV(s,energyIntervals,temp,geoParam)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*(geoParam[1]^geoParam[2])/2
            x = x + ener^2 *exp(big(s[i])-ener/temp)
        end
        x=x/Partition(s,energyIntervals,temp,geoParam)
        return (x-DOSEnergy(s,energyIntervals,temp,geoParam)^2)/temp^2
    end

end
