module StatEnsemble
    #using Auxiliar,GraphTheory
    include("auxiliar.jl")
    include("graphTheory.jl")
    include("cyclesAta.jl")
    function Energy(latt,param,neigLatt;printLog=false)
        e=0
        s=size(latt)
        dim=length(s)
        B=param[3]
        J=param[2]
        for pos in CartesianIndices(size(latt))
            e=e-B*latt[pos]-J*latt[pos]*Auxiliar.NeighborSum(latt,neigLatt,pos)*1/2
        end
        return e
    end

    function PenalizedEnergy(latt,param,neigLatt;printLog=false)
        e=0
        s=size(latt)
        dim=length(s)
        B=param[3]
        J=param[2]
        for pos in CartesianIndices(size(latt))
            e=e-(B+J*Auxiliar.NeighborSum(latt,neigLatt,pos)*1/2)*latt[pos]
        end
        cycles=GraphTheory.Cycles(latt,neigLatt,printLog=printLog)
        for pos in cycles
            e=e+param[7]*latt[pos]
        end
        return e
    end

    function PenalizedEnergy2(latt,param,neigLatt;printLog=false)
        e=0
        s=size(latt)
        dim=length(s)
        B=param[3]
        J=param[2]
        for pos in CartesianIndices(size(latt))
            e=e-(B+J*Auxiliar.NeighborSum(latt,neigLatt,pos)*1/2)*latt[pos]
        end
        #println("hay energia penalizada")
        cycles=CyclesAta.ciclos2(latt)
        if sum(cycles) > 0
            println("bien")
        end
        e=e+param[7]*sum(cycles)
        return e
    end


    function ProbCanonical(latt,pos,param,enerFunc,neigLatt,temp)
        e1=enerFunc(latt,param,neigLatt)
        e2=enerFunc(Auxiliar.ChangeSpin(latt,pos),param,neigLatt)
        return exp((-e2+e1)/temp)
    end

    function MagArray(energyIntervals,enerFunc,param;printLog=false)
        m=ones(param[1],param[1])
        neigLatt=Auxiliar.NeighborIndexLattice(m,Auxiliar.SquareLatticeNeighborsIndex)
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

    function Partition(s,energyIntervals,temp,param)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*param[1]^2/2
            x=x+exp(big(s[i])-ener/temp)
        end
        return x
    end

    function DOSMag(s,energyIntervals,mag,enerFunc,temp,param;printLog=false)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*param[1]^2/2
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
        return abs(x)/(Partition(s,energyIntervals,temp,param))
    end

    function DOSMagNoSave(s,energyIntervals,enerFunc,temp,param;printLog=false)
        m=MagArray(energyIntervals,enerFunc,param)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*param[1]^2/2
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
        return x/(Partition(s,energyIntervals,temp,param))
    end

    function DOSEnergy(s,energyIntervals,temp,param)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*param[1]^2/2
            x=x+ener*exp(big(s[i])-ener/temp)
        end
        return x/(Partition(s,energyIntervals,temp,param))
    end
        
    function DOSCV(s,energyIntervals,temp,param)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*param[1]^2/2
            x = x + ener^2 *exp(big(s[i])-ener/temp)
        end
        x=x/Partition(s,energyIntervals,temp,param)
        return (x-DOSEnergy(s,energyIntervals,temp,param)^2)/temp^2
    end

end
