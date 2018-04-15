include("auxiliar.jl")
include("graphTheory.jl")
module StatEnsemble
    using Auxiliar,GraphTheory
    #funcion modulo auxiliar. Igual que el módulo normal, pero tal que mod(n,n)=n




    #función para calcular la energía. Para mayor generalidad, esta función depende a su vez de una función de vecinos, denotada por la variable "func". Esta función nos debe de regresar la suma de los spines vecinos a un punto en un arreglo.

    function Energy(latt,param,neigLatt;printLog=false)
        e=0
        s=size(latt)
        dim=length(s)
        B=param[3]
        J=param[2]
        for pos in CartesianRange(size(latt))
            e=e-B*latt[pos]-J*latt[pos]*Auxiliar.NeighborSum(latt,neigLatt,pos)*1/2
        end
        return e
    end

    function PenalizedEnergy(latt,param,func,pen;printLog=false) 
        e=0
        s=size(latt)
        dim=length(s)
        B=param[3]
        J=param[2]
        if dim==1
            for i in 1:s[1]
                e = e - B*latt[i] - J*func(latt,[i])*latt[i]*(1/2)
            end
        elseif dim==2
            for i in 1:s[1]
                for j in 1:s[2]
                    e = e - B*latt[i,j] - J*func(latt,[i,j])*latt[i,j]*(1/2)
                end
            end
            l=GraphTheory.SearchAllCycles1(latt)
            for pos in l 
                e=e+pen*Auxiliar.GetValue(latt,pos)
            end
        elseif dim==2
            for i in 1:s[1]
                for j in 1:s[2]
                    for k in 1:s[3]
                        e = e - B*latt[i,j,k] - J*func(latt,[i,j,k])*latt[i,j,k]*(1/2)
                    end
                end
            end
        else
            error("dimension not soported")
        end
        return e
    end


    function ProbCanonical(latt,pos,param,neigLatt,temp)
        aux=param[2]*Auxiliar.NeighborSum(latt,neigLatt,pos)+param[3]
        beta=1/temp
        return exp(-2*latt[pos]*aux*beta)
    end

    function MagArray(energyIntervals,param;printLog=false)
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
                e=Energy(m,param,neigLatt)/(param[1]^2)
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

    function DOSMag(s,energyIntervals,temp,param;printLog=false)
        m=MagArray(energyIntervals,param)
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

    function WangLandauValues(s,energyIntervals,tempArray,param)
        x=[]
        y=[]
        for temp in tempArray
            push!(x,DOSEnergy(s,energyIntervals,temp,param))
            push!(y,DOSCV(s,energyIntervals,temp,param))
        end
        return (x,y)
    end
end
