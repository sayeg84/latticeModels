include("auxiliar.jl")
module StatEnsemble
export Energy, ProbCanonical
    using Auxiliar
    #funcion modulo auxiliar. Igual que el módulo normal, pero tal que mod(n,n)=n




    #función para calcular la energía. Para mayor generalidad, esta función depende a su vez de una función de vecinos, denotada por la variable "func". Esta función nos debe de regresar la suma de los spines vecinos a un punto en un arreglo.
    function Energy(latt,param,func;test=false)
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


    function ProbCanonical(latt,pos,param,func,temp)
        aux=param[2]*func(latt,pos)+param[3]
        dim=length(size(latt))
        beta=1/temp
        if dim==1
            return exp(-2*latt[pos[1]]*aux*beta)
        elseif dim==2
            return exp(-2*latt[pos[1],pos[2]]*aux*beta)
        elseif dim==3
            return exp(-2*latt[pos[1],pos[2],pos[3]]*aux*beta)
        else
            error("dimension not supported")
        end
    end

    function MagArray(energyIntervals,param;test=false)
        m=ones(param[1],param[1])
        aux=[]
        for i in 1:(length(energyIntervals)-1)
            push!(aux,[])
        end

        for i in 1:param[1]
            k=Auxiliar.Modl(i,2)
            for j in k:2:param[1]
                mag=abs(sum(m))
                e=Energy(m,param,Auxiliar.SquareLatticeNeighbors)/(param[1]^2)
                pos=Auxiliar.SearchSortedMod(energyIntervals,e)
                push!(aux[pos],mag)
                if test
                    println(m)
                    println(mag)
                    println(e)
                    println(pos)
                end
                Auxiliar.ChangeSpin!(m,[i,j])
            end
        end
        fin=[mean(aux[i]) for i in 1:length(aux)]
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

    function DOSMag(s,energyIntervals,temp,param)
        m=MagArray(energyIntervals,param)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])*param[1]^2/2
            x=x+m[i]*exp(big(s[i])-ener/temp)
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
#=
x=[1 -1 -1 1 -1; -1 1 -1 -1 1; 1 1 1 1 -1; -1 1 -1 -1 1; -1 1 1 -1 1]
println(StatEnsemble.Energy(x,[10,1,0,100,10],StatEnsemble.Auxiliar.SquareLatticeNeighbors))
=#