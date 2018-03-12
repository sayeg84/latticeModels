module StatEnsemble
export Modl, SquareLatticeNeighbors, Energy, ProbCanonical


    #funcion modulo auxiliar. Igual que el módulo normal, pero tal que mod(n,n)=n
    function Modl(a::Int64,b::Int64)
        x=mod(a,b)
        if x==0
            return b
        else
            return x
        end
    end


    #función de vecinos inmediatos utilizada para la función de peso del algoritmo metrópolis
    function SquareLatticeNeighbors(latt,pos;test=false)
        s=size(latt)
        dim=length(s)
        if dim==1
            v1=latt[Modl(pos[1]-1,s[1])]
            v2=latt[Modl(pos[1]-2,s[1])]
            if test
                println([v1 0 v2])
            end
            return v1+v2
        elseif dim==2
            v1 = latt[Modl(pos[1]-1,s[1]), pos[2]]
            v2 = latt[Modl(pos[1]+1,s[1]), pos[2]]
            v3 = latt[pos[1],Modl(pos[2]-1,s[2])]
            v4 = latt[pos[1],Modl(pos[2]+1,s[2])]
            if test
                println([0 v1 0; v3 0 v4 ; 0 v2 0])
            end
            return v1+v2+v3+v4
        elseif dim==3
            v1 = latt[Modl(pos[1]-1,s[1]), pos[2], pos[3]]
            v2 = latt[Modl(pos[1]+1,s[1]), pos[2], pos[3]]
            v3 = latt[pos[1], Modl(pos[2]-1,s[2]), pos[3]]
            v4 = latt[pos[1], Modl(pos[2]+1,s[2]), pos[3]]
            v5 = latt[pos[1], pos[2], Modl(pos[3]-1,s[3])]
            v6 = latt[pos[1], pos[2], Modl(pos[3]+1,s[3])]
            if test
                println([0 v1 0; v3 0 v4 ; 0 v2 0])
                println([v5 0 v6])
            end
            return v1+v2+v3+v4+v5+v6
        else
            error("dimension not soported")
        end
    end


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

    #=
    Probabilidad de cambio para el algoritmo Metrópolis.

    Definida como el cociente de la densidad de probabilidad de estados
    del espacio muestra entre dos estados deseados. Para mayor información,
    consultar la referencia bibliográfica
    =#

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

    function Partition(s,energyIntervals,temp,param)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])/2
            x=x+exp(big(s[i])-ener*param[1]^2/temp)
        end
        return x
    end

    function DOSEnergy(s,energyIntervals,temp,param)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])/2
            x=x+ener*exp(big(s[i])-ener*param[1]^2/temp)
        end
        return x/(Partition(s,energyIntervals,temp,param)*param[1])
    end
        
    function DOSCV(s,energyIntervals,temp,param)
        x=0
        for i in 1:length(s)
            ener=(energyIntervals[i]+energyIntervals[i+1])/2
            x = x + ener^2 *exp(big(s[i])-ener*param[1]^2/temp)
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
println(StatEnsemble.Energy(x,[10,1,0,100,10],StatEnsemble.SquareLatticeNeighbors))
=#