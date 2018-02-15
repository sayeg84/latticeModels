include("statEnsemble.jl")
module Algorithms
    using StatEnsemble
    #=
    Implementación del algoritmo metrópolis para crear N matrices de Nx*Ny 
    usando el algoritmo Metrópolis. Suponemos los valores de spin son 1 y -1.
    =#
    function ChangeSpin!(latt,pos)
        dim=length(size(latt))
        if dim==1
            latt[pos[1]]=-1*latt[pos[1]]
        elseif dim==2
            latt[pos[1],pos[2]]=-1*latt[pos[1],pos[2]]
        elseif dim==3
            latt[pos[1],pos[2],pos[3]]=-1*latt[pos[1],pos[2],pos[3]]
        end
    end
    function Metropolis(temp,param,initLatt;test=false)
        matArr=[]
        mat=initLatt
        push!(matArr,copy(mat))
        for i in 1:(param[4]+1)
            pos=rand(1:param[1],length(size(mat)))
            res=StatEnsemble.ProbCanonical(mat,
                                        pos,
                                        param,
                                        StatEnsemble.SquareLatticeNeighbors,
                                        temp)
            t=rand()
            if test
                print("iter = ")
                println(i-1)
                println("lattice")
                println(mat)
                print("Postion = ")
                println(pos)
                print("Probability = ")
                println(res)
                print("Random number = ")
                println(t)
                println()
                println()
            end
            if res>t
                ChangeSpin!(mat,pos)
            end
            if mod(i,param[5])==1
                push!(matArr,copy(mat))
            end
        end
        return matArr
    end

    # =============================================
    # Wang-Landau beggining
    function isFlat(hist)
        min=minimum(hist)
        max=maximum(hist)
        avg=mean(hist)
        if (max-avg)/avg > 0.2 || (avg-min)/avg > 0.2
            return false
        else
            return true
        end
    end
    function WangLandau(param,func,intervals::Int64;test=false)
        hist=zeros(intervals)
        dos=zeros(intervals)
        energyIntervals=linspace(-2,0,intervals+1)
        modfact=exp(1)
        for i in 1:param[4]+1
            pos=rand(1:param[1],1:param[1])
        end
    end
end