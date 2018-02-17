include("statEnsemble.jl")
module Algorithms
    using StatEnsemble
    #=
    Implementación del algoritmo metrópolis para crear N matrices de Nx*Ny 
    usando el algoritmo Metrópolis. Suponemos los valores de spin son 1 y -1.
    =#
    function SetValue!(mat,pos,val)
        dim=length(size(mat))
        if dim!=length(pos)
            error("Dimensions must match")
        elseif eltype(mat)!=typeof(val) && eltype(mat)==Any
            error("Types must match")
        elseif dim==1
            mat[pos[1]]=val
        elseif dim==2
            mat[pos[1],pos[2]]=val
        elseif dim==3
            mat[pos[1],pos[2],pos[3]]=val
        else
            error("Dimension not supported")
        end

    end
    function GetValue(mat,pos)
        dim=length(size(mat))
        if dim!=length(pos)
            error("Dimensions must match")
        elseif dim==1
            return mat[pos[1]]
        elseif dim==2
            return mat[pos[1],pos[2]]
        elseif dim==3
            return mat[pos[1],pos[2],pos[3]]
        else
            error("Dimension not supported")
        end
    end
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
    function Metropolis(temp::Float64,param,initLatt;test=false)
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
        if min >= avg*0.5
            return true
        else
            return false
        end
    end
    function WangLandau(param,initLatt,intervals::Int64;test=false)
        hist=zeros(intervals)
        dos=ones(intervals)
        energyIntervals=linspace(-2,2,intervals+1)
        modfact=exp(1)
        latt=copy(initLatt)
        n=0
        while (modfact>=1e-5+1)
            pos=rand(1:param[1],length(size(latt)))
            energyBefore=StatEnsemble.Energy(latt,param,StatEnsemble.SquareLatticeNeighbors)
            energyAfter=energyBefore+2*GetValue(latt,pos)*(param[2]*StatEnsemble.SquareLatticeNeighbors(latt,pos)+param[3])
            energyBefore=energyBefore/(param[1])^2
            energyAfter=energyAfter/(param[1])^2
            p1=searchsortedlast(energyIntervals,energyBefore)
            p2=searchsortedlast(energyIntervals,energyAfter)
            η=dos[p1]/dos[p2]
            tes=rand()
            if test
                println("latt")
                println(latt)
                println("pos")
                println(pos)
                println("hist")
                println(hist)
                println("DOS")
                println(dos)
                
                println("before")
                println(energyBefore)
                println("after")
                println(energyAfter)
                println("probab")
                println(η)
                println("rand")
                println(tes)
                println("avg")
                println(mean(hist))
                println("min")
                println(minimum(hist))
            end
            if tes < η
                ChangeSpin!(latt,pos) 
            end
            dos[p1]=dos[p1]*modfact
            hist[p1]=hist[p1]+1
            if isFlat(hist)
                modfact=sqrt(modfact)
                println(n)
                println("change")
                println(modfact)
            end
            if n==10^4
                println("exceded tolerance")
                break
            end
            n=n+1
        end
        print("Iterations: ")
        println(n)
        return (energyIntervals,dos,hist)
    end
end