include("statEnsemble.jl")
include("auxiliar.jl")
module Algorithms
    using StatEnsemble, Auxiliar
    #=
    Implementación del algoritmo metrópolis para crear N matrices de Nx*Ny 
    usando el algoritmo Metrópolis. Suponemos los valores de spin son 1 y -1.
    =#
    function Metropolis(temp::Float64,param,initLatt;test=false)
        matArr=[]
        mat=initLatt
        push!(matArr,copy(mat))
        for i in 1:(param[4]+1)
            pos=rand(1:param[1],length(size(mat)))
            res=StatEnsemble.ProbCanonical(mat,
                                        pos,
                                        param,
                                        Auxiliar.SquareLatticeNeighbors,
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
                Auxiliar.ChangeSpin!(mat,pos)
            end
            if mod(i,param[5])==1
                push!(matArr,copy(mat))
            end
        end
        return matArr
    end

    # =============================================
    # Wang-Landau beggining
    function isFlat(hist;per=0.5)
        min=minimum(hist)
        max=maximum(hist)
        avg=mean(hist)
        if min >= avg*per
            return true
        else
            return false
        end
    end
    function WangLandau(param,initLatt;test=false)
        last=[]
        hist=zeros(param[6])
        s=zeros(param[6])
        energyIntervals=collect(linspace(-2,2,param[6]+1))
        modfact=1
        latt=copy(initLatt)
        n=0
        while (modfact>=1e-5)
            pos=rand(1:param[1],length(size(latt)))
            energyBefore=StatEnsemble.Energy(latt,param,Auxiliar.SquareLatticeNeighbors)
            energyAfter=energyBefore+2*Auxiliar.GetValue(latt,pos)*(param[2]*Auxiliar.SquareLatticeNeighbors(latt,pos)+param[3])
            energyBefore=energyBefore/(param[1])^2
            energyAfter=energyAfter/(param[1])^2
            p1=Auxiliar.SearchSortedMod(energyIntervals,energyBefore)
            p2=Auxiliar.SearchSortedMod(energyIntervals,energyAfter)
            tes=rand()
            last=[]
            if test
                println("latt")
                println(latt)
                println("pos")
                println(pos)
                println("hist")
                println(hist)
                println("entropy")
                println(s)
                
                println("before")
                println(energyBefore)
                println("bin")
                println(p1)
                println("after")
                println(energyAfter)
                println("bin")
                println(p2)
                println("probab")
                println(exp(s[p1]-s[p2]))
                println("rand")
                println(tes)
                println("avg")
                println(mean(hist))
                println("min")
                println(minimum(hist))
            end
            a=exp(s[p1])
            b=exp(s[p2])
            η=exp(big(s[p1]-s[p2]))
            if tes < η
                Auxiliar.ChangeSpin!(latt,pos) 
            end
            s[p1]=s[p1]+modfact
            hist[p1]=hist[p1]+1
            if isFlat(hist)
                modfact=modfact*1/2
                last=copy(hist)
                hist=zeros(param[6])
                println(n)
                println("change")
                println(modfact)
            end
            if mod(n,10^5) == 0 && ~test
                println(n/10^5)
            end
            if n==10^9
                println("Exceded tolerance")
                return (energyIntervals,s,last)
            end
            n=n+1
        end
        print("Iterations: ")
        println(n)
        return (energyIntervals,s,last)
    end
    function WangLandauSimple(param,initLatt;test=false)
        last=[]
        hist=zeros(param[6])
        s=zeros(param[6])
        energyIntervals=collect(linspace(-2,0,param[6]+1))
        modfact=1
        latt=copy(initLatt)
        n=0
        while (modfact>=1e-5)
            pos=rand(1:param[1],length(size(latt)))
            energyBefore=StatEnsemble.Energy(latt,param,Auxiliar.SquareLatticeNeighbors)
            energyAfter=energyBefore+2*Auxiliar.GetValue(latt,pos)*(param[2]*Auxiliar.SquareLatticeNeighbors(latt,pos)+param[3])
            energyBefore=energyBefore/((param[1])^2*param[2])
            energyAfter=energyAfter/((param[1])^2*param[2])
            p1=Auxiliar.SearchSortedMod(energyIntervals,energyBefore)
            p2=Auxiliar.SearchSortedMod(energyIntervals,energyAfter)
            tes=rand()
            last=[]
            if test
                println("latt")
                println(latt)
                println("pos")
                println(pos)
                println("hist")
                println(hist)
                println("entropy")
                println(s)
                
                println("before")
                println(energyBefore)
                println("bin")
                println(p1)
                println("after")
                println(energyAfter)
                println("bin")
                println(p2)
                println("probab")
                println(exp(s[p1]-s[p2]))
                println("rand")
                println(tes)
                println("avg")
                println(mean(hist))
                println("min")
                println(minimum(hist))
            end
            a=exp(s[p1])
            b=exp(s[p2])
            η=exp(big(s[p1]-s[p2]))
            if  tes < η && energyAfter <=0
                Auxiliar.ChangeSpin!(latt,pos) 
            end
            s[p1]=s[p1]+modfact
            hist[p1]=hist[p1]+1
            if isFlat(hist)
                modfact=modfact*1/2
                last=copy(hist)
                hist=zeros(param[6])
                println("change")
                println()
                print("iter: ")
                println(n)
                print("modfact: ")
                println(modfact)
                println()
            end
            if mod(n,10^5) == 0 && ~test
                print("iter: ")
                println(n)
                println()
            end
            if n==10^10
                println("Exceded tolerance")
                append!(energyIntervals,[abs(energyIntervals[i]) for i in (param[6]):-1:1])
                Auxiliar.MirrorList!(s)
                Auxiliar.MirrorList!(last)
                return (energyIntervals,s,last)
            end
            n=n+1
        end
        print("Iterations: ")
        println(n)
        append!(energyIntervals,[abs(energyIntervals[i]) for i in (param[6]):-1:1])
        Auxiliar.MirrorList!(s)
        Auxiliar.MirrorList!(last)
        return (energyIntervals,s,last)
    end
end