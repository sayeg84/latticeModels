include("statEnsemble.jl")
include("auxiliar.jl")
module Algorithms
    using StatEnsemble, Auxiliar
    #=
    Metropolis algorithm implementation
    =#
    function Metropolis(temp::Float64,param,enerFunc,initLatt,neigLatt;printLog=false)
        matArr=[]
        mat=initLatt
        push!(matArr,copy(mat))
        for i in 1:(param[4]+1)
            pos=Auxiliar.RandomPosition(mat)
            res=StatEnsemble.ProbCanonical(mat,pos,param,enerFunc,neigLatt,temp)
            t=rand()
            if printLog
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
                mat[pos]=-1*mat[pos]
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
    function isFlatCount(hist;n=300)
        if maximum(hist)>=n
            return true
        else
            return false
        end
    end
    function WangLandau(param,initLatt,neigLatt;printLog=false)
        last=[]
        globalHist=zeros(param[6])
        mag=zeros(param[6])
        hist=zeros(param[6])
        s=zeros(param[6])
        energyIntervals=collect(linspace(-2,2,param[6]+1))
        modfact=1
        latt=copy(initLatt)
        n=0
        while (modfact>=1e-5)
            pos=Auxiliar.RandomPosition(latt)
            energyBefore=StatEnsemble.Energy(latt,param,neigLatt)
            energyAfter=energyBefore+2*latt[pos]*(param[2]*Auxiliar.NeighborSum(latt,neigLatt,pos)+param[3])
            energyBefore=energyBefore/(param[1])^2
            energyAfter=energyAfter/(param[1])^2
            p1=Auxiliar.SearchSortedMod(energyIntervals,energyBefore)
            p2=Auxiliar.SearchSortedMod(energyIntervals,energyAfter)
            tes=rand()
            last=[]
            if printLog
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
            globalHist[p1]+=1
            mag[p1]=(mag[p1]*(globalHist[p1]-1)+abs(sum(latt))/(param[1]^2))*1/(globalHist[p1])
            η=exp(big(s[p1]-s[p2]))
            if tes < η
                latt[pos]=-1*latt[pos]
            end
            globalHist[p1]+=1
            s[p1]+=modfact
            hist[p1]+=1
            if isFlat(hist)
                modfact=modfact*1/2
                last=copy(hist)
                hist=zeros(param[6])
                println(n)
                println("change")
                println(modfact)
            end
            if mod(n,10^5) == 0 && ~printLog
                print("iter: ")
                println(n)
                println()
            end
            if n==10^9
                println("Exceded tolerance")
                s=[s[i]-s[1]+log(2) for i in 1:length(s)]
                mag=mag*(param[1]^2)
                return (energyIntervals,s,mag,last,n)
            end
            n=n+1
        end
        print("Iterations: ")
        println(n)
        s=[s[i]-s[1]+log(2) for i in 1:length(s)]
        mag=mag*(param[1]^2)
        return (energyIntervals,s,mag,last,n)
    end

    function WangLandauSimple(param,initLatt,neigLatt;printLog=false)
        last=[]
        globalHist=zeros(param[6])
        mag=zeros(param[6])
        hist=zeros(param[6])
        s=zeros(param[6])
        energyIntervals=collect(linspace(-2,0,param[6]+1))
        modfact=1
        latt=copy(initLatt)
        n=0
        while (modfact>=1e-5)
            pos=Auxiliar.RandomPosition(latt)
            energyBefore=StatEnsemble.Energy(latt,param,neigLatt)
            energyAfter=energyBefore+2*latt[pos]*(param[2]*Auxiliar.NeighborSum(latt,neigLatt,pos)+param[3])
            energyBefore=energyBefore/(param[1])^2
            energyAfter=energyAfter/(param[1])^2
            p1=Auxiliar.SearchSortedMod(energyIntervals,energyBefore)
            p2=Auxiliar.SearchSortedMod(energyIntervals,energyAfter)
            tes=rand()
            last=[]
            if printLog
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
            globalHist[p1]+=1
            mag[p1]=(mag[p1]*(globalHist[p1]-1)+abs(sum(latt))/(param[1]^2))*1/(globalHist[p1])
            η=exp(big(s[p1]-s[p2]))
            if  tes < η && energyAfter <=0
                latt[pos]=-1*latt[pos]
            end
            s[p1]+=modfact
            hist[p1]+=1
            if isFlat(hist)
                modfact=modfact*1/2
                last=copy(hist)
                hist=zeros(param[6])
                println(n)
                println("change")
                println(modfact)
            end
            if mod(n,10^5) == 0 && ~printLog
                print("iter: ")
                println(n)
                println()
            end
            if n==10^9
                println("Exceded tolerance")
                append!(energyIntervals,[abs(energyIntervals[i]) for i in (param[6]):-1:1])
                Auxiliar.MirrorList!(s)
                #normalization
                s=[s[i]-s[1]+log(2) for i in 1:length(s)]
                Auxiliar.MirrorList!(mag)
                Auxiliar.MirrorList!(last)
                mag=mag*(param[1]^2)
                return (energyIntervals,s,mag,last,n)
            end
            n=n+1
        end
        print("Iterations: ")
        println(n)
        append!(energyIntervals,[abs(energyIntervals[i]) for i in (param[6]):-1:1])
        Auxiliar.MirrorList!(s)
        #normalization
        s=[s[i]-s[1]+log(2) for i in 1:length(s)]
        Auxiliar.MirrorList!(mag)
        Auxiliar.MirrorList!(last)
        mag=mag*(param[1]^2)
        return (energyIntervals,s,mag,last,n)
    end

    function WangLandauCycle(param,initLatt,neigLatt;printLog=false)
        last=[]
        globalHist=zeros(param[6])
        mag=zeros(param[6])
        hist=zeros(param[6])
        s=zeros(param[6])
        energyIntervals=collect(linspace(-2,2,param[6]+1))
        modfact=1
        latt=copy(initLatt)
        n=0
        while (modfact>=1e-5)
            pos=Auxiliar.RandomPosition(latt)
            if printLog
                println(latt)
            end
            try energyBefore=StatEnsemble.PenalizedEnergy(latt,param,neigLatt,printLog=printLog)
            catch LoadError
                #println("antes")
                #println(latt)
                #println(energyBefore)
                #error("Wrong lattice")
            end
            try energyAfter=StatEnsemble.PenalizedEnergy(Auxiliar.ChangeSpin(latt,pos),param,neigLatt,printLog=printLog)
            catch LoadError 
                #println("después")
                #println(Auxiliar.ChangeSpin(latt,pos))
                #println(energyAfter)
                #error("Wrong lattice")
            end
            energyBefore=energyBefore/(param[1])^2
            energyAfter=energyAfter/(param[1])^2
            p1=Auxiliar.SearchSortedMod(energyIntervals,energyBefore)
            p2=Auxiliar.SearchSortedMod(energyIntervals,energyAfter)
            tes=rand()
            last=[]
            if printLog
                println("iter")
                println(n)
                println("latt")
                println(latt)
                println("pos")
                println(pos)
                println("neigs")
                println(neigLatt[pos])
                println()
                println()
                println(neigLatt)
                println()
                println()
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
            globalHist[p1]+=1
            mag[p1]=(mag[p1]*(globalHist[p1]-1)+abs(sum(latt))/(param[1]^2))*1/(globalHist[p1])
            η=exp(big(s[p1]-s[p2]))
            if  tes < η
                latt[pos]=-1*latt[pos]
            end
            s[p1]+=modfact
            hist[p1]+=1
            if isFlat(hist)
                modfact=modfact*1/2
                last=copy(hist)
                hist=zeros(param[6])
                println(n)
                println("change")
                println(modfact)
            end
            if mod(n,10^5) == 0 && ~printLog
                print("iter: ")
                println(n)
                println()
            end
            if n==10^9
                println("Exceded tolerance")
                #normalization
                s=[s[i]-s[1] for i in 1:length(s)]
                mag=mag*(param[1]^2)
                return (energyIntervals,s,mag,last,n)
            end
            n=n+1
        end
        print("Iterations: ")
        println(n)
        #normalization
        s=[s[i]-s[1] for i in 1:length(s)]
        mag=mag*(param[1]^2)
        return (energyIntervals,s,mag,last,n)
    end
end