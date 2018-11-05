
module Algorithms
    include("statEnsemble.jl")
    include("auxiliar.jl")
    include("geometry.jl")
    using Statistics
    #=
    Metropolis algorithm implementation
    =#


    function Metropolis(simulParam,algoParam,initLatt,neigLatt,model;printLog=false)
        
        if model=="normal"
            enerFunc=StatEnsemble.NormalEnergy
            nrml=true
        else
            enerFunc=StatEnsemble.PenalizedEnergy2
            nrml=false
        end
        matArr=[]
        mat=initLatt
        push!(matArr,copy(mat))
        for i in 0:(algoParam[1])
            pos=Auxiliar.RandomPosition(mat)
            res=0.0
            try res=StatEnsemble.Prob(mat,neigLatt,pos,simulParam,enerFunc,nrml)
            catch LoadError 
                res=1
            end
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
                mat=Auxiliar.ChangeSpin(mat,pos,normal=nrml)
            end
            if mod(i,algoParam[2])==0
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
        avg=Statistics.mean(hist)
        if min >= avg*per
            return true
        else
            return false
        end
    end
    function isFlatCount(hist;n=1000)
        if maximum(hist)>=n
            return true
        else
            return false
        end
    end

    function WangLandau(simulParam,algoParam,geoParam,initLatt,neigLatt;printLog=false)
        last=[]
        Nbins=Int64(algoParam[1])
        globalHist=zeros(Nbins)
        mag=zeros(Nbins)
        hist=zeros(Nbins)
        s=zeros(Nbins)
        energyIntervals=collect(range(-2,stop=2,length=Int64(Nbins+1)))
        modfact=algoParam[3]
        latt=copy(initLatt)
        n=0
        while (modfact>=1e-5)
            pos=Auxiliar.RandomPosition(latt)
            energyBefore=StatEnsemble.NormalEnergy(latt,neigLatt,simulParam)
            energyAfter=energyBefore+2*latt[pos]*(simulParam[2]*Geometry.NeighborSum(latt,neigLatt,pos)+simulParam[1])
            energyBefore=energyBefore/((geoParam[1])^(geoParam[2]))
            energyAfter=energyAfter/((geoParam[1])^(geoParam[2]))
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
                println(Statistics.mean(hist))
                println("min")
                println(minimum(hist))
            end
            globalHist[p1]+=1
            mag[p1]=(mag[p1]*(globalHist[p1]-1)+abs(sum(latt))/(geoParam[1]^geoParam[2]))*1/(globalHist[p1])
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
                hist=zeros(Nbins)
                println(n)
                println("change")
                println(modfact)
            end
            if mod(n,10^5) == 0 && ~printLog
                print("iter: ")
                println(n)
                println()
            end
            if n==algoParam[4]
                println("Exceded tolerance")
                s=[s[i]-s[1]+log(2) for i in 1:length(s)]
                mag=mag*(geoParam[1]^geoParam[2])
                return (energyIntervals,s,mag,last,n)
            end
            n=n+1
        end
        print("Iterations: ")
        println(n)
        s=[s[i]-s[1]+log(2) for i in 1:length(s)]
        mag=mag*(geoParam[1]^geoParam[2])
        return (energyIntervals,s,mag,last,n)
    end

    function WangLandauCycle(simulParam,algoParam,geoParam,initLatt,neigLatt;printLog=false)
        last=[]
        Nbins=Int64(algoParam[1])
        globalHist=zeros(Nbins)
        mag=zeros(Nbins)
        hist=zeros(Nbins)
        s=zeros(Nbins)
        energyIntervals=collect(range(-2,stop=2,length=Int64(algoParam[1]+1)))
        modfact=1
        latt=copy(initLatt)
        n=0
        while (modfact>=1e-5)
            pos=Auxiliar.RandomPosition(latt)
            if printLog
                println(latt)
            end
            try energyBefore=StatEnsemble.PenalizedEnergy2(latt,neigLatt,simulParam,printLog=printLog)
            catch LoadError
                #println("antes")
                #println(latt)
                #println(energyBefore)
                #error("Wrong lattice")
            end
            try energyAfter=StatEnsemble.PenalizedEnergy2(Auxiliar.ChangeSpin(latt,pos,normal=false),param,neigLatt,printLog=printLog)
            catch LoadError 
                #println("después")
                #println(Auxiliar.ChangeSpin(latt,pos))
                #println(energyAfter)
                #error("Wrong lattice")
            end
            energyBefore=energyBefore/(geoParam[1]^geoParam[2])
            energyAfter=energyAfter/(geoParam[1]^geoParam[2])
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
                println(Statistics.mean(hist))
                println("min")
                println(minimum(hist))
            end
            globalHist[p1]+=1
            mag[p1]=(mag[p1]*(globalHist[p1]-1)+abs(sum(latt))/(geoParam[1]^geoParam[2]))*1/(globalHist[p1])
            η=exp(big(s[p1]-s[p2]))
            if  tes < η
                latt=Auxiliar.ChangeSpin(latt,pos,normal=false)
            end
            s[p1]+=modfact
            hist[p1]+=1
            if isFlat(hist)
                modfact=modfact*1/2
                last=copy(hist)
                hist=zeros(Nbins)
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
                mag=mag*(geoParam[1]^geoParam[2])
                return (energyIntervals,s,mag,last,n)
            end
            n=n+1
        end
        print("Iterations: ")
        println(n)
        #normalization
        s=[s[i]-s[1] for i in 1:length(s)]
        mag=mag*(geoParam[1]^geoParam[2])
        return (energyIntervals,s,mag,last,n)
    end

end