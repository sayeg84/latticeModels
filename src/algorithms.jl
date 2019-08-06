
include("statEnsemble.jl")

module Algorithms

    import ..SpinLattice
    import ..LatticeGas
    import ..IsingModel
    import ..ChangeSpin
    import ..ChangeSpin!
    import ..N
    import ..StatEnsemble.Prob
    import ..StatEnsemble.ProbOptimal
    import ..StatEnsemble.NormalEnergy
    import ..RandomPosition



    using Statistics

    """
    Metropolis(x,enerFunc,simulParam::AbstractArray,algoParam::AbstractArray)

    Using x as initial array, makes a metropolis-based simulation using enerFunc, simulParam, and algoParam as parameters. 
    The acceptance rate for spin changes is calculated by calculating explicitly the energy, which makes it slow.

    """
    #@show StatEnsemble.Prob
    function Metropolis(initSys,enerFunc,simulParam::AbstractArray,algoParam::AbstractArray)
        sys = copy(initSys)
        changes = zeros(Int32,algoParam[1])
        for i in 1:algoParam[1]
            pos = RandomPosition(sys)
            acep = Prob(sys,pos,enerFunc,simulParam)
            prob = rand()
            if prob<acep
                ChangeSpin!(sys,pos)
                changes[i] = pos
            else 
                changes[i] = 0
            end
        end
        return changes, sys
    end

    """
    
    MetropolisOptimal(x,enerFunc,simulParam,algoParam)
    
    Using x as initial array, makes a metropolis-based simulation using enerFunc, simulParam, and algoParam as parameters. 
    The acceptance rate for spin changes is calculated by calculating explicitly the energy, which makes it slow.
    
    """
    function MetropolisFast(initSys,enerFunc,simulParam::AbstractArray,algoParam::AbstractArray)
        sys = copy(initSys)
        changes = zeros(Int32,algoParam[1])
        ener = enerFunc(sys,simulParam)
        for i in 1:algoParam[1]
            pos = RandomPosition(initSys)
            newEner = enerFunc(ChangeSpin(sys,pos),simulParam)
            acep = exp(big((-newEner+ener)/simulParam[4]))
            prob = rand()
            if prob<acep
                ChangeSpin!(sys,pos)
                ener = newEner
                changes[i] = pos
            else 
                changes[i] = 0
            end
        end
        return changes, sys
    end


    """
    MetropolisOptimal(x,enerFunc,simulParam,algoParam)
    
    Using x as initial array, makes a metropolis-based simulation using enerFunc, simulParam, and algoParam as parameters. 
    The acceptance rate for spin changes is calculated by calculating explicitly the energy, which makes it slow.
    """
    function MetropolisOptimal(initSys,enerFunc,simulParam,algoParam)
        sys = copy(initSys)
        changes = zeros(Int32,algoParam[1])
        ener = enerFunc(initSys,simulParam)
        for i in 1:algoParam[1]
            pos = RandomPosition(sys)
            acep = ProbOptimal(sys,pos,enerFunc,simulParam)
            prob = rand()
            if prob<acep
                ChangeSpin!(sys,pos)
                ener = ener - log(acep)*simulParam[4]
                changes[i] = pos
            else 
                changes[i] = 0
            end

        end
        return changes, sys
    end
    """
    isFlat(hist, per)

    Function to find out if an histogram is flat dependenig if the minimum is larger than a percentage of the average
    """
    function isFlat(hist;per=0.5)
        min = minimum(hist)
        max = maximum(hist)
        avg = Statistics.mean(hist)
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

    """
        SearchSortedMod(x,a)

        Modified version of searchsortedlast to perform binary search in array x in search of largest index i such that x[i] < a and returning n-1 when  a > x[end]
    """
    function SearchSortedMod(x,a;printLog=false)
        b=searchsortedlast(x,a)
        if printLog
            println("busq")
            println(x)
            println(a)
            println(b)
            println()
        end
            if b==length(x)
            return length(x)-1
        else
            return b
        end
    end

    """
        Index(x,a)
        function that searches for the index of item a in array x. Returns -1 if a is not in x
    """
    function Index(array,element)
        b=find(x -> x==element,array)
        if ~isempty(b)
            return b[1]
        else
            return -1
        end
    end

    function WangLandau(simulParam,algoParam,metaParam,init;printLog=false)
        last = []
        Nbins = Int64(algoParam[1])
        globalHist = zeros(Nbins)
        mag = zeros(Nbins)
        hist = zeros(Nbins)
        s = zeros(Nbins)
        energyIntervals = collect(range(-2,stop=2,length=Int64(Nbins+1)))
        modfact = algoParam[3]
        sys = copy(init)
        n = 0
        while (modfact>=1e-5)
            pos = RandomPosition(sys)
            energyBefore = NormalEnergy(sys,simulParam)
            energyAfter = energyBefore + 2*sys.linearLatt[pos]*(simulParam[2]*sys.linearNeigSumLatt[pos]+simulParam[1])
            energyBefore = energyBefore/(N(sys))
            energyAfter = energyAfter/(N(sys))
            p1 = SearchSortedMod(energyIntervals,energyBefore)
            p2 = SearchSortedMod(energyIntervals,energyAfter)
            tes = rand()
            last = []
            if printLog
                println("latt")
                println(sys.linearLatt)
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
                println(exp(big(s[p1]-s[p2])))
                println("rand")
                println(tes)
                println("avg")
                println(Statistics.mean(hist))
                println("min")
                println(minimum(hist))
            end
            globalHist[p1] += 1
            mag[p1] = (mag[p1]*(globalHist[p1]-1)+abs(sum(sys.linearLatt))/(N(sys)))*1/(globalHist[p1])
            if p2 == 0
                    println("latt")
                    println(sys.linearLatt)
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
                    println(exp(big(s[p1]-s[p2])))
                    println("rand")
                    println(tes)
                    println("avg")
                    println(Statistics.mean(hist))
                    println("min")
                    println(minimum(hist))
                
            end
            η = exp(big(s[p1]-s[p2]))
            if tes < η
                ChangeSpin!(sys,pos)
            end
            globalHist[p1] += 1
            s[p1] += modfact
            hist[p1] += 1
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
            if n == algoParam[4]
                println("Exceded tolerance")
                s = [s[i]-s[1]+log(2) for i in 1:length(s)]
                mag = mag*(N(sys))
                return (energyIntervals,s,mag,last,n)
            end
            n=n+1
        end
        print("Iterations: ")
        println(n)
        s = [s[i]-s[1]+log(2) for i in 1:length(s)]
        mag = mag*(N(sys))
        return (energyIntervals,s,mag,last,n)
    end

    function WangLandauOptimal(simulParam,algoParam,metaParam,init;printLog=false)
        last = []
        Nbins = Int64(algoParam[1])
        globalHist = zeros(Nbins)
        mag = zeros(Nbins)
        hist = zeros(Nbins)
        s = zeros(Nbins)
        energyIntervals = collect(range(-2,stop=2,length=Int64(Nbins+1)))
        modfact = algoParam[3]
        sys = copy(init)
        n = 0
        energyBefore = NormalEnergy(sys,simulParam)
        while (modfact>=1e-5)
            pos = RandomPosition(sys)
            energyAfter = energyBefore + 2*sys.linearLatt[pos]*(simulParam[2]*sys.linearNeigSumLatt[pos]+simulParam[1])
            energyBefore = energyBefore/(N(sys))
            energyAfter = energyAfter/(N(sys))
            p1 = SearchSortedMod(energyIntervals,energyBefore)
            p2 = SearchSortedMod(energyIntervals,energyAfter)
            tes = rand()
            last = []
            if printLog
                println("latt")
                println(sys.linearLatt)
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
                println(exp(big(s[p1]-s[p2])))
                println("rand")
                println(tes)
                println("avg")
                println(Statistics.mean(hist))
                println("min")
                println(minimum(hist))
            end
            globalHist[p1] += 1
            mag[p1] = (mag[p1]*(globalHist[p1]-1)+abs(sum(sys.linearLatt))/(N(sys)))*1/(globalHist[p1])
            if p2 == 0
                    println("latt")
                    println(sys.linearLatt)
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
                    println(exp(big(s[p1]-s[p2])))
                    println("rand")
                    println(tes)
                    println("avg")
                    println(Statistics.mean(hist))
                    println("min")
                    println(minimum(hist))
                
            end
            η = exp(big(s[p1]-s[p2]))
            if tes < η
                ChangeSpin!(sys,pos)
                energyBefore = energyAfter
            end
            energyBefore *= N(sys)
            globalHist[p1] += 1
            s[p1] += modfact
            hist[p1] += 1
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
            if n == algoParam[4]
                println("Exceded tolerance")
                s = [s[i]-s[1]+log(2) for i in 1:length(s)]
                mag = mag*(N(sys))
                return (energyIntervals,s,mag,last,n)
            end
            n=n+1
        end
        print("Iterations: ")
        println(n)
        s = [s[i]-s[1]+log(2) for i in 1:length(s)]
        mag = mag*(N(sys))
        return (energyIntervals,s,mag,last,n)
    end
end
    #=
    function Metropolis(simulParam,algoParam,init,neigLatt,model;printLog=false)
        
        if model=="normal"
            enerFunc=StatEnsemble.NormalEnergy
            nrml=true
        else
            enerFunc=StatEnsemble.PenalizedEnergy3
            nrml=false
        end
        matArr=[]
        mat=initLatt
        push!(matArr,copy(mat))
        for i in 1:(algoParam[1])
            pos=Auxiliar.RandomPosition(mat)
            res=0.0
            res=StatEnsemble.Prob(mat,neigLatt,pos,simulParam,enerFunc,nrml)
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
                mat = Auxiliar.ChangeSpin(mat,pos,normal=nrml)
            end
            if mod(i,algoParam[2])==0
                push!(matArr,copy(mat))
                println(i)
            end
        end
        return matArr
    end
    
    # =============================================
    # Wang-Landau beggining
    

    function WangLandauCycle(simulParam,algoParam,metaParam,initLatt,neigLatt;printLog=false)
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
            energyBefore=energyBefore/(metaParam[1]^metaParam[2])
            energyAfter=energyAfter/(metaParam[1]^metaParam[2])
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
            mag[p1]=(mag[p1]*(globalHist[p1]-1)+abs(sum(latt))/(metaParam[1]^metaParam[2]))*1/(globalHist[p1])
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
                mag=mag*(metaParam[1]^metaParam[2])
                return (energyIntervals,s,mag,last,n)
            end
            n=n+1
        end
        print("Iterations: ")
        println(n)
        #normalization
        s=[s[i]-s[1] for i in 1:length(s)]
        mag=mag*(metaParam[1]^metaParam[2])
        return (energyIntervals,s,mag,last,n)
    end

end
=#
