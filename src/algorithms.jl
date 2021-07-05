
include("statEnsemble.jl")

module Algorithms

    import ..SpinLattice
    import ..LatticeGas
    import ..IsingModel
    import ..ChangeSpin
    import ..ChangeSpin!
    import ..N
    import ..StatEnsemble.LogProb
    import ..StatEnsemble.LogProbAta
    import ..StatEnsemble.LogProbOptimal
    import ..StatEnsemble.NormalEnergy
    import ..StatEnsemble.Magnetization
    import ..StatEnsemble.MagnetizationUpdate
    import ..RandomPosition



    using Statistics

    """
    Metropolis(initSys,enerFunc,simulParam,algoParam)

    Using initSys as initial array, makes a metropolis-based simulation using enerFunc, simulParam, and algoParam as parameters. 
    The acceptance rate for spin changes is calculated by calculating explicitly the energy, which makes it slow.

    """
    #@show StatEnsemble.LogProb
    function Metropolis(initSys,enerFunc,simulParam,algoParam)
        sys = copy(initSys)
        changes = zeros(Int32,algoParam[1])
        eners = Array{Float32,1}(undef,algoParam[1]+1)
        mags = Array{Float32,1}(undef,algoParam[1]+1)
        ener = enerFunc(sys,simulParam)
        mag = Magnetization(sys,absolute=false)*N(sys)
        eners[1] = ener
        mags[1] = abs(mag)
        for i in 1:algoParam[1]
            pos = RandomPosition(sys)
            acep = LogProb(sys,pos,enerFunc,simulParam)
            prob = rand()
            if prob<acep
                ener = newEner
                mag +=  MagnetizationUpdate(sys,pos)
                ChangeSpin!(sys,pos)
                changes[i] = pos
            else 
                changes[i] = 0
            end
            eners[i+1] = ener
            mags[i+1] = mag
        end
        return changes, sys, mags, eners
    end

    """
    MetropolisNewAta(initSys,enerFunc,simulParam,algoParam)

    Using initSys as initial array, makes a metropolis-based simulation using enerFunc, simulParam, and algoParam as parameters. The probability transition is calculated from `ProbNewAta`, so that it runs faster. 

    """
    function MetropolisNewAta(initSys,enerFunc,simulParam,algoParam)
        sys = copy(initSys)
        changes = zeros(Int32,algoParam[1])
        eners = Array{Float32,1}(undef,algoParam[1]+1)
        mags = Array{Float32,1}(undef,algoParam[1]+1)
        ener = enerFunc(sys,simulParam)
        mag = Magnetization(sys,absolute=false)*N(sys)
        eners[1] = ener
        mags[1] = abs(mag)
        for i in 1:algoParam[1]
            pos = RandomPosition(sys)
            acep = LogProbAta(ener,sys,pos,simulParam)
            prob = log(rand())
            if prob<acep
                ChangeSpin!(sys,pos)
                ener = ener - log(acep)*simulParam[4]
                mag +=  MagnetizationUpdate(sys,pos)
                changes[i] = pos
            else 
                changes[i] = 0
            end
            eners[i+1] = ener
            mags[i+1] = mag
        end
        return changes, sys
    end

    """
    
    MetropolisFast(initSys,enerFunc,simulParam,algoParam)
    
    Using initSys as initial array, makes a metropolis-based simulation using enerFunc, simulParam, and algoParam as parameters. 
    The acceptance rate for spin changes is calculated by calculating explicitly the energy, which makes it slow.
    
    """
    function MetropolisFast(initSys,enerFunc,simulParam,algoParam)
        sys = copy(initSys)
        changes = zeros(Int32,algoParam[1])
        eners = Array{Float32,1}(undef,algoParam[1]+1)
        mags = Array{Float32,1}(undef,algoParam[1]+1)
        ener = enerFunc(sys,simulParam)
        mag = Magnetization(sys,absolute=false)*N(sys)
        eners[1] = ener
        mags[1] = abs(mag)
        for i in 1:algoParam[1]
            pos = RandomPosition(initSys)
            newEner = enerFunc(ChangeSpin(sys,pos),simulParam)
            # using logarithms to make numbers smaller
            acep = (-newEner+ener)/simulParam[4]
            prob = log(rand())
            if prob < acep
                ener = newEner
                mag +=  MagnetizationUpdate(sys,pos)
                ChangeSpin!(sys,pos)
                changes[i] = pos
            else 
                changes[i] = 0
            end
            eners[i+1] = ener
            mags[i+1] = mag
        end
        return changes, sys , mags, eners
    end


    """
    MetropolisOptimal(initSys,enerFunc,simulParam,algoParam)
    
    Using initSys as initial array, makes a metropolis-based simulation using enerFunc, simulParam, and algoParam as parameters. 
    The acceptance rate for spin changes is calculated by calculating explicitly the energy, which makes it slow.
    """
    function MetropolisOptimal(initSys,enerFunc,simulParam,algoParam)
        sys = copy(initSys)
        changes = zeros(Int32,algoParam[1])
        ener = enerFunc(initSys,simulParam)
        for i in 1:algoParam[1]
            pos = RandomPosition(sys)
            acep = LogProbOptimal(sys,pos,enerFunc,simulParam)
            prob = log(rand())
            if prob < acep
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
    function EnergyBounds(sys::IsingModel,simulParam)
        n = N(sys)
        min = 0
        max = max = N(sys)*(simulParam[1]+Order(sys)*(simulParam[2]+simulParam[3]))
    end
    function EnergyBounds(sys::LatticeGas,simulParam)
        n = N(sys)
        min = N(sys)*(-simulParam[1]-Order(sys)*simulParam[2])
        max = N(sys)*(+simulParam[1]+Order(sys)*(simulParam[2]+simulParam[3]))
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
            energyAfter = energyBefore + 2*sys.sites[pos]*(simulParam[2]*sys.neigSum[pos]+simulParam[1])
            energyBefore = energyBefore/(N(sys))
            energyAfter = energyAfter/(N(sys))
            p1 = SearchSortedMod(energyIntervals,energyBefore)
            p2 = SearchSortedMod(energyIntervals,energyAfter)
            tes = rand()
            last = []
            if printLog
                println("latt")
                println(sys.sites)
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
            mag[p1] = (mag[p1]*(globalHist[p1]-1)+abs(sum(sys.sites))/(N(sys)))*1/(globalHist[p1])
            if p2 == 0
                    println("latt")
                    println(sys.sites)
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
            energyAfter = energyBefore + 2*sys.sites[pos]*(simulParam[2]*sys.neigSum[pos]+simulParam[1])
            energyBefore = energyBefore/(N(sys))
            energyAfter = energyAfter/(N(sys))
            p1 = SearchSortedMod(energyIntervals,energyBefore)
            p2 = SearchSortedMod(energyIntervals,energyAfter)
            tes = rand()
            last = []
            if printLog
                println("latt")
                println(sys.sites)
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
            mag[p1] = (mag[p1]*(globalHist[p1]-1)+abs(sum(sys.sites))/(N(sys)))*1/(globalHist[p1])
            if p2 == 0
                    println("latt")
                    println(sys.sites)
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