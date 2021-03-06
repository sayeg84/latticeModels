include("structs.jl")

module InOut
    import ..EdgList
    import ..AdjMat
    import ..SpinLattice
    import ..LatticeGas
    import ..IsingModel
    import ..ChangeSpin!
    using DelimitedFiles, Statistics, Dates
    dest = Dates.format(Dates.now(),"dd-mm-Y_HH-MM")

    """
        ParseArray(array::String)

        Converts an string of form `"[a,b,..,c]"` where a,b,c are Ints or floats to an array of floats.
    """
    function ParseArray(array::String)
        a=split(array,",")
        #remove bracket from first object
        a[1]=a[1][2:end]
        #remove bracket from last object
        a[end]=a[end][1:(end-1)]
        a=[parse(Float64,x) for x in a]
        return a
    end

    function StringToTuple(str::AbstractString)
        nums = split(str,',')
        nums = [parse(Int64,num) for num in nums]
        return Tuple(nums)
    end

    """
        Folders()

        Function to return array with strings of folder names in current working directory
    """
    function Folders()
        a=[file for file in readdir() if isdir(file) && ~(file[1]=='.')]
        return a
    end

    """
        MakeAndEnterDirectories()

        Makes directories to save output of a single current simulation.
    """
    function MakeAndEnterDirectories()
        cd("..")
        if ~(isdir("outputs"))
            mkdir("outputs")
        end
        cd("outputs")
        if ~(isdir(dest))
            try mkdir(dest) catch SystemError end
        end
        cd(dest)
    end



    """
        ExitDirectories()

        Exits from the directory of a simulation.
    """
    function ExitDirectories()
        cd("..")
        cd("..")
        cd("src")
    end

    function WriteSimulParamTable(simulParam)
        open("simulationParameters.csv","w") do io 
            write(io,"B,$(simulParam[1])\n")
            write(io,"J,$(simulParam[2])\n")
            write(io,"C,$(simulParam[3])\n")
            write(io,"kT,$(simulParam[4])\n")
        end
    end

    function WriteMetaParamTable(metaParam)
        open("metaParameters.csv","w") do io 
            write(io,"N,$(metaParam[1])\n")
            write(io,"Periodic,$(metaParam[2])\n")
            write(io,"Geometry,$(metaParam[3])\n")
            write(io,"EnergyFunction,$(metaParam[4])\n")
            write(io,"Model,$(metaParam[5])\n")
        end
    end

    function WriteAlgoParamTable(algoParam,algo)
        open("algorithmParameters.csv","w") do io
            if algo=="metropolis" 
                write(io,"steps,$(algoParam[1])\n")
                write(io,"averages,$(algoParam[2])\n")
            else
                write(io,"Nbins,$(algoParam[1])\n")
                write(io,"Flatness percentage,$(algoParam[2])\n")
                write(io,"Change factor,$(algoParam[3])\n")
                write(io,"Maximum steps,$(algoParam[4])\n")
            end
        end
    end

    function WriteAdjMat(sys::IsingModel)
        M = AdjMat(sys.edgList)
        open("adjMat.csv","w") do io
            DelimitedFiles.writedlm(io,M,',')
        end
    end


    function WriteDOSTable(s,mag,energyIntervals)
        #write("DOS.csv")
        open("DOS.csv","w") do io 
            write(io," E_i, E_f, S, DOS, mag \n  ")
            for i in 1:length(s)
                write(io,"$(energyIntervals[i]), $(energyIntervals[i+1]), $(s[i]), $(round(exp(big(s[i])),digits=5)), $(mag[i]) \n")
            end
        end
    end


    function ReadSimulParamTable()
        table = DelimitedFiles.readdlm("simulationParameters.csv",',')
        table = Array{Float64,1}(table[:,end])
        return table
    end

    function ReadMetaParamTable()
        table = readlines("metaParameters.csv")
        res = []
        for (i,row) in enumerate(table)
            arr = split(row,',')
            if arr[1]=="N"
                tuple = join(arr[2:end],',')
                # remove parenthesis
                tuple = tuple[2:(end-1)]
                push!(res,StringToTuple(tuple))
            else
                push!(res,arr[2])
            end
        end
        return res
    end

    function ReadAlgoParamTable()
        table = DelimitedFiles.readdlm("algorithmParameters.csv",',')
        table = Array(table[:,end])
        return table
    end

    function ReadAdjMat()
        M = DelimitedFiles.readdlm("adjMat.csv",',',Int64)
        return M
    end

    function ReadSimulParamDict()
        table = DelimitedFiles.readdlm("simulParamDict.csv",',',skipstart=1)
        return table
    end



    function ReadDOSTable()
        x=DelimitedFiles.readdlm("DOS.csv",',')
        energyIntervals = Array{Float64,1}(x[2:end,1]) 
        push!(energyIntervals,Float64(x[end,2]))
        s = Array{Float64,1}(x[2:end,3])
        mag = Array{Float64,1}(x[2:end,end])
        return (energyIntervals, s, mag)
    end

    """
        MetropolisSystemsOut(sysArray,param,name)

        Outputs array of matrices sysArray from Metropolis-based simulation with params `param` into a folder with name `name`
    """

    function MetropolisSystemsOut(sysArray,algoParam)
        if ~(isdir("systems"))
            mkdir("systems")
        end
        cd("systems")
        for i in 1:length(sysArray)
            n=Int64((i-1)*algoParam[2])
            name="$(n).csv"
            open(name,"w") do io
                DelimitedFiles.writedlm(io,sysArray[i].sites)
            end
        end
        cd("..")
    end

    #=
    function MetropolisAllOut(initSys,changes,algoParam)
        open("initial.csv","w") do io
            DelimitedFiles.writedlm(io,initSys.sites,',')
        end
        open("changes.csv","w") do io
            for i in 1:algoParam[1]
                write(io,"$(changes[i])\n")
            end
        end
    end

    =#

    function MetropolisAllOut(initSys,res,name)
        open(string(name,"_initial.csv"),"w") do io
            DelimitedFiles.writedlm(io,initSys.sites,',')
        end
        open(string(name,"_changes.csv"),"w") do io
            DelimitedFiles.writedlm(io,res[1],',')
        end
        open(string(name,"_mag.csv"),"w") do io
            DelimitedFiles.writedlm(io,res[3],',')
        end
        open(string(name,"_ener.csv"),"w") do io
            DelimitedFiles.writedlm(io,res[4],',')
        end
    end

    function WriteSimulParamDict(simulParamDict)
        open("simulParamDict.csv","w") do io
            write(io,"index,B,J,C,kT \n  ")
            for pair in simulParamDict
                write(io,string("$(pair.second),",join([string(v) for v in pair.first],",")," \n"))
            end
        end
    end


    function ReadSingleSimul(name)
        original = pwd()
        metaParam = InOut.ReadMetaParamTable()
        algoParam = InOut.ReadAlgoParamTable()
        adjMat = InOut.ReadAdjMat()  
        cd(name)   
        print("Reading ")
        println(name)
        simulParam = InOut.ReadSimulParamTable()
        sys = DelimitedFiles.readdlm("initial.csv",',',Int8)
        sys = reshape(sys,length(sys))
        sys = getfield(Main,Symbol(metaParam[5]))(sys,adjMat)
        changes = DelimitedFiles.readdlm("changes.csv",',',Int32)
        changes = reshape(changes,length(changes))
        cd(original)
        return sys, changes , metaParam, simulParam , algoParam, adjMat
    end

    function ReadSingleSimul(simulParamIndex,iterIndex,metaParam,adjMat)
        # original = pwd()
        # metaParam = InOut.ReadMetaParamTable()
        # algoParam = InOut.ReadAlgoParamTable()
        # adjMat = InOut.ReadAdjMat()
        # paramArray = InOut.ReadParamDict()
        # cd(name)   
        print("Reading ")
        println("$(simulParamIndex),$(iterIndex)")
        sys = DelimitedFiles.readdlm(join([simulParamIndex,iterIndex,"initial.csv"],"_"),',',Int8)
        sys = reshape(sys,length(sys))
        println(metaParam)
        sys = getfield(Main,Symbol(metaParam[5]))(sys,adjMat)
        changes = DelimitedFiles.readdlm(join([simulParamIndex,iterIndex,"changes.csv"],"_"),',',Int32)
        changes = reshape(changes,length(changes))
        M = DelimitedFiles.readdlm(join([simulParamIndex,iterIndex,"mag.csv"],"_"),',',Float64)
        E = DelimitedFiles.readdlm(join([simulParamIndex,iterIndex,"ener.csv"],"_"),',',Float64)
        return sys, changes, M, E
    end
    
    #=

    """
        MetropolisIn(name)

        Read all the matrices of a set of simulations, along with its parameters. 
    """
    function MetropolisIn(name)
        #getting list of directories
        original=pwd()
        cd(name)
        dirs=Folders()
        simulParamArray=[]
        matricesArray=[]
        for folder in dirs
            matrices=[]
            cd(folder)
            cd("matrices")
            simulParam=folder
            #getting list of files
            matricesPaths=[x for x in readdir() if isfile(x) && split(x,".")[end]=="csv"]
            sort!(matricesPaths, by =mat->parse(Int64,split(mat,".")[1]) )
            for mat in matricesPaths
                push!(matrices,DelimitedFiles.readdlm(mat, Float64))
            end
            push!(simulParamArray,simulParam)
            push!(matricesArray,matrices)
            cd("..")
            cd("..")
        end
        cd(original)
        #saving the steps simulation as the last one 
        return Dict(zip(simulParamArray,matricesArray))
        
    end

    """
        MetropolisMeanIn(name)

        Read all the matrices of a set of simulations, along with its parameters. Retrun list of means with the same parameters
    """
    function MetropolisMeanIn(name)
        #getting list of directories
        original=pwd()
        cd(name)
        dirs=Folders()
        sims=Set([split(d,"_")[1] for d in dirs])
        simulParamArray = []
        matricesArray = []
        for simulParam in sims
            print("Reading ")
            println(simulParam)
            newDirs  = [d for d in dirs if split(d,"_")[1] == simulParam]
            simMat= [] 
            @time for folder in newDirs
                matrices=[]
                cd(folder)
                cd("matrices")
                simulParam=folder
                #getting list of files
                matricesPaths=[x for x in readdir() if isfile(x) && split(x,".")[end]=="csv"]
                for mat in matricesPaths
                    push!(matrices,DelimitedFiles.readdlm(mat, Float16))
                end
                push!(simMat,matrices)
                cd("..")
                cd("..")
            end
            push!(simulParamArray,simulParam)
            push!(matricesArray,Statistics.mean(simMat))
        end
        cd(original)
        #saving the steps simulation as the last one 
        return Dict(zip(simulParamArray,matricesArray))
    end
=#

end