include("structs.jl")

module InOut
    import ..EdgList
    import ..AdjMat
    import ..SpinLattice
    import ..LatticeGas
    import ..IsingModel
    import ..ChangeSpin!
    using DelimitedFiles, Statistics, Dates
    dest = Dates.format(Dates.now(),"dd-mm-Y HH:MM")

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

    """
        Folders()

        Function to return array with strings of folder names in current working directory
    """
    function Folders()
        a=[x for x in readdir() if isdir(x) && ~(x[1]=='.')]
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
            write(io,"Dimension,$(metaParam[2])\n")
            write(io,"Geometry,$(metaParam[3])\n")
            write(io,"Energy Function,$(metaParam[4])\n")
            write(io,"N,$(metaParam[5])\n")
        end
    end

    function WriteAlgoParamTable(algoParam,algo)
        open("algorithmParameters.csv","w") do io
            if algo=="metropolis" 
                write(io,"steps,$(algoParam[1])\n")
                write(io,"frecuency,$(algoParam[2])\n")
                write(io,"averages,$(algoParam[3])\n")
            else
                write(io,"Nbins,$(algoParam[1])\n")
                write(io,"Flatness percentage,$(algoParam[2])\n")
                write(io,"Change factor,$(algoParam[3])\n")
                write(io,"Maximum steps,$(algoParam[4])\n")
            end
        end
    end

    function WriteAdjMat(x::IsingModel)
        M = AdjMat(x.linearNeigLatt)
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
        x=DelimitedFiles.readdlm("simulationParameters.csv",',')
        x=Array{Float64,1}(x[:,end])
        return x
    end

    function ReadMetaParamTable()
        x=DelimitedFiles.readdlm("metaParameters.csv",',')
        x=Array(x[:,end])
        for i in 1:length(x)
            if typeof(x[i]) == SubString{String}
                x[i] = replace(x[i]," " =>"")
            end
        end
        return x
    end

    function ReadAlgoParamTable()
        x=DelimitedFiles.readdlm("algorithmParameters.csv",',')
        x=Array(x[:,end])
        return x
    end

    function ReadAdjMat()
        M = DelimitedFiles.readdlm("adjMat.csv",',',Int64)
        return M
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
        MetropolisSystemsOut(X,param,name)

        Outputs array of matrices X from Metropolis-based simulation with params `param` into a folder with name `name`
    """

    function MetropolisSystemsOut(X,algoParam)
        if ~(isdir("systems"))
            mkdir("systems")
        end
        cd("systems")
        for i in 1:length(X)
            n=Int64((i-1)*algoParam[2])
            name="$(n).csv"
            open(name,"w") do io
                DelimitedFiles.writedlm(io,X[i].linearLatt)
            end
        end
        cd("..")
    end

    function MetropolisAllOut(initSys,changes,algoParam)
        open("initial.csv","w") do io
            DelimitedFiles.writedlm(io,initSys.linearLatt,',')
        end
        open("changes.csv","w") do io
            for i in 1:algoParam[1]
                write(io,"$(changes[i])\n")
            end
        end
    end

    """
        ReadSingleSimul(name)

        Function to read matrices from a single folder.
    """
    function ReadSingleSimul(name)
        original=pwd()
        cd(name)
        print("Reading ")
        println(name)
        simulParam=InOut.ReadSimulParamTable()
        metaParam=InOut.ReadMetaParamTable()
        algoParam=InOut.ReadAlgoParamTable()
        adjMat=InOut.ReadAdjMat()
        systems=Array{IsingModel,1}()
        cd("systems")
        #getting list of files
        #
        #
        systemsPaths = [x for x in readdir() if isfile(x) &&split(x,".")[end] == "csv"]
        sort!(systemsPaths, by = sys -> parse(Int64,split(sys,".")[1]) )
        for sys in systemsPaths
            sys = DelimitedFiles.readdlm(sys,',', Int64)
            tu = (Int64(sqrt(length(sys))),Int64(sqrt(length(sys))))
            x = SpinLattice(reshape(sys,length(sys)),adjMat,tu)
            push!(systems,x)
        end
        cd(original)
        return systems , metaParam, simulParam , algoParam, adjMat
    end
    function ReadSingleSimul(name,frequency)
        original = pwd()
        cd(name)
        print("Reading ")
        println(name)
        simulParam = InOut.ReadSimulParamTable()
        metaParam = InOut.ReadMetaParamTable()
        algoParam = InOut.ReadAlgoParamTable()
        adjMat = InOut.ReadAdjMat()
        sys = DelimitedFiles.readdlm("initial.csv",',',Int8)
        sys = reshape(sys,length(sys))
        sys = getfield(Main,Symbol(metaParam[5]))(sys,adjMat,(metaParam[1],metaParam[2]))
        changes = DelimitedFiles.readdlm("changes.csv",',',Int32)
        changes = reshape(changes,length(changes))
        systems = Array{IsingModel,1}()
        push!(systems,copy(sys))
        for i in 1:algoParam[1]
            if changes[i] != 0
                ChangeSpin!(sys,changes[i])
            end
            if mod(i,frequency) == 0
                push!(systems,copy(sys))
            end
        end
        cd(original)
        return systems , metaParam, simulParam , algoParam, adjMat
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