module InOut
    using Plots, Rsvg, Dates, DelimitedFiles, Statistics
    gr()
    dest=Dates.format(Dates.now(),"dd-mm-Y HH:MM")

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
            mkdir(dest)
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


    """
        MetropolisOut(X,param,name)

        Outputs array of matrices X from Metropolis-based simulation with params `param` into a folder with name `name`
    """
    function MetropolisOut(X,algoParam)
        if ~(isdir("matrices"))
            mkdir("matrices")
        end
        cd("matrices")
        for i in 1:length(X)
            n=Int64((i-1)*algoParam[2])
            name="$(n).csv"
            open(name,"w") do f
                DelimitedFiles.writedlm(f,X[i])
            end
        end
        cd("..")
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
        geoParam=InOut.ReadGeoParamTable()
        algoParam=InOut.ReadAlgoParamTable()
        matrices=[]
        cd("matrices")
        #getting list of files
        #
        #
        matricesPaths = [x for x in readdir() if isfile(x) && split(x,".")[end] == "csv"]
        sort!(matricesPaths, by =mat->parse(Int64,split(mat,".")[1]) )
        for mat in matricesPaths
            push!(matrices,DelimitedFiles.readdlm(mat, Float64))
        end
        cd(original)
        return matrices , geoParam, simulParam  , algoParam
    end

    """
        ReadSingleMeanSimul(name)

        Function to read mean of all matrices with same simulation parameters.
    """
    function ReadSingleMeanSimul(simulParam)
        original=pwd()
        dirs=Folders()
        print("Reading ")
        println(simulParam)
        newDirs  = [d for d in dirs if split(d,"_")[1] == simulParam]
        simMat= [] 
        for folder in newDirs
            matrices=[]
            cd(folder)
            cd("matrices")
            simulParam=folder
            #getting list of files
            matricesPaths=[x for x in readdir() if isfile(x) && split(x,".")[end]=="csv"]
            sort!(matricesPaths, by =mat->parse(Int64,split(mat,".")[1]) )
            for mat in matricesPaths
                push!(matrices,DelimitedFiles.readdlm(mat, Float16))
            end
            push!(simMat,matrices)
            cd("..")
            cd("..")
        end
        cd(original)
        return mean(simMat)
    end

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

    function WriteSimulParamTable(simulParam)
        open("simulationParameters.csv","w") do f 
            write(f,"B,$(simulParam[1]) \n")
            write(f,"J,$(simulParam[2]) \n")
            write(f,"C,$(simulParam[3]) \n")
            write(f,"kT,$(simulParam[4]) \n")
        end
    end

    function WriteGeoParamTable(geoParam)
        open("geometryParameters.csv","w") do f 
            write(f,"N,$(geoParam[1]) \n")
            write(f,"Dimension,$(geoParam[2]) \n")
            write(f,"Geometry,$(geoParam[3]) \n")
        end
    end

    function WriteAlgoParamTable(algoParam,algo)
        open("algorithmParameters.csv","w") do f
            if algo=="metropolis" 
                write(f,"steps,$(algoParam[1]) \n")
                write(f,"frecuency,$(algoParam[2]) \n")
                write(f,"averages,$(algoParam[3]) \n")
            else
                write(f,"Nbins,$(algoParam[1]) \n")
                write(f,"Flatness percentage,$(algoParam[2]) \n")
                write(f,"Change factor,$(algoParam[3]) \n")
                write(f,"Maximum steps,$(algoParam[4]) \n")
            end
        end
    end

    function ReadSimulParamTable()
        x=DelimitedFiles.readdlm("simulationParameters.csv",',')
        x=Array{Float64,1}(x[:,end])
        return x
    end

    function ReadGeoParamTable()
        x=DelimitedFiles.readdlm("geometryParameters.csv",',')
        x=Array(x[:,end])
        for i in 1:length(x)
            if typeof(x[i]) == SubString{String}
                x[i] = split(x[i]," ")[2]
            end
        end
        return x
    end

    function ReadAlgoParamTable()
        x=DelimitedFiles.readdlm("algorithmParameters.csv",',')
        x=Array(x[:,end])
        return x
    end


    function MakeParamsTable(param,temp,geo)
        cd("..")
        cd("outputs")
        cd(dest)
        open("params.txt","w") do f
            write(f,"Modelo de ising. Fecha: $(dest) \n")
            write(f," ,Tiempo de ejecución: $(round(time,digits=2)) \n")
            write(f,"n: param[1]")
            write(f,"B: param[3]")
            write(f,"J: param[2]")
            write(f,"C: param[7]")
            write(f,"KT: $(temp)")
            writef(f,"Geometría: $(geo)")
        end
        cd("..")
        cd("..")
        cd("src")
    end

    function MakePlots(t,a1,a2,a3,param)
        cd("..")
        cd("outputs")
        cd(dest)
        Plots.scatter(t,a1,
            title="Magnetizacion",
            xlabel="Temperatura (KB=1)",
            ylabel="Magnetizacion absoluta promedio por spin",
            label="simulación con N=$(param[1]^2)"
        )
        #Plots.plot!(T,[magOnsager(t) for t in T],label="Solución teórica")
        #Plots.plot!([2.3,2.3],[0,1],linestyle=:dash,label="Temperatura crítica     ")
        Plots.savefig("mag.png")
        #display(p)


        scatter(t,a2,
            title="Energia",
            xlabel="Temperatura (KB=1)",
            ylabel="Energía promedio por spin",
            label="simulación con N=$(param[1]^2)"
        )

        #Plots.plot!(T,[energiaOnsager(t) for t in T],label="Solución teórica")
        #Plots.plot!([2.3,2.3],[-2,0],linestyle=:dash,label="Temperatura crítica     ")
        Plots.savefig("ener.png")
        #display(p)



        scatter(t,a3,
            ylim=(0.0,maximum(a3[5:end])),
            title="Capacidad calorífica",
            xlabel="Temperatura (KB=1)",
            ylabel="Capacidad calorífica por spin",
            label="simulación con N=$(param[1]^2)"
        )
        #Plots.plot!(T,[cvOnsager(t) for t in T],label="Solución teórica")
        #Plots.plot!([2.3,2.3],[0,2],ylim=(0,2),linestyle=:dash,label="Temperatura crítica      ")
        Plots.savefig("cv.png")
        cd("..")
        cd("..")
        cd("src")
    end

    function MakeTable(t,a1,a2,a3,param,time)
        cd("..")
        cd("outputs")
        cd(dest)
        #write("datos.csv")
        open("datos.csv","w") do f 
            write(f,"Modelo de ising,Fecha, $(dest) \n")
            write(f," ,Tiempo de ejecución, $(round(time,digits=2)) \n")
            write(f," Parámetros \n")
            write(f," n, J, B, pasos, frecuencia, bins, C \n ")
            write(f," $(param[1]), $(param[2]), $(param[3]), $(param[4]), $(param[5]), $(param[6]), $(param[7]) \n ")
            #write(f," $n, $J, $B, $(pasos), $(freq) \n ")
            write(f,"\n")
            write(f,"T,M,E,CV \n")
            for i in 1:length(t)
                write(f,"$(t[i]), $(a1[i]),$(a2[i]), $(a3[i]) \n")
            end
        end
        cd("..")
        cd("..")
        cd("src")
    end


    function WriteDOSTable(s,mag,energyIntervals)
        #write("DOS.csv")
        open("DOS.csv","w") do f 
            write(f," E_i, E_f, S, DOS, mag \n  ")
            for i in 1:length(s)
                write(f,"$(energyIntervals[i]), $(energyIntervals[i+1]), $(s[i]), $(round(exp(big(s[i])),digits=5)), $(mag[i]) \n")
            end
        end
    end


    function ReadDOSTable()
        x=DelimitedFiles.readdlm("DOS.csv",',')
        energyIntervals = Array{Float64,1}(x[2:end,1]) 
        push!(energyIntervals,Float64(x[end,2]))
        s = Array{Float64,1}(x[2:end,3])
        mag = Array{Float64,1}(x[2:end,end])
        return (energyIntervals, s, mag)
    end
end