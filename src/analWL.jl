println()
println("Wang-Landau Data analyzer")
println()

println()
println("Importing libraries")
println()

include("statEnsemble.jl")
include("inOut.jl")

using DelimitedFiles, ArgParse, Statistics

original=pwd()

println()
println("Parsing arguments")
println()

function ParseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-D", "--Dirname" 
            help = "Directory path to analyse set"
            arg_type = String
            required=true
    end
    return parse_args(s)
end
parsedArgs = ParseCommandline()



function MacroscopicVariables(name,model,temparray)
    M=[]
    E=[]
    CV=[]
    Xi=[]
    #getting list of directories
    original=pwd()
    cd(name)
    cd("Single")
    energyIntervals, s, mag = InOut.ReadDOSTable()
    metaParam = InOut.ReadMetaParamTable()
    for temp in tempArray
        push!(M,StatEnsemble.DOSMagnetization(s,energyIntervals,mag,temp,metaParam))
        push!(E,StatEnsemble.DOSInternalEnergy(s,energyIntervals,temp,metaParam))
        push!(CV,StatEnsemble.DOSHeatCapacity(s,energyIntervals,temp,metaParam))
        push!(Xi,StatEnsemble.DOSMagneticSucep(s,energyIntervals,mag,temp,metaParam))
    end
    cd(original)
    return ( M , E , CV , Xi )
end


function MacroscopicTables(name,tempArray,vals)
    original = pwd()
    cd(name)
    cd("Single")
    open("Macroscopic.csv","w") do f
        write(f,"kT,M,E,CV,Xi\n")
        for i in 1:length(tempArray)
            write(f," $(tempArray[i]) , $(vals[1][i]), $(vals[2][i]), $(vals[3][i]), $(vals[4][i]) \n")
        end
    end
    cd(original)
end


println()
println("Making Analysis")
println()


model = "normal"
tempArray = range(0.1,5,length=50)
vals = MacroscopicVariables(parsedArgs["Dirname"],model,tempArray)
MacroscopicTables(parsedArgs["Dirname"],tempArray,vals)




