println("Wang-Landau ising model")
println()
println()
println("Importing libraries")
println()

include("inOut.jl")
include("algorithms.jl")
include("statEnsemble.jl")
include("auxiliar.jl")
include("geometry.jl")
#using Algorithms,InOut,StatEnsemble, ArgParse
using ArgParse, Statistics, Dates

println()
println("Parsing Arguments")
println()
function ParseCommandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-N", "--Nlatt" 
            help = "Lattice size"
            arg_type = Int64
            default = 10
        "-D", "--dim" 
            help = "Dimension"
            arg_type = Int64
            default = 2
        "-G", "--Geometry" 
            help = "Lattice geometry"
            arg_type = String
            default = "square"
        "-B", "--Bfield"
            help = "Magnetic field"
            arg_type = Float64
            default = 0.0
         "-J", "--Jconst"
            help = "Coupling constant of Ising model"
            arg_type = Float64
            default = 1.0
        "-C", "--Cconst"
            help = "Cycle constant"
            arg_type = Float64
            default = 1.0
        "-T", "--temp"
            help = "Temperature"
            arg_type = Float64
            default = 5.0
        "-P", "--flatness"
            help = "logarithm base 10 of saving frecuency. Must be less than steps"
            arg_type = Float64
            default = 4.0
        "-F", "--Ffactor"
            help = "Modification factor for histogram"
            arg_type = Float64
            default = 1.0
        "-M", "--MaxSteps"
            help = "Logarithm base 10 of the maximum number of steps "
            arg_type = Float64
            default = 9.0
    end
    return parse_args(s)
end

parsedArgs = ParseCommandline()



simulParam=Array{Float64,1}([
    parsedArgs["Bfield"],
    parsedArgs["Jconst"],
    parsedArgs["Cconst"],
    parsedArgs["temp"]
])

geoParam=[
    parsedArgs["Nlatt"],
    parsedArgs["dim"],
    parsedArgs["Geometry"]
]

algoParam=Array{Float64,1}([
    ceil((geoParam[1])*(geoParam[1]/2 - 1)),
    parsedArgs["flatness"],
    parsedArgs["Ffactor"],
    10^parsedArgs["MaxSteps"]
])

println()
println("Initializing Lattice")
println()
# Initializing first lattice

latt , neigLatt =  Geometry.BuildLattices(geoParam)

println()
println("Making simulation")
println()
println()

X=Algorithms.WangLandau(simulParam,algoParam,geoParam,latt,neigLatt,printLog=false)
energyIntervals=X[1]
s=X[2]
tempArray=range(0.1,stop=5,length=50)
ener=[]
cv=[]
mag=[]
for temp in tempArray
    push!(ener,StatEnsemble.DOSEnergy(s,energyIntervals,temp,geoParam)/((geoParam[1]^geoParam[2])))
    push!(cv,StatEnsemble.DOSCV(s,energyIntervals,temp,geoParam)/(geoParam[1]^geoParam[2]))
    push!(mag,StatEnsemble.DOSMag(s,energyIntervals,X[3],StatEnsemble.NormalEnergy,temp,geoParam)/(geoParam[1]^geoParam[2]))
end
X[5]


println()
println("Writing Output")
println()


InOut.MakeAndEnterDirectories()
name="Single"
mkdir(name)
cd(name)
InOut.WriteDOSTable(s,energyIntervals)
InOut.WriteAlgoParamTable(algoParam,"WL")
InOut.WriteSimulParamTable(simulParam)
InOut.WriteGeoParamTable(geoParam)
InOut.ExitDirectories()