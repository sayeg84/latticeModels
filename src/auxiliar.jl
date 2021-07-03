using ArgParse, Statistics, Dates

include("lattices.jl")
include("algorithms.jl")
include("inOut.jl")


function setParseCommandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "-N", "--Nlatt" 
            help = "Lattice shape. Must be compatible with Network shape. Is ignored if network is custom."
            arg_type = String
            default = "10,10"
        "-P", "--Periodic" 
            help = "Boolean to cerify if lattice has periodic boundary conditions. Must be compatible with Network shape. Is ignored if network is custom."
            arg_type = Bool
            default = true
        "-L", "--Network" 
            help = "Network for running the model. Must be direction of a CSV file containing an adjacency matrix or the name of a compatible lattice-creator function"
            arg_type = String
            default = "square"
        "-E", "--EnerFunc" 
            help = "Energy function"
            arg_type = String
            default = "NormalEnergy"
        "-M", "--Model" 
            help = "Lattice geometry"
            arg_type = String
            default = "SpinLattice"
        "-S", "--Sweeps"
            help = "logarithm base 10 of total sweeps"
            arg_type = Float64
            default = 3.0
        "-A", "--Systems"
            help = "Number of independent systems simulated for same simulation parameters"
            arg_type = Int64
            default = 20
        "-B", "--Barray"
            help = "String describing the begining, end and length of the Barray. Reverse order (begining > end) is supported"
            arg_type = String
            default = "-3.5,0.0,21"
        "-J", "--Jarray"
            help = "String describing the begining, end and length of the Jarray. Reverse order (begining > end) is supported"
            arg_type = String
            default = "2.0,2.0,1"
        "-C", "--Carray"
            help = "String describing the begining, end and length of the Carray. Reverse order (begining > end) is supported"
            arg_type = String
            default = "0.2,1.8,21"
        "-T", "--kTarray"
            help = "String describing the begining, end and length of the kTarray. Reverse order (begining > end) is supported"
            arg_type = String
            default = "0.5,0.5,1"
        "-O", "--order"
            help = "String describing the order of parameters to retrive. String must be a permutation of B,J,C,kT separated by commas"
            arg_type = String
            default = "B,J,C,kT"
    end
    return parse_args(s)
end

function StringToTuple(str::AbstractString)
    nums = split(str,',')
    nums = [parse(Int64,num) for num in nums]
    return Tuple(nums)
end

function TupleToRange(str)
    numbers = split(str,',')
    if length(numbers) != 3
        error("String $(str) non parsable")
    else
        tup = (parse(Float64,numbers[1]),parse(Float64,numbers[2]),parse(Int64,numbers[3]))
        if tup[3] == 1
            return tup[1]
        else
            return range(tup[1], stop = tup[2], length = tup[3])
        end
    end
end

function getParams(parsedArgs)
    metaParam=(
        StringToTuple(parsedArgs["Nlatt"]),
        parsedArgs["Periodic"],
        parsedArgs["Network"],
        parsedArgs["EnerFunc"],
        parsedArgs["Model"]
    )

    algoParam=NTuple{2,Int64}((
        floor(10^parsedArgs["Sweeps"])*prod(metaParam[1]),
        parsedArgs["Systems"]
    ))

    return metaParam, algoParam
end

function getSimulParams(parsedArgs)
    arrays = []
    varsOrder = split(parsedArgs["order"],',')
    for var in varsOrder
        println(var)
        array = TupleToRange(parsedArgs[string(var,"array")])
        println(array)
        push!(arrays,array)
    end
    orderDict = Dict(varsOrder[i] => i for i in 1:length(varsOrder))
    params1 = [simulParam for simulParam in Iterators.product(arrays...)]
    params1 = reshape(params1,length(params1))
    # combination of parameters in desired order
    params1 = [simulParam[[orderDict["B"],orderDict["J"],orderDict["C"],orderDict["kT"]]] for simulParam in params1]
    # adding number of iterations
    params2 = [(par...,iter) for iter in 1:algoParam[2] for par in params1]
    simulParamDict = Dict(params1[i] => i for i in 1:length(params1))
    return params2, simulParamDict
end
