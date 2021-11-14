using MappedArrays

include("lattices.jl")

function EdgList(adjMat::Array{T,2}) where {T<:Integer}
    n=size(adjMat)[1]
    A=Array{Array{Int64,1},1}()
    for i in 1:n
        B=Array{Int64,1}()
        for j in 1:n
            if adjMat[i,j] == 1
                push!(B,j)
            end
        end
        push!(A,B)
    end
    return A
end

function AdjMat(edgList::Array{Array{T,1},1}) where {T<:Integer}  
    M = zeros(Int8,length(edgList),length(edgList))
    for i in 1:length(edgList)
        M[i,edgList[i]] = ones(length(edgList[i]))
    end
    return M
end


abstract type AbstractSystem end

"""
    SpinLattice

    Construct to store an ising model of -1  and +1 spins. It has the following fields:

    * sites::Array{Int8,1}: the linear list of +1,-1 spins
    * edgList::Array{Array{Int32,1},1}: the edge list of the network where the model is defined
    * neigSum::ReadonlyMappedArray: an array that has the sum of the neighboring spins for a given spin. It updates automatically when the sites change.
    * shape::NTuple: a tuple that represents the real shape of the `sites` array. Its value defaults to `(length(sites))` for non-lattice-inspired networks.
    * n: number of elements
"""
struct SpinLattice <: AbstractSystem

    sites::Array{Int8,1}
    graph::Array{Array{Int32,1},1}
    neigSum::ReadonlyMappedArray
    shape::NTuple
    n::Int32
    """

    SpinLattice(name::AbstractString,shape;periodic::Bool=false,random::Bool=false)

    Create a lattice-inspired system with lattice type `name` and shape `shape`. 
    
    `random` indicates wheter to initialize the lattice with random values on the sites or wheter to initialize then at `1`.

    `periodic` indicates wheter the lattice should have periodic or closed boundary conditions.
    """
    function SpinLattice(name::AbstractString,shape;periodic::Bool=false,random::Bool=false)
        if random
            latt = rand(Array{Int8}([-1,1]),shape)                
        else
            latt = ones(Int8,shape)
        end
        sites = reshape(latt,prod(shape))
        edgList = makeLattice(name,shape,periodic=periodic)[2]
        neigSum = mappedarray(x -> sum(sites[x]),edgList)
        new(sites,edgList,neigSum,shape,prod(shape))
    end
    
    """
    SpinLattice(adjMat::Array{T,2};random::Bool=true) where {T<:Integer}

    Create a SpinLattice from the adyacency matrix of a graph. For now, directed graphs are not supported, so `adjMat` must be symmetri.
    """
    function SpinLattice(adjMat::Array{T,2};random::Bool=true) where {T<:Integer}
        if size(adjMat)[1] != size(adjMat)[2]
            return error("Lattice and adjacency matrix must match")
        end
        if random
            sites = rand(Array{Int8}([-1,1]),size(adjMat)[1])
        else
            sites = ones(Int8,size(adjMat)[1])
        end
        edgList = EdgList(adjMat)
        neigSum = mappedarray(x -> sum(sites[x]),edgList)
        new(sites,edgList,neigSum,(size(adjMat)[1],),size(adjMat)[1])
    end
    
    """
    SpinLattice(latt,adjMat::Array{T,2}) where {T<:Integer}

    Create SpinLattice from values in `latt` and with network given by `adjMat`.
    """
    function SpinLattice(latt,adjMat::Array{T,2}) where {T<:Integer}
        if length(latt) != size(adjMat)[1] 
            return error("Lattice and adjacency matrix must match")
        elseif size(adjMat)[1] != size(adjMat)[2]
            return error("Lattice and adjacency matrix must match")
        end
        sites = copy(latt)
        edgList = EdgList(adjMat)
        neigSum = mappedarray(x -> sum(sites[x]),edgList)
        new(sites,edgList,neigSum,(size(adjMat)[1],),size(adjMat)[1])
    end

    """
    SpinLattice(latt,adjMat,shape)

    Create SpinLattice from values in `latt`, network given by `adjMat` and shape given by `shape`
    """
    function SpinLattice(latt,adjMat,shape)
        if length(latt) != size(adjMat)[1] 
            return error("Lattice and adjacency matrix must match")
        elseif size(adjMat)[1] != size(adjMat)[2]
            return error("Lattice and adjacency matrix must match")
        end
        sites = copy(latt)
        edgList = EdgList(adjMat)
        neigSum = mappedarray(x -> sum(sites[x]),edgList)
        new(sites,edgList,neigSum,shape,prod(shape))
    end
end


"""
    LatticeGas

    Construct to store a lattice gas model of 0 and 1 entries. Fields and constructors are equivalent to the `SpinLattice`. Refer to its documentation for more info.
"""
struct LatticeGas<:AbstractSystem

    sites::Array{Int8,1}
    edgList::Array{Array{Int32,1},1}
    neigSum::ReadonlyMappedArray
    shape::NTuple
    n::Int32

     function LatticeGas(name::AbstractString,shape;periodic::Bool=false,random::Bool=false)
        if random
            latt = rand(Array{Int8}([0,1]),shape)                
        else
            latt = ones(shape)
        end
        sites = reshape(latt,prod(shape))
        edgList = makeLattice(name,shape,periodic=periodic)[2]
        neigSum = mappedarray(x -> sum(sites[x]),edgList)
        new(sites,edgList,neigSum,shape,prod(shape))
    end
    
    function LatticeGas(adjMat::Array{T,2};random::Bool=true) where {T<:Integer}
        if size(adjMat)[1] != size(adjMat)[2]
            return error("Lattice and adjacency matrix must match")
        end
        sites = rand(Array{Int8}([0,1]),size(adjMat)[1])
        edgList = EdgList(adjMat)
        neigSum = mappedarray(x -> sum(sites[x]),edgList)
        new(sites,edgList,neigSum,(size(adjMat)[1],),size(adjMat)[1])
    end
    
    function LatticeGas(latt,adjMat::Array{T,2}) where {T<:Integer}
        if length(latt) != size(adjMat)[1] 
            return error("Lattice and adjacency matrix must match")
        elseif size(adjMat)[1] != size(adjMat)[2]
            return error("Lattice and adjacency matrix must match")
        end
        sites = copy(latt)
        edgList = EdgList(adjMat)
        neigSum = mappedarray(x -> sum(sites[x]),edgList)
        new(sites,edgList,neigSum,(size(adjMat)[1],),size(adjMat)[1])
    end
    
    function LatticeGas(latt,adjMat,shape::Tuple)
        if length(latt) != size(adjMat)[1] 
            return error("Lattice and adjacency matrix must match")
        elseif size(adjMat)[1] != size(adjMat)[2]
            return error("Lattice and adjacency matrix must match")
        end
        sites = copy(latt)
        edgList = EdgList(adjMat)
        neigSum = mappedarray(x -> sum(sites[x]),edgList)
        new(sites,edgList,neigSum,shape,prod(shape))
    end
end

#=

#############################

This was an attempt to use a Graphs.jl object instead of an edge list 
to store the network information. Although it seemed like a good idea,
the recomputation of the edg list any time we needed it made everything too slow 

#############################

using Graphs


struct LatticeGas2<:AbstractSystem

    sites::Array{Int8,1}
    graph::AbstractGraph
    neigSum::ReadonlyMappedArray
    shape::NTuple
    n::Int32

     function LatticeGas2(name::AbstractString,shape;periodic::Bool=false,random::Bool=false)
        if random
            latt = rand(Array{Int8}([0,1]),shape)                
        else
            latt = ones(shape)
        end
        sites = reshape(latt,prod(shape))
        edgList = makeLattice(name,shape,periodic=periodic)[2]
        neigSum = mappedarray(x -> sum(sites[x]),edgList)
        new(sites,SimpleGraph(AdjMat(edgList)),neigSum,shape,prod(shape))
    end
    
    function LatticeGas2(adjMat::Array{T,2};random::Bool=true) where {T<:Integer}
        if size(adjMat)[1] != size(adjMat)[2]
            return error("Lattice and adjacency matrix must match")
        end
        sites = rand(Array{Int8}([0,1]),size(adjMat)[1])
        G = SimpleGraph(adjMat)
        neigSum = mappedarray(x -> sum(sites[x]),G.fadjlist)
        new(sites,G,neigSum,(size(adjMat)[1],),size(adjMat)[1])
    end
    
    function LatticeGas2(latt,adjMat::Array{T,2}) where {T<:Integer}
        if length(latt) != size(adjMat)[1] 
            return error("Lattice and adjacency matrix must match")
        elseif size(adjMat)[1] != size(adjMat)[2]
            return error("Lattice and adjacency matrix must match")
        end
        sites = copy(latt)
        G = SimpleGraph(adjMat)
        neigSum = mappedarray(x -> sum(sites[x]),G.fadjlist)
        new(sites,G,neigSum,(size(adjMat)[1],),size(adjMat)[1])
    end
    
    function LatticeGas2(latt,adjMat,shape::Tuple)
        if length(latt) != size(adjMat)[1] 
            return error("Lattice and adjacency matrix must match")
        elseif size(adjMat)[1] != size(adjMat)[2]
            return error("Lattice and adjacency matrix must match")
        end
        sites = copy(latt)
        G = SimpleGraph(adjMat)
        neigSum = mappedarray(x -> sum(sites[x]),G.fadjlist)
        new(sites,G,neigSum,shape,prod(shape))
    end
end

function Base.sizeof(sys::LatticeGas2)
    return sizeof(sys.sites) + sizeof(sys.graph.fadjlist) + sizeof(sys.neigSum) + sizeof(sys.shape) 
end


=#

function Base.copy(sys::AbstractSystem)
    return typeof(sys)(deepcopy(sys.sites),copy(AdjMat(sys.edgList)),sys.shape)
end



function Base.sizeof(sys::AbstractSystem)
    return sizeof(sys.sites) + sizeof(sys.edgList) + sizeof(sys.neigSum) + sizeof(sys.shape) 
end

function N(sys::AbstractSystem)
    return sys.n
end

function Order(sys::AbstractSystem)
    return sum([length(y) for y in sys.edgList])/2
end

function RandomPosition(sys::AbstractSystem)
    return rand(1:N(sys))
end

function NeigborSum(sys::AbstractSystem,pos::Integer)
    return sum(sys.sites[sys.edgList[pos]])
end

function MutateSingle(sys::SpinLattice,pos::Integer)
    return -1*y.sites[pos]
end

function ChangeSpin(sys::SpinLattice,pos::Integer)
    y = copy(sys)
    y.sites[pos] = -1*y.sites[pos]
    return y
end

function ChangeSpin(sys::LatticeGas,pos::Integer)
    y = copy(sys)
    y.sites[pos] = 1-y.sites[pos]
    return y
end

function ChangeSpin!(sys::SpinLattice,pos::Integer)
    sys.sites[pos] = -1*sys.sites[pos]
end

function ChangeSpin!(sys::LatticeGas,pos::Integer)
    sys.sites[pos] = 1-sys.sites[pos]
end



if abspath(PROGRAM_FILE) == @__FILE__
    #   Testing
    a = SpinLattice("square",(3,3),periodic=true)
    println(a.sites)
    println(a.edgList[1])
    ChangeSpin!(a,2)
    println(a.sites)
    println(a.sites)
    b = copy(a)
    println("copiado")
    println(sum(a.sites))
    println(sum(b.sites))


    n=10
    adjMat = rand([0,1],(n,n))
    for i in 1:n
        adjMat[i,i] = 0
    end
    @show adjMat

    a = SpinLattice(adjMat)
    println(a.sites)
    println(a.edgList[1])

    b=rand([-1,1],n)
    @show b
    a = SpinLattice(b,adjMat)
    println(a.sites)
    println(a.edgList[1])

    n=20
    adjMat = rand([0,1],(n,n))
    for i in 1:n
        adjMat[i,i] = 0
        for j in i+1:n
            if adjMat[i,j]==1
                adjMat[j,i]=1
            end
        end
    end
    @show adjMat

    a = LatticeGas(adjMat)
    println(a.sites)
    println(a.edgList[1])

    @show N(a)
    @show RandomPosition(a)
    @show a.neigSum[1]
end


