    using MappedArrays


    function EdgList(adjMat)
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

    function AdjMat(edgList)
        M = zeros(Int8,length(edgList),length(edgList))
        for i in 1:length(edgList)
            M[i,edgList[i]] = ones(length(edgList[i]))
        end
        return M
    end

    """
        SpinLattice

        Construct to store an ising model of -1  and +1 spins
    """


    struct SpinLattice
        latt::AbstractArray
        neigLatt::AbstractArray 
        neigSumLatt::ReadonlyMappedArray
        linearLatt::Array{Int8,1}
        linearNeigLatt::Array{Array{Int64,1},1}
        linearNeigSumLatt::ReadonlyMappedArray
        
        
        function SpinLattice(neigFunc::Function,l::Int64,d::Int64)
            sizes = fill(l,d)
            sizes = tuple(sizes...)
            latt = rand(Array{Int8}([-1,1]),sizes)
            neigLatt = Array{Array{CartesianIndex{length(sizes)},1},length(sizes)}(undef,sizes)
            for pos in CartesianIndices(sizes)
                neigLatt[pos] = neigFunc(latt,pos)
            end
            neigSumLatt = mappedarray(x -> sum(latt[x]),neigLatt) 
            linearLatt = reshape(latt,prod(sizes))
            linearNeigLatt = Array{Array{Int64,1},length(sizes)}(undef,sizes)
            aux = LinearIndices(sizes)
            for pos in CartesianIndices(sizes)
                linearNeigLatt[pos] = aux[neigLatt[pos]]
            end
            linearNeigLatt = reshape(linearNeigLatt,prod(sizes))
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)
            new(latt,neigLatt,neigSumLatt,linearLatt,linearNeigLatt,linearNeigSumLatt)
        end
        
        function SpinLattice(adjMat::Array{Int64,2})
            if size(adjMat)[1] != size(adjMat)[2]
                return error("Lattice and adjacency matrix must match")
            end
            linearLatt = rand(Array{Int8}([-1,1]),size(adjMat)[1])
            linearNeigLatt = EdgList(adjMat)
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)
            latt = linearLatt
            neigLatt = linearNeigLatt
            neigSumLatt = mappedarray(x -> sum(latt[x]),neigLatt) 
            new(latt,neigLatt,neigSumLatt,linearLatt,linearNeigLatt,linearNeigSumLatt)
        end
        
        function SpinLattice(latt::Array{Int64,1},adjMat::Array{Int64,2})
            if length(latt) != size(adjMat)[1] 
                return error("Lattice and adjacency matrix must match")
            elseif size(adjMat)[1] != size(adjMat)[2]
                return error("Lattice and adjacency matrix must match")
            end
            linearLatt = copy(latt)
            linearNeigLatt = EdgList(adjMat)
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)
            latt = linearLatt
            neigLatt = linearNeigLatt
            neigSumLatt = mappedarray(x -> sum(latt[x]),neigLatt) 
            new(latt,neigLatt,neigSumLatt,linearLatt,linearNeigLatt,linearNeigSumLatt)
        end
    
        function SpinLattice(latt::AbstractArray,neigLatt::AbstractArray)
            sizes=size(latt)
            neigSumLatt = mappedarray(x -> sum(latt[x]),neigLatt) 
            linearLatt = reshape(latt,prod(sizes))
            linearNeigLatt = Array{Array{Int64,1},length(sizes)}(undef,sizes)
            aux = LinearIndices(sizes)
            for pos in CartesianIndices(sizes)
                linearNeigLatt[pos] = aux[neigLatt[pos]]
            end
            linearNeigLatt = reshape(linearNeigLatt,prod(sizes))
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)
            new(latt,neigLatt,neigSumLatt,linearLatt,linearNeigLatt,linearNeigSumLatt)
        end
        
    end
    """
        LatticeGas

        Construct to store a lattice gas model of 0 and 1 entries
    """
    struct LatticeGas
        latt::AbstractArray
        neigLatt::AbstractArray 
        neigSumLatt::AbstractArray
        linearLatt::Array{Int8,1}
        linearNeigLatt::Array{Array{Int64,1},1}
        linearNeigSumLatt::Array{Int8,1}

        function LatticeGas(neigFunc::Function,l::Int64,d::Int64)
            sizes = fill(l,d)
            sizes = tuple(sizes...)
            latt = rand(Array{Int8}([0,1]),sizes)
            neigLatt = Array{Array{CartesianIndex{length(sizes)},1},length(sizes)}(undef,sizes)
            for pos in CartesianIndices(sizes)
                neigLatt[pos] = neigFunc(latt,pos)
            end
            neigSumLatt = mappedarray(x -> sum(latt[x]),neigLatt) 
            linearLatt = reshape(latt,prod(sizes))
            linearNeigLatt = Array{Array{Int64,1},length(sizes)}(undef,sizes)
            aux = LinearIndices(sizes)
            for pos in CartesianIndices(sizes)
                linearNeigLatt[pos] = aux[neigLatt[pos]]
            end
            linearNeigLatt = reshape(linearNeigLatt,prod(sizes))
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)
            new(latt,neigLatt,neigSumLatt,linearLatt,linearNeigLatt,linearNeigSumLatt)
        end
        
        function LatticeGas(adjMat::Array{Int64,2})
            if size(adjMat)[1] != size(adjMat)[2]
                return error("Lattice and adjacency matrix must match")
            end
            linearLatt = rand(Array{Int8}([0,1]),size(adjMat)[1])
            linearNeigLatt = EdgList(adjMat)
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)
            latt = linearLatt
            neigLatt = linearNeigLatt
            neigSumLatt = mappedarray(x -> sum(latt[x]),neigLatt) 
            new(latt,neigLatt,neigSumLatt,linearLatt,linearNeigLatt,linearNeigSumLatt)
        end
        

        function LatticeGas(latt::Array{Int64,1},adjMat::Array{Int64,2})
            if length(latt) != size(adjMat)[1] 
                return error("Lattice and adjacency matrix must match")
            elseif size(adjMat)[1] != size(adjMat)[2]
                return error("Lattice and adjacency matrix must match")
            end
            linearLatt = copy(latt)
            linearNeigLatt = EdgList(adjMat)
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)
            latt = linearLatt
            neigLatt = linearNeigLatt
            neigSumLatt = mappedarray(x -> sum(latt[x]),neigLatt) 
            new(latt,neigLatt,neigSumLatt,linearLatt,linearNeigLatt,linearNeigSumLatt)
        end
        
        function LatticeGas(str::String,latt::AbstractArray,neigLatt::AbstractArray)
            sizes=size(latt)
            neigSumLatt = mappedarray(x -> sum(latt[x]),neigLatt) 
            linearLatt = reshape(latt,prod(sizes))
            linearNeigLatt = Array{Array{Int64,1},length(sizes)}(undef,sizes)
            aux = LinearIndices(sizes)
            for pos in CartesianIndices(sizes)
                linearNeigLatt[pos] = aux[neigLatt[pos]]
            end
            linearNeigLatt = reshape(linearNeigLatt,prod(sizes))
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)
            new(latt,neigLatt,neigSumLatt,linearLatt,linearNeigLatt,linearNeigSumLatt)
        end


    end
    # Auxiliary struct to define functions for both systems

    IsingModel = Union{LatticeGas, SpinLattice}

    function Size(x::IsingModel)
        return size(x.latt)
    end

    function Dim(x::IsingModel)
        return length(Size(x))
    end

    function N(x::IsingModel)
        return length(x.linearLatt)
    end

    #function RandomPosition(x::IsingModel)
    function RandomPosition(x::IsingModel)
        return rand(1:N(x))
    end

    function NeigborSum(x::IsingModel,pos::Int64)
        return sum(x.linearLatt[x.linearNeigLatt[pos]])
    end

    function ChangeSpin(x::SpinLattice,pos::Int64)
        y=copy(x)
        y.linearLatt[pos] = -1*y.linearLatt[pos]
        return y
    end

    function ChangeSpin(y::LatticeGas,pos::Int64)
        y=copy(x)
        y.linearLatt[pos] = 1-y.linearLatt[pos]
        return y
    end

    function ChangeSpin!(x::SpinLattice,pos::Int64)
        x.linearLatt[pos] = -1*x.linearLatt[pos]
    end

    function ChangeSpin!(x::LatticeGas,pos::Int64)
        x.linearLatt[pos] = 1-x.linearLatt[pos]
    end

    function Base.copy(x::SpinLattice)
        return SpinLattice(deepcopy(x.linearLatt),deepcopy(x.linearNeigLatt))
    end

    function Base.copy(x::LatticeGas)
        return LatticeGas(deepcopy(x.linearLatt),deepcopy(x.linearNeigLatt))
    end

    #   Testing
    include("lattices.jl")
    #=
    
    a = SpinLattice(Lattices.PeriodicSquareLatticeNeighbors,3,2)
    println(a.latt)
    println(a.linearLatt)
    println(a.neigLatt[1,1])
    println(a.latt[a.neigLatt[1,1]])
    println(a.neigSumLatt[1,1])
    println(a.linearNeigLatt[1])
    ChangeSpin!(a,2)
    println(a.neigSumLatt[1,1])
    println(a.latt)
    println(a.linearLatt)
    println(a.linearLatt)
    b = copy(a)
    println("copiado")
    println(sum(a.linearLatt))
    println(sum(b.linearLatt))
=#

    a = SpinLattice(Lattices.PeriodicSquareLatticeNeighbors,100,2)
    println(sizeof(a))
    #InteractiveUtils.varinfo()
    
    #=
    n=10
    adjMat = rand([0,1],(n,n))
    for i in 1:n
        adjMat[i,i] = 0
    end
    @show adjMat

    a = SpinLattice(adjMat)
    println(a.latt)
    println(a.linearLatt)
    println(a.neigLatt[1,1])
    println(a.linearNeigLatt[1])

    b=rand([-1,1],n)
    @show b
    a = SpinLattice(b,adjMat)
    println(a.latt)
    println(a.linearLatt)
    println(a.neigLatt[1,1])
    println(a.linearNeigLatt[1])

    n=20
    adjMat = rand([0,1],(n,n))
    for i in 1:n
        adjMat[i,i] = 0
    end
    @show adjMat

    a = LatticeGas(adjMat)
    println(a.latt)
    println(a.linearLatt)
    println(a.neigLatt[1,1])
    println(a.linearNeigLatt[1])

    @show Size(a)
    @show Dim(a)
    @show N(a)
    @show RandomPosition(a)

    =#