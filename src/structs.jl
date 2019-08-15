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
        linearLatt::Array{Int8,1}
        linearNeigLatt::Array{Array{Int32,1},1}
        linearNeigSumLatt::ReadonlyMappedArray
        shape::Tuple{Int16,Int16}
        
        function SpinLattice(neigFunc::Function,l::Integer,d::Integer)
            sizes = fill(l,d)
            sizes = tuple(sizes...)
            latt = rand(Array{Int8}([-1,1]),sizes)
            neigLatt = Array{Array{CartesianIndex{length(sizes)},1},length(sizes)}(undef,sizes)
            for pos in CartesianIndices(sizes)
                neigLatt[pos] = neigFunc(latt,pos)
            end
            neigSumLatt = mappedarray(x -> sum(latt[x]),neigLatt) 
            linearLatt = reshape(latt,prod(sizes))
            linearNeigLatt = Array{Array{Int32,1},length(sizes)}(undef,sizes)
            aux = LinearIndices(sizes)
            for pos in CartesianIndices(sizes)
                linearNeigLatt[pos] = aux[neigLatt[pos]]
            end
            linearNeigLatt = reshape(linearNeigLatt,prod(sizes))
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)

            new(linearLatt,linearNeigLatt,linearNeigSumLatt,(l,d))
        end
        
        function SpinLattice(adjMat)
            if size(adjMat)[1] != size(adjMat)[2]
                return error("Lattice and adjacency matrix must match")
            end
            linearLatt = rand(Array{Int8}([-1,1]),size(adjMat)[1])
            linearNeigLatt = EdgList(adjMat)
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)
            new(linearLatt,linearNeigLatt,linearNeigSumLatt,(size(adjMat)[1],1))
        end
        
        function SpinLattice(latt,adjMat)
            if length(latt) != size(adjMat)[1] 
                return error("Lattice and adjacency matrix must match")
            elseif size(adjMat)[1] != size(adjMat)[2]
                return error("Lattice and adjacency matrix must match")
            end
            linearLatt = copy(latt)
            linearNeigLatt = EdgList(adjMat)
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)
            new(linearLatt,linearNeigLatt,linearNeigSumLatt,(size(adjMat)[1],1))
        end
    
        function SpinLattice(latt,adjMat,shape)
            if length(latt) != size(adjMat)[1] 
                return error("Lattice and adjacency matrix must match")
            elseif size(adjMat)[1] != size(adjMat)[2]
                return error("Lattice and adjacency matrix must match")
            end
            linearLatt = copy(latt)
            linearNeigLatt = EdgList(adjMat)
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)
            new(linearLatt,linearNeigLatt,linearNeigSumLatt,shape)
        end
    end
    """
        LatticeGas

        Construct to store a lattice gas model of 0 and 1 entries
    """
    struct LatticeGas
        linearLatt::Array{Int8,1}
        linearNeigLatt::Array{Array{Int32,1},1}
        linearNeigSumLatt::ReadonlyMappedArray
        shape::Tuple{Int16,Int16}

        function LatticeGas(neigFunc::Function,l::Integer,d::Integer)
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
            new(linearLatt,linearNeigLatt,linearNeigSumLatt,(l,d))
        end
        
        function LatticeGas(adjMat)
            if size(adjMat)[1] != size(adjMat)[2]
                return error("Lattice and adjacency matrix must match")
            end
            linearLatt = rand(Array{Int8}([0,1]),size(adjMat)[1])
            linearNeigLatt = EdgList(adjMat)
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)
            new(linearLatt,linearNeigLatt,linearNeigSumLatt,(size(adjMat)[1],1))
        end
        

        function LatticeGas(latt,adjMat)
            if length(latt) != size(adjMat)[1] 
                return error("Lattice and adjacency matrix must match")
            elseif size(adjMat)[1] != size(adjMat)[2]
                return error("Lattice and adjacency matrix must match")
            end
            linearLatt = copy(latt)
            linearNeigLatt = EdgList(adjMat)
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)
            new(linearLatt,linearNeigLatt,linearNeigSumLatt,(size(adjMat)[1],1))
        end
        
        function LatticeGas(latt,adjMat,shape::Tuple)
            if length(latt) != size(adjMat)[1] 
                return error("Lattice and adjacency matrix must match")
            elseif size(adjMat)[1] != size(adjMat)[2]
                return error("Lattice and adjacency matrix must match")
            end
            linearLatt = copy(latt)
            linearNeigLatt = EdgList(adjMat)
            linearNeigSumLatt = mappedarray(x -> sum(linearLatt[x]),linearNeigLatt)
            new(linearLatt,linearNeigLatt,linearNeigSumLatt,shape)
        end


    end
    # Auxiliary struct to define functions for both systems

    IsingModel = Union{LatticeGas, SpinLattice}


    function Dim(x::IsingModel)
        return x.shape[2]
    end

    function N(x::IsingModel)
        return length(x.linearLatt)
    end
    function Order(x::IsingModel)
        return sum([length(y) for y in x.linearNeigLatt])/2
    end

    #function RandomPosition(x::IsingModel)
    function RandomPosition(x::IsingModel)
        return rand(1:N(x))
    end

    function NeigborSum(x::IsingModel,pos::Integer)
        return sum(x.linearLatt[x.linearNeigLatt[pos]])
    end

    function ChangeSpin(x::SpinLattice,pos::Integer)
        y=copy(x)
        y.linearLatt[pos] = -1*y.linearLatt[pos]
        return y
    end

    function ChangeSpin(x::LatticeGas,pos::Integer)
        y=copy(x)
        y.linearLatt[pos] = 1-y.linearLatt[pos]
        return y
    end

    function ChangeSpin!(x::SpinLattice,pos::Integer)
        x.linearLatt[pos] = -1*x.linearLatt[pos]
    end

    function ChangeSpin!(x::LatticeGas,pos::Integer)
        x.linearLatt[pos] = 1-x.linearLatt[pos]
    end

    function Base.copy(x::SpinLattice)
        return SpinLattice(deepcopy(x.linearLatt),deepcopy(AdjMat(x.linearNeigLatt)),x.shape)
    end

    function Base.copy(x::LatticeGas)
        return LatticeGas(deepcopy(x.linearLatt),deepcopy(AdjMat(x.linearNeigLatt)),x.shape)
    end
    function Base.sizeof(x::IsingModel)
        return sizeof(x.linearLatt) + sizeof(x.linearNeigLatt) + sizeof(x.linearNeigSumLatt) + sizeof(x.shape)
    end
    

    #   Testing
    #=
    include("lattices.jl")
    
    
    a = SpinLattice(Lattices.PeriodicSquareLatticeNeighbors,3,2)
    println(a.linearLatt)
    println(a.linearNeigLatt[1])
    ChangeSpin!(a,2)
    println(a.linearLatt)
    println(a.linearLatt)
    b = copy(a)
    println("copiado")
    println(sum(a.linearLatt))
    println(sum(b.linearLatt))
    

    n=10
    adjMat = rand([0,1],(n,n))
    for i in 1:n
        adjMat[i,i] = 0
    end
    @show adjMat

    a = SpinLattice(adjMat)
    println(a.linearLatt)
    println(a.linearNeigLatt[1])

    b=rand([-1,1],n)
    @show b
    a = SpinLattice(b,adjMat)
    println(a.linearLatt)
    println(a.linearNeigLatt[1])

    n=20
    adjMat = rand([0,1],(n,n))
    for i in 1:n
        adjMat[i,i] = 0
    end
    @show adjMat

    a = LatticeGas(adjMat)
    println(a.linearLatt)
    println(a.linearNeigLatt[1])

    @show Dim(a)
    @show N(a)
    @show RandomPosition(a)
    @show a.linearNeigSumLatt[1]
=#

  