module Geometry
    function Modl(a::Int64,b::Int64)
        x=mod(a,b)
        if x==0 
            return b
        else
            return x
        end
    end
    ####################################################################################################################
    #Legacy code
    
    #square lattice neighbors sum
    #=
    function SquareLatticeNeighbors(latt,pos;printLog=false)
        s=size(latt)
        dim=length(s)
        if dim==1
            v1=latt[Modl(pos[1]-1,s[1])]
            v2=latt[Modl(pos[1]+1,s[1])]
            if printLog
                println([v1 0 v2])
            end
            return v1+v2
        elseif dim==2
            v1 = latt[Modl(pos[1]-1,s[1]), pos[2]]
            v2 = latt[Modl(pos[1]+1,s[1]), pos[2]]
            v3 = latt[pos[1],Modl(pos[2]-1,s[2])]
            v4 = latt[pos[1],Modl(pos[2]+1,s[2])]
            if printLog
                println([0 v1 0; v3 0 v4 ; 0 v2 0])
            end
            return v1+v2+v3+v4
        elseif dim==3
            v1 = latt[Modl(pos[1]-1,s[1]), pos[2], pos[3]]
            v2 = latt[Modl(pos[1]+1,s[1]), pos[2], pos[3]]
            v3 = latt[pos[1], Modl(pos[2]-1,s[2]), pos[3]]
            v4 = latt[pos[1], Modl(pos[2]+1,s[2]), pos[3]]
            v5 = latt[pos[1], pos[2], Modl(pos[3]-1,s[3])]
            v6 = latt[pos[1], pos[2], Modl(pos[3]+1,s[3])]
            if printLog
                println([0 v1 0; v3 0 v4 ; 0 v2 0])
                println([v5 0 v6])
            end
            return v1+v2+v3+v4+v5+v6
        else
            error("dimension not soported")
        end
    end
    =#

    #=
    function SquareLatticeNeighborsArray(latt,pos;printLog=false)
        s=size(latt)
        dim=length(s)
        if dim==1
            v1=latt[Modl(pos[1]-1,s[1])]
            v2=latt[Modl(pos[1]+1,s[1])]
            if printLog
                println([v1 0 v2])
            end
            return [v1,v2]
        elseif dim==2
            v1 = latt[Modl(pos[1]-1,s[1]), pos[2]]
            v2 = latt[Modl(pos[1]+1,s[1]), pos[2]]
            v3 = latt[pos[1],Modl(pos[2]-1,s[2])]
            v4 = latt[pos[1],Modl(pos[2]+1,s[2])]
            if printLog
                println([0 v1 0; v3 0 v4 ; 0 v2 0])
            end
            return [v1,v2,v3,v4]
        elseif dim==3
            v1 = latt[Modl(pos[1]-1,s[1]), pos[2], pos[3]]
            v2 = latt[Modl(pos[1]+1,s[1]), pos[2], pos[3]]
            v3 = latt[pos[1], Modl(pos[2]-1,s[2]), pos[3]]
            v4 = latt[pos[1], Modl(pos[2]+1,s[2]), pos[3]]
            v5 = latt[pos[1], pos[2], Modl(pos[3]-1,s[3])]
            v6 = latt[pos[1], pos[2], Modl(pos[3]+1,s[3])]
            if printLog
                println([0 v1 0; v3 0 v4 ; 0 v2 0])
                println([v5 0 v6])
            end
            return [v1,v2,v3,v4,v5,v6]
        else
            error("dimension not soported")
        end
    end
    =#
    ####################################################################################################
    """
        SquareLatticeNeighbors(latt,pos)

        For a given position pos in a lattice latt, it returns an array of neighbor positions (represented as cartesian indexes) for an square lattice geometry. Supports dimensions up to 3.
    """
    function SquareLatticeNeighbors(latt,pos;printLog=false)
        s=size(latt)
        dim=length(s)
        if dim==1
            v1=CartesianIndex(Modl(pos[1]-1,s[1]))
            v2=CartesianIndex(Modl(pos[1]+1,s[1]))
            if printLog
                println([v1 0 v2])
            end
            return [v1,v2]
        elseif dim==2
            v1 = CartesianIndex(Modl(pos[1]-1,s[1]), pos[2])
            v2 = CartesianIndex(Modl(pos[1]+1,s[1]), pos[2])
            v3 = CartesianIndex(pos[1],Modl(pos[2]-1,s[2]))
            v4 = CartesianIndex(pos[1],Modl(pos[2]+1,s[2]))
            if printLog
                println([0 v1 0; v3 0 v4 ; 0 v2 0])
            end
            return [v1,v2,v3,v4]
        elseif dim==3
            v1 = CartesianIndex(Modl(pos[1]-1,s[1]), pos[2], pos[3])
            v2 = CartesianIndex(Modl(pos[1]+1,s[1]), pos[2], pos[3])
            v3 = CartesianIndex(pos[1], Modl(pos[2]-1,s[2]), pos[3])
            v4 = CartesianIndex(pos[1], Modl(pos[2]+1,s[2]), pos[3])
            v5 = CartesianIndex(pos[1], pos[2], Modl(pos[3]-1,s[3]))
            v6 = CartesianIndex(pos[1], pos[2], Modl(pos[3]+1,s[3]))
            if printLog
                println([0 v1 0; v3 0 v4 ; 0 v2 0])
                println([v5 0 v6])
            end
            return [v1,v2,v3,v4,v5,v6]
        else
            error("dimension not soported")
        end
    end

    """
        TriangularLatticeNeighbors(latt,pos)

        For a given position pos in a lattice latt, it returns an array of neighbor positions (represented as cartesian indexes) for an triangular lattice geometry. Supports dimensions up to 3.
    """
    function TriangularLatticeNeighbors(latt,pos;printLog=false)
        s=size(latt)
        dim=length(s)
        if dim==1
            v1=CartesianIndex(Modl(pos[1]-1,s[1]))
            v2=CartesianIndex(Modl(pos[1]+1,s[1]))
            if printLog
                println([v1 0 v2])
            end
            return [v1,v2]
        elseif dim==2
            v1 = CartesianIndex(Modl(pos[1]-1,s[1]), pos[2])
            v2 = CartesianIndex(Modl(pos[1]+1,s[1]), pos[2])
            v3 = CartesianIndex(pos[1],Modl(pos[2]-1,s[2]))
            v4 = CartesianIndex(pos[1],Modl(pos[2]+1,s[2]))
            v5 = CartesianIndex(Modl(pos[1]-1,s[1]),Modl(pos[2]+1,s[2]))
            v6 = CartesianIndex(Modl(pos[1]+1,s[1]),Modl(pos[2]+1,s[2]))
            if printLog
                println([0 v1 v5; v3 0 v4 ; 0 v2 v6])
            end
            return [v1,v2,v3,v4,v5,v6]
        elseif dim==3           
            #suposing that triangular 3D is the hexagonal close packing 
            #normal 2d neighbors
            v1 = CartesianIndex(Modl(pos[1]-1,s[1]), pos[2],pos[3])
            v2 = CartesianIndex(Modl(pos[1]+1,s[1]), pos[2],pos[3])
            v3 = CartesianIndex(pos[1],Modl(pos[2]-1,s[2]),pos[3])
            v4 = CartesianIndex(pos[1],Modl(pos[2]+1,s[2]),pos[3])
            v5 = CartesianIndex(Modl(pos[1]-1,s[1]),Modl(pos[2]+1,s[2]),pos[3])
            v6 = CartesianIndex(Modl(pos[1]+1,s[1]),Modl(pos[2]+1,s[2]),pos[3])
            #upper neighbors
            v7 = CartesianIndex(pos[1],pos[2],Modl(pos[3]+1,s[3]))
            v8 = CartesianIndex(Modl(pos[1]-1,s[1]),pos[2],Modl(pos[3]+1,s[3]))
            v9 = CartesianIndex(Modl(pos[1]-1,s[1]),Modl(pos[2]+1,s[2]),Modl(pos[3]+1,s[3]))
            #lower neighbors
            v10 = CartesianIndex(pos[1],pos[2],Modl(pos[3]-1,s[3]))
            v11 = CartesianIndex(Modl(pos[1]-1,s[1]),pos[2],Modl(pos[3]-1,s[3]))
            v12 = CartesianIndex(Modl(pos[1]-1,s[1]),Modl(pos[2]+1,s[2]),Modl(pos[3]-1,s[3]))
            if printLog
                println([0 v1 v5; v3 0 v4 ; 0 v2 v6])
                println([0 v8 v9; 0 v7 0 ; 0 0 0 ])
                println([0 v11 v12; 0 v10 0 ; 0 0 0 ])
            end
            return [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12]
        else
            error("dimension not soported")
        end
    end
    
    """
    HexagonalLatticeNeighbors(latt,pos)

    For a given position pos in a lattice latt, it returns an array of neighbor positions (represented as cartesian indexes) for an hexagonal lattice geometry. Supports dimensions up to 2.

    """
    function HexagonalLatticeNeighbors(latt,pos;printLog=false)
        s=size(latt)
        dim=length(s)
        if dim==1
            v1=CartesianIndex(Modl(pos[1]-1,s[1]))
            v2=CartesianIndex(Modl(pos[1]+1,s[1]))
            if printLog
                println([v1 0 v2])
            end
            return [v1,v2]
        elseif dim==2
            v1 = CartesianIndex(Modl(pos[1]-1,s[1]), pos[2])
            v2 = CartesianIndex(Modl(pos[1]+1,s[1]), pos[2])
            if mod(pos[1],2)==1
                if mod(pos[2],2)==1
                    v3=CartesianIndex(pos[1],Modl(pos[2]+1,s[2]))
                else
                    v3=CartesianIndex(pos[1],Modl(pos[2]-1,s[2]))
                end
            else
                if mod(pos[2],2)==0
                    v3=CartesianIndex(pos[1],Modl(pos[2]+1,s[2]))
                else
                    v3=CartesianIndex(pos[1],Modl(pos[2]-1,s[2]))
                end
            end

            if printLog
                println([0 v1 v5; v3 0 v4 ; 0 v2 v6])
            end
            return [v1,v2,v3]
            #=
        elseif dim==3           
            #suposing that triangular 3D is the hexagonal close packing 
            #normal 2d neighbors
            v1 = CartesianIndex(Modl(pos[1]-1,s[1]), pos[2],pos[3])
            v2 = CartesianIndex(Modl(pos[1]+1,s[1]), pos[2],pos[3])
            v3 = CartesianIndex(pos[1],Modl(pos[2]-1,s[2]),pos[3])
            v4 = CartesianIndex(pos[1],Modl(pos[2]+1,s[2]),pos[3])
            v5 = CartesianIndex(Modl(pos[1]-1,s[1]),Modl(pos[2]+1,s[2]),pos[3])
            v6 = CartesianIndex(Modl(pos[1]+1,s[1]),Modl(pos[2]+1,s[2]),pos[3])
            #upper neighbors
            v7 = CartesianIndex(pos[1],pos[2],Modl(pos[3]+1,s[3]))
            v8 = CartesianIndex(Modl(pos[1]-1,s[1]),pos[2],Modl(pos[3]+1,s[3]))
            v9 = CartesianIndex(Modl(pos[1]-1,s[1]),Modl(pos[2]+1,s[2]),Modl(pos[3]+1,s[3]))
            #lower neighbors
            v10 = CartesianIndex(pos[1],pos[2],Modl(pos[3]-1,s[3]))
            v11 = CartesianIndex(Modl(pos[1]-1,s[1]),pos[2],Modl(pos[3]-1,s[3]))
            v12 = CartesianIndex(Modl(pos[1]-1,s[1]),Modl(pos[2]+1,s[2]),Modl(pos[3]-1,s[3]))
            if printLog
                println([0 v1 v5; v3 0 v4 ; 0 v2 v6])
                println([0 v8 v9; 0 v7 0 ; 0 0 0 ])
                println([0 v11 v12; 0 v10 0 ; 0 0 0 ])
            end
            return [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12]
            =#
        else
            error("dimension not soported")
        end
    end

    """
        Index`Lattice(latt,func)

        Creates a lattice of the same size than latt that in each site i has an array of the indexes j_1 , ... , j_n of the neighbors of i under the neighbor function `func`
    """
    function IndexLattice(latt,func;printLog=false)
        s=size(latt)
        fin=Array{Array{CartesianIndex{length(s)},1},length(s)}(undef,s)
        for pos in CartesianIndices(s)
            fin[pos]=func(latt,pos)
        end
        return fin
    end

    """
        ValueLattice(latt,func)

        Creates a lattice of the same size than latt that in each site i has an array of the values  latt[j_1] , ... , latt[j_n] of the neighbors j_1 , ... , j_n of under the neighbor function `func`
    """
    function ValueLattice(latt,func)
        s=size(latt)
        fin=Array{Array{Int8,1},length(s)}(s)
        for pos in CartesianIndices(s)
            idxs=func(latt,pos)
            for i in idxs
                push!(fin[pos],latt[i])
            end
        end
        return fin
    end
    """
    SumLattice(latt,func)

        Creates a lattice of the same size than latt that in each site i has the value sum_{<i,j>} latt[j] =  latt[j_1] + ... + latt[j_n] of the neighbors j_1 , ... , j_n of under the neighbor function `func`.
    """
    function SumLattice(latt,func)
        s=size(latt)
        fin=Array{Int8,length(s)}(s)
        for pos in CartesianIndices(s)
            idxs=func(latt,pos)
            for i in idxs
                fin[pos]+=latt[i]
            end
        end
        return fin
    end
    
    """
        NeighborSum(latt,neigLatt,pos)

        For position `pos` it returns the sum of its neighbors values, that is sum_{j in neigLatt[pos]} latt[j] .
    """
    function NeighborSum(latt,neigLatt,pos)
        neigVal=0
        for neigPos in neigLatt[pos]
            neigVal+=latt[neigPos]
        end
        return neigVal
    end

    """
        BuildLattices(N,dim,g)

        Builds a N-side lattice of dimension `dim` and its neighbor index lattice with geometry `g`
    """
    function BuildLattices(geoParam,model)
        N=geoParam[1]
        dim=geoParam[2]
        g=geoParam[3]
        
        if ~(g in ["square","triangle","hexagonal"])
            error("Geometry string $g is not supported")
        end
        if model=="normal"
            pick=Array{Int8,1}([1,-1])
        else
            pick=Array{Int8,1}([1,0])
        end
        if dim==1
            latt=rand(pick,N)
        elseif dim==2
            latt=rand(pick,N,N)
        elseif dim==3
            latt=rand(pick,N,N,N)
        else
            error("Dimension not supported")
        end
        if g=="hexagonal"
            neigLatt=IndexLattice(latt,HexagonalLatticeNeighbors)
        elseif g=="triangle"
            neigLatt=IndexLattice(latt,TriangularLatticeNeighbors)
        else
            neigLatt=IndexLattice(latt,SquareLatticeNeighbors)
        end
        return latt , neigLatt
    end
end
