module Lattices
    #using Plots
    #gr()
    """
        PeriodicSquareLatticeNeighbors(latt,pos)

        For a given position pos in a lattice latt, it returns an array of neighbor positions (represented as cartesian indexes) for an square lattice geometry. Supports dimensions up to 3.
    """
    function PeriodicSquareLatticeNeighbors(latt,pos;printLog=false)
        sizes=size(latt)
        dim=length(sizes)
        if dim==1
            v1=CartesianIndex(mod1(pos[1]-1,sizes[1]))
            v2=CartesianIndex(mod1(pos[1]+1,sizes[1]))
            if printLog
                println([v1 0 v2])
            end
            return [v1,v2]
        elseif dim==2
            v1 = CartesianIndex(mod1(pos[1]-1,sizes[1]), pos[2])
            v2 = CartesianIndex(mod1(pos[1]+1,sizes[1]), pos[2])
            v3 = CartesianIndex(pos[1],mod1(pos[2]-1,sizes[2]))
            v4 = CartesianIndex(pos[1],mod1(pos[2]+1,sizes[2]))
            if printLog
                println([0 v1 0; v3 0 v4 ; 0 v2 0])
            end
            return [v1,v2,v3,v4]
        elseif dim==3
            v1 = CartesianIndex(mod1(pos[1]-1,sizes[1]), pos[2], pos[3])
            v2 = CartesianIndex(mod1(pos[1]+1,sizes[1]), pos[2], pos[3])
            v3 = CartesianIndex(pos[1], mod1(pos[2]-1,sizes[2]), pos[3])
            v4 = CartesianIndex(pos[1], mod1(pos[2]+1,sizes[2]), pos[3])
            v5 = CartesianIndex(pos[1], pos[2], mod1(pos[3]-1,sizes[3]))
            v6 = CartesianIndex(pos[1], pos[2], mod1(pos[3]+1,sizes[3]))
            if printLog
                println([0 v1 0; v3 0 v4 ; 0 v2 0])
                println([v5 0 v6])
            end
            return [v1,v2,v3,v4,v5,v6]
        else
            error("dimension not soported")
        end
    end

    
    function SquareLatticeNeighbors(latt,pos;printLog=false)
        sizes=size(latt)
        dim=length(sizes)
        if dim==1
            neigs = Array{CartesianIndex{1},1}()
            if pos[1] > 1
                push!(neigs,CartesianIndex(pos[1]-1))
            end
            if pos[1] < sizes[1]
                push!(neigs,CartesianIndex(pos[1]+1))
            end
            if printLog
                println(neigs)
            end
            return neigs
        elseif dim==2
            neigs = Array{CartesianIndex{2},1}()
            if pos[1] > 1
                push!(neigs,CartesianIndex(pos[1]-1, pos[2]))
            end
            if pos[1] < sizes[1]
                push!(neigs,CartesianIndex(pos[1]+1, pos[2]))
            end
            if pos[2] > 1
                push!(neigs,CartesianIndex(pos[1], pos[2]-1))
            end
            if pos[2] < sizes[2]
                push!(neigs,CartesianIndex(pos[1], pos[2]+1))
            end
            if printLog
                println(neigs)
            end
            return neigs
        elseif dim==3
            neigs = Array{CartesianIndex{3},1}()
            if pos[1] > 1
                push!(neigs,CartesianIndex(pos[1]-1, pos[2], pos[3]))
            end
            if pos[1] < sizes[1]
                push!(neigs,CartesianIndex(pos[1]+1, pos[2], pos[3]))
            end
            if pos[2] > 1
                push!(neigs,CartesianIndex(pos[1], pos[2]-1, pos[3]))
            end                
            if pos[2] < sizes[2]
                push!(neigs,CartesianIndex(pos[1], pos[2]+1, pos[3]))
            end 
            if pos[3] > 1
                push!(neigs,CartesianIndex(pos[1], pos[2], pos[3]-1)) 
            end
            if pos[3] < sizes[3]
                push!(neigs,CartesianIndex(pos[1], pos[2], pos[3]+1))
            end
            if printLog
                println(neigs)
            end
            return neigs
        else
            error("dimension not soported")
        end
    end

    
    """
        PeriodicTriangularLatticeNeighbors(latt,pos)

        For a given position pos in a lattice latt, it returns an array of neighbor positions (represented as cartesian indexes) for an triangular lattice geometry. Supports dimensions up to 3.
    """
    function PeriodicTriangularLatticeNeighbors(latt,pos;printLog=false)
        sizes=size(latt)
        dim=length(sizes)
        if dim==1
            v1=CartesianIndex(mod1(pos[1]-1,sizes[1]))
            v2=CartesianIndex(mod1(pos[1]+1,sizes[1]))
            if printLog
                println([v1 0 v2])
            end
            return [v1,v2]
        elseif dim==2
            v1 = CartesianIndex(mod1(pos[1]-1,sizes[1]), pos[2])
            v2 = CartesianIndex(mod1(pos[1]+1,sizes[1]), pos[2])
            v3 = CartesianIndex(pos[1],mod1(pos[2]-1,sizes[2]))
            v4 = CartesianIndex(pos[1],mod1(pos[2]+1,sizes[2]))
            if mod(pos[2],2) == 0
                   v5 = CartesianIndex(mod1(pos[1]+1,sizes[1]),mod1(pos[2]-1,sizes[2]))
                   v6 = CartesianIndex(mod1(pos[1]+1,sizes[1]),mod1(pos[2]+1,sizes[2]))
            else
                   v5 = CartesianIndex(mod1(pos[1]-1,sizes[1]),mod1(pos[2]-1,sizes[2]))
                   v6 = CartesianIndex(mod1(pos[1]-1,sizes[1]),mod1(pos[2]+1,sizes[2]))
            end
            if printLog
                println([0 v1 v5; v3 0 v4 ; 0 v2 v6])
            end
            return [v1,v2,v3,v4,v5,v6]
        elseif dim==3           
            #suposing that triangular 3D is the FCC close packing 
            #normal 2d neighbors
            v1 = CartesianIndex(mod1(pos[1]-1,sizes[1]), pos[2],pos[3])
            v2 = CartesianIndex(mod1(pos[1]+1,sizes[1]), pos[2],pos[3])
            v3 = CartesianIndex(pos[1],mod1(pos[2]-1,sizes[2]),pos[3])
            v4 = CartesianIndex(pos[1],mod1(pos[2]+1,sizes[2]),pos[3])
            if mod(pos[2],2) == 0
                v5 = CartesianIndex(mod1(pos[1]+1,sizes[1]),mod1(pos[2]-1,sizes[2]),pos[3])
                v6 = CartesianIndex(mod1(pos[1]+1,sizes[1]),mod1(pos[2]+1,sizes[2]),pos[3])
            else
                v5 = CartesianIndex(mod1(pos[1]-1,sizes[1]),mod1(pos[2]-1,sizes[2]),pos[3])
                v6 = CartesianIndex(mod1(pos[1]-1,sizes[1]),mod1(pos[2]+1,sizes[2]),pos[3])
            end
            #upper neighbors
            v7 = CartesianIndex(pos[1],pos[2],mod1(pos[3]+1,sizes[3]))
            v8 = CartesianIndex(mod1(pos[1]-1,sizes[1]),pos[2],mod1(pos[3]+1,sizes[3]))
            v9 = CartesianIndex(mod1(pos[1]-1,sizes[1]),mod1(pos[2]+1,sizes[2]),mod1(pos[3]+1,sizes[3]))
            #lower neighbors
            v10 = CartesianIndex(pos[1],pos[2],mod1(pos[3]-1,sizes[3]))
            v11 = CartesianIndex(mod1(pos[1]-1,sizes[1]),pos[2],mod1(pos[3]-1,sizes[3]))
            v12 = CartesianIndex(mod1(pos[1]-1,sizes[1]),mod1(pos[2]+1,sizes[2]),mod1(pos[3]-1,sizes[3]))
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
        TriangularLatticeNeighbors(latt,pos)

        For a given position pos in a lattice latt, it returns an array of neighbor positions (represented as cartesian indexes) for an triangular lattice geometry. Supports dimensions up to 3.
    """
    function TriangularLatticeNeighbors(latt,pos;printLog=false)
        sizes=size(latt)
        dim=length(sizes)
        if dim==1
            neigs = Array{CartesianIndex{1},1}()
            if pos[1] > 1
                push!(neigs,CartesianIndex(pos[1]-1))
            end
            if pos[1] < sizes[1]
                push!(neigs,CartesianIndex(pos[1]+1))
            end
            if printLog
                println(neigs)
            end
            return neigs
        elseif dim==2
            neigs = Array{CartesianIndex{2},1}()
            if mod(pos[2],2) == 0
                if pos[1] < sizes[1] && pos[2] > 1
                    push!(neigs,CartesianIndex(pos[1]+1,pos[2]-1))
                end
                if pos[1] < sizes[1] && pos[2] < sizes[2]
                    push!(neigs,CartesianIndex(pos[1]+1,pos[2]+1))
                end
            else
                if pos[1] > 1 && pos[2] > 1
                    push!(neigs,CartesianIndex(pos[1]-1,pos[2]-1))
                end
                if pos[1] > 1 && pos[2] < sizes[2]
                    push!(neigs,CartesianIndex(pos[1]-1,pos[2]+1))
                end
            end
            if pos[1] > 1
                push!(neigs,CartesianIndex(pos[1]-1, pos[2]))
            end
            if pos[2] > 1 
                push!(neigs,CartesianIndex(pos[1],pos[2]-1))
            end
            if pos[1] < sizes[1]
                push!(neigs,CartesianIndex(pos[1]+1, pos[2]))
            end
            if pos[2] < sizes[2]
                push!(neigs,CartesianIndex(pos[1],pos[2]+1))
            end
            if printLog
                println(neigs)
            end
            return neigs
        elseif dim==3           
            neigs = Array{CartesianIndex{3},1}()
            #suposing that triangular 3D is the FCC close packing 
            #normal 2d neighbors
            if mod(pos[2],2) == 0
                if pos[1] < sizes[1] && pos[2] > 1
                    push!(neigs,CartesianIndex(pos[1]+1,pos[2]-1))
                end
                if pos[1] < sizes[1] && pos[2] < sizes[2]
                    push!(neigs,CartesianIndex(pos[1]+1,pos[2]+1))
                end
            else
                if pos[1] > 1 && pos[2] > 1
                    push!(neigs,CartesianIndex(pos[1]-1,pos[2]-1))
                end
                if pos[1] > 1 && pos[2] < sizes[2]
                    push!(neigs,CartesianIndex(pos[1]-1,pos[2]+1))
                end
            end
            if pos[1] > 1
                push!(neigs,CartesianIndex(pos[1]-1, pos[2],pos[3]))
            end
            if pos[1] < sizes[1]
                push!(neigs,CartesianIndex(pos[1]+1, pos[2],pos[3]))
            end
            if pos[2] > 1
                push!(neigs,CartesianIndex(pos[1],pos[2]-1,pos[3]))
            end
            if pos[2] < sizes[2]
                push!(neigs,CartesianIndex(pos[1],pos[2]+1,pos[3]))
            end
            #upper neighbors
            if pos[3] < sizes[3]
                push!(neigs,CartesianIndex(pos[1],pos[2],pos[3]+1))
            end
            if pos[1] > 1 && pos[3] < sizes[3]
                push!(neigs,CartesianIndex(pos[1]-1,pos[2],pos[3]+1))
            end
            if pos[1] > 1 && pos[2] < sizes[2] && pos[3] < sizes[3]
                push!(neigs,CartesianIndex(pos[1]-1,pos[2]+1,pos[3]+1))    
            end
            #lower neighbors
            if pos[3] > 1
                push!(neigs, CartesianIndex(pos[1],pos[2],pos[3]-1))
            end
            if pos[1] > 1 && pos[3] > 1
                push!(neigs, CartesianIndex(pos[1]-1,pos[2],pos[3]-1))
            end
            if pos[1] > 1 && pos[2] < sizes[2] && pos[3] > 1
                push!(neigs, CartesianIndex(pos[1]-1,pos[2]+1,pos[3]-1))
            end
            if printLog
                println(neigs)
            end
            return neigs
        else
            error("dimension not soported")
        end
    end
    
    """
    PeriodicHexagonalLatticeNeighbors(latt,pos)

    For a given position pos in a lattice latt, it returns an array of neighbor positions (represented as cartesian indexes) for an hexagonal (honeycomb) lattice geometry. Supports dimensions up to 3.

    """
    function PeriodicHexagonalLatticeNeighbors(latt,pos;printLog=false)
        sizes=size(latt)
        dim=length(sizes)
        if dim==1
            v1=CartesianIndex(mod1(pos[1]-1,sizes[1]))
            v2=CartesianIndex(mod1(pos[1]+1,sizes[1]))
            if printLog
                println([v1 0 v2])
            end
            return [v1,v2]
        elseif dim==2
            v1 = CartesianIndex(pos[1], mod1(pos[2]-1,sizes[2]))
            v2 = CartesianIndex(pos[1], mod1(pos[2]+1,sizes[2]))
            if mod(pos[1],2)==1
                if mod(pos[2],2)==1
                    v3=CartesianIndex(mod1(pos[1]+1,sizes[1]),pos[2])
                else
                    v3=CartesianIndex(mod1(pos[1]-1,sizes[1]),pos[2])
                end
            else
                if mod(pos[2],2)==1
                    v3=CartesianIndex(mod1(pos[1]-1,sizes[1]),pos[2])
                else
                    v3=CartesianIndex(mod1(pos[1]+1,sizes[1]),pos[2])
                end
            end

            if printLog
                println([v1,v2,v3])
            end
            return [v1,v2,v3]
        elseif dim==3           
            #suposing that hexagonal 3D is layered 2D hexagonal 
            #normal 2d neighbors
            v1 = CartesianIndex(pos[1], mod1(pos[2]-1,sizes[2]), pos[3])
            v2 = CartesianIndex(pos[1], mod1(pos[2]+1,sizes[2]), pos[3])
            if mod(pos[1],2)==1
                if mod(pos[2],2)==1
                    v3=CartesianIndex(mod1(pos[1]+1,sizes[1]),pos[2],pos[3])
                else
                    v3=CartesianIndex(mod1(pos[1]-1,sizes[1]),pos[2],pos[3])
                end
            else
                if mod(pos[2],2)==0
                    v3=CartesianIndex(mod1(pos[1]-1,sizes[1]),pos[2],pos[3])
                else
                    v3=CartesianIndex(mod1(pos[1]+1,sizes[1]),pos[2],pos[3])
                end
            end
            v4 = CartesianIndex(pos[1], pos[2], mod1(pos[3]-1,sizes[3]))
            v5 = CartesianIndex(pos[1], pos[2], mod1(pos[3]+1,sizes[3]))
            if printLog
                println([v1,v2,v3,v4,v5])
            end
            return [v1,v2,v3,v4,v5]
        else
            error("dimension not soported")
        end
    end
    
    """
    HexagonalLatticeNeighbors(latt,pos)

    For a given position pos in a lattice latt, it returns an array of neighbor positions (represented as cartesian indexes) for an hexagonal (honeycomb) lattice geometry. Supports dimensions up to 3.

    """
    function HexagonalLatticeNeighbors(latt,pos;printLog=false)
        sizes=size(latt)
        dim=length(sizes)
        if dim==1
            neigs = Array{CartesianIndex{1},1}()
            if pos[1] > 1
                push!(neigs,CartesianIndex(pos[1]-1))
            end
            if pos[1] < sizes[1]
                push!(neigs,CartesianIndex(pos[1]+1))
            end
            if printLog
                println(neigs)
            end
            return neigs
        elseif dim==2
            neigs = Array{CartesianIndex{2},1}()
            if pos[2] > 1
                push!(neigs,CartesianIndex(pos[1], pos[2]-1))
            end
            if pos[2] < sizes[2]
                push!(neigs,CartesianIndex(pos[1], pos[2]+1))
            end
            if mod(pos[1],2)==1
                if mod(pos[2],2) == 1 && pos[1] < sizes[1]
                    push!(neigs,CartesianIndex(pos[1]+1,pos[2]))
                elseif pos[1] > 1
                    push!(neigs,CartesianIndex(pos[1]-1,pos[2]))
                end
            
            else
                if mod(pos[2],2) == 1 && pos[1] > 1
                    push!(neigs,CartesianIndex(pos[1]-1,pos[2]))
                elseif pos[1] < sizes[1]
                    push!(neigs,CartesianIndex(pos[1]+1,pos[2]))
                end
            end

            if printLog
                println(neigs)
            end
            return neigs
        elseif dim==3           
            neigs = Array{CartesianIndex{3},1}()
            #suposing that hexagonal 3D is layered 2D hexagonal 
            #normal 2d neighbors
            if pos[2] > 2
                push!(neigs,CartesianIndex(pos[1], pos[2]-1, pos[3]))
            end
            if pos[2] < sizes[2]
                push!(neigs,CartesianIndex(pos[1], pos[2]+1, pos[3]))
            end
            if mod(pos[1],2)==1
                if mod(pos[2],2)==1 && pos[1] < sizes[1]
                    push!(neigs,CartesianIndex(pos[1]+1,pos[2],pos[3]))
                elseif pos[1] > 2
                    push!(neigs,CartesianIndex(pos[1]-1,pos[2],pos[3]))
                end
            else
                if mod(pos[2],2)==0 && pos[1] > 2
                    push!(neigs,CartesianIndex(pos[1]-1,pos[2],pos[3]))
                elseif pos[1] < sizes[1]
                    push!(neigs,CartesianIndex(pos[1]+1,pos[2],pos[3]))
                end
            end
            if pos[3] > 1
                v4 = CartesianIndex(pos[1], pos[2], mod1(pos[3]-1,sizes[3]))
            end
            if pos[3] < sizes[3]
                v5 = CartesianIndex(pos[1], pos[2], mod1(pos[3]+1,sizes[3]))
            end
            if printLog
                println(neigs)
            end
            return neigs
        else
            error("dimension not soported")
        end
    end
    
    """
        IndexLattice(latt,func)

        Creates a lattice of the same size than latt that in each site i has an array of the indexes j_1 , ... , j_n of the neighbors of i under the neighbor function `func`
    """
    function IndexLattice(sizes,func;printLog=false)
        fin=Array{Array{CartesianIndex{length(sizes)},1},length(sizes)}(undef,sizes)
        for pos in CartesianIndices(sizes)
            fin[pos]=func(latt,pos)
        end
        return fin
    end

    """
    SquarePositionLattice(size)

    Creates a lattice of the real positions of the points of a square Lattice

    """
    function SquarePositionLattice(sizes)
        d = length(sizes)
        if d == 1
            return collect(1:sizes[1])
        elseif d == 2
            unitVectors = [[1,0],[0,1]]
            posLatt = Array{Array{Float32,1},2}(undef,sizes)
            for pos in CartesianIndices(sizes)
                posLatt[pos] = pos[1]*unitVectors[1] + pos[2]*unitVectors[2]
            end
            return posLatt
        elseif d == 3
            unitVectors = [[1,0,0],[0,1,0],[0,0,1]]
            posLatt = Array{Array{Float32,1},3}(undef,sizes)
            for pos in CartesianIndices(sizes)
                posLatt[pos] = pos[1]*unitVectors[1] + pos[2]*unitVectors[2] + pos[3]*unitVectors[3]
            end
            return posLatt
        else
            error("Dimension not supported")
        end
    end

    """
    TriangularPositionLattice(size)

    Creates a lattice of the real positions of the points of a triangular lattice

    """
    function TriangularPositionLattice(sizes)
        d = length(sizes)
        if d == 1
            return collect(1:sizes[1])
        elseif d == 2
            unitVectors = [[1,0],[cos(pi/3),sin(pi/3)]]
            posLatt = Array{Array{Float32,1},2}(undef,sizes)
            for pos in CartesianIndices(sizes)
                posLatt[pos] = (pos[1]-div(pos[2]-1,2))*unitVectors[1] + pos[2]*unitVectors[2]
                #posLatt[pos] = (pos[1])*unitVectors[1] + pos[2]*unitVectors[2]            
            end
            return posLatt
        elseif d == 3
            v1  = [1,0,0]
            v2 = [cos(pi/3),sin(pi/3),0]
            v3 = [(1+cos(pi/3))/3,sin(pi/3)/3,sqrt(2/3)]
            unitVectors = [v1,v2,v3]
            posLatt = Array{Array{Float32,1},3}(undef,sizes)
            for pos in CartesianIndices(sizes)
                posLatt[pos] = (pos[1]-div(pos[2]-1,2) - div(pos[3]-1,2) )*unitVectors[1] + (pos[2] - div(pos[3]-1,2))*unitVectors[2] + pos[3]*unitVectors[3]
            end
            return posLatt
        else
            error("Dimension not supported")
        end
    end
    """
    TwoDTriangularToHex(pos)

    Auxiliar function to help with the building of hexagonal lattices. As the hexagonal (honeycomb) lattice
    """
    function TwoDTriangularToHex(pos)
        if length(pos) != 2
            error("Dimension must be 2")
        end
        if (mod(pos[1],3) == 1 && mod(pos[2],2) == 1 ) || ( mod(pos[1],3) == 2 && mod(pos[2],2) == 0 )
            return false
        else
            return true
        end
    end
    """
    HexagonalPositionLattice(size)

    Creates a lattice of the real positions of the points of an hexagonal lattice

    """
    function HexagonalPositionLattice(sizes)
        d = length(sizes)
        if d == 1
            newSizes = Tuple(Int16(sizes[1]*3/2))
            posLatt = TriangularPositionLattice(newSizes)
            posLatt = [posLatt[i] for i in 1:length(posLatt) if mod(i,3) !=1]
            return posLatt
        elseif d == 2
            
            newSizes = Tuple((Int16(sizes[1]*3/2),Int16(sizes[2])))
            @show 
            auxLatt = TriangularPositionLattice(newSizes)
            posLatt = Array{Array{Float32,1},2}(undef,sizes)
            posLatt = [auxLatt[indx] for indx in CartesianIndices(size(auxLatt)) if TwoDTriangularToHex(indx)]
            posLatt = reshape(posLatt,sizes)
            return posLatt
        elseif d == 3
            basePosLatt =  HexagonalPositionLattice(sizes[1:2])
            posLatt = Array{Array{Float32,1},3}(undef,sizes)
            for i in 1:sizes[3]
                posLatt[:,:,i] = [vcat(x,[i]) for x in basePosLatt]
            end
            return posLatt
        else
            error("Dimension not supported")
        end  
    end
    #=
    function PlotLattice(sizes,neigFunc,posFunc)
        println("Loading Packages")
        println("Done")
        posLatt = posFunc(sizes)
        #@show posLatt
        d = length(sizes)
        p = Plots.plot(leg = false,aspect_ratio=:equal)
        if d == 1
            for indx in CartesianIndices(sizes)
                pos = posLatt[indx]
                Plots.scatter!([0],[pos],color = :blue)
                neigs = neigFunc(posLatt,indx)
                for neig in neigs
                    Plots.plot!([0,0],[pos,posLatt[neig]],color = :black)
                end
            end
        elseif d == 2
            for indx in CartesianIndices(sizes)
                pos = posLatt[indx]
                Plots.scatter!([pos[1]],[pos[2]],color = :blue)
                neigs = neigFunc(posLatt,indx)
                for neig in neigs
                    Plots.plot!([pos[1],posLatt[neig][1]],[pos[2],posLatt[neig][2]],color = :black)
                end
            end
        elseif d == 3
            for indx in CartesianIndices(sizes)
                pos = posLatt[indx]
                Plots.scatter!([pos[1]],[pos[2]],[pos[3]],color = :blue)
                neigs = neigFunc(posLatt,indx)
                for neig in neigs
                    @show neig
                    Plots.plot!([pos[1],posLatt[neig][1]],[pos[2],posLatt[neig][2]],[pos[3],posLatt[neig][3]],color = :black)
                end
            end
        end
        Plots.savefig("../lattice.png")
        #display(p)
    end
    =#
end
    #=
    """
        SquareLatticeNeighbors(latt,pos)

        For a given position pos in a lattice latt, it returns an array of neighbor positions (represented as cartesian indexes) for an square lattice geometry. Supports dimensions up to 2.
    """
    function SquareLatticeNeighbors(latt,pos;printLog=false)
        sizes=size(latt)
        dim=length(sizes)
        if dim==1
            if pos[1] == 1
                return [CartesianIndex(2)]
            elseif pos[1] == sizes[1]
                return [CartesianIndex(sizes[1]-1)]
            else
                v1=CartesianIndex(mod1(pos[1]-1,sizes[1]))
                v2=CartesianIndex(mod1(pos[1]+1,sizes[1]))
                if printLog
                    println([v1 0 v2])
                end
                return [v1,v2]
            end
        elseif dim==2
            #four corners
            if pos[1] == 1 && pos[2] == 1
                v1 = CartesianIndex(1,2)
                v2 = CartesianIndex(2,1)
                return [v1,v2]
            elseif pos[1] == sizes[1] && pos[2] == 1
                v1 = CartesianIndex(sizes[1],2)
                v2 = CartesianIndex(sizes[1]-1,1)
                return [v1,v2]
            elseif pos[1] == 1 && pos[2] == sizes[2]
                v1 = CartesianIndex(2,sizes[2])
                v2 = CartesianIndex(1,sizes[2]-1)
                return [v1,v2]
            elseif pos[1] == sizes[1] && pos[2] == sizes[2]
                v1 = CartesianIndex(sizes[1]-1,sizes[2])
                v2 = CartesianIndex(sizes[1],sizes[2]-1)
                return [v1,v2]
            #four walls
            elseif pos[1] == 1
                v1 = CartesianIndex(1,pos[2]-1)
                v2 = CartesianIndex(1,pos[2]+1)
                v3 = CartesianIndex(2,pos[2])
                return [v1,v2,v3]
            elseif pos[2] == 1
                v1 = CartesianIndex(pos[1]-1,1)
                v2 = CartesianIndex(pos[1]+1,1)
                v3 = CartesianIndex(pos[1],2)
                return [v1,v2,v3]
            elseif pos[1] == sizes[1]
                v1 = CartesianIndex(sizes[1],pos[2]-1)
                v2 = CartesianIndex(sizes[1],pos[2]+1)
                v3 = CartesianIndex(sizes[1]-1,pos[2])
                return [v1,v2,v3]
            elseif pos[2] == sizes[2]
                v1 = CartesianIndex(pos[1]-1,sizes[2])
                v2 = CartesianIndex(pos[1]+1,sizes[2])
                v3 = CartesianIndex(pos[1],sizes[2]-1)
                return [v1,v2,v3]
            #inside the lattice
            else
                v1 = CartesianIndex(pos[1]-1, pos[2])
                v2 = CartesianIndex(pos[1]+1, pos[2])
                v3 = CartesianIndex(pos[1],pos[2]-1)
                v4 = CartesianIndex(pos[1],pos[2]+1)
                return [v1,v2,v3,v4]
            end
        
        elseif dim==3
            #eight corners
            if pos[1] == 1 && pos[2] == 1 && pos[3] == 1
                v1 = CartesianIndex(2,1,1)
                v2 =  CartesianIndex(1,2,1)
                v3 =  CartesianIndex(1,1,2)
                return [v1,v2,v3]
            elseif pos[1] == 1 && pos[2] == 1 && pos[3] == sizes[3]
                v1 = CartesianIndex(2,1,1)
                v2 =  CartesianIndex(1,2,1)
                v3 =  CartesianIndex(1,1,sizes[3]-1)
                return [v1,v2,v3]
            elseif pos[1] == 1 && pos[2] == sizes[2] && pos[3] == 1
                v1 = CartesianIndex(2,1,1)
                v2 =  CartesianIndex(1,sizes[2]-1,1)
                v3 =  CartesianIndex(1,1,2)
                return [v1,v2,v3]
            elseif pos[1] == 1 && pos[2] == sizes[2] && pos[3] == sizes[3]
                v1 = CartesianIndex(2,1,1)
                v2 =  CartesianIndex(1,sizes[2]-1,1)
                v3 =  CartesianIndex(1,1,sizes[3]-1)
                return [v1,v2,v3]
            elseif pos[1] == sizes[1] && pos[2] == 1 && pos[3] == 1
                v1 = CartesianIndex(sizes[1]-1,1,1)
                v2 =  CartesianIndex(1,2,1)
                v3 =  CartesianIndex(1,1,2)
                return [v1,v2,v3]
            elseif pos[1] == sizes[1] && pos[2] == 1 && pos[3] == sizes[3]
                v1 = CartesianIndex(sizes[1]-1,1,1)
                v2 =  CartesianIndex(1,2,1)
                v3 =  CartesianIndex(1,1,sizes[3]-1)
                return [v1,v2,v3]
            elseif pos[1] == sizes[1] && pos[2] == sizes[2] && pos[3] == 1
                v1 = CartesianIndex(sizes[1]-1,1,1)
                v2 =  CartesianIndex(1,sizes[2]-1,1)
                v3 =  CartesianIndex(1,1,2)
                return [v1,v2,v3]
            elseif pos[1] == sizes[1] && pos[2] == sizes[2] && pos[3] == sizes[3]
                v1 = CartesianIndex(sizes[1]-1,1,1)
                v2 =  CartesianIndex(1,sizes[2]-1,1)
                v3 =  CartesianIndex(1,1,sizes[3]-1)
                return [v1,v2,v3]
            #twelve edges
            elseif pos[1] == 1 && pos[2] == 1
                v1 = CartesianIndex(2,1,1)
                v2 =  CartesianIndex(1,2,1)
                v3 =  CartesianIndex(1,1,pos[3]-1)
                v4 =  CartesianIndex(1,1,pos[3]+1)
                return [v1,v2,v3,v4]
            elseif pos[2] == 1 && pos[3] == 1
                v1 = CartesianIndex(pos[1]+1,1,1)
                v2 =  CartesianIndex(pos[1]-1,1,1)
                v3 =  CartesianIndex(1,2,1)
                v4 =  CartesianIndex(1,1,2)
                return [v1,v2,v3,v4]
            elseif pos[1] == 1 && pos[3] == 1
                v1 = CartesianIndex(1,pos[2]+1,1)
                v2 =  CartesianIndex(1,pos[2]-1,1)
                v3 =  CartesianIndex(1,1,1)
                v4 =  CartesianIndex(1,1,1)
                return [v1,v2,v3,v4]
            elseif pos[1] == 1 && pos[3] == sizes[3]
                v1 = CartesianIndex(1,1,sizes[3]-1)
                v2 =  CartesianIndex(2,1,1)
                v3 =  CartesianIndex(1,pos[2]+1,1)
                v4 =  CartesianIndex(1,pos[2]-1,1)
                return [v1,v2,v3,v4]
            elseif pos[2] == 1 && pos[3] == sizes[3]
                v1 = CartesianIndex(1,1,sizes[3]-1)
                v2 =  CartesianIndex(1,2,1)
                v3 =  CartesianIndex(pos[1]+1,1,1)
                v4 =  CartesianIndex(pos[1]-1,1,1)
                return [v1,v2,v3,v4]
            elseif pos[1] == 1 && pos[2] == sizes[2]
                v1 = CartesianIndex(1,sizes[2]-1,1)
                v2 =  CartesianIndex(2,1,1)
                v3 =  CartesianIndex(1,1,pos[3]+1)
                v4 =  CartesianIndex(1,1,pos[3]-1)
                return [v1,v2,v3,v4]
            elseif pos[2] == sizes[2] && pos[3] == 1
                v1 = CartesianIndex(1,sizes[2]-1,1)
                v2 =  CartesianIndex(1,1,2)
                v3 =  CartesianIndex(1,pos[2]+1,1)
                v4 =  CartesianIndex(1,pos[2]-1,1)
                return [v1,v2,v3,v4]
            elseif pos[1] == sizes[1] && pos[2] == 1
                v1 = CartesianIndex(1,2,1)
                v2 =  CartesianIndex(sizes[1]-1,1,1)
                v3 =  CartesianIndex(1,1,pos[3]+1)
                v4 =  CartesianIndex(1,1,pos[3]-1)
                return [v1,v2,v3,v4]
            elseif pos[1] == sizes[1] && pos[3] == 1
            elseif pos[1] == sizes[1] && pos[2] == sizes[2]
            elseif pos[1] == sizes[1] && pos[3] == sizes[3]
            elseif pos[2] == sizes[2] && pos[3] == sizes[3]
            #six faces
            elseif pos[1] == 1 
            elseif pos[2] == 1 
            elseif pos[3] == 1 
            elseif pos[1] == sizes[1] 
            elseif pos[2] == sizes[2] 
            elseif pos[3] == sizes[3] 
            else
                v1 = CartesianIndex(mod1(pos[1]-1,sizes[1]), pos[2], pos[3])
                v2 = CartesianIndex(mod1(pos[1]+1,sizes[1]), pos[2], pos[3])
                v3 = CartesianIndex(pos[1], mod1(pos[2]-1,sizes[2]), pos[3])
                v4 = CartesianIndex(pos[1], mod1(pos[2]+1,sizes[2]), pos[3])
                v5 = CartesianIndex(pos[1], pos[2], mod1(pos[3]-1,sizes[3]))
                v6 = CartesianIndex(pos[1], pos[2], mod1(pos[3]+1,sizes[3]))
                if printLog
                    println([0 v1 0; v3 0 v4 ; 0 v2 0])
                    println([v5 0 v6])
                end
                return [v1,v2,v3,v4,v5,v6]
            end
        else
            error("dimension not soported")
        end
    end
    =#