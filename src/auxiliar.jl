module Auxiliar
    
    function Modl(a::Int64,b::Int64)
        x=mod(a,b)
        if x==0
            return b
        else
            return x
        end
    end
    function MirrorList!(l;printLog=false)
        n=length(l)
        aux=[l[i] for i in n:-1:1]
        append!(l,aux)
    end
    function MeanMod(x)
        if length(x)>0
            return mean(x)
        else
            return 0
        end
    end

    function SetValue!(mat,pos,val)
        dim=length(size(mat))
        if dim!=length(pos)
            error("Dimensions must match")
        elseif eltype(mat)!=typeof(val) && eltype(mat)==Any
            error("Types must match")
        elseif dim==1
            mat[pos[1]]=val
        elseif dim==2
            mat[pos[1],pos[2]]=val
        elseif dim==3
            mat[pos[1],pos[2],pos[3]]=val
        else
            error("Dimension not supported")
        end

    end

    function GetValue(mat,pos)
        dim=length(size(mat))
        if dim!=length(pos)
            error("Dimensions must match")
        elseif dim==1
            return mat[pos[1]]
        elseif dim==2
            return mat[pos[1],pos[2]]
        elseif dim==3
            return mat[pos[1],pos[2],pos[3]]
        else
            error("Dimension not supported")
        end
    end

    function ChangeSpin(latt,pos)
        m=copy(latt)
        m[pos]=-1*m[pos]
        return m
    end
    function ChangeSpin!(latt,pos)
        dim=length(size(latt))
        if dim==1
            latt[pos[1]]=-1*latt[pos[1]]
        elseif dim==2
            latt[pos[1],pos[2]]=-1*latt[pos[1],pos[2]]
        elseif dim==3
            latt[pos[1],pos[2],pos[3]]=-1*latt[pos[1],pos[2],pos[3]]
        end
    end
    function RandomPosition(latt)
        s=size(latt)
        dim=length(s)
        if dim==1
            return CartesianIndex{1}(rand(1:s[1]))
        elseif dim==2
            return CartesianIndex{2}(rand(1:s[1]),rand(1:s[2]))
        elseif dim==3
            return CartesianIndex{3}(rand(1:s[1]),rand(1:s[2]),rand(1:s[3]))
        else
            error("Dimension not supported")
        end
    
    end
    function SearchSortedMod(x,a;printLog=false)
        b=searchsortedlast(x,a)
        if printLog
            println("busq")
            println(x)
            println(a)
            println(b)
            println()
        end
            if b==length(x)
            return length(x)-1
        else
            return b
        end
    end
    #función de vecinos inmediatos utilizada para la función de peso del algoritmo metrópolis
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

    function SquareLatticeNeighborsValue(latt,pos;printLog=false)
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

    function SquareLatticeNeighborsIndex(latt,pos;printLog=false)
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

    function NeighborIndexLattice(latt,func;printLog=false)
        s=size(latt)
        fin=Array{Array{CartesianIndex{length(s)},1},length(s)}(s)
        for pos in CartesianRange(s)
            fin[pos]=func(latt,pos)
        end
        return fin
    end
    function NeighborSum(latt,neigLatt,pos)
        neigVal=0
        for neigPos in neigLatt[pos]
            neigVal=neigVal+latt[neigPos]
        end
        return neigVal
    end
    function Index(array,element)
        b=find(x -> x==element,array)
        if ~isempty(b)
            return b[1]
        else
            return -1
        end
    end
    
end