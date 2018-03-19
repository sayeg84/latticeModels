module Auxiliar
    
    function Modl(a::Int64,b::Int64)
        x=mod(a,b)
        if x==0
            return b
        else
            return x
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
    
    function SearchSortedMod(x,a)
        b=searchsortedlast(x,a)
        if b==length(x)
            return length(x)-1
        else
            return b
        end
    end
    #función de vecinos inmediatos utilizada para la función de peso del algoritmo metrópolis
    function SquareLatticeNeighbors(latt,pos;test=false)
        s=size(latt)
        dim=length(s)
        if dim==1
            v1=latt[Modl(pos[1]-1,s[1])]
            v2=latt[Modl(pos[1]-2,s[1])]
            if test
                println([v1 0 v2])
            end
            return v1+v2
        elseif dim==2
            v1 = latt[Modl(pos[1]-1,s[1]), pos[2]]
            v2 = latt[Modl(pos[1]+1,s[1]), pos[2]]
            v3 = latt[pos[1],Modl(pos[2]-1,s[2])]
            v4 = latt[pos[1],Modl(pos[2]+1,s[2])]
            if test
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
            if test
                println([0 v1 0; v3 0 v4 ; 0 v2 0])
                println([v5 0 v6])
            end
            return v1+v2+v3+v4+v5+v6
        else
            error("dimension not soported")
        end
    end
end