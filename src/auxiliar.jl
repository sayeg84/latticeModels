module Auxiliar
    using Statistics

    """
        Modl(a,b)

        Function equivalent to the traditional modulo n ( % n , mod(x,n)) operator but returning values between 1 and n.
    """
    function Modl(a::Int64,b::Int64)
        x=mod(a,b)
        if x==0 
            return b
        else
            return x
        end
    end

    """
        MirrorList!(l)

        Liven a list [a_1,...a_n] it "mirrors" it converting it in [a_1,...,a_n,a_n,a_{n-1},...,a_1]  

    """
    function MirrorList!(l;printLog=false)
        n=length(l)
        aux=[l[i] for i in n:-1:1]
        append!(l,aux)
    end

    """
        MeanMod(x)

        Modified mean function that returns 0 for empty lists
    """
    function MeanMod(x)
        if length(x)>0
            return Statistics.mean(x)
        else
            return 0.0
        end
    end



    ########################################################################
    #legacy functions. Not used anymore
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

    ########################################################################

    """
        ChangeSpin(lattPos)

        Wrapper that changes the spin in that given position. It can alter between the normal ising model and our modified version
    """
    function ChangeSpin(latt,pos;normal=true)
        m=copy(latt)
        if normal
            m[pos]=-1*m[pos]
        else
            if m[pos]==0
                m[pos]=1
            else
                m[pos]=0
            end
        end
        return m
    end
    #############################################################################################
    #=Legacy function
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
    =#

    """
        RandomPosition(latt)

        Wrapper for returning a CartesianIndex corresponding to a position in the n-dimensional array latt
    """
    function RandomPosition(latt)
        return rand(CartesianIndices(size(latt)))
    end

    """
        SearchSortedMod(x,a)

        Modified version of searchsortedlast to perform binary search in array x in search of largest index i such that x[i] < a and returning n-1 when  a > x[end]
    """
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
    
    """
        Index(x,a)
        function that searches for the index of item a in array x. Returns -1 if a is not in x
    """
    function Index(array,element)
        b=find(x -> x==element,array)
        if ~isempty(b)
            return b[1]
        else
            return -1
        end
    end
end