module CyclesRefactored

    #include("geometry.jl")
    #include("cyclesAta.jl")

    """
        Vectorize(mat) 
    
        Converts an n-dimensional array `mat` in a vector.
    """
    function Vectorize(mat)
        return reshape(copy(mat),length(mat))
    end

    """
        Matrize(vec,dimsVec) 
    
        Converts a vector of length `prod(dimsVec)` into an length  (dimsVec)-array with dimensions dimsVec[1],..,dimsVec[n]
    """
    function Matrize(vec,dimsVec)
        if prod(dimsVec) != length(vec)
            error("dimentions must match")
        end
        dimsVec=Tuple(dimsVec)
        return reshape(copy(vec),dimsVec)
    end

    """
        MakeLists(latt,neigLatt)

        Converts latt and neigLatt into one-dimensional arrays
    """
    function MakeLists(latt, neigLatt)
        newLatt = Array{Int8,1}(Vectorize(latt))
        aux1 = LinearIndices(latt)
        aux2 = Vectorize(neigLatt)
        newNeigLatt = Array{Array{Int16,1},1}(length(newLatt))
        for i in 1:length(newLatt)
            newNeigLatt[i] = [aux1[pos] for pos in aux2[i]]
        end
        return newLatt, newNeigLatt
    end

    """
        SquareLattice(n,m) 

        generates a list of neighbors of a square n x m network with    periodic boundary conditions
    """
    function SquareLatticePer(n,m)  
        k=n*m
        ve=Array[]
        for j in 1:k
            ve1 = [j%n+1+÷(j-1, n)*n; mod1(j-1,n)+÷(j-1,n)*n;mod1(j-n,n*m);mod1(j+n,n*m)]
            push!(ve,ve1)
        end
        return ve
    end
    """
        SquareLatticeFix(n,m) 

        Generates a list of neighbors for a square n x m network with   fixed boundary conditions"""
    function SquareLatticeFix(n,m) 
        k=n*m
        ve=Array[]
        push!(ve,[2,1+n])
        for j in 2:k-1
            if j == n
                push!(ve,[n-1,2n])
            elseif j == n*m-n+1
                push!(ve,[n*m-n+2,n*m-2n+1])
            elseif j< n
                push!(ve, [j-1,j+1,j+n])
            elseif j> n*m-n+1
                push!(ve, [j-1,j+1,j-n])
            elseif mod(j,n) == 1 && j<n*m-n+1
                push!(ve, [j+1, j+n, j-n])
            elseif mod(j,n) == 0 && j>n
                push!(ve, [j-1, j+n, j-n])
            else
                ve1 = [j%n+1+÷(j-1, n)*n; mod1(j-1,n)+÷(j-1,n)*n;mod1(j-n,n*m);mod1(j+n,n*m)]
                push!(ve,ve1)
            end
        end
        push!(ve,[n*m-1,n*(m-1)])
        return ve
    end
    """
        CycNeig(latt,neigLatt)

        Generates a list of neighbors from a lattice `latt` and neigbor connectivity `neigLatt`. For example, CycNeig(latt,neigLatt)[i] will be an array that cointains all j's such that latt[j] = 1 and j is in neigLatt[i]
    """
    function CycNeig(latt,neigLatt)  
        newNeigLatt = Array[]
        n = length(latt)
        for i in 1:n
            ve = []
            if latt[i] == 1
                for j in neigLatt[i]
                    if latt[j] == 1
                        push!(ve,j)
                    end
                end
            end
            push!(newNeigLatt,ve)
        end
        return newNeigLatt
    end
    
    """
        RemoveTails!(latt,neigLatt,n,m)

        Removes recursively all elements with connectivity less or equal than 1 from the lattice.
    """
    function RemoveTails!(latt,neigLatt,n,m)
        tot = n*m
        σ2 = copy(latt)
        bool1 = true
        σn = zeros(tot)
        while bool1 
            tailsTest = 1
            for i in 1:tot
                if latt[i] != 0
                    links = sum(latt[neigLatt[i]])
                    if links <= 1  
                        latt[i] = 0
                        σn[i] = 1
                    end
                    tailsTest *= (1-links)
                    if abs(tailsTest) > 2
                        tailsTest = 1
                    end
                end
            end
            if tailsTest == 0
                bool1 = true
            else
                bool1 = false
            end
        end
        return latt, σn
    end
    function quita_pelos!(σ,vecino,n1,m)
        n = n1*m
        σ2 = copy(σ)
        test = true
        σn = zeros(n)
        while test 
            test2 = 1
            for i in 1:n
                if σ[i] != 0
                    enlaces = sum(σ[vecino[i]])
                    if enlaces <= 1  
                        σ[i] = 0
                        σn[i] = 1
                    end
                    test2 *=(1-enlaces)
                    if abs(test2) >2
                        test2 = 1
                    end
                end
            end
            if test2 == 0
                test = true
            else
                test = false
            end
        end
        return σ, σn
    end
    """
        ListOddDeg(neigLatt)

        Returns a list of all elements in `neigLatt` with   odd degree higher than 3
    """
    function ListOddDeg(neigLatt)
        n = length(neigLatt)
        res = []
        for i in 1:n
            if length(neigLatt[i]) == 3
                push!(res, i)
            end
        end
        return res
    end
    """
    Walk(neigLatt)

    Walks the graph defined by `neigLatt` starting from a vertex with odd degree, at random until it finds a cycle. If it bumps into another point with odd degree, it makes special considereations.
    """
    function Walk(neigLatt)
        list1 = ListOddDeg(neigLatt)
        if length(list1) == 0
            return [], neigLatt, true, false
        end
        newNeigLatt = Array[]
        for i in 1: length(neigLatt)
            push!(newNeigLatt, copy(neigLatt[i]))
        end
        i = rand(list1)
        j = rand(neigLatt[i])
        #points that are being traveled
        list2 = [i] 
        #points of list2 with degree higher than two
        list3 = [1] 
        # list of elements in a cycle. These are the lements that are elimitaed in the links of their neighbors
        list4 = [] 
        #iterating while searching for a cycle
        while j ∉ list2 
            push!(list2,j)
            n = 1
            if length(neigLatt[j]) > 2
                push!(list3, length(list2))
            end
            k = neigLatt[j][n]
            #making sure that it does not come from where it goes
            while k == i  # tiene que elegir un camino que no sea regresar. 
                n += 1
                if n > length(neigLatt[j])
                    @show "error", j, n, length(neigLatt[j]), list4, newNeigLatt, true, true
                    error("no new point to move")
                end
                k = neigLatt[j][n]
            end
            i = j
            j = k
        end
        # boolean to see if there is an point with only one neighbor
        bool = false  
        k = list2[list3[end]]
        #case 1: cycle is formed without any point with conectivity larger than two in between
        if length(list3) == 1 
            s = findall((in)(list2[end]),newNeigLatt[j])
            deleteat!(newNeigLatt[j],s[1])
            s = findall((in)(list2[2]),newNeigLatt[j])
            deleteat!(newNeigLatt[j],s[1])
            if length(newNeigLatt[j]) == 1
                bool = true
            end
            for i in 2:length(list2)
                push!(list4,list2[i])
                newNeigLatt[list2[i]] = []
            end
            push!(list4,list2[1])
            return list4, newNeigLatt, bool, true
        end
        # case 2: cycle is formed with the last element in list3
        if list3[end] == length(list2) 
            s = findall((in)(list2[end]),newNeigLatt[j])
            deleteat!(newNeigLatt[j],s[1])
            s = findall((in)(j),newNeigLatt[k])
            deleteat!(newNeigLatt[k],s[1])
            if length(newNeigLatt[j]) == 1
                bool = true
            end
            x = findall((in)(j),list2)
            for i in x[1]:length(list2)
                push!(list4,list2[i])
            end
            return list4, newNeigLatt, bool, true
        else  # case 3: any other case
            s = findall((in)(list2[end]),newNeigLatt[j])
            deleteat!(newNeigLatt[j],s[1])
            s = findall((in)(list2[list3[end]+1]),newNeigLatt[k])
            deleteat!(newNeigLatt[k],s[1])
            x = findall((in)(j),list2)
            for i in list3[end]+1:length(list2)
                newNeigLatt[list2[i]] = []
            end
            if length(newNeigLatt[j]) == 1
                bool = true
            end
            if length(newNeigLatt[k]) == 1
                bool = true
            end
            for i in x[1]:length(list2)
                push!(list4,list2[i])
            end
            return list4, newNeigLatt, bool, true
        end
    end

    "sigma(vecino) generates an array `latt` a from neighbor lattice `neigLatt`"
    function NewLatt(neigLatt)
        latt = Int[]
        n = length(neigLatt)
        for i in 1:n
            x = 0
            if length(neigLatt[i])>1
                x = 1
            end
            push!(latt,x)
        end
        return latt
    end

    "AddFromList!(list, latt) 
    
    adds to `latt` the elements from `list`.
    "
    function AddFromList!(list, latt)
        for i in list
            latt[i] = 1
        end
        return latt
    end  
 
  
    "FindCycles(latt,neigLatt,n1,m)
    
    Returns three arrays, one with cycles, the original array and the diference between the two
    "
    function FindCycles(latt,neigLatt,n,m)
        σo = copy(latt)
        latt, σn = RemoveTails!(latt,neigLatt,n,m)
        newNeigLatt = CycNeig(latt,neigLatt)
        list, newNeigLatt, test, test2 = Walk(newNeigLatt)
        σ2 = zeros(Int,n*m)   
        cont = 0
        while test2   
            latt = NewLatt(newNeigLatt)
            σ2 = AddFromList!(list, σ2)
            if test
                cont += 1
                σcop = copy(latt)
                latt, σn1 = RemoveTails!(latt,neigLatt,n,m)
                σn += σn1
                newNeigLatt = CycNeig(latt,neigLatt)
            end
            list, newNeigLatt, test, test2 = Walk(newNeigLatt)
        end
        latt, σn1 = RemoveTails!(latt,neigLatt,n,m)
        σn += σn1
        σc = latt+σ2
        for i in 1:length(σc)
            if σc[i]>1
                σc[i] =1
            end
        end
        return σo, σc, σo-σc
    end
    "Cycles(latt)
    
    Returns matrix of cycles for a square lattice with periodic boundary conditions
    "
    function Cycles(latt)
        n=size(latt)[1]
        ##println("hay ciclos")
        a=Vectorize(latt); 
        b=SquareLatticePer(n,n)
        vec=FindCycles(a,b,n,n)[2]
        return Matrize(vec,(n,n))
    end
    #=
    function TestSimilarity(nTests=100)

        c = 0
        t1 = 0
        t2 = 0
        for i in 1:nTests
            println(i)
            latt , neigLatt =  Geometry.BuildLattices( [100,2,"square"], "cycle")
            a,t1aux =  (@timed Cycles(latt) )[1:2]
            b,t2aux = (@timed CyclesAta.ciclos2(latt) )[1:2]
            t1 += t1aux
            t2 += t2aux
            if a != b
                c += 1
                @show a
                @show b
                println("fail")
            end
        end
        println("Performed $(nTests) tests, $(c) failed")
        println(t1)
        println(t2)
    end
    
    TestSimilarity(2)
    TestSimilarity()
    =#
end

