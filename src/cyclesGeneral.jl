module CyclesGeneral

    "vectorize(mat) converts an n-dimensional array mat in a vector"
    function vectorize(mat)
        return reshape(mat,length(mat))
    end

    "matrize(vec,dimsVec) converts a vector of length prod(dimsVec) into an length(dimsVec)-array with dimensions dimsVec[1],..,dimsVec[n]"
    function matrize(vec,dimsVec)
        if prod(dimsVec) != length(vec)
            error("dimentions must match")
        end
        dimsVec=Tuple(dimsVec)
        return reshape(vec,dimsVec)
    end

    "vecinos(n,m) genera una lista de vecinos de una red cuadrada de nxm, con condiciones periódicas a la frontera"
    function vecinos(n,m)  
        k=n*m
        ve=Array[]
        for j in 1:k
            ve1 = [j%n+1+÷(j-1, n)*n; mod1(j-1,n)+÷(j-1,n)*n;mod1(j-n,n*m);mod1(j+n,n*m)]
            push!(ve,ve1)
        end
        return ve
    end

    "vecinos2(n,m) genera una lista de vecinos de una red cuadrada de nxm, con condiciones fijas a la frontera"
    function vecinos2(n,m) 
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

    "dibuja_sigma(σ,n,m,estilo) Dibuja un mapa de calor de σ, en un red de nxm. Estilo puede ser igual a false o :ijulia "
    function dibuja_sigma(σ,n,m,estilo) 
        σ =reshape(σ,n,m)
        heatmap!(1:n,1:m,σ,aspect_ratio=:equal,  show = estilo)
    end

    "vecinos_s(σ,vecino) Genera una lista de vecinos a partir de un arreglo σ y la conectividad de la red donde se encuentra σ. Ejemplo, v = vecinos_s(σ,vecino), entonces, v[i] será un arreglo que contenga todos j's tales que σ[j] = 1, y j un elemento de vecino[i]"
    function vecinos_s(σ,vecino)  
        vecino2 = Array[]
        n = length(σ)
        for i in 1:n
            ve = []
            if σ[i] == 1
                for j in vecino[i]
                    if σ[j] == 1
                        push!(ve,j)
                    end
                end
            end
            push!(vecino2,ve)
        end
        return vecino2
    end

    "CycleNeigs(latt,neigLatt)

    returns lattice of neighbors of elements that are ones and their neighbors that are one as well. For elements that are zero it returns nothing
    "
    function CycleNeigs(latt,neigLatt)
        res=copy(neigLatt)
        for pos in CartesianIndices(size(neigLatt))
            if latt[pos]==1
                res[pos] = [p for p in neigLatt[pos] if latt[p]==1]
            else
                res[pos]=[] 
            end
        end
        return res
    end

    function RemoveTails(latt,neigLatt)
        res = copy(latt)
        tails = zeros(Int8,size(latt))
        bool = true
        while bool
            flag=1
            for pos in CartesianIndices(size(latt))
                if res[pos]!=0
                    deg = sum(res[neigLatt[pos]])
                    if deg <=1
                        res[pos] = 0
                        tails[pos] = 1
                    end
                    flag *= (1-deg)
                    if abs(flag) > 2
                        flag = 1
                    end
                end
            end 
            if flag == 0
                bool = true
            else
                bool = false
            end
        end
        return res , tails
    end


    "quita_pelos!(σ,vecino,n,m) quita del arreglo σ, los elementos cuya conectividad sea menor que 2. Esto lo hace de forma recurrente, hasta que ningún elemento de σ tenga una conectividad menor que 2. n y m representan los tamaños de la caja donde se encuentra el arreglo σ. vecino es un arreglo que contiene la conectividad de la cuadrada. Regresa los valores de σ y un arreglo con los elementos que se quitaron (los pelos). También hay la versión con quita_pelos!(σ,vecino, vecino2,n,m), en cuyo caso tomará en cuenta los enlaces."
    function quita_pelos!(σ,vecino,n1,m)
        n = n1*m
        σ2 = copy(σ)
        test = true
        contador = 0
        σn = zeros(n)
        while test 
            contador +=1
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

    function quita_pelos!(σ,vecino, vecino2,n1,m)
        n = n1*m
        σ2 = copy(σ)
        test = true
        contador = 0
        σn = zeros(n)
        while test
            contador +=1
            test2 = 1
            for i in 1:n
                if σ[i] != 0
                    enlaces = sum(σ[vecino[i]])
                    enlaces2 = length(vecino2[i])
                    if enlaces2 == 1
                        @show enlaces, enlaces2
                        enlaces == enlaces2
                    end
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
                    
    "Lista(vecino), usa un arreglo de conectividad para hacer una lista de los elementos que tienen conectividad igual a 3"
    function Lista(vecino)
        n = length(vecino)
        lista = []
        for i in 1:n
            if length(vecino[i]) == 3
                push!(lista, i)
            end
        end
        return lista
    end

    function listDeg3(latt,neigLatt)
        res = []
        for pos in CartesianIndices(size(latt))
            if length(neigLatt[pos])==3
                push!(res,pos)
            end
        end
        return res
    end

    "recorrer(vecino): dado una arreglo de conectividad (donde todos los elementos tienen conectividad mayor a 1), elige aleatoriamente entre los puntos que tenga conectividad igual a 3, para recorrer el grafo de forma aleatoria, hasta que encuentre un ciclo. Si no hay ningún punto con conectivdad 3, regresa un arreglo vacío, la misma lista de vecinos, true, false. En otro caso, regresa una lista de elementos que están seguro en un ciclo, el arreglo vecino2 descartando los elementos de la lista, un test que es true si algún elemento de vecino2 tiene medida 1 y false si no, y un valor false, que significa que sí hubo puntos con conectividad 3."
    function recorrer(vecino)
        lista = Lista(vecino)
        if length(lista) == 0
            return [], vecino, true, false
        end
        vecino2 = Array[]
        for i in 1: length(vecino)
            push!(vecino2, copy(vecino[i]))
        end
        i = rand(lista)
        j = rand(vecino[i])
        lista2 = [i] # Son los puntos que se están recorriendo 
        lista3 = [1] # Son las posiciones en lista2, donde la conectividad es mayor a 2.
        lista4 = [] # Es la lista de elementos dentro de un ciclo. Estos son los elementos que se eliminan en los enlaces de vecino
        while j ∉ lista2 #hace lo siguiente hasta que encuentre un ciclo. 
            push!(lista2,j)
            n = 1
            if length(vecino[j]) > 2
                push!(lista3, length(lista2))
            end
            k = vecino[j][n]
            while k == i  # tiene que elegir un camino que no sea regresar. 
                n += 1
                if n > length(vecino[j])
                    @show "error", j, n, length(vecino[j]), lista4, vecino2, true, true
                end
                k = vecino[j][n]
            end
            i = j
            j = k
        end
        test = false  #Nos dice si hay un vecino2[i] que tenga sólo un elemento.
        k = lista2[lista3[end]]
        if length(lista3) == 1 # caso 1: El ciclo se forma sin ningún punto con conectividad mayor a 2 de pormedio
            s = findall((in)(lista2[end]),vecino2[j])
            deleteat!(vecino2[j],s[1])
            s = findall((in)(lista2[2]),vecino2[j])
            deleteat!(vecino2[j],s[1])
            if length(vecino2[j]) == 1
                test = true
            end
            for i in 2:length(lista2)
                push!(lista4,lista2[i])
                vecino2[lista2[i]] = []
            end
            push!(lista4,lista2[1])
            return lista4, vecino2, test, true
        end
        if lista3[end] == length(lista2) # caso 2: El ciclo se forma con el último elemento de la lista3. 
            s = findall((in)(lista2[end]),vecino2[j])
            deleteat!(vecino2[j],s[1])
            s = findall((in)(j),vecino2[k])
            deleteat!(vecino2[k],s[1])
            if length(vecino2[j]) == 1
                test = true
            end
            x = findall((in)(j),lista2)
            for i in x[1]:length(lista2)
                push!(lista4,lista2[i])
            end
            return lista4, vecino2, test, true
        else  # caso 3: cualquier otro caso. 
            s = findall((in)(lista2[end]),vecino2[j])
            deleteat!(vecino2[j],s[1])
            s = findall((in)(lista2[lista3[end]+1]),vecino2[k])
            deleteat!(vecino2[k],s[1])
            x = findall((in)(j),lista2)
            for i in lista3[end]+1:length(lista2)
                vecino2[lista2[i]] = []
            end
            if length(vecino2[j]) == 1
                test = true
            end
            if length(vecino2[k]) == 1
                test = true
            end
            for i in x[1]:length(lista2)
                push!(lista4,lista2[i])
            end
            return lista4, vecino2, test, true
        end
    end


    function Walk(latt,neigLatt)
        #println("aca")
        #list of poins with degree 3 to start walking the graph
        list1 = listDeg3(latt,neigLatt)
        if length(list1) == 0
            #println("sali trivialmente")
            return [], neigLatt, true, false
        end
        newNeigLatt = deepcopy(neigLatt)
        #select point to start traveling graph
        pos1 = rand(list1)
        pos2 = rand(neigLatt[pos1])
        # points that are being traveled
        list2 = [pos1]
        # indexes of points in list2 where degree is higher than 2
        list3 = [1]
        # list of elements in cycle
        list4=[]
        while ~(pos2 in list2)
            #println("ciclo uno")
            push!(list2,pos2)
            n = 1
            if length(neigLatt[pos2]) > 2
                push!(list3, length(list2))
            end
            pos3 = neigLatt[pos2][n]
            #choose a different neighbor to move than where you come from
            while pos3 == pos1
                n += 1
                if n > length(neigLatt[pos2])
                    @show pos2, n, length(neigLatt[pos2]), list4, newNeigLatt, true, true
                    error("nowehere to move")
                end
                pos3 = neigLatt[pos2][n]
            end
            pos1 = pos2
            pos2 = pos3
        end
        #println("sali ciclo 1")
        # tells if in newNeigLatt there is an element with only one neighbor
        oneNeig = false 
        #take the last element with degree = 3 from list
        pos3 = list2[list3[end]]
        #Case 1 : A cycle is found with no elements with degree 3 in between:
        if length(list3)==1
            #println("caso 1")
            #erase links between first and second
            s = findall((in)([list2[end]]),newNeigLatt[pos2])
            deleteat!(newNeigLatt[pos2],s[1])
            s = findall((in)([list2[2]]),newNeigLatt[pos2])
            deleteat!(newNeigLatt[pos2],s[1])
            if length(newNeigLatt[pos2]) == 1
                oneNeig = true
            end
            #erase links of the elements with even connectivity on cycle and put them in cycle list
            for i in 2:length(list2)
                push!(list4,list2[i])
                newNeigLatt[list2[i]] = []
            end
            #add starting point
            push!(list4,list2[1])
            #println("sali")
            return list4, newNeigLatt, oneNeig, true
        end
        # Case 2: Cycle is formed with the last element of list3 
        if list3[end] == length(list2)
            #println("caso 2")
            #deleting links between last element and cycle-beggining element
            s = findall((in)([list2[end]]),newNeigLatt[pos2])
            deleteat!(newNeigLatt[pos2],s[1])
            #println("es esto")
            s = findall((in)([pos2]),newNeigLatt[pos3])
            #println("no es esto")
            deleteat!(newNeigLatt[pos3],s[1])
            if length(newNeigLatt[pos2]) == 1
                oneNeig = true
            end
            #erase links of the elements with even connectivity in cycle and put them in cycle list
            x = findall((in)([pos2]),list2)
            for i in x[1]:length(list2)
                push!(list4,list2[i])
            end
            #println("sali")
            return list4, newNeigLatt, oneNeig, true
        # Case 3: else
        else
            #println("caso 3")
            #deleting links between last element and cycle-beggining element
            #println("es esto")
            s = findall((in)([list2[end]]),newNeigLatt[pos2])
            #println("no es esto")
            deleteat!(newNeigLatt[pos2],s[1])
            s = findall((in)([list2[list3[end]+1]]),newNeigLatt[pos3])
            deleteat!(newNeigLatt[pos3],s[1])
            #println("es esto")
            x = findall((in)([pos2]),list2)
            #println("no es esto")
            #erase links of the elements with even connectivity in cycle and put them in cycle list
            x = findall((in)([pos2]),list2)
            #removing elements in cycle
            for i in list3[end]+1:length(list2)
                newNeigLatt[list2[i]] = []
            end
            if length(newNeigLatt[pos2]) == 1
                oneNeig = true
            end
            if length(newNeigLatt[pos3]) == 1
                oneNeig = true
            end
            #adding elements in cycle
            for i in x[1]:length(list2)
                push!(list4,list2[i])
            end
            #println("sali")
            return list4, newNeigLatt, oneNeig, true
        end
    end

    "sigma(vecino) genera un arreglo σ a partir del arreglo de conectividad vecino."
    function sigma(vecino)
        σ = Int[]
        n = length(vecino)
        for i in 1:n
            x = 0
            if length(vecino[i])>1
                x = 1
            end
            push!(σ,x)
        end
        return σ
    end

    function makeLatticeFromNeigbors(neigLatt)
        latt = zeros(Int8,size(neigLatt))
        for pos in CartesianIndices(size(neigLatt))
            if length(neigLatt[pos]) > 1
                latt[pos] = 1
            end
        end
        return latt
    end

    "sigma2(lista, σ) agrega a σ los elementos del arreglo lista"
    function sigma2(lista, σ)
        for i in lista
            σ[i] = 1
        end
        return σ
    end  
    
    function addFromList!(latt,posList)
        for pos in posList
            latt[pos] = 1
        end
    end

    "ciclos(σ,vecino,n,m) regresa 3 arreglos, el que corresponde a los ciclos, el arreglo original y la diferencia entre los ciclos y el original"    
    function ciclos(σ,vecino,n1,m)
        σo = copy(σ)
        σ, σn = quita_pelos!(σ,vecino,n1,m)
        vecino2 = vecinos_s(σ,vecino)
        lista, vecino2, test, test2 = recorrer(vecino2)
        σ2 = zeros(Int,n1*m)   
        cont = 0
        while test2   
            σ = sigma(vecino2)
            σ2 = sigma2(lista, σ2)
            if test
                cont += 1
                σcop = copy(σ)
                σ, σn1 = quita_pelos!(σ,vecino, vecino2,n1,m)
                σn += σn1
                vecino2 = vecinos_s(σ,vecino)
            end
            
            lista, vecino2, test, test2 = recorrer(vecino2)
        end
        σ, σn1 = quita_pelos!(σ,vecino,n1,m)
        σn += σn1
        σc = σ+σ2
        for i in 1:length(σc)
            if σc[i]>1
                σc[i] =1
            end
        end
        return σo, σc, σo-σc
    end

    function cycles(latt,neigLatt)
        original = deepcopy(latt)
        sigma, tails = RemoveTails(latt,neigLatt)
        cycNeigLatt = CycleNeigs(sigma,neigLatt)
        list, cycNeigLatt, bool1, bool2 = Walk(sigma, cycNeigLatt)
        res = zeros(Int8,size(latt))
        cont=0
        while bool2
            sigma = makeLatticeFromNeigbors(cycNeigLatt)
            addFromList!(res,list)
            if bool1
                cont += 1
                cop = deepcopy(sigma)
                sigma , sigman1 = RemoveTails(sigma,neigLatt)
                tails += sigman1
                cycNeigLatt = CycleNeigs(sigma,neigLatt)
            end
            list, cycNeigLatt, bool1, bool2 = Walk(sigma,cycNeigLatt)
        end
        sigma, sigman1 = RemoveTails(sigma, neigLatt)
        tails += sigman1
        cyc = sigma + res
        for pos in CartesianIndices(size(cyc))
            if cyc[pos]>1
                cyc[pos] =1
            end
        end
        return original, cyc, original - cyc
    end
    function ciclos2(latt)
        n=size(latt)[1]
        ##println("hay ciclos")
        a=vectorize(latt)
        b=vecinos(n,n)
        vec=ciclos(a,b,n,n)[2]
        return matrize(vec,(n,n))
    end
    
    #=
    include("geometry.jl")
    latt , neigLatt =  Geometry.BuildLattices( [4,2,"square"], "cycle")
    
    cycles(latt,neigLatt)
    ciclos2(latt)

    latt , neigLatt =  Geometry.BuildLattices( [100,2,"square"], "cycle")
    
    @time cycles(latt,neigLatt)
    @time ciclos2(latt)
    =#
end

