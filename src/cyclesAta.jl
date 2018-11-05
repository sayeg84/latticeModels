module CyclesAta

    "vectorizar(mat) convierte un arreglo n-dimensional mat en un vector de n*m"
    function vectorizar(mat)
        res=[]
        for pos in CartesianRange(size(mat))
            push!(res,mat[pos])
        end
        return res
    end
    "matrizar convierte un vector de longitud n*m en una matriz de n*m"
    function matrizar(vec,n,m)
        if n*m != length(vec)
            error("las dimensiones deben ser iguales")
        end
        res=zeros(typeof(vec[1]),n,m)
        for i in 1:n
            for j in 1:m
                index= (i-1)*n +j
                res[i,j]=vec[index]
            end
        end
        return res
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

    "quita_pelos!(σ,vecino,n,m) quita del arreglo σ, los elementos cuya conectividad sea menor que 2. Esto lo hace de forma recurrente, hasta que ningún el emento de σ tenga una conectividad menor que 2. n y m representan los tamaños de la caja donde se encuentra el arreglo σ. vecino es un arreglo que contiene la conectividad de la cuadrada. Regresa los valores de σ y un arreglo con los elementos que se quitaron (los pelos). También hay la versión con quita_pelos!(σ,vecino, vecino2,n,m), en cuyo caso tomará en cuenta los enlaces."
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
                    @show "error", j, n, length(vecino[j])
                    lista4, vecino2, true, true
                end
                k = vecino[j][n]
            end
            i = j
            j = k
        end
        test = false  #Nos dice si hay un vecino2[i] que tenga sólo un elemento.
        k = lista2[lista3[end]]
        if length(lista3) == 1 # caso 1: El ciclo se forma sin ningún punto con conectividad mayor a 2 de pormedio
            s = findin(vecino2[j],lista2[end])
            deleteat!(vecino2[j],s[1])
            s = findin(vecino2[j],lista2[2])
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
            s = findin(vecino2[j],lista2[end])
            deleteat!(vecino2[j],s[1])
            s = findin(vecino2[k],j)
            deleteat!(vecino2[k],s[1])
            if length(vecino2[j]) == 1
                test = true
            end
            x = findin(lista2,j)
            for i in x[1]:length(lista2)
                push!(lista4,lista2[i])
            end
            return lista4, vecino2, test, true
        else  # caso 3: cualquier otro caso. 
            s = findin(vecino2[j],lista2[end])
            deleteat!(vecino2[j],s[1])
            s = findin(vecino2[k],lista2[lista3[end]+1])
            deleteat!(vecino2[k],s[1])
            x = findin(lista2,j)
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

    "sigma2(lista, σ) agrega a σ los elementos del arreglo lista"
    function sigma2(lista, σ)
        n = length(σ)
        for i in lista
            σ[i] = 1
        end
        return σ
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
    function ciclosprueba()
        println("si se puede")
    end
    function ciclos2(latt)
        n=size(latt)[1]
        #println("hay ciclos")

        a=vectorizar(latt)
        b=vecinos(n,n)
        vec=ciclos(a,b,n,n)[2]
        return matrizar(vec,n,n)
    end
    #=
    #gr()
    n= 10
    vecino = vecinos2(n,n)
    σ = rand([0,1],n^2)
    #quita_pelos!(σ,vecino, n,n)
    @time σo,σc, σn = ciclos(σ,vecino,n,n)
    @time σo,σc, σn = ciclos(σ,vecino,n,n)
    @time σo,σc, σn = ciclos(σ,vecino,n,n)
    =#
end
