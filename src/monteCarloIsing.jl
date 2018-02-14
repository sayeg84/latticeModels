
using Plots, Elliptic, Rsvg
plotlyjs()
init=Dates.time()
#funcion modulo auxiliar. Igual que el módulo normal, pero tal que mod(n,n)=n
println("Iniciando")

function modl(a::Int64,b::Int64)
    x=mod(a,b)
    if x==0
        return b
    else
        return x
    end
end
#función de vecinos inmediatos utilizada para la función de peso del algoritmo metrópolis
function vecinos(M,P;test=false)
    s=size(M)
    v1 = M[modl(P[1]-1,s[1]), P[2]]
    v2 = M[modl(P[1]+1,s[1]), P[2]]
    v3 = M[P[1],modl(P[2]-1,s[2])]
    v4 = M[P[1],modl(P[2]+1,s[2])]
    if test
        println([0 v1 0; v3 0 v4 ; 0 v2 0])
    end
    return v1+v2+v3+v4
end
#función auxiliar para calcular la energía
function energia(M)
    e=0
    for i in 1:size(M)[1]
        for j in 1:size(M)[2]
            #=
            al contabilizar todos los vecinos para cada punto, debemos de dividir entre dos pues sumamos cada
            par de vecinos inmediatos dos veces
            =#
            e=e-B*M[i,j]-J*vecinos(M,[i,j])*M[i,j]*(1/2)
        end
    end
    return e
end



#funciones auxiliares de la solución teórica de Onsager para los parámetros del modelo
function magOnsager(a)
    if a<2.268
        return (1-(sinh(2*J/a)^(-4)))^(1/8)
    else
        return 0
    end
end
function energiaOnsager(a)
    x=2J/a
    k1=2/(cosh(x)*coth(x))
    k2=2(tanh(x)^2)-1
    k3=Elliptic.K(k1)
    return -J*coth(x)*(1+2/pi*k2*k3)
end
function cvOnsager(a)
    x=2J/a
    k1=2/(cosh(x)*coth(x))
    k2=2(tanh(x)^2)-1
    k3=Elliptic.K(k1)
    k4=Elliptic.E(k1)
    k5=(x*coth(x))^2/(2*pi)
    return k5*(2*k3-2*k4-(1-k2)*(pi/2 + k2*k3)) 
end

#=
Probabilidad de cambio para el algoritmo Metrópolis.

Definida como el cociente de la densidad de probabilidad de estados
del espacio muestra entre dos estados deseados. Para mayor información,
consultar la referencia bibliográfica
=#
function r(M,P)
    aux=J*vecinos(M,P)+B
    return exp(-2*M[P[1],P[2]]*aux*β)
end
#=
Implementación del algoritmo metrópolis para crear N matrices de Nx*Ny 
usando el algoritmo Metrópolis. Suponemos los valores de spin son 1 y -1.
=#
function metro(N=100, Nx=10, Ny=10,f=10;test=false)
    X=[]
    #calculamos la energía como parte de la simulación
    Y=rand([-1,1],Nx,Ny)
    push!(X,Y)
    for i in 2:N
        res=0
        P=[rand(1:Nx),rand(1:Ny)]
        res=r(Y,P)
        t=rand()
        if test
            print("iter = ")
            println(i-1)
            println("lattice")
            println(Y)
            print("P = ")
            println(P)
            print("r = ")
            println(res)
            print("tes = ")
            println(t)
            println()
            println()
        end
        if res>t
  
            Y[P[1],P[2]]=-1*Y[P[1],P[2]]
        end
        if mod(i,f)==1
            push!(X,copy(Y))
        end
    end
    return X
end

#=
Elección de parámetros importantes.

Por simplicicidad, tomamos que la Constante de Boltzmann tiene valor de 1
=#
#=
T=0.1
β=1/T
J=1
B=0
X=metro(22,5,5,test=true)
=#

n = parse(Int64,ARGS[1])
J = 1
B = 0
β = 0
aux = parse(Int64,ARGS[2])
pasos = 10^aux + 1
aux = parse(Int64,ARGS[3])
freq = 10^aux
T = linspace(0.1,5,50)
M = []
E = []
EE = []
println()
println("Simulando")

for t in T
    β=1/t
    temp1=[]
    temp2=[]
    temp3=[]
    for i in 1:5
        X=metro(pasos,n,n,freq)
        inicio=convert(Int64,floor(length(X)*3/4))
        mag=[abs(sum(X[i])/n^2) for i in inicio:length(X)]
        push!(temp1,mean(mag))
        en=[energia(x) for x in X]
        temp4=[en[i] for i in inicio:length(X)]
        push!(temp2,mean(temp4))
        push!(temp3,mean([x^2 for x in temp4]))
        print(t)
        print("-")
        println(i)
    end
    push!(M,mean(temp1))
    push!(E,mean(temp2))
    push!(EE,mean(temp3))
end
CV=[(EE[i]-E[i]^2)/(T[i]^2) for i in 1:length(T)]


#escribimos los datos

println("Escribiendo resultados")





if ~(isdir("../outputs"))
    mkdir("../outputs")
end
dest=Dates.format(Dates.now(),"HH:MM:SS dd-mm-Y")
mkdir("../outputs/$(dest)")


T=linspace(0.1,5,50)
p1=scatter(T,M,
    title="Magnetizacion",
    xlabel="Temperatura (KB=1)",
    ylabel="Magnetizacion absoluta promedio por spin",
    label="simulación con N=10000"
)
#plot!(T,[magOnsager(t) for t in T],label="Solución teórica")
#plot!([2.3,2.3],[0,1],linestyle=:dash,label="Temperatura crítica")
savefig("../outputs/$(dest)/mag.png")
#display(p)


p=scatter(T,E/n^2,
    title="Energia",
    xlabel="Temperatura (KB=1)",
    ylabel="Energía promedio por spin",
    label="simulación con N=10000"
)

#plot!(T,[energiaOnsager(t) for t in T],label="Solución teórica")
#plot!([2.3,2.3],[-2,0],linestyle=:dash,label="Temperatura crítica")
savefig("../outputs/$(dest)/ener.png")
#display(p)



p=scatter([T[i] for i in 10:length(T)],[CV[i]/n^2 for i in 10:length(CV)],
    title="Capacidad calorífica",
    xlabel="Temperatura (KB=1)",
    ylabel="Capacidad calorífica por spin",
    label="simulación con N=10000"
)
#plot!(T,[cvOnsager(t) for t in T],label="Solución teórica")
#plot!([2.3,2.3],[0,2],linestyle=:dash,label="Temperatura crítica")
savefig("../outputs/$(dest)/cv.png")


write("../outputs/$(dest)/datos.csv")
open("../outputs/$(dest)/datos.csv","w") do f 
        write(f,"Modelo de ising,Fecha, $(dest) \n")
        write(f," ,Tiempo de ejecución, $(round(Dates.time()-init,2)) \n")
        write(f," Parámetros \n")
        write(f," n, J, B, pasos, frecuencia \n ")
        write(f," $n, $J, $B, $(pasos), $(freq) \n ")
        write(f,"\n")
        write(f,"T,M,E,CV \n")
        for i in 1:length(T)
            write(f,"$(T[i]), $(M[i]),$(E[i]), $(CV[i]) \n")
        end
    end