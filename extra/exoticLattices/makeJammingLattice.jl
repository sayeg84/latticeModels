using ArgParse
using CSV
using DataFrames
using Dates
using DelimitedFiles
using LinearAlgebra
using LsqFit
using Plots

function ParseCommandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "-N","--Ndiscs"
            help = "Number of discs in the simulation"
            arg_type = Int64
            default = 100
        "-p","--phi"
            help = "Packaging Fraction"
            arg_type = Float64
            default = 0.82
        "-S","--save"
            help = "Save prefix for filenames"
            arg_type = String
            default = "jam"
    end
    return parse_args(s)
end

const parsed_args = ParseCommandline()

gr()

#maximo tiempo simulacion
global tinf = 1e6;

function modulo(Numero::Int64, Mod::Int64)
    if Numero%Mod == 0
        Mod
    else
        Numero%Mod
    end
end

function contar_valor_en_arreglo(Arreglo::Array{Float64, 1}, valor::Real)
    N = length(Arreglo)
    contador = 0
    for i in 1:N
        if Arreglo[i] == valor
            contador += 1
        end
    end
   contador
end

mutable struct esferas
    N::Int64 #numero de discos
    Ntipos::Int64
    posiciones::Array{Float64, 2}
    velocidades::Array{Float64, 2}
    radios::Array{Float64,1}
    CR::Array{Float64,1} #crecimiento de radio
    tiempos::Array{Float64,1}
    masas::Array{Float64,1}
    vecinos::Array{Array{Int64,1},1}
    numeroVecinos::Array{Int64,1}
    numeroColisiones::Array{Int64,1}
    transferenciasPeriodicas::Array{Int64, 2}
    caja::Array{Float64,1} #caja cuadrada del espacio de contencion
    function esferas(N::Int64, Caja::Array{Float64,1}, cr::Float64, ntipos::Int64, vmax::Float64 = 3.0)
        posiciones = zeros(N, 2)
        radios = zeros(N)
        crecimientos = zeros(N)
        velocidades = zeros(N, 2)
        masas = ones(N)
        tiempos = zeros(N)
        numeroVecinos = zeros(Int64, N)
        numeroColisiones = zeros(Int64, N)
        transferenciasPeriodicas = zeros(Int64, N, 2)
        # tamaño del dominio
        L = abs(Caja[2] - Caja[1])
        for i in 1:N
            posiciones[i,:] = [rand()*L, rand()*L] 
            if ntipos == 1
                crecimientos[i] = cr
            elseif ntipos == 2 
                if rand() < 2/3
                    crecimientos[i] = cr
                else
                    crecimientos[i] = cr*1.4
                    masas[i] = 1.4^2
                end
            else
                error("Solo se aceptan sistemas con 1 o 2 crecimientos")
            end
            
            #Velocidades
            r = rand()*vmax
            ϕ = rand()*2*pi
            # generador de velocidad con distribución circular
            velocidades[i,:] = [r*cos(ϕ), r*sin(ϕ)]
        end
        
        vecinos = Array{Array{Int64,1},1}(undef, N)
        for i in 1:N
            vecinos[i] = []
        end
        new(N, ntipos, posiciones, velocidades, radios, crecimientos, tiempos, masas, vecinos, numeroVecinos, numeroColisiones, transferenciasPeriodicas, Caja)
    end
end

function circulos(Esferas::esferas)
    #rdp = radio particulas
    θ = 0:0.05*pi:2pi
    nθ = length(θ)
    circulo = zeros(nθ, 2)
    
    # generador de círculos alrededor de la partícula
    for j in 1:Esferas.N
        circulo[1:nθ, 1] = collect(Esferas.posiciones[j,1] .+ Esferas.radios[j]*cos.(θ'))'
        circulo[1:nθ, 2] = collect(Esferas.posiciones[j,2] .+ Esferas.radios[j]*sin.(θ'))'
        plot!(circulo[:,1],circulo[:,2], c=:black, fillrange=0, fillcolor=RGB(0.4 + Esferas.radios[j], 0.4 + Esferas.radios[j]*3, 0.3))
    end
    plot!()
end

function arraydearrayint(m::Int64)
    str = "["
    M = (m+2)^2 - 1
    for i in 1:M
        str *= "Array{Int64}(undef, 0), "
    end
    str *= "Array{Int64}(undef, 0)]"
    ex = Meta.parse(str)
    prov = Meta.eval(ex)
    reshape(prov, m+2, m+2)
end

function arraydearraytup(m::Int64)
    M = (m+2)^2 - 1
    str = "["
    for i in 1:M
        str *= "Array{Array{Int64,1},1}(), "
    end
    str *= "Array{Array{Int64,1},1}()]"
    ex = Meta.parse(str)
    prov = Meta.eval(ex)
    reshape(prov, m+2, m+2)
end

function celdas_vecinas_iniciales(m::Int64)
    cvi = arraydearraytup(m)
    m1 = m + 1

    for i in 2:m1
        for j in 2:m1
            for h in -1:0, l in -1:1
                if h == 0 && l == -1  else push!(cvi[i,j], [i+h,j+l]) end
            end
        end
    end
    cvi
end

function celdas_vecinas(m::Int64)
    cv = arraydearraytup(m)
    m1 = m + 1

    for i in 2:m1
        for j in 2:m1
            for h in -1:1, l in -1:1
                push!(cv[i,j], [i+h,j+l])
            end
        end
    end
    cv
end


function celdas_imagen(m::Int64)
    #ci celdas para esferas imagen
    ci = Array{Array{Int64,1},1}(undef, 0)
    
    for i in 1:(m+2)
        for j in 1:(m+2)
            if i == 1 
                push!(ci, [i,j])
            elseif j == 1 && i != 1 
                push!(ci, [i,j])
            elseif i == m + 2 
                push!(ci, [i,j])
            elseif j == m + 2 && i != m + 2
                push!(ci, [i,j])
            end
        end
    end
    ci
end

function celdas_frontera(m::Int64)
#cf celdas frontera
    cf = Array{Array{Int64,1},1}(undef, 0)
    
    for i in 2:(m+1)
        for j in 2:(m+1)
            if i == 2
                push!(cf, [i,j])
            elseif j == 2 && i != 2
                push!(cf, [i,j])
            elseif i == m + 1
                push!(cf, [i,j])
            elseif j == m + 1 && i != m + 1
                push!(cf, [i,j])
            end
        end
    end
    cf
end

function celdas_pared(m::Int64)
#cf celdas frontera
    cp = Array{Array{Int64,1},1}(undef, 0)
    
    for i in 2:(m+1)
        for j in 2:(m+1)
            if i == 2
                push!(cp, [i,j])
            elseif j == 2 && i != 2
                push!(cp, [i,j])
            elseif i == m + 1
                push!(cp, [i,j])
            elseif j == m + 1 && i != m + 1
                push!(cp, [i,j])
            end
        end
    end
    cp
end

mutable struct celdas
    m::Int64 #numero de celdas por lado sin contar celdas de esferas imagen
    L::Float64
    residentes::Array{Array{Int64,1},2}
    presenteEnCelda::Array{Array{Int64,1},1}
    vecinas::Array{Array{Array{Int64,1},1},2}
    vecinasIniciales::Array{Array{Array{Int64,1},1},2}
    fronteraPeriod::Array{Array{Int64,1},1}
    imagen::Array{Array{Int64,1},1}
    pared::Array{Array{Int64,1},1}
    LC::Array{Float64,1} #limites celdas x y
    function celdas(Esferas::esferas)
        unitParticle = 0.0
        if Esferas.Ntipos == 1
            unitParticle = Esferas.N
            L = 1 / floor(sqrt(unitParticle - unitParticle/10) )
        else
            for n1 in 1:Esferas.N
                unitParticle += Esferas.CR[n1]/maximum(Esferas.CR)
            end
            L = 1 / floor(sqrt(unitParticle - unitParticle/4) )
        end
        M = convert(Int64, ceil(abs(Esferas.caja[2] - Esferas.caja[1]) / L)) # numero de celdas por dimension (sin contar fantasma)
    
        LC = collect(Esferas.caja[1]:L:Esferas.caja[2]) #inicio:L:final
        
        preResidentes = arraydearrayint(M-2)
        Residentes = arraydearrayint(M)
        PresenteEnCelda = Array{Array{Int64,1},1}(undef, Esferas.N)
        
        for n1 in 1:Esferas.N
            i::Int = cld(Esferas.posiciones[n1,1], L)
            j::Int = cld(Esferas.posiciones[n1,2], L)
            push!(preResidentes[i,j], n1)
            PresenteEnCelda[n1] = [i+1, j+1]
        end
        
        Residentes[2:(M+1),2:(M+1)] = preResidentes
        
        cf = celdas_frontera(M)
        ci = celdas_imagen(M)
        cp = celdas_pared(M)
        
        new(M, L, Residentes, PresenteEnCelda, celdas_vecinas(M), celdas_vecinas_iniciales(M), cf, ci, cp, LC)
    end
end

function esferasCeldasImagen(Esferas::esferas, Celdas::celdas) 
    N = Esferas.N
    m1 = Celdas.m + 1 #celdas frontera
    m2 = Celdas.m + 2 #celdas fantasma

    for C in Celdas.fronteraPeriod
        for i in 1:2
            cf = copy(C)
            if C[i] == 2
                cf[i] = m2
                Celdas.residentes[cf[1],cf[2]] = Celdas.residentes[C[1],C[2]]
            elseif C[i] == m1
                cf[i] = 1
                Celdas.residentes[cf[1],cf[2]] = Celdas.residentes[C[1],C[2]]
            end
        end
    end
    
    #vertices
    for i in [2,m1], j in [2,m1]
        if i == 2
            nuevai = m2
        elseif i == m1
            nuevai = 1
        end
        if j == 2
            nuevaj = m2
        elseif j == m1
            nuevaj = 1
        end
        Celdas.residentes[nuevai,nuevaj] = Celdas.residentes[i,j]
    end
    
    Celdas
end

function posicionEsferasImagen(Esferas::esferas, Celdas::celdas)
    NI = 0
    for C in Celdas.imagen
        NI += length(Celdas.residentes[C[1], C[2]])
    end
    
    posicionesImagen = zeros(NI, 2)
    residentesImagen = zeros(Int64, NI)
    PresenteEnCelda = Array{Array{Int64,1},1}(undef, NI)
    NIcont = 0
    
    L = abs(Esferas.caja[2] - Esferas.caja[1])
            
    for C in Celdas.imagen
        for ni in Celdas.residentes[C[1], C[2]]
            NIcont += 1
            residentesImagen[NIcont] = ni
            PresenteEnCelda[NIcont] = C
            posicionesImagen[NIcont, :] = Esferas.posiciones[ni,:]
            for i in 1:2
                if C[i] == 1 
                    posicionesImagen[NIcont, i] -= L
                elseif C[i] == Celdas.m + 2 
                    posicionesImagen[NIcont, i] += L
                end
            end
        end
    end

    NI, posicionesImagen, residentesImagen, PresenteEnCelda
end

mutable struct esferasImagen
    NI::Int64
    posiciones::Array{Float64, 2}
    residentes::Array{Int64, 1}
    presenteEnCelda::Array{Array{Int64,1},1}
    function esferasImagen(Esferas::esferas, Celdas::celdas)
        Ni, posiciones, residentes, presenteEnCelda = posicionEsferasImagen(Esferas, Celdas)
        new(Ni, posiciones, residentes, presenteEnCelda)
    end
end

function circulosImagen(Esferas::esferas, EsferasImagen::esferasImagen)
    #rdp = radio particulas
    θ = 0:0.05*pi:2pi
    nθ = length(θ)
    circulo = zeros(nθ, 2)
    
    # generador de círculo alrededor de la partícula
    for j in 1:EsferasImagen.NI
        circulo[1:nθ, 1] = collect(EsferasImagen.posiciones[j,1] .+ Esferas.radios[EsferasImagen.residentes[j]]*cos.(θ'))'
        circulo[1:nθ, 2] = collect(EsferasImagen.posiciones[j,2] .+ Esferas.radios[EsferasImagen.residentes[j]]*sin.(θ'))'
        plot!(circulo[:,1],circulo[:,2], c=:black, fillrange=0, fillcolor=RGB(0.4 + Esferas.radios[EsferasImagen.residentes[j]], 0.4 + Esferas.radios[EsferasImagen.residentes[j]]*3, 0.7))
    end
    plot!()
end

function findindex(esfera2I::Int64, Celda::Array{Int64, 1}, Celdas::celdas, EsferasImagen::esferasImagen)
    index = 0
    for i in 1:EsferasImagen.NI
        if Celda == EsferasImagen.presenteEnCelda[i] && EsferasImagen.residentes[i] == esfera2I 
            index = i
        end
    end
    index
end

function Base.findfirst(X::Array{Float64,N} where N, a::Real)
    for i in 1:length(X)
        if X[i] == a
            return i
        end
    end
    return 0
end

function tiempo_para_colision(n1::Int64, n2::Int64, Esferas::esferas, Celdas::celdas, gtime = 0.0)
    p1 = Esferas.posiciones[n1,:] + Esferas.velocidades[n1,:].*(gtime - Esferas.tiempos[n1])
    p2 = Esferas.posiciones[n2,:] + Esferas.velocidades[n2,:].*(gtime - Esferas.tiempos[n2])
    v1 = Esferas.velocidades[n1,:]
    v2 = Esferas.velocidades[n2,:]
 
    Δv = v1 - v2
    Δp = p1 - p2
    
    σ = Esferas.CR[n1] + Esferas.CR[n2]
    if σ == 0.0
        ρ = Esferas.radios[n1] + Esferas.radios[n2]
    else
        ρ = gtime*σ 
    end
    a, b, c = (norm(Δv)^2 - σ^2, Δv⋅Δp - σ*ρ, norm(Δp)^2 - ρ^2)
    if c < -1e-4*ρ
#         @show "n", n1, n2, a, b, c
        if 2*Esferas.radios[n1] > Celdas.L || 2*Esferas.radios[n2] > Celdas.L
            error("Alguna esfera ha crecido mas que la celda")
        end
#         error("Traslape encontrado, $c")
    end
    
    t = quadratic(a, b, c)
end

function tiempo_para_colision_imagen(n1::Int64, n2::Int64, nEI::Int64, Esferas::esferas, EsferasImagen::esferasImagen, Celdas::celdas, gtime = 0.0)
    p1 = Esferas.posiciones[n1,:] + Esferas.velocidades[n1,:].*(gtime - Esferas.tiempos[n1])
    p2 = EsferasImagen.posiciones[nEI,:] + Esferas.velocidades[n2,:].*(gtime - Esferas.tiempos[n2])
    v1 = Esferas.velocidades[n1,:]
    v2 = Esferas.velocidades[n2,:]

    Δv = v1 - v2
    Δp = p1 - p2
    
    σ = Esferas.CR[n1] + Esferas.CR[n2]
    if σ == 0.0
        ρ = Esferas.radios[n1] + Esferas.radios[n2]
    else
        ρ = gtime*σ 
    end
        
    a, b, c = (norm(Δv)^2 - σ^2, Δv⋅Δp - σ*ρ, norm(Δp)^2 - ρ^2)
    if c < -1e-4*ρ
#         @show "f", n1, n2, a, b, c
        if 2*Esferas.radios[n1] > Celdas.L || 2*Esferas.radios[n2] > Celdas.L
            error("Alguna esfera ha crecido mas que la celda")
        end
#         error("Traslape encontrado, $c")
    end
    t = quadratic(a, b, c)
end

function quadratic(a::Float64, b::Float64, c::Float64)
    d2 = b^2 - a*c
    t = tinf
    if c <= 0.0 && b < 0.0
        t = 0.0
    elseif d2 >= -10*eps()
        if (d2 < 0.0); d2 = 0.0; end
        if b < 0.0
            t = c / (-b + sqrt(d2))
        elseif a < 0.0 && b > 0.0
            t = -(b + sqrt(d2))/a
        end
    else 
        t = tinf
    end
    t
end

function tiempo_colision_con_pared(n1::Int64, Esferas::esferas, gtime = 0.0)
    tsP = zeros(2)
    cps = zeros(Int64, 2)
    pn1 = Esferas.posiciones[n1,:] + Esferas.velocidades[n1,:].*(gtime - Esferas.tiempos[n1])
    rnow = Esferas.radios[n1] + (gtime - Esferas.tiempos[n1])*Esferas.CR[n1] 
    for i in 1:2
        if Esferas.velocidades[n1,i] + Esferas.CR[n1] > 0.0
            cps[i] = -1 - 2*(i-1)
            tsP[i] = (Esferas.caja[2] - rnow - pn1[i]  ) / (Esferas.velocidades[n1,i] + Esferas.CR[n1])
        elseif Esferas.velocidades[n1,i] - Esferas.CR[n1] < 0.0  
            cps[i] = -2 - 2*(i-1)
            tsP[i] = (Esferas.caja[1] + rnow - pn1[i] ) / (Esferas.velocidades[n1,i] - Esferas.CR[n1])
        else
            cps[i] = 0
            tsP[i] = tinf
        end
    end
    
    t = minimum(tsP)
    cp = cps[findfirst(tsP, t)]
    
    t, cp
end

function tiempos_iniciales_para_colision(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen)
    ts = zeros(Esferas.N)
    vs = zeros(Int64, Esferas.N)
    
    for n1 in 1:Esferas.N
        tsn1 = Array{Float64, 1}(undef, 0)
        vsn1 = Array{Int64, 1}(undef, 0)
        i, j = Celdas.presenteEnCelda[n1]
        
        #checa esferas en las celdas vecinas iniciales
        for C in Celdas.vecinasIniciales[i,j]
            for n2 in Celdas.residentes[C[1],C[2]]
                push!(vsn1, n2)
                if C ∈ Celdas.imagen
                    nEI = findindex(n2, C, Celdas, EsferasImagen)
                    t = tiempo_para_colision_imagen(n1, n2, nEI, Esferas, EsferasImagen, Celdas)
                elseif n1 != n2
                    t = tiempo_para_colision(n1, n2, Esferas, Celdas)
                else
                    t = tinf
                end
                push!(tsn1, t)
            end
        end
        
        #checa colision con paredes si esta presente en celda de frontera
#         if [i,j] ∈ Celdas.pared
#             t, cp = tiempo_colision_con_pared(n1, Esferas)
#             push!(tsn1, t)
#             push!(vsn1, cp)
#         end
        
        ts[n1] = minimum(tsn1)
        vs[n1] = vsn1[findfirst(tsn1, ts[n1])]
        
    end
    ts, vs
end

function tiempo_cambio_de_celda(n1::Int64, Esferas::esferas, Celdas::celdas, gtime = 0.0)
    indices = Celdas.presenteEnCelda[n1]
    tsc = zeros(2)
    transc = zeros(Int64, 2)
    
    pn1 = Esferas.posiciones[n1,:] + Esferas.velocidades[n1,:].*(gtime - Esferas.tiempos[n1])
    
    for i in 1:2
        if Esferas.velocidades[n1,i] > 0.0
            transc[i] = 1 + 2*(i-1)
            tsc[i] = (Celdas.LC[indices[i]] - pn1[i]) / Esferas.velocidades[n1,i]
        elseif Esferas.velocidades[n1,i] < 0.0 
            transc[i] = 2 + 2*(i-1)
            tsc[i] = (Celdas.LC[indices[i]-1] - pn1[i]) / Esferas.velocidades[n1,i]
        else
            transc[i] = 0
            tsc[i] = tinf
        end
    end
    t = minimum(tsc)
    t, transc[findfirst(tsc, t)]
end

function tiempos_iniciales_cambio_de_celda(Esferas::esferas, Celdas::celdas)
    ts = zeros(Esferas.N)
    transferencia = zeros(Int64, Esferas.N)
    
    for n1 in 1:Esferas.N
        ts[n1], transferencia[n1] = tiempo_cambio_de_celda(n1, Esferas, Celdas)
    end
    ts, transferencia
end

mutable struct eventos
    tiempo::Float64
    tiemposCyT::Array{Float64, 2}
    eventoCyT::Array{Int64, 2}
    function eventos(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen)
        ts = zeros(Esferas.N, 2)
        ev = zeros(Int64, Esferas.N, 2)
        ts[:,1], ev[:,1] = tiempos_iniciales_para_colision(Esferas, Celdas, EsferasImagen)
        ts[:,2], ev[:,2] = tiempos_iniciales_cambio_de_celda(Esferas, Celdas)
        new(0.0, ts, ev)
    end
end

function mover_esfera(n1::Int64, Esferas::esferas, Eventos::eventos)
    #movimiento y crecimiento
    Δt = (Eventos.tiempo - Esferas.tiempos[n1])
    Esferas.posiciones[n1,:] = Esferas.posiciones[n1,:] + Δt.*Esferas.velocidades[n1,:]
    if Esferas.CR[n1] != 0.0
        Esferas.radios[n1] = Eventos.tiempo*Esferas.CR[n1]
    end 
    Esferas.tiempos[n1] = Eventos.tiempo
    Esferas
end

function mover_esferas(Esferas::esferas, Eventos::eventos)
    for n1 in 1:Esferas.N
        Esferas = mover_esfera(n1, Esferas, Eventos)
    end
    Esferas
end

function tiempo_nuevo_para_colision(n1::Int64, Celdaold::Array{Int64, 1}, Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos)
    tsn1 = Array{Float64, 1}(undef, 0)
    vsn1 = Array{Int64, 1}(undef, 0)
    
    push!(vsn1, 0)
    push!(tsn1, tinf)
    
    i, j = Celdas.presenteEnCelda[n1]
    
    seguir = false
    if [i,j] ∈ Celdas.fronteraPeriod
        seguir = true
    end
    
    nuevasNoVecinas = Celdas.vecinas[i,j] ∩ Celdas.vecinas[Celdaold[1],Celdaold[2]]
    
    for C in Celdas.vecinas[i,j]
        if seguir || C ∉ nuevasNoVecinas
            for n2 in Celdas.residentes[C[1],C[2]]
                push!(vsn1, n2)
                if C ∈ Celdas.imagen
                    nEI = findindex(n2, C, Celdas, EsferasImagen)
                    t = tiempo_para_colision_imagen(n1, n2, nEI, Esferas, EsferasImagen, Celdas, Eventos.tiempo)
                elseif n1 != n2
                    t = tiempo_para_colision(n1, n2, Esferas, Celdas, Eventos.tiempo)
                else 
                    t = tinf
                end
                push!(tsn1, t)
            end
        end
    end
    
#     if [i,j] ∈ Celdas.pared
#         t, cp = tiempo_colision_con_pared(n1, Esferas)
#         push!(tsn1, t)
#         push!(vsn1, cp)
#     end  
    
    tmin = minimum(tsn1)
    index = findfirst(tsn1, tmin)
    
    tc = tmin + Eventos.tiempo
    
    if vsn1[index] != 0 && Eventos.tiemposCyT[n1, 1] > tc 
        Eventos.tiemposCyT[n1, 1] = tc
        Eventos.eventoCyT[n1, 1] = vsn1[index]
    end

    Eventos
end

function tiempos_nuevos_para_colision(n1s::Array{Int64,1}, Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos)
    for n1 in n1s
        tsn1 = Array{Float64, 1}(undef, 0)
        vsn1 = Array{Int64, 1}(undef, 0)

        i, j = Celdas.presenteEnCelda[n1]
        
        for C in Celdas.vecinas[i,j]
            for n2 in Celdas.residentes[C[1],C[2]]
                push!(vsn1, n2)
                if C ∈ Celdas.imagen && n2 ∉ n1s
                    nEI = findindex(n2, C, Celdas, EsferasImagen)
                    t = tiempo_para_colision_imagen(n1, n2, nEI, Esferas, EsferasImagen, Celdas, Eventos.tiempo)
                elseif n2 ∉ n1s
                    t = tiempo_para_colision(n1, n2, Esferas, Celdas, Eventos.tiempo)
                else 
                    t = tinf
                end
                push!(tsn1, t)
            end
        end

#         if [i,j] ∈ Celdas.pared
#             t, cp = tiempo_colision_con_pared(n1, Esferas)
#             push!(tsn1, t)
#             push!(vsn1, cp)
#         end
        
    tmin = minimum(tsn1)
    index = findfirst(tsn1, tmin)
        
    Eventos.tiemposCyT[n1, 1] = tmin + Eventos.tiempo
    Eventos.eventoCyT[n1, 1] = vsn1[index]

    end
    
    Eventos
end

function tiempo_nuevo_cambio_de_celda(n1::Int64, Esferas::esferas, Celdas::celdas, Eventos::eventos)
    t, trans = tiempo_cambio_de_celda(n1, Esferas, Celdas, Eventos.tiempo)
    
    Eventos.tiemposCyT[n1, 2] = t + Eventos.tiempo
    Eventos.eventoCyT[n1, 2] = trans
    
    Eventos
end

function recalcular_tiempo_para_colision(n1s::Array{Int64,1}, Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos)
    for n1 in n1s
        i, j = Celdas.presenteEnCelda[n1]
        for C in Celdas.vecinas[i,j]
            for n2 in Celdas.residentes[C[1],C[2]]
                if Eventos.eventoCyT[n2,1] ∈ n1s
                    Eventos = tiempos_nuevos_para_colision([n2], Esferas, Celdas, EsferasImagen, Eventos)
                end
            end
        end
    end 
    Eventos
end

function colision(n1::Int64, n2::Int64, Esferas::esferas, Celdas::celdas)
    v1 = Esferas.velocidades[n1,:] + Esferas.CR[n1].*[sign(Esferas.velocidades[n1,1]), sign(Esferas.velocidades[n1,2])]
    v2 = Esferas.velocidades[n2,:] + Esferas.CR[n2].*[sign(Esferas.velocidades[n2,1]), sign(Esferas.velocidades[n2,2])]

    Δp = Esferas.posiciones[n1,:] - Esferas.posiciones[n2,:]
    Δv = v1 - v2
    M = ((2*Esferas.masas[n1]*Esferas.masas[n2])/(Esferas.masas[n1]+Esferas.masas[n2]))
    
#     σ = Esferas.radios[n1] + Esferas.radios[n2]
#     if norm(Δp) > σ + 1e-6
#         @show norm(Δp), σ
#         error("No en contacto")
#     end
    
    J =  M*(Δv⋅Δp) / (norm(Δp))^2
    Jxy = Δp * J
    
    Esferas.velocidades[n1,:] = v1 - (Jxy/Esferas.masas[n1])
    Esferas.velocidades[n2,:] = v2 + (Jxy/Esferas.masas[n2])
       
#     @show Esferas.masas[n1]*(v1'*v1) + Esferas.masas[n2]*(v2'*v2)
#     @show Esferas.masas[n1]*(Esferas.velocidades[n1,:]'*Esferas.velocidades[n1,:]) + Esferas.masas[n2]*(Esferas.velocidades[n2,:]'*Esferas.velocidades[n2,:])
    
    Esferas                           
end

function colisionImagen(n1::Int64, n2::Int64, nEI::Int64, Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen)
    v1 = Esferas.velocidades[n1,:] + Esferas.CR[n1].*[sign(Esferas.velocidades[n1,1]), sign(Esferas.velocidades[n1,2])]
    v2 = Esferas.velocidades[n2,:] + Esferas.CR[n2].*[sign(Esferas.velocidades[n2,1]), sign(Esferas.velocidades[n2,2])]

    
    Δp = Esferas.posiciones[n1,:] - EsferasImagen.posiciones[nEI,:]
    Δv = v1 - v2
    M = ((2*Esferas.masas[n1]*Esferas.masas[n2])/(Esferas.masas[n1]+Esferas.masas[n2]))
    
#     σ = Esferas.radios[n1] + Esferas.radios[n2]
#     if norm(Δp) > σ + 1e-6
#         @show norm(Δp), σ
#         error("No en contacto")
#     end
    
    J =  M*(Δv⋅Δp) / (norm(Δp))^2
    Jxy = Δp * J
    
    Esferas.velocidades[n1,:] = v1 - (Jxy/Esferas.masas[n1])
    Esferas.velocidades[n2,:] = v2 + (Jxy/Esferas.masas[n2])
       
#     @show Esferas.masas[n1]*(v1'*v1) + Esferas.masas[n2]*(v2'*v2)
#     @show Esferas.masas[n1]*(Esferas.velocidades[n1,:]'*Esferas.velocidades[n1,:]) + Esferas.masas[n2]*(Esferas.velocidades[n2,:]'*Esferas.velocidades[n2,:])
        
    Esferas                           
end

function colision_con_pared(n1::Int64, colconpared::Int64, Esferas::esferas)
    if colconpared == -1  || colconpared == -2
        Esferas.velocidades[n1,1] = - ( Esferas.velocidades[n1,1] + sign(Esferas.velocidades[n1,1])*Esferas.CR[n1])
    elseif colconpared == -3  || colconpared == -4
        Esferas.velocidades[n1,2] = - ( Esferas.velocidades[n1,2] + sign(Esferas.velocidades[n1,2])*Esferas.CR[n1])
    end
    Esferas
end

function transferenciaCelda(n1::Int64, transferencia::Int64, Esferas::esferas, Celdas::celdas)
    i, j = Celdas.presenteEnCelda[n1]
    indiceEnCelda = 0
    
    for l in 1:length(Celdas.residentes[i,j])
        if Celdas.residentes[i,j][l] == n1
            indiceEnCelda = l
        end
    end
    
    m = Celdas.m
    
    Ts = [[i+1,j], [i-1,j], [i,j+1], [i,j-1]]
    TsCpprov = [[m+2,j], [1,j], [i,m+2], [i,1]]
    TsCp = [[2,j], [m+1,j], [i,2], [i,m+1]]
    
    L = abs(Esferas.caja[2] - Esferas.caja[1])
    
    correcionPosiciones = [[-L,0],[L,0],[0,-L],[0,L]]
    
    for T in 1:4
        if transferencia == T
            if Ts[T] != TsCpprov[T]
                ni, nj = Ts[T]
            else #condiciones periodicas
                ni, nj = TsCp[T]
                Esferas.posiciones[n1,:] = Esferas.posiciones[n1,:] + correcionPosiciones[T]
                Esferas.transferenciasPeriodicas[n1, :] += -correcionPosiciones[T]
            end
            push!(Celdas.residentes[ni,nj], Celdas.residentes[i,j][indiceEnCelda])
            deleteat!(Celdas.residentes[i,j], indiceEnCelda)
            Celdas.presenteEnCelda[n1] = [ni,nj]
        end
    end
    
    Esferas, Celdas
end

function agregar_vecinos(n1::Int64, n2::Int64, Esferas::esferas)
    Esferas.numeroColisiones[n1] += 1
    Esferas.numeroColisiones[n2] += 1
    vecino = false
    for n in Esferas.vecinos[n1]
        if n == n2
            vecino = true
        end
    end
    if !vecino 
        push!(Esferas.vecinos[n1], n2)
        push!(Esferas.vecinos[n2], n1)
        Esferas.numeroVecinos[n1] += 1
        Esferas.numeroVecinos[n2] += 1
    end
    Esferas
end

function n2esImagen(n1::Int64, n2::Int64, Celdas::celdas)
    n2esimagen = false
    C = [0,0]
    if Celdas.presenteEnCelda[n1] ∈ Celdas.fronteraPeriod
        i, j = Celdas.presenteEnCelda[n1]
        for c ∈ Celdas.vecinas[i,j]
            if c ∈ Celdas.imagen && n2 ∈ Celdas.residentes[c[1],c[2]]
                n2esimagen = true
                C = c
            end
        end
    end
    n2esimagen, C
end

function condiciones_iniciales(N::Int64, ntipos::Int64, Crecimiento::Float64, vmax::Float64, L::Array{Float64,1})
    Esferas = esferas(N, L, Crecimiento, ntipos, vmax)
    Celdas = celdas(Esferas)
    Celdas = esferasCeldasImagen(Esferas, Celdas)
    EsferasImagen = esferasImagen(Esferas, Celdas)
    Eventos = eventos(Esferas, Celdas, EsferasImagen)
    Esferas, Celdas, EsferasImagen, Eventos
end

function proximo_evento(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos)
    indiceTiempoMin = findfirst(Eventos.tiemposCyT, minimum(Eventos.tiemposCyT))
    tiempoMin = Eventos.tiemposCyT[indiceTiempoMin]
        
    Eventos.tiempo = tiempoMin
    
    n1, j = modulo(indiceTiempoMin, Esferas.N), ceil(Int64, indiceTiempoMin/Esferas.N)
    n2 = Eventos.eventoCyT[n1, j]

    
    if j == 1 #colision
        if n2 > 0 #con esfera
            Esferas = mover_esfera(n1, Esferas, Eventos)
            Esferas = mover_esfera(n2, Esferas, Eventos)
            EsferasImagen = esferasImagen(Esferas, Celdas)
            n2esimagen, C = n2esImagen(n1, n2, Celdas)
            if n2esimagen
                nEI = findindex(n2, C, Celdas, EsferasImagen)
                Esferas = colisionImagen(n1, n2, nEI, Esferas, Celdas, EsferasImagen)
            else
                Esferas = colision(n1, n2, Esferas, Celdas)
            end
            Esferas = agregar_vecinos(n1, n2, Esferas)
            
            Eventos = tiempo_nuevo_cambio_de_celda(n1, Esferas, Celdas, Eventos)
            Eventos = tiempo_nuevo_cambio_de_celda(n2, Esferas, Celdas, Eventos)
            Eventos = tiempos_nuevos_para_colision([n1,n2], Esferas, Celdas, EsferasImagen, Eventos)
            Eventos = recalcular_tiempo_para_colision([n1,n2], Esferas, Celdas, EsferasImagen, Eventos)
        else #paredes
            Esferas = colision_con_pared(n1, n2, Esferas)
            
            Eventos = tiempo_nuevo_cambio_de_celda(n1, Esferas, Celdas, Eventos)
            Eventos = tiempos_nuevos_para_colision([n1], Esferas, Celdas, EsferasImagen, Eventos)
        end
    else #tranferencia
        Esferas = mover_esfera(n1, Esferas, Eventos)
        Celdaold = Celdas.presenteEnCelda[n1]
        Esferas, Celdas = transferenciaCelda(n1, n2, Esferas, Celdas)
        Celdas = esferasCeldasImagen(Esferas, Celdas)
        EsferasImagen = esferasImagen(Esferas, Celdas)
        
        Eventos = tiempo_nuevo_cambio_de_celda(n1, Esferas, Celdas, Eventos)
        Eventos = tiempo_nuevo_para_colision(n1, Celdaold, Esferas, Celdas, EsferasImagen, Eventos)
    end
    Esferas, Celdas, EsferasImagen, Eventos
end

function packingFraction(Esferas::esferas, Eventos::eventos)
    phi = 0.0
    if Esferas.Ntipos == 1
        if Esferas.CR[1] == 0.0
            phi = Esferas.N*pi*(Esferas.radios[1])^2
        else
            phi = Esferas.N*pi*(Eventos.tiempo*Esferas.CR[1])^2
        end
    elseif Esferas.Ntipos == 2
        if Esferas.CR[1] == 0.0
            r1 = minimum(Esferas.radios)
            r2 = maximum(Esferas.radios)
            N1 = contar_valor_en_arreglo(Esferas.radios, r1)
            N2 = Esferas.N - N1
            phi = N1*pi*(r1)^2 + N2*pi*(r2)^2
        else
            cr = minimum(Esferas.CR)
            N1 = contar_valor_en_arreglo(Esferas.CR, cr)
            N2 = Esferas.N - N1
            phi = N1*pi*(Eventos.tiempo*cr)^2 + N2*pi*(1.4*Eventos.tiempo*cr)^2
        end   
    else
        for n1 in 1:Esferas.N
            phi += (Esferas.radios[n1] + (Eventos.tiempo - Esferas.tiempos[n1])*Esferas.CR[n1] )^2
        end
        π * phi# / 1^ 2
    end
end

function g(Esferas::esferas, rs::Array{Float64,1}, radio::Float64, maximo::Int64 = 20)
    g = zeros(length(rs))
    dr = rs[2] - rs[1]
    Lcaja = Esferas.caja[2] - Esferas.caja[1]
    for n1 in 1:Esferas.N-1
        for n2 in (n1+1):Esferas.N
            X = Esferas.posiciones[n1,:] - Esferas.posiciones[n2,:]
            for i in [1, 2]
                if abs(X[i]) > Lcaja / 2
                    X[i] = X[i] - sign(X[i])*Lcaja
                end
            end
            if norm(X) < maximo*radio
                j::Int64 = ceil(norm(X)/dr)
                g[j] += 2
            end
        end
    end
    for k in 1:length(rs)
        g[k] = g[k] / (Esferas.N * 2π * rs[k] )
    end
    g[1] = 0.0 #valor grande (autocorrelacion)
    g
end

function gbinaria(Esferas::esferas, rs::Array{Float64,1}, radios::Array{Float64, 1}, maximo::Int64 = 20)
    g = zeros(length(rs))
    g11 = zeros(length(rs))
    g12 = zeros(length(rs))
    g22 = zeros(length(rs))
    dr = rs[2] - rs[1]
    Lcaja = Esferas.caja[2] - Esferas.caja[1]
    for n1 in 1:Esferas.N-1
        for n2 in (n1+1):Esferas.N
            X = Esferas.posiciones[n1,:] - Esferas.posiciones[n2,:]
            for i in [1, 2]
                if abs(X[i]) > Lcaja / 2
                    X[i] = X[i] - sign(X[i])*Lcaja
                end
            end
            if norm(X) < maximo*radios[2]
                j::Int64 = ceil(norm(X)/dr)
                g[j] += 2
                if Esferas.radios[n1] == radios[1] && Esferas.radios[n2] == radios[1]
                    g11[j] += 2
                elseif Esferas.radios[n1] == radios[1] && Esferas.radios[n2] == radios[2]
                    g12[j] += 2
                elseif Esferas.radios[n1] == radios[2] && Esferas.radios[n2] == radios[1]
                    g12[j] += 2
                elseif Esferas.radios[n1] == radios[2] && Esferas.radios[n2] == radios[2]
                    g22[j] += 2
                end
            end
        end
    end
    for k in 1:length(rs)
        g[k] = g[k] / (Esferas.N * 2π * rs[k] )
        g11[k] = g11[k] / (Esferas.N * 2π * rs[k] )
        g12[k] = g12[k] / (Esferas.N * 2π * rs[k] )
        g22[k] = g22[k] / (Esferas.N * 2π * rs[k] )
    end
    g[1], g11[1], g12[1], g22[1] = (0.0, 0.0, 0.0, 0.0)
    g, g11, g12, g22
end

function angulo(n1::Int64, n2::Int64, Esferas::esferas)
    X = Esferas.posiciones[n1,:] - Esferas.posiciones[n2,:]
    Lcaja = Esferas.caja[2] - Esferas.caja[1]
    for i in [1,2]
        if abs(X[i]) > Lcaja / 2
            X[i] = X[i] - sign(X[i])*Lcaja
        end
    end
    θ = atan(X[2], X[1])
end

function anguloD(n1::Int64, n2::Int64, Esferas::esferas)
    X = Esferas.posiciones[n1,:] - Esferas.posiciones[n2,:]
    Lcaja = Esferas.caja[2] - Esferas.caja[1]
    for i in [1,2]
        if abs(X[i]) > Lcaja / 2
            X[i] = X[i] - sign(X[i])*Lcaja
        end
    end
    
    if norm(X) > 4*minimum(Esferas.radios)
        return nothing
    end
    θ = atan(X[2], X[1])
end

function psi6(Esferas::esferas) 
    PO6 = 0
    for n1 in 1:Esferas.N
        psi = 0
        for n2 in Esferas.vecinos[n1]
            psi += exp(im*6*angulo(n1, n2, Esferas))
        end
        if Esferas.numeroVecinos[n1] > 2
            psi /= Esferas.numeroVecinos[n1]
        else
            psi = 0
        end
        PO6 += norm(psi)
    end
    PO6 /= Esferas.N
    PO6 
end

function psi6D(Esferas::esferas) 
    PO6 = 0
    for n1 in 1:Esferas.N
        psi = 0
        for n2 in Esferas.vecinos[n1]
            θ = anguloD(n1, n2, Esferas)
            if typeof(θ) == Float64
                psi += exp(im*6*θ)
            else 
                psi += 0
            end
        end
        if Esferas.numeroVecinos[n1] > 2
            psi /= Esferas.numeroVecinos[n1]
        else
            psi = 0
        end
        PO6 += norm(psi)
    end
    PO6 /= Esferas.N
    PO6
end

function psi(Esferas::esferas, orden::Float64) 
    PO6 = 0
    for n1 in 1:Esferas.N
        psi = 0
        for n2 in Esferas.vecinos[n1]
            psi += exp(im*orden*angulo(n1, n2, Esferas))
        end
        if Esferas.numeroVecinos[n1] > 2
            psi /= Esferas.numeroVecinos[n1]
        else
            psi = 0
        end
        PO6 += norm(psi)
    end
    PO6 /= Esferas.N
    PO6
end

function psiD(Esferas::esferas, orden::Float64) 
    PO6 = 0
    for n1 in 1:Esferas.N
        psi = 0
        for n2 in Esferas.vecinos[n1]
            θ = anguloD(n1, n2, Esferas)
            if typeof(θ) == Float64
                psi += exp(im*orden*θ)
            else 
                psi += 0
            end
        end
        if Esferas.numeroVecinos[n1] > 2
            psi /= Esferas.numeroVecinos[n1]
        else
            psi = 0
        end
        PO6 += norm(psi)
    end
    PO6 /= Esferas.N
    PO6
end

function meanSquareDisplacement(N::Int64, x0::Array{Float64, 2}, xf::Array{Float64, 2}, periodicTransfer::Array{Int64, 2})
    MSD = 0
    for i in 1:N
        MSD += (norm(xf[i,:] + periodicTransfer[i,:] - x0[i,:]))^2
    end
    MSD /= N
    MSD
end 

function grafica(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos)
    ti = Eventos.tiempo
    l = abs(Esferas.caja[2] - Esferas.caja[1])
    Lc = Celdas.L
    phi = packingFraction(Esferas, Eventos)
    p = round(phi, digits = 3)
    time = round(ti, digits = 4)
    plot(title = "Time = $time s, phi = $p", 
        xlims = (-Lc, l+Lc), ylims = (-Lc, l+Lc), 
        
#         xlims = (0, l), ylims = (0, l), 
    xticks = -Lc:Lc:(l+Lc), yticks = -Lc:Lc:(l+Lc), aspect_ratio=:equal, leg=false)

    circulos(Esferas)
    circulosImagen(Esferas, EsferasImagen)
    
    plot!([0,l],[l,l], c=:orange) #arriba
    plot!([l,l],[l,0], c=:orange) #derecha
    plot!([0,0],[0,l], c=:orange) #izquierda
    plot!([l,0],[0,0], c=:orange) #abajo
end

function graficaST(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos)
    ti = Eventos.tiempo
    l = abs(Esferas.caja[2] - Esferas.caja[1])
    Lc = Celdas.L
    phi = packingFraction(Esferas, Eventos)
    ps6 = psi6(Esferas)
    ps6D = psi6D(Esferas)
    
    p = round(phi, digits = 3)
    ps = round(ps6, digits = 3)
    psD = round(ps6D, digits = 3)
    
    
    time = round(ti, digits = 4)
    plot(title = "phi = $p , psi6 = $ps , $psD", 
        xlims = (-Lc, l+Lc), ylims = (-Lc, l+Lc), 
        
#         xlims = (0, l), ylims = (0, l), 
        aspect_ratio=:equal, leg=false)

    circulos(Esferas)
    circulosImagen(Esferas, EsferasImagen)
    
    plot!([0,l],[l,l], c=:orange) #arriba
    plot!([l,l],[l,0], c=:orange) #derecha
    plot!([0,0],[0,l], c=:orange) #izquierda
    plot!([l,0],[0,0], c=:orange) #abajo
end

function rapidez(Esferas::esferas)
    vels = zeros(Esferas.N)
    for n1 in 1:Esferas.N
        vels[n1] = sqrt(Esferas.velocidades[n1,:]⋅Esferas.velocidades[n1,:])
    end
    vels
end

function energia(Esferas::esferas)
    Energia  = 0
    for n1 in 1:Esferas.N
        Energia += Esferas.masas[n1]*(Esferas.velocidades[n1,:]⋅Esferas.velocidades[n1,:]) / 2
    end
    Energia
end

function crecer(N::Int64, crecimiento::Float64, n::Int64, v::Float64=4.0)
    caja::Array{Float64,1} = [0.0, 1.0]
    Esferas, Celdas, EsferasImagen, Eventos = condiciones_iniciales(N, crecimiento, v, caja)
    prints = 1
    for i in 1:n
        if prints/3 < i/n && prints <= 3
            println("Avance de crecimiento a $prints / 3")
            prints += 1
        end
        Esferas, Celdas, EsferasImagen, Eventos = proximo_evento(Esferas, Celdas, EsferasImagen, Eventos)
    end
    Esferas, Celdas, EsferasImagen, Eventos
end

function reiniciar(Esferas::esferas)
    for i in 1:Esferas.N
        Esferas.numeroVecinos[i] = 0
        Esferas.numeroColisiones[i] = 0
        Esferas.vecinos[i] = []
        Esferas.transferenciasPeriodicas[i, :] = [0, 0]
    end
    Esferas
end

function parar_crecimiento(Esferas::esferas, Eventos::eventos)
    Esferas = reiniciar(Esferas)
    if Esferas.Ntipos == 1
        Esferas.radios = fill(Eventos.tiempo*maximum(Esferas.CR), Esferas.N)
    else
        for i in 1:Esferas.N
            Esferas.radios[i] = Eventos.tiempo*Esferas.CR[i]
        end
    end
    Esferas.CR = zeros(Esferas.N)
    Esferas
end

function CrecerDiscosAPhi(N::Int64, tipo::Int64, crecimiento::Float64, ϕ::Float64, v::Float64=5.0)
    caja::Array{Float64,1} = [0.0, 1.0]
    Esferas, Celdas, EsferasImagen, Eventos = condiciones_iniciales(N, tipo, crecimiento, v, caja)
    tPhi = 0
    if tipo == 1
        tPhi = (ϕ / (N*pi*crecimiento^2))^(1/2)
    else 
        N1 = contar_valor_en_arreglo(Esferas.CR, crecimiento)
        N2 = N - N1
        tPhi = (ϕ / (pi*N1*crecimiento^2 + pi*N2*(1.4*crecimiento)^2))^(1/2)
    end
    prints = 1
    while Eventos.tiempo < tPhi
        if prints*tPhi/5 < Eventos.tiempo && prints <= 5
            println("Avance de crecimiento a $prints / 5")
            prints += 1
        end
        Esferas, Celdas, EsferasImagen, Eventos = proximo_evento(Esferas, Celdas, EsferasImagen, Eventos)
    end
    
    Esferas = parar_crecimiento(Esferas, Eventos)

    Esferas, Celdas, EsferasImagen, Eventos
end

function evolucionarNumeroPasos(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos, pasos::Int64)
    for i in 1:pasos
        Esferas, Celdas, EsferasImagen, Eventos = proximo_evento(Esferas, Celdas, EsferasImagen, Eventos)
    end
    Esferas, Celdas, EsferasImagen, Eventos
end

function evolucionarFlu(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos, a = 1/4)
    # tiempo dependiente de phi o con nppp? o nuevos datos en esferas para tiempo del vecino, ir quitando los
    # hace 3 segundos dependiente a phi o con nppp?
    np::Int64 = a*Esferas.N^2
    for i in 1:np
        if i == floor(np/2)
            println("Avance de dinamica a la mitad")
        end
        Esferas, Celdas, EsferasImagen, Eventos = proximo_evento(Esferas, Celdas, EsferasImagen, Eventos)
    end
    Esferas, Celdas, EsferasImagen, Eventos
end

function histograma(Arreglo::Array{Float64, 1}, nbins = 100)
    delta = maximum(Arreglo) / nbins
    
    centros = zeros(nbins)
    conteos = zeros(nbins)
    
    start = 0
    for k in 1:nbins
        stop = start + delta
        conteos[k] = count(i -> start <= i < stop, Arreglo)
        if delta == 1
            centros[k] = start #+ delta/2.
        else
            centros[k] = start + delta/2
        end
        start = stop
   end
        
   conteos, centros 
  end

  function histogramaMono(delta::Int64, Esferas::esferas)
    #tipos == 1
    nbins::Int = Esferas.N/delta
    conteos = zeros(nbins)
    centros = zeros(nbins)
    
    start = 0
    for k in 1:nbins
        stop = start + delta
        conteos[k] = count(i -> start <= i < stop, Esferas.numeroVecinos)
        if delta == 1
            centros[k] = start #+ delta/2.
        else
            centros[k] = start + delta/2
        end
        start = stop
   end
        
   conteos, centros
  end

  function histogramaBi(delta::Int64, Esferas::esferas)
    nbins::Int = Esferas.N/delta
    conteos = zeros(nbins)
    conteosChicas = zeros(nbins)
    conteosGrandes = zeros(nbins)
    centros = zeros(nbins)
    
    start = 0
    for k in 1:nbins
        stop = start + delta
        conteos[k] = count(i -> start <= i < stop, Esferas.numeroVecinos)
        cCh = 0
        cGr = 0
        rCh = minimum(Esferas.radios)
        rGr = maximum(Esferas.radios)
        for i in 1:Esferas.N
            if start <= Esferas.numeroVecinos[i] < stop && Esferas.radios[i] == rCh
                cCh += 1
            elseif start <= Esferas.numeroVecinos[i] < stop && Esferas.radios[i] == rGr
                cGr += 1
            end
        end
        
        conteosChicas[k] = cCh
        conteosGrandes[k] = cGr
        conteos[k] == cCh + cGr ? nothing : error("Suma de chicas y grandes no es la total")
        
        if delta == 1
            centros[k] = start #+ delta/2.
        else
            centros[k] = start + delta/2
        end
        start = stop
   end
        
   conteos, conteosChicas, conteosGrandes, centros
end

function gaussDis(Nesferas::Array{Float64,1} , Nvecinos::Array{Float64,1}, promedio::Float64)
    @.  model(x, p) = 1 / (p[1]*sqrt(2*pi))*exp(- (x - p[2])^2 / 2*p[1]^2)
    p0 = [1.0, promedio]       # adivinanza inicial
    fit = curve_fit(model, Nvecinos, Nesferas, p0)
    p0 = fit.param
    p0
end

function date()
    Date = string(Dates.today())*"-"*string(Dates.hour(now()))*"h-"*string(Dates.minute(now()))*"m"
end

function datosNombreArchivo(N::Int64, dimension::Int64, tipo::Int64, phi::Float64, crecimiento::Float64)
    DatosNombreArchivo = "$dimension"*"D-"*"N_$N"*"-t_$tipo"*"-phi_$phi"*"-CR_$crecimiento"*"_"*date()
end

function medirEnergia(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos, phi::Float64, crecimiento::Float64, modNombre::String="")
    # n es numero de muestras
    N = Esferas.N
    tipo = Esferas.Ntipos
    dimension = length(Esferas.posiciones[1,:])
    
    DatosNombreArchivo = datosNombreArchivo(N, dimension, tipo, phi, crecimiento)
    NombreArchivoEnergia = "Energy"*modNombre*DatosNombreArchivo*".csv"
    
    Energia = energia(Esferas)
    VelAvg = sum(rapidez(Esferas))/ N
    
    CSV.write(NombreArchivoEnergia, DataFrame(energia = Energia, velAvg = VelAvg))
end

function medirVariacionEnergia(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos, phi::Float64, crecimiento::Float64, n::Int64)
    # n es numero de muestras
    N = Esferas.N
    tipo = Esferas.Ntipos
    dimension = length(Esferas.posiciones[1,:])
    
    DatosNombreArchivo = datosNombreArchivo(N, dimension, tipo, phi, crecimiento)
    NombreArchivoEnergia = "Energy"*DatosNombreArchivo*".csv"
    
    Energias = zeros(n)
    Energias[1] = energia(Esferas)
    Tiempos = zeros(n)
    Tiempos[1] = Eventos.tiempo
    prints = 1
    for i in 1:(n)-1
        if prints / 5 < i / (n)
            println("Avance de calculo de Energia a $prints / 5 ")
            prints += 1
        end
        Esferas, Celdas, EsferasImagen, Eventos = evolucionarNumeroPasos(Esferas, Celdas, EsferasImagen, Eventos, Esferas.N)
        Energias[i+1] = energia(Esferas)
        Tiempos[i+1] = Eventos.tiempo
    end
    
    CSV.write(NombreArchivoEnergia, DataFrames(energia = Energias, tiempo = Tiempos))
    
    plot(Tiempos, Energias)
    plot!(label = "phi = $phi", xlabel = "tiempo", ylabel = "Energia")
    savefig("E-"*DatosNombreArchivo*".png")
end

function medirDistribucionVecinosDinamicos(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos, phi::Float64, crecimiento::Float64, modNombre::String="")
    N = Esferas.N
    tipo = Esferas.Ntipos
    dimension = length(Esferas.posiciones[1,:])
    
    DatosNombreArchivo = datosNombreArchivo(N, dimension, tipo, phi, crecimiento)
    NombreArchivoHist = "Hist"*modNombre*DatosNombreArchivo*".csv"
    NombreArchivoAvg = "Avg"*modNombre*DatosNombreArchivo*".csv"
    NombreArchivoFitG = "FitG"*modNombre*DatosNombreArchivo*".csv"
    
    if tipo == 1 
        nesferas, nvecinos = histogramaMono(1, Esferas)
        AvgVecinos = sum(Esferas.numeroVecinos)/Esferas.N
        AvgColisiones = sum(Esferas.numeroColisiones)/Esferas.N
        FitGauss = gaussDis(nesferas, nvecinos, AvgVecinos)
        #Genera archivo de histograma
        CSV.write(NombreArchivoHist, DataFrame(nvecinos = nvecinos, nesferas = nesferas))
    elseif tipo == 2 
        nesferas, nesferasCh, nesferasGr, nvecinos = histogramaBi(1, Esferas)
        AvgVecinos = sum(Esferas.numeroVecinos)/Esferas.N
        AvgColisiones = sum(Esferas.numeroColisiones)/Esferas.N
        FitGauss = gaussDis(nesferas, nvecinos, AvgVecinos)
        CSV.write(NombreArchivoHist, DataFrame(nvecinos = nvecinos, nesferas = nesferas, nesferasCh = nesferasCh, nesferasGr = nesferasGr))
    else    
        error("Aun no se implementa bien polidisperso 3D")
    end
    
    #Genera archivos con promedio de vecinos y ajuste de gauss
    CSV.write(NombreArchivoAvg, DataFrame(phi = phi, AvgVecinos = AvgVecinos, AvgColisiones = AvgColisiones))
    CSV.write(NombreArchivoFitG, DataFrame(phi = phi, FitGauss = FitGauss))
        
end


function medirParametroOrdenOrientacional(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos, phi::Float64, crecimiento::Float64, modNombre::String="")
    N = Esferas.N
    tipo = Esferas.Ntipos
    dimension = length(Esferas.posiciones[1,:])
    
    DatosNombreArchivo = datosNombreArchivo(N, dimension, tipo, phi, crecimiento)
    NombreArchivoPO = "PO"*modNombre*DatosNombreArchivo*".csv"

    if dimension == 2
        PO = psi6(Esferas)
        CSV.write(NombreArchivoPO, DataFrame(phi = phi, psi6 = PO))
        
        NombreArchivoPOX = "POX"*modNombre*DatosNombreArchivo*".csv"
        O = 0:0.01:15
        POX = zeros(length(O))
        for i in 1:length(O)
            POX[i] = psi(Esferas, O[i])
        end
        CSV.write(NombreArchivoPOX, DataFrame(orden = O, psix = POX))
        plot(title = "Con $tipo radios, $N discos",
            xlabel = "Orden",  ylabel = "psi 6", xtick=0:1:15, label = "phi = $phi")
        plot!(O, POX)
        savefig("pox-"*modNombre*DatosNombreArchivo*".png")
    elseif dimension == 3 
    PO1 = Q6(Esferas)
    PO2 = Q6V2(Esferas)
    CSV.write(NombreArchivoPO, DataFrame(phi = phi, Q6 = PO1, Q6V2 = PO2))
    end
end

function medirDistribucionRadial(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos, phi::Float64, crecimiento::Float64, modNombre::String="", maximo::Int64 = 20)
    N = Esferas.N
    tipo = Esferas.Ntipos
    dimension = length(Esferas.posiciones[1,:])
    
    DatosNombreArchivo = datosNombreArchivo(N, dimension, tipo, phi, crecimiento)
    NombreArchivoGrS = "GrS"*modNombre*DatosNombreArchivo*".csv" 
    NombreArchivoGr = "Gr"*modNombre*DatosNombreArchivo*".csv" 

    radioMax = Esferas.radios[1]
    
    println("Empezando calculo de g(r)")

    dr = 1 / (2*N)
    rs = collect(0:dr:maximo*radioMax)
    G = g(Esferas, rs, radioMax, maximo) 
    CSV.write(NombreArchivoGrS, DataFrame(distancia = rs./radioMax, gr = G))

    prints = 1
    for i in 1:(N/2)-1
        if prints / 5 <= i / ((N/2) - 1)
            println("Avance de g(r) a $prints / 5 ")
            prints += 1
        end
        Esferas, Celdas, EsferasImagen, Eventos = evolucionarNumeroPasos(Esferas, Celdas, EsferasImagen, Eventos, N)
        Esferas = mover_esferas(Esferas, Eventos)
        EsferasImagen = esferasImagen(Esferas, Celdas)
        G += g(Esferas, rs, radioMax, maximo)
    end
    G /= (N/2)
    rs = rs ./ radioMax
    
    CSV.write(NombreArchivoGr, DataFrame(distancia = rs, gr = G))
    
    plot(rs, G)
    plot!(xlim=(0,10), label = "phi = $phi", xlabel = "r", ylabel = "g(r)")
    savefig("gr"*modNombre*"-"*DatosNombreArchivo*".png")
end

function medirDistribucionRadialBinario(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos, phi::Float64, crecimiento::Float64,  modNombre::String="", maximo::Int64 = 20)
    N = Esferas.N
    tipo = Esferas.Ntipos
    dimension = length(Esferas.posiciones[1,:])
    
    DatosNombreArchivo = datosNombreArchivo(N, dimension, tipo, phi, crecimiento)
    NombreArchivoGrS = "GrS"*modNombre*DatosNombreArchivo*".csv" 
    NombreArchivoGr = "Gr"*modNombre*DatosNombreArchivo*".csv" 
    
    dr = 1 / (2*N)
    radios = [minimum(Esferas.radios), maximum(Esferas.radios)]
    rs = collect(0:dr:maximo*radios[2])
    G, G11, G12, G22 = gbinaria(Esferas, rs, radios, maximo) 
    CSV.write(NombreArchivoGrS, DataFrame(distancia = rs1, g = G, g11 = G11, g12 = G12, g22 = G22))
       
    prints = 1
    for i in 1:(N/2)-1
        if prints / 5 <= i / ((N/2) - 1)
            println("Avance de g(r) a $prints / 5 ")
            prints += 1
        end
        Esferas, Celdas, EsferasImagen, Eventos = evolucionarNumeroPasos(Esferas, Celdas, EsferasImagen, Eventos, N)
        Esferas = mover_esferas(Esferas, Eventos)
        EsferasImagen = esferasImagen(Esferas, Celdas)
        
        Gs, G11s, G12s, G22s = gbinaria(Esferas, rs, radios, maximo)
        G += Gs
        G11 += G11s
        G12 += G12s
        G22 += G22s
    end
    G /= (N/2)
    G11 /= (N/2)
    G12 /= (N/2)
    G22 /= (N/2)
    rs1 = rs ./ radios[1]
    rs2 = rs ./ radios[2]
    
    CSV.write(NombreArchivoGr, DataFrame(distancia = rs1, g = G, g11 = G11, g12 = G12, g22 = G22))
    
    plot(rs1, G, label = "g")
    plot!(rs1, G11, label = "g11")
    plot!(rs1, G12, label = "g12")
    plot!(rs1, G22, label = "g22")
    plot!(xlim=(0,10), title = "phi = $phi", xlabel = "r", ylabel = "g(r)")
    savefig("gr"*modNombre*"-"*DatosNombreArchivo*".png")
end

function medirDifusion(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos, phi::Float64, crecimiento::Float64)
    N = Esferas.N
    tipo = Esferas.Ntipos
    dimension = length(Esferas.posiciones[1,:])
    DatosNombreArchivo = datosNombreArchivo(N, dimension, tipo, phi, crecimiento)
    NombreArchivoD = "Dif"*DatosNombreArchivo*".csv"
    NombreArchivoPs = "Pos12enT"*DatosNombreArchivo*".csv"

    posicionesIniciales = copy(Esferas.posiciones)
    tiempoInicial = Eventos.tiempo
    
    MSDS = Float64[]
    tiempos = Float64[]
    colisiones = Float64[]
    # posP1 = Esferas.posiciones[1,:]'
    # posP2 = Esferas.posiciones[2,:]'
    # tsP12 = [Eventos.tiempo]

    prints = 1
    for i in 1:8N
        if prints / 5 < i / 8N
            println("Avance de calculo de difusion a $prints / 5 ")
            prints += 1
        end
        # for j in 1:4
            Esferas, Celdas, EsferasImagen, Eventos = evolucionarNumeroPasos(Esferas, Celdas, EsferasImagen, Eventos, convert(Int64, N))
            Esferas = mover_esferas(Esferas, Eventos)
            EsferasImagen = esferasImagen(Esferas, Celdas)

        #     posP1 = vcat(posP1, Esferas.posiciones[1,:]' +  periodicTransfer[1,:]')
        #     posP2 = vcat(posP2, Esferas.posiciones[2,:]' +  periodicTransfer[2,:]')
        
        #     tsP12 = push!(tsP12, Eventos.tiempo)
        # end
        posicionesFinales = copy(Esferas.posiciones)
        dt = Eventos.tiempo - tiempoInicial
        MSD = meanSquareDisplacement(N, posicionesIniciales, posicionesFinales, Esferas.transferenciasPeriodicas)
        
        push!(MSDS, MSD)
        push!(tiempos, dt)
        push!(colisiones, sum(Esferas.numeroColisiones)/Esferas.N)
    end
    CSV.write(NombreArchivoD, DataFrame(MSD = MSDS, tiempo = tiempos, colisiones = colisiones ))
    # MSDS, tiempos
end

function medirLineasDeTiempo(Esferas::esferas, Celdas::celdas, EsferasImagen::esferasImagen, Eventos::eventos, phi::Float64, crecimiento::Float64, n::Int64)
    #     n es numero de muestras de la serie
    N = Esferas.N
    tipo = Esferas.Ntipos
    dimension = length(Esferas.posiciones[1,:])
    DatosNombreArchivo = datosNombreArchivo(N, dimension, tipo, phi, crecimiento)
    NombreArchivoLT = "LT"*DatosNombreArchivo*".csv"

    AvgVecinos = Float64[]
    AvanceTiempo = Float64[]
    NumeroPasos = Int64[]
    
    prints = 1
    pasos = convert(Int64, N / 100)
    for i in 1:n
        if prints / 5 < i / n 
            println("Avance de linea de tiempo a $prints / 5 ")
            prints += 1
        end
        Esferas, Celdas, EsferasImagen, Eventos = evolucionarNumeroPasos(Esferas, Celdas, EsferasImagen, Eventos, 100)
        push!(AvgVecinos, sum(Esferas.numeroVecinos)/Esferas.N)
        push!(AvanceTiempo, Eventos.tiempo)
        push!(NumeroPasos, pasos*i)
    end
    CSV.write(NombreArchivoLT, DataFrame(Repeticion = collect(1:n), AvgVecino = AvgVecinos,  numeroPasos = NumeroPasos,  avanceTiempo = AvanceTiempo))
end

function medirCasiTodo(N::Int64, tipo::Int64, crecimiento::Float64, phi::Float64)
    dimension = 2
    DatosNombreArchivo = datosNombreArchivo(N, dimension, tipo, phi, crecimiento)
    
    println("Creciendo particulas con phi = $phi , con $tipo radios en $dimension D")
    Esferas, Celdas, EsferasImagen, Eventos = CrecerDiscosAPhi(N, tipo,  crecimiento, phi)
    println("Termalizando")
    Esferas, Celdas, EsferasImagen, Eventos = evolucionarFlu(Esferas, Celdas, EsferasImagen, Eventos, 0.1)
    (1500^2)*0.1
    
    Esferas = reiniciar(Esferas)
    medirEnergia(Esferas, Celdas, EsferasImagen, Eventos, phi, crecimiento, "I")
    println("Haciendo dinamica corta con phi = $phi , con $tipo radios  en $dimension D")
    Esferas, Celdas, EsferasImagen, Eventos = evolucionarFlu(Esferas, Celdas, EsferasImagen, Eventos, 0.2)
    println("Obteniendo distribucion de vecinos dinamicos tiempo cortito")
    medirDistribucionVecinosDinamicos(Esferas, Celdas, EsferasImagen, Eventos, phi, crecimiento, "MC")
    
    println("Haciendo dinamica normal con phi = $phi , con $tipo radios  en $dimension D")
    Esferas, Celdas, EsferasImagen, Eventos = evolucionarFlu(Esferas, Celdas, EsferasImagen, Eventos, 0.3)
    println("Obteniendo distribucion de vecinos dinamicos tiempo corto")
    medirDistribucionVecinosDinamicos(Esferas, Celdas, EsferasImagen, Eventos, phi, crecimiento, "C")

    println("Haciendo dinamica normal con phi = $phi , con $tipo radios  en $dimension D")
    Esferas, Celdas, EsferasImagen, Eventos = evolucionarFlu(Esferas, Celdas, EsferasImagen, Eventos, 0.5)
    println("Obteniendo distribucion de vecinos dinamicos tiempo normal")
    medirDistribucionVecinosDinamicos(Esferas, Celdas, EsferasImagen, Eventos, phi, crecimiento)

    Esferas = mover_esferas(Esferas, Eventos)
    EsferasImagen = esferasImagen(Esferas, Celdas)

    graficaST(Esferas, Celdas, EsferasImagen, Eventos)
    savefig(DatosNombreArchivo*".png")
    savefig(DatosNombreArchivo*".pdf")
    
    println("Obteniendo parametro de orden")
    medirParametroOrdenOrientacional(Esferas, Celdas, EsferasImagen, Eventos, phi, crecimiento, "I")
    println("Obteniendo g(r)")
    if tipo == 1 
        medirDistribucionRadial(Esferas, Celdas, EsferasImagen, Eventos, phi, crecimiento, "I")
    elseif tipo == 2
        medirDistribucionRadialBinario(Esferas, Celdas, EsferasImagen, Eventos, phi, crecimiento, "I")
    end
    
    Esferas = reiniciar(Esferas)
    println("Midiendo difusion")
    medirDifusion(Esferas, Celdas, EsferasImagen, Eventos, phi, crecimiento)
    
    Esferas = reiniciar(Esferas)
    println("Haciendo dinamica larga con phi = $phi , con $tipo radios  en $dimension D")
    Esferas, Celdas, EsferasImagen, Eventos = evolucionarFlu(Esferas, Celdas, EsferasImagen, Eventos, 2)
    println("Obteniendo distribucion de vecinos dinamicos tiempo largo")
    medirDistribucionVecinosDinamicos(Esferas, Celdas, EsferasImagen, Eventos, phi, crecimiento, "L")
    medirEnergia(Esferas, Celdas, EsferasImagen, Eventos, phi, crecimiento, "F")
    println("Obteniendo parametro de orden")
    medirParametroOrdenOrientacional(Esferas, Celdas, EsferasImagen, Eventos, phi, crecimiento, "F")
    println("Obteniendo g(r)")
    if tipo == 1 
        medirDistribucionRadial(Esferas, Celdas, EsferasImagen, Eventos, phi, crecimiento, "F")
    elseif tipo == 2
        medirDistribucionRadialBinario(Esferas, Celdas, EsferasImagen, Eventos, phi, crecimiento, "F")
    end
    # println("Obteniendo lineas de tiempo")
    # medirLineasDeTiempo(Esferas, Celdas, EsferasImagen, Eventos, phi, crecimiento, 3_000_000)
    println("Termino")
end

function medicionesAldo(N::Int64, tipo::Int64, crecimiento::Float64, phi::Float64)
    dimension = 2
    DatosNombreArchivo = datosNombreArchivo(N, dimension, tipo, phi, crecimiento)
    
    println("Creciendo particulas con phi = $phi , con $tipo radios en $dimension D")
    Esferas, Celdas, EsferasImagen, Eventos = CrecerDiscosAPhi(N, tipo,  crecimiento, phi)
    println("Termalizando")
    Esferas, Celdas, EsferasImagen, Eventos = evolucionarFlu(Esferas, Celdas, EsferasImagen, Eventos, 0.1)
    Esferas = reiniciar(Esferas)
    
    println("Haciendo dinamica corta con phi = $phi , con $tipo radios  en $dimension D")
    Esferas, Celdas, EsferasImagen, Eventos = evolucionarFlu(Esferas, Celdas, EsferasImagen, Eventos, 0.2)

    Esferas = mover_esferas(Esferas, Eventos)
    EsferasImagen = esferasImagen(Esferas, Celdas)

    # graficaST(Esferas, Celdas, EsferasImagen, Eventos)
    # savefig(DatosNombreArchivo*".png")
    # savefig(DatosNombreArchivo*".pdf")

    println("Termino")

    return Esferas.vecinos, Esferas.posiciones, Esferas.radios
end

                # N 2radios cr  phi
neigs, pos, rs = medicionesAldo(parsed_args["Ndiscs"], 2, 8e-3, parsed_args["phi"])

function adjMat(neighList)
    n = length(neighList)
    mat = zeros(Int16,n,n)
    for i in 1:n
        mat[neighList[i],i] .= 1
    end
    return mat
end

DelimitedFiles.writedlm(string("../networks/","mat_",parsed_args["save"],".csv"),adjMat(neigs),',')
DelimitedFiles.writedlm(string("../networks/","pos_",parsed_args["save"],".csv"),pos,',')
DelimitedFiles.writedlm(string("../networks/","rs_",parsed_args["save"],".csv"),rs,',')

#open("mat.csv","w+") do io
#    DelimitedFiles.writedlm(adjMat(neigs),)
#end