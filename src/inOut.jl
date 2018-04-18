module InOut
export MakeDirectories, MakePlots, MakeTable
    using Plots, Rsvg
    plotlyjs()
    dest=Dates.format(Dates.now(),"dd-mm-Y HH:MM:SS")
    function MakeDirectories()
        cd("..")
        if ~(isdir("outputs"))
            mkdir("outputs")
        end
        cd("outputs")
        mkdir(dest)
        cd("..")
        cd("src")
    end
    function MakePlots(t,a1,a2,a3,param)
        cd("..")
        cd("outputs")
        cd(dest)
        Plots.scatter(t,a1,
            title="Magnetizacion",
            xlabel="Temperatura (KB=1)",
            ylabel="Magnetizacion absoluta promedio por spin",
            label="simulación con N=$(param[1]^2)"
        )
        #Plots.plot!(T,[magOnsager(t) for t in T],label="Solución teórica")
        Plots.plot!([2.3,2.3],[0,1],linestyle=:dash,label="Temperatura crítica     ")
        Plots.savefig("mag.png")
        #display(p)


        scatter(t,a2,
            title="Energia",
            xlabel="Temperatura (KB=1)",
            ylabel="Energía promedio por spin",
            label="simulación con N=$(param[1]^2)"
        )

        #Plots.plot!(T,[energiaOnsager(t) for t in T],label="Solución teórica")
        Plots.plot!([2.3,2.3],[-2,0],linestyle=:dash,label="Temperatura crítica     ")
        Plots.savefig("ener.png")
        #display(p)



        scatter(t,a3,
            ylim=(0.0,maximum(a3[5:end])),
            title="Capacidad calorífica",
            xlabel="Temperatura (KB=1)",
            ylabel="Capacidad calorífica por spin",
            label="simulación con N=$(param[1]^2)"
        )
        #Plots.plot!(T,[cvOnsager(t) for t in T],label="Solución teórica")
        Plots.plot!([2.3,2.3],[0,2],linestyle=:dash,label="Temperatura crítica      ")
        Plots.savefig("cv.png")
        cd("..")
        cd("..")
        cd("src")
    end

    function MakeTable(t,a1,a2,a3,param,time)
        cd("..")
        cd("outputs")
        cd(dest)
        write("datos.csv")
        open("datos.csv","w") do f 
            write(f,"Modelo de ising,Fecha, $(dest) \n")
            write(f," ,Tiempo de ejecución, $(round(time,2)) \n")
            write(f," Parámetros \n")
            write(f," n, J, B, pasos, frecuencia, bins \n ")
            write(f," $(param[1]), $(param[2]), $(param[3]), $(param[4]), $(param[5]), $(param[6]) \n ")
            #write(f," $n, $J, $B, $(pasos), $(freq) \n ")
            write(f,"\n")
            write(f,"T,M,E,CV \n")
            for i in 1:length(t)
                write(f,"$(t[i]), $(a1[i]),$(a2[i]), $(a3[i]) \n")
            end
        end
        cd("..")
        cd("..")
        cd("src")
    end
    function MakeDOSTable(s,energyIntervals,param,time)
        cd("..")
        cd("outputs")
        cd(dest)
        write("DOS.csv")
        open("DOS.csv","w") do f 
            write(f,"Densidad de estados")
            write(f," ,Tiempo de ejecución, $(round(time,2)) \n")
            write(f," Parámetros \n")
            write(f," n, J, B, pasos, frecuencia, bins, penalización \n ")
            write(f," $(param[1]), $(param[2]), $(param[3]), $(param[4]), $(param[5]), $(param[6]), $(param[7]) \n ")
            write(f," E_i, E_f, S, DOS \n  ")
            for i in 1:length(s)
                write(f,"$(energyIntervals[i]), $(energyIntervals[i+1]), $(s[i]), $(round(exp(big(s[i])),5)) \n")
            end
        end
        cd("..")
        cd("..")
        cd("src")
    end
end