sizes=(8 10 12 16 20 26 32)
mus=(-0.5  -2.0 -10.0 )
temps=(0.2 0.5 0.8 1.1 1.4 1.7 2.0 2.3 2.6 2.9)
                                                                                                                                               
for si in "${sizes[@]}"  ; do
    for mu in "${mus[@]}" ; do
       julia runTempHystSet.jl -N $si -A 10 -M LatticeGas -E  PenalizedEnergy -S 2.5 -T 2.3 -B $mu -C 0.9
       sleep 60
    done
    #for te in "${temps[@]}"  ; do
    #    julia runMuHystSet.jl -N $si -A 10 -M LatticeGas -E  PenalizedEnergy -S 2.5 -T $te -C 0.9
    #    sleep 60
    #done
done

for folder in ../outputs/*/ ; do
    echo $folder
    if [[ -d "$folder/endothermic" ]]
    then
        echo "temp cycle"
        julia analMet.jl --path "$folder/endothermic" --ncut 0.5
        julia analMet.jl --path "$folder/exothermic" --ncut 0.5
       python3 tempHystPlots.py --path "$folder"
    elif [[ -d "$folder/mu-increasing" ]]
    then
        echo "mu cycle"
        julia analMet.jl --path "$folder/mu-increasing" --ncut 0.5
        julia analMet.jl --path "$folder/mu-decreasing" --ncut 0.5
       python3 mu-cHystPlots.py --path "$folder"
    else
        echo "normal"
        julia analMet.jl --path "$folder" --ncut 0.5
        julia analMet.jl --path "$folder" --ncut 0.5
    fi
done
tar cvfz simuls.tar.gz ../outputs





#for folder in ../../simulacionesValiosas/ultimas/*/ ; do
#    echo $folder
#    if [[ -d "$folder/endothermic" ]]
#    then
#        
#        python3 tempHystPlots.py --path "$folder"
#    else [[ -d "$folder/mu-increasing" ]]
#        python3 mu-cHystPlots.py --path "$folder"
#    fi
#done