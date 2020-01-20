sizes=(8 10 12 16 20 26 32)
mus=(-0.5  -2.0 -10.0 )
#sizes=(16)                                                                                                                                                                            
for si in "${sizes[@]}"  ; do
    for mu in "${mus[@]}" ; do
       julia runTempHystSet.jl -N $si -A 10 -M LatticeGas -E PenalizedEnergy -S 2.5 -T 2.3 -B $mu -C 0.9
       sleep 60
    done
done
for folder in ../outputs/*/ ; do
    echo $folder
    if [[ -d "$folder/endothermic" ]]
    then
        echo "temp cycle"
        julia analMet.jl --path "$folder/endothermic" --ncut 0.5
        julia analMet.jl --path "$folder/exothermic" --ncut 0.5
    elif [[ -d "$folder/mu-increasing" ]]
    then
        echo "mu cycle"
        julia analMet.jl --path "$folder/mu-increasing" --ncut 0.5
        julia analMet.jl --path "$folder/mu-decreasing" --ncut 0.5
    else
        echo "normal"
        julia analMet.jl --path "$folder" --ncut 0.5
        julia analMet.jl --path "$folder" --ncut 0.5
    fi
done
tar cvfz simuls.tar.gz ../outputs




