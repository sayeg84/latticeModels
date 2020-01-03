

sizes=(8 10 12 16 20 26 32)
#sizes=(16)
for si in "${sizes[@]}"  ; do
    julia runTempHystSet.jl -N $si -A 15 -M LatticeGas -E PenalizedEnergy -S 2.5 
    sleep 60
done
for folder in ../outputs/*/ ; do
    echo $folder
    julia analMet.jl --path "$folder/endothermic" --ncut 0.5
    julia analMet.jl --path "$folder/exothermic" --ncut 0.5
done

