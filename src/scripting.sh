#sizes=(8 10 12 16 20 26 32)
sizes=(50 64 80 100)
#mus=(-0.5  -2.0 -10.0 )
#temps=(0.2 0.5 0.8 1.1 1.4 1.7 2.0 2.3 2.6 2.9)



#for si in "${sizes[@]}"  ; do
#    for mu in "${mus[@]}" ; do
#       julia runTempHystSet.jl -N $si -A 10 -M LatticeGas -E  PenalizedEnergy -S 2.5 -T 2.3 -B $mu -C 0.9
#       sleep 60
#    done
#    for te in "${temps[@]}"  ; do
#        julia runMuHystSet.jl -N $si -A 10 -M LatticeGas -E  PenalizedEnergy -S 2.5 -T $te -C 0.9
#        sleep 60
#    done
#done

julia runParSet.jl --Barray "0,-3.5,21" --Carray "0.9,0.9,1" --Jarray "2,2,1" --kTarray "0.5,0.5,1" --order "B,kT,C,J" -N 32 -S 2.0 -M LatticeGas -E PenalizedEnergy 

#julia runParSet.jl --Barray "-3.5,0,21" --Carray "0.9,0.9,1" --Jarray "2,2,1" --kTarray "0.5,0.5,1" --order "B,kT,C,J" -N 32 -S 2.0 -M LatticeGas -E PenalizedEnergy 

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
    fi
done
tar cvfz simuls.tar.gz ../outputs
