

#sizes=(8 10 12 16 20 26 32 40 50 64)
sizes=(16 20)
for si in "${sizes[@]}"  ; do
    julia runAtaIdea.jl -N $si -A 15 
    sleep 60
done
for folder in ../outputs/*/ ; do
    echo $folder
    julia analysisMetropolis.jl --path "$folder" --ncut 0.5
done

