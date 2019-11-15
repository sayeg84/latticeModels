

sizes=(8 10 12 16 20 26 32 40 50 64)
for si in "${sizes[@]}"  ; do
    julia runParallelSet.jl -N $si -A 3 -S 3.5
    sleep 60
done
for folder in ../outputs/*/ ; do
    echo $folder
    julia analysisMetropolis.jl --path "$folder" --ncut 0.75
done

