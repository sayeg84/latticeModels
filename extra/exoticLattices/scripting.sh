# building quasiperiodic lattices
symetries=(7 10 13 23 31 53 73 89 101 113)
iterations=(3 3 3 4 4 5 7 8 10 11)
echo "quasiperiodic lattices"
for ((i=0; i < ${#symetries[@]}; i ++ )); do
    echo $i ;
    julia makeQPLattice.jl --Nsides ${symetries[i]} --Iterations ${iterations[i]} --save "QP$i"
done

# building jamming/random lattices
echo "jamming lattices"
sizes=(200 300 400 500 600 700 800 900 1000)
for ((i=0; i < 10; i ++ )); do
    echo $i ; 
    julia makeJammingLattice.jl --Ndiscs ${sizes[i]} --save "jam$i" 
done

