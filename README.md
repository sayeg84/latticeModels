# Lattice Models

Code for making Ising Model-like simulations using Metropolis-Hastings Algorithm and Wang-Landau algorithm.

## Instalation

Code for simulations is written in Julia and plots and data analysis is made on Python. Dependencies can be checked on the |`.toml` files for Julia (which can also be used to build the project) and in `REQUIREpython` for Python.

For instalation, run `install.sh` file on a terminal. To test instalation, run `test.sh`

## Documentation

All of the documentation regarding the functions can be checked in the sources files inside `src` directory.

## Running

The scripts that run simulations are the `src/run.jl` and `src/runDistributed.jl` files. By default, `src/runDistributed.jl` tries to use all of the cores available in the machine to parallelize the simulation. For more information regarding their use, run them with the `--help` flag.

## Outputs

All of the simulations are stored in a folder with the name of the date and time they were performed. Inside the each folder, there will be the following files:

|name|descrition|
|:-:|:-:|
|`adjMat.csv`|Adjacency matrix of the network of the simulation|
|`algoParam.csv`|Parameters of the algorithm used (number of steps, name)|
|`metaParam.csv`|Parameters of the network, model  the energy function|
|`simulParamDict.csv`|Table with all the values of the energy function parameters that were simulated and the order in which they were performed. Order is important for hystheresis analysis|

The files beggining with `i_j` represent the `j` independent iteration of the simulation ran with the energy parameters of the `i+1` csv of `simulParamDict.csv` . Their meaning is the following.

|name|descrition|
|:-:|:-:|
|`i_j_initial.csv`|Linearized array of the sites values at the beggining of each simulation.|
|`i_j_changes.csv`|Linearized array representing the changes at each MC step. If `k=0`, no change was made at that step. If `k!=0`, the site `k` was changed in that moment|
|`i_j_mag.csv`,`i_j_ener.csv`|Energy and magnetization values in each time step|

### Warning

**If the `outputs` folder hasn't been created, if is possible that the first run of `src/runDistributed.jl` will throw an error**.
