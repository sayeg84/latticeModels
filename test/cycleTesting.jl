#=
    m1=[
    0 0 0 0 1 1 1 1 1 ;
    0 0 0 0 1 0 0 0 1 ;
    0 1 1 1 1 0 1 0 1 ;
    0 0 0 0 0 0 0 0 1 ;
    0 0 0 0 0 0 0 0 0 ;
    ]
    m2=[
    0 0 0 0 1 1 1 0 0 ;
    0 1 1 1 1 0 1 0 0 ;
    0 0 0 0 1 1 1 1 1 ;
    0 0 0 0 0 0 0 0 1 ;
    0 0 0 0 0 0 0 0 0 ;
    ]
    m3=[
    0 0 0 0 1 1 1 0 0 ;
    0 2 2 2 1 0 1 0 0 ;
    0 0 0 0 2 1 1 2 2 ;
    0 0 0 0 0 0 0 0 2 ;
    0 0 0 0 0 0 0 0 0 ;
    ]
    
    m4=[
    0 0 0 0 1 1 1 0 0 ;
    0 2 2 2 1 0 1 0 0 ;
    0 0 0 0 1 1 1 2 2 ;
    0 0 0 0 1 0 1 0 2 ;
    0 0 0 0 1 1 1 0 0 ;
    ]
 
    #l=RemoveTails(m2)
    #x=WalkPath([3,5],[2,5],m3)

    
    m3=[
    0 0 0 0 1 1 1 1 0 0 ;
    0 1 1 1 1 0 0 1 0 0 ;
    0 0 0 0 1 1 1 1 1 1 ;
    0 0 0 0 1 0 0 1 0 1 ;
    0 0 0 0 1 1 1 1 0 1 ;
    0 0 0 0 0 0 0 0 0 0 ;
    ]
   
    m3=[
        0 1 1 1 1 1 1 1 1 0 0 ;
        0 0 0 0 1 0 0 0 1 0 0 ;
        0 0 1 1 1 1 1 1 1 1 1 ;
        0 0 1 0 1 0 1 0 1 0 1 ;
        0 0 1 0 1 1 1 1 1 0 1 ;
        0 0 1 0 0 0 1 0 0 0 0 ;
        0 0 1 1 1 1 1 0 0 0 0 ;
        0 0 0 0 0 0 0 0 0 0 0 ;
        ]
   

m3= [
    0 0 0 0 0 0 0 0 0 0 0 0 ;
    0 1 1 1 0 0 1 1 1 0 0 0 ;
    0 1 0 1 1 1 1 0 1 0 0 0 ;
    0 1 1 1 0 0 1 1 1 0 0 0 ;
    0 0 0 0 0 0 0 0 0 0 0 0 ;
    
    ]

    #neigLatt=Auxiliar.NeighborIndexLattice(m3,Auxiliar.SquareLatticeNeighborsIndex)
    #println(RemoveOthers(m3,-1))
    #println(RemoveTails(RemoveOthers(m3,-1),neigLatt))
    neigLatt=Auxiliar.NeighborIndexLattice(m3,Auxiliar.SquareLatticeNeighborsIndex)
    @time a=Cycles(m3,neigLatt,printLog=false)
    println(a)

    for i in 1:2000
         @time x=Cycles(m3,neigLatt,printLog=false)
         if x!=a 
            println("fail")
         end
     end
  
   =#

    #@time x=WalkComplicatedPath(m3,pos,neigLatt)
    #println(x[1])
    #println(x[2][CartesianIndex((3,6))])
    #@time y=WalkComplicatedPath(m3,pos,neigLatt)
    #println(y[1])
    #println(y[2][CartesianIndex((3,6))])
    #@time z=WalkComplicatedPath(m3,pos,neigLatt)
    #println(z[1])
    #println(z[2][CartesianIndex((3,6))])
    #@time x=WalkSimplePath(m3,pos,neigLatt)
    #println(x)
    #pos=CartesianIndex((3,5))
    #@time x=WalkPath(m3,pos,neigLatt)
    #println(x)
    #=
    g=Erdos.Graph(ConstructAdjacencyMatrix(m))
    x=Erdos.conected_components(g)
    =#