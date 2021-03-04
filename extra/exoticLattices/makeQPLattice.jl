include("basicOperations.jl")
include("dualNetwork.jl")
include("mainClusterIteratedAlgorithm.jl")
include("neighborhoods.jl")
include("separateXY.jl")
include("voronoi.jl")
include("voronoiAlgorithm.jl")

using Plots, DelimitedFiles, ArgParse

function ParseCommandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "-N","--Nsides"
            help = "Symmetry of the quasiperiodic lattice"
            arg_type = Int64
            default = 10
        "-E","--Nerror"
            help = "Margin of error in the integer numbers generated as the projection of the arbitrary point into the star vectors"
            arg_type = Int64
            default = 2
        "-I","--Iterations"
            help = "Number of iterations for the area algorithm. Higher values for higher symmetry is required. Value too large will cause error"
            arg_type = Int64
            default = 2
        "-A","--Alpha"
            help = "Constant for algorithm"
            arg_type = Float64
            default = 0.2
        "-S","--save"
            help = "Save prefix for filenames"
            arg_type = String
            default = "qp"
    end
    return parse_args(s)
end

const parsed_args = ParseCommandline()

function qpLattice(Nsides,Nerror;Iterations=5,Alpha=0.2)
    #Nsides = 10; #Symmetry of the quasiperiodic lattice
    Star_Vectors = [[BigFloat(1),0]]; #Array wich will contain the Star Vectors
    for i in 1:(Nsides-1)
        push!(Star_Vectors, [cos((2*i)*pi/Nsides), sin((2*i)*pi/Nsides)]); #Complete the Star_Vector Array
    end
    Alphas_Array = fill(Alpha, Nsides); #Array of the alphas constant
    Average_Distance_Stripes = fill(Nsides/2, Nsides); #Array with the average distance between stripes
    #Nerror = 2; #Margin of error in the integer numbers generated as the projection of the arbitrary point into the star vectors
    SL = 1e6; #Size of a half side of the square in which the algorithm generate a random point inside it

    #Let's generate the arbitrary point around which we will generate the neighborhood
    APoint = Float64[]; #An Float64 array that will held the coordinates of the arbitrary point

    #Generate two random numbers that will determine in which cuadrant will be the arbitrary point
    x = rand();
    y = rand();

    if (x > 0.5) && (y > 0.5)
        APoint = [rand()*SL, rand()*SL];
    elseif (x > 0.5) && (y < 0.5)
        APoint = [rand()*SL, -rand()*SL];
    elseif (x < 0.5) && (y > 0.5)
        APoint = [-rand()*SL, rand()*SL];
    elseif (x < 0.5) && (y < 0.5)
        APoint = [-rand()*SL, -rand()*SL];
    end

    #Let's get the coordinates of the vertices of the polygons, also get the coordinates of the point around which
    #the neighborhood was generated
    Dual_Points = region_Local_Voronoi(Nerror, Average_Distance_Stripes, Star_Vectors, Alphas_Array, APoint);

    #Let's get the coordinates as tuples of the centroids and the dictionary that relates the centroid's coordinates 
    #with the polygons vertices' coordinates of the polygon that generate the centroid.
    Centroids, Centroids_Dictionary = centroides(Dual_Points);

    #Iterations = 5; #Number of iterations of the areas algorithm desired.
    Bounded_Area = 1.2; #The value for the area of the polygons that will be a discriminator value in the areas algorithm

    #Let's get the centroids that remains after the iterations of the areas algorithm
    Inside_Clusters_Centroids = centroides_Area_Acotada_Iterada(Centroids, Bounded_Area, Iterations);

    #Let's get the X and Y coordinates of the vertices of the retained polygons in the quasiperiodic lattice
    return centroides_A_Vertices(Inside_Clusters_Centroids, Centroids_Dictionary);

end

#Let's visualize the new neighborhood of the quasiperiodic lattice
#Beware: Plotting this image consume a lot of RAM.
#=
plot()
for i in 1:4:length(X)
    plot!([X[i],X[i+1],X[i+2],X[i+3],X[i]],[Y[i],Y[i+1],Y[i+2],Y[i+3],Y[i]], markersize = 0.2, key = false, aspect_ratio=:equal, grid = false, color =:black, xaxis = nothing, yaxis = nothing)
end
scatter!([APoint[1]], [APoint[2]], color = "red", xaxis = nothing, yaxis = nothing) #Graph the arbitrary point as a red circle
#plot!(xlimit = [APoint[1]-13, APoint[1]+13], ylimit = [APoint[2]-13, APoint[2]+13])
=#

function adjMat(X,Y)
    pairs = [(X[i],Y[i]) for i in 1:length(X)]
    unPairs = unique(pairs)
    n = length(unPairs)
    vals = 1:n
    d = Dict(unPairs[i]=>i for i in 1:n)
    mat = zeros(Int8,n,n)
    for i in 1:4:length(pairs)
        for j in 1:4
            mat[d[pairs[i+(j-1)]],d[pairs[i+(j%4)]]] = 1
            # making symetric conection
            mat[d[pairs[i+(j%4)]],d[pairs[i+(j-1)]]] = 1
        end
    end
    # positions post processing
    # converting to a matrix
    pos = transpose(hcat([collect(p) for p in unPairs]...))
    pos = [Float64(a) for a in pos]
    # centering positions
    centroid = transpose(sum(pos,dims=1)/size(pos)[1])
    for i in 1:size(pos)[1]
        pos[i,:] -= centroid
    end
    return pos, mat
end

function neighList(mat)
    list = []
    for i in 1:size(mat)[1]
        push!(list,[j for j in 1:size(mat)[2] if mat[i,j]!=0])
    end
    return list
end

X,Y = qpLattice(parsed_args["Nsides"],parsed_args["Nerror"];Iterations=parsed_args["Iterations"],Alpha=parsed_args["Alpha"])
pos,mat = adjMat(X,Y)
neig = neighList(mat);
#=
p = scatter(pos[:,1],pos[:,2])
for i in 1:size(mat)[1]
    for j in neig[i]
        plot!([pos[i,1],pos[j,1]],[pos[i,2],pos[j,2]],c=:black,lw=0.4,alpha=0.4,legend=false,aspect_ratio=1)
    end
end
p
=#
writedlm(string("../networks/","mat_",parsed_args["save"],".csv"),mat,',')
writedlm(string("../networks/","pos_",parsed_args["save"],".csv"),pos,',')


