# CartesianIndex are array-like objects that allow us to 
# index a multidimensional array with a single object.
# We use the CI abreviation to ease writting
const CI = CartesianIndex


"""
inBounds(shape,index)::Bool

Takes a tuple `shape` and a CartesianIndex `index`.

Returns true if the cartesian index is in bound for an array with shape `shape`
"""
function inBounds(shape,index)::Bool
    dim = length(shape)
    flag = true
    i = 1
    while flag && i <= dim
        flag = flag && (1 <= index[i] <= shape[i] )
        i = i+1
    end
    return flag
end


"""
periodicBoundary(shape,index)

Returns a CartesianIndex with same indexes than `index` 
but after begin applied the `mod1` function in order to keep it
inside the boundaries of an array with shape `shape`
"""
function periodicBoundary(shape,index)
    d = length(shape)
    newIndex = [mod1(index[i],shape[i]) for i in 1:d]
    return CI(newIndex...)
end

# all the functions that continue take a CartesianIndex as argument
# and return its neighbors or its cartesian coordinates 
# in different lattice configurations


# linear chanin
function linearNeighbors(pos)
    n1 = CI(pos[1]-1)
    n2 = CI(pos[1]+1)
    return [n1,n2]
end

function linearCoordinates(pos)
    return [pos[1]]
end

# square 2D lattice
function squareNeighbors(pos)
    n1 = CI(pos[1]+1,pos[2])
    n2 = CI(pos[1]-1,pos[2])
    n3 = CI(pos[1],pos[2]+1)
    n4 = CI(pos[1],pos[2]-1)
    return [n1,n2,n3,n4]
end

function squareCoordinates(pos)
    return [pos[2],pos[1]]
end

# cubic 3D lattice
function cubicNeighbors(pos)
    n1 = CI(pos[1]+1,pos[2],pos[3])
    n2 = CI(pos[1]-1,pos[2],pos[3])
    n3 = CI(pos[1],pos[2]+1,pos[3])
    n4 = CI(pos[1],pos[2]-1,pos[3])
    n5 = CI(pos[1],pos[2],pos[3]+1)
    n6 = CI(pos[1],pos[2],pos[3]-1)
    return [n1,n2,n3,n4,n5,n6]
end

function cubicCoordinates(pos)
    return [pos[2],pos[1],pos[3]]
end

# triangular 2D lattice
function triangularNeighbors(pos)
    # normal square neighbors
    n1 = CI(pos[1]+1,pos[2])
    n2 = CI(pos[1]-1,pos[2])
    n3 = CI(pos[1],pos[2]+1)
    n4 = CI(pos[1],pos[2]-1)
    # diagonal neighbors
    n5 = CI(pos[1]+1,pos[2]-1)
    n6 = CI(pos[1]-1,pos[2]+1)
    return [n1,n2,n3,n4,n5,n6]
end

function triangularCoordinates(pos)
    unitVectors = [[1/2,sqrt(3)/2],[1,0]]
    #posLatt[pos] = (pos[1])*unitVectors[1] + pos[2]*unitVectors[2]            
    return  pos[1]*unitVectors[1] + pos[2]*unitVectors[2]
end


# triangular 2D lattice but centered so that it doesnt look
# like a slided square
function centeredTriangularNeighbors(pos)
    # normal square neighbors
    n1 = CI(pos[1]+1,pos[2])
    n2 = CI(pos[1]-1,pos[2])
    n3 = CI(pos[1],pos[2]+1)
    n4 = CI(pos[1],pos[2]-1)
    # diagonal neighbors
    # even line
    if mod(pos[1],2) == 0
        n5 = CI(pos[1]+1,pos[2]+1)
        n6 = CI(pos[1]-1,pos[2]+1)
    # odd line
    else
        n5 = CI(pos[1]+1,pos[2]-1)
        n6 = CI(pos[1]-1,pos[2]-1)
    end
    return [n1,n2,n3,n4,n5,n6]
end

function centeredTriangularCoordinates(pos)
    unitVectors = [[1/2,sqrt(3)/2],[1,0]]
    #posLatt[pos] = (pos[1])*unitVectors[1] + pos[2]*unitVectors[2]            
    return pos[1]*unitVectors[1] + (pos[2]-div(pos[1]-1,2))*unitVectors[2]
end


# layered triangular 3D lattice 
function hexagonalNeighbors(pos)
    # normal square neighbors
    n1 = CI(pos[1]+1,pos[2],pos[3])
    n2 = CI(pos[1]-1,pos[2],pos[3])
    n3 = CI(pos[1],pos[2]+1,pos[3])
    n4 = CI(pos[1],pos[2]-1,pos[3])
    # diagonal neighbors
    n5 = CI(pos[1]+1,pos[2]-1,pos[3])
    n6 = CI(pos[1]-1,pos[2]+1,pos[3])
    # upper neighbors
    n7 = CI(pos[1],pos[2],pos[3]+1)
    n8 = CI(pos[1],pos[2],pos[3]-1)
    return [n1,n2,n3,n4,n5,n6,n7,n8]
end

function hexagonalCoordinates(pos)
    unitVectors = [[1/2,sqrt(3)/2,0],[1,0,0],[0,0,1]]
    #posLatt[pos] = (pos[1])*unitVectors[1] + pos[2]*unitVectors[2]            
    return  pos[1]*unitVectors[1] + pos[2]*unitVectors[2] + pos[3]*unitVectors[3]
end

# layered centeredTriangular 3D lattice 
function centeredHexagonalNeighbors(pos)
    n1 = CI(pos[1]+1,pos[2],pos[3])
    n2 = CI(pos[1]-1,pos[2],pos[3])
    n3 = CI(pos[1],pos[2]+1,pos[3])
    n4 = CI(pos[1],pos[2]-1,pos[3])
    # diagonal neighbors
    # even line
    if mod(pos[1],2) == 0
        n5 = CI(pos[1]+1,pos[2]+1,pos[3])
        n6 = CI(pos[1]-1,pos[2]+1,pos[3])
    # odd line
    else
        n5 = CI(pos[1]+1,pos[2]-1,pos[3])
        n6 = CI(pos[1]-1,pos[2]-1,pos[3])
    end
    # upper neighbors
    n7 = CI(pos[1],pos[2],pos[3]+1)
    n8 = CI(pos[1],pos[2],pos[3]-1)
    return [n1,n2,n3,n4,n5,n6,n7,n8]
end

function centeredHexagonalCoordinates(pos)
    unitVectors = [[1/2,sqrt(3)/2,0],[1,0,0],[0,0,1]]
    #posLatt[pos] = (pos[1])*unitVectors[1] + pos[2]*unitVectors[2]            
    return  pos[1]*unitVectors[1] + (pos[2]-div(pos[1]-1,2))*unitVectors[2] + pos[3]*unitVectors[3]
end

# fcc 3D lattice
function fccNeighbors(pos)
    # first the 6 triangular neighbors on the 2D layer
    # normal square neighbors
    n1 = CI(pos[1]+1,pos[2],pos[3])
    n2 = CI(pos[1]-1,pos[2],pos[3])
    n3 = CI(pos[1],pos[2]+1,pos[3])
    n4 = CI(pos[1],pos[2]-1,pos[3])
    # diagonal neighbors
    n5 = CI(pos[1]+1,pos[2]-1,pos[3])
    n6 = CI(pos[1]-1,pos[2]+1,pos[3])
    # 3 lower neighbors
    n7 = CI(pos[1],pos[2],pos[3]-1)
    n8 = CI(pos[1]+1,pos[2],pos[3]-1)
    n9 = CI(pos[1],pos[2]+1,pos[3]-1)
    # 3 upper neighbors
    n10 = CI(pos[1],pos[2],pos[3]+1)
    n11 = CI(pos[1]-1,pos[2],pos[3]+1)
    n12 = CI(pos[1],pos[2]-1,pos[3]+1)
    return [n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12]
end

function fccCoordinates(pos)
    r3 = sqrt(3)
    r2 = sqrt(2)
    unitVectors = [[1/2,r3/2,0],[1,0,0],[1/2,1/(2*r3),r2/r3]]
    #posLatt[pos] = (pos[1])*unitVectors[1] + pos[2]*unitVectors[2]            
    return  pos[1]*unitVectors[1] + pos[2]*unitVectors[2] + pos[3]*unitVectors[3]
end

# fcc 3D lattice but centered
function centeredFccNeighbors(pos)
    # first the 6 triangular neighbors on the 2D layer
    n1 = CI(pos[1]+1,pos[2],pos[3])
    n2 = CI(pos[1]-1,pos[2],pos[3])
    n3 = CI(pos[1],pos[2]+1,pos[3])
    n4 = CI(pos[1],pos[2]-1,pos[3])
    # even line
    if mod(pos[1],2) == 0
        n5 = CI(pos[1]+1,pos[2]+1,pos[3])
        n6 = CI(pos[1]-1,pos[2]+1,pos[3])
        # second layer
        if mod(pos[3],3) == 2
            #lower neighbors
            n7 = CI(pos[1],pos[2],pos[3]-1)
            n8 = CI(pos[1],pos[2]+1,pos[3]-1)
            n9 = CI(pos[1]+1,pos[2],pos[3]-1)
            # 3 upper neighbors
            n10 = CI(pos[1],pos[2],pos[3]+1)
            n11 = CI(pos[1],pos[2]+1,pos[3]+1)
            n12 = CI(pos[1]+1,pos[2]+1,pos[3]+1)
        # fisrt layer
        elseif mod(pos[3],3) == 1
            #lower neighbors
            n7 = CI(pos[1],pos[2],pos[3]-1)
            n8 = CI(pos[1]-1,pos[2],pos[3]-1)
            n9 = CI(pos[1]-1,pos[2]+1,pos[3]-1)
            # 3 upper neighbors
            n10 = CI(pos[1],pos[2],pos[3]+1)
            n11 = CI(pos[1],pos[2]-1,pos[3]+1)
            n12 = CI(pos[1]-1,pos[2],pos[3]+1)
        #third layer
        else
            #lower neighbors
            n7 = CI(pos[1],pos[2],pos[3]-1)
            n8 = CI(pos[1],pos[2]-1,pos[3]-1)
            n9 = CI(pos[1]+1,pos[2],pos[3]-1)
            # 3 upper neighbors
            n10 = CI(pos[1],pos[2],pos[3]+1)
            n11 = CI(pos[1]+1,pos[2],pos[3]+1)
            n12 = CI(pos[1]+1,pos[2]+1,pos[3]+1)
        end
    # odd line
    else
        # diagonal neighbors
        n5 = CI(pos[1]+1,pos[2]-1,pos[3])
        n6 = CI(pos[1]-1,pos[2]-1,pos[3])
        # second layer
        if mod(pos[3],3) == 2
            #lower neighbors
            n7 = CI(pos[1],pos[2],pos[3]-1)
            n8 = CI(pos[1],pos[2]+1,pos[3]-1)
            n9 = CI(pos[1]+1,pos[2],pos[3]-1)
            # 3 upper neighbors
            n10 = CI(pos[1],pos[2],pos[3]+1)
            n11 = CI(pos[1],pos[2]+1,pos[3]+1)
            n12 = CI(pos[1]-1,pos[2],pos[3]+1)
        # fisrt layer
        elseif mod(pos[3],3) == 1
            #lower neighbors
            n7 = CI(pos[1],pos[2],pos[3]-1)
            n8 = CI(pos[1]-1,pos[2],pos[3]-1)
            n9 = CI(pos[1]-1,pos[2]-1,pos[3]-1)
            # 3 upper neighbors
            n10 = CI(pos[1],pos[2],pos[3]+1)
            n11 = CI(pos[1],pos[2]-1,pos[3]+1)
            n12 = CI(pos[1]-1,pos[2]-1,pos[3]+1)
        #third layer
        else
            #lower neighbors
            n7 = CI(pos[1],pos[2],pos[3]-1)
            n8 = CI(pos[1],pos[2]-1,pos[3]-1)
            n9 = CI(pos[1]+1,pos[2]-1,pos[3]-1)
            # 3 upper neighbors
            n10 = CI(pos[1],pos[2],pos[3]+1)
            n11 = CI(pos[1]+1,pos[2],pos[3]+1)
            n12 = CI(pos[1]+1,pos[2]-1,pos[3]+1)
        end
    end
    return [n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12]
end

function centeredFccCoordinates(pos)
    r3 = sqrt(3)
    r2 = sqrt(2)
    unitVectors = [[1/2,r3/2,0],[1,0,0],[1/2,1/(2*r3),r2/r3]]
    r = Int(mod(pos[3],3)==0)
    s = div(pos[3]-1,3)
    t = div(pos[1]-1,2)
    #posLatt[pos] = (pos[1])*unitVectors[1] + pos[2]*unitVectors[2]            
    return   (pos[1]-s)*unitVectors[1] + (pos[2]-t-r-s)*unitVectors[2] + pos[3]*unitVectors[3]
end

# hpc 3D lattice
function hcpNeighbors(pos)
    # first the 6 triangular neighbors on the 2D layer
    # normal square neighbors
    n1 = CI(pos[1]+1,pos[2],pos[3])
    n2 = CI(pos[1]-1,pos[2],pos[3])
    n3 = CI(pos[1],pos[2]+1,pos[3])
    n4 = CI(pos[1],pos[2]-1,pos[3])
    # diagonal neighbors
    n5 = CI(pos[1]+1,pos[2]-1,pos[3])
    n6 = CI(pos[1]-1,pos[2]+1,pos[3])
    # 3 lower neighbors
    n7 = CI(pos[1],pos[2],pos[3]-1)
    n8 = CI(pos[1]+1,pos[2],pos[3]-1)
    n9 = CI(pos[1],pos[2]+1,pos[3]-1)
    # 3 upper neighbors
    n10 = CI(pos[1]-1,pos[2]-1,pos[3]+1)
    n11 = CI(pos[1]-1,pos[2],pos[3]+1)
    n12 = CI(pos[1],pos[2]-1,pos[3]+1)
    return [n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12]
end

function hcpCoordinates(pos)
    r3 = sqrt(3)
    r2 = sqrt(2)
    unitVectors = [[1/2,r3/2,0],[1,0,0],[1/2,1/(2*r3),r2/r3]]
    #posLatt[pos] = (pos[1])*unitVectors[1] + pos[2]*unitVectors[2]            
    return  pos[1]*unitVectors[1] + pos[2]*unitVectors[2] + (pos[3])*unitVectors[3] + (div(pos[3]-1,2))*[1/2,1/(2*sqrt(3)),0]
end

# hpc 3D lattice but centered
function centeredHcpNeighbors(pos)
    # normal centered triangular neighbors
    n1 = CI(pos[1]+1,pos[2],pos[3])
    n2 = CI(pos[1]-1,pos[2],pos[3])
    n3 = CI(pos[1],pos[2]+1,pos[3])
    n4 = CI(pos[1],pos[2]-1,pos[3])
    # second layer
    if mod(pos[3],2) == 0
        # even line
        if mod(pos[1],2) == 0
            # diagonal neighbors
            n5 = CI(pos[1]+1,pos[2]+1,pos[3])
            n6 = CI(pos[1]-1,pos[2]+1,pos[3])
            # lower neighbors
            n7 = CI(pos[1],pos[2],pos[3]-1)
            n8 = CI(pos[1],pos[2]+1,pos[3]-1)
            n9 = CI(pos[1]+1,pos[2]+1,pos[3]-1)
            # upper neighbors
            n10 = CI(pos[1],pos[2],pos[3]+1)
            n11 = CI(pos[1],pos[2]+1,pos[3]+1)
            n12 = CI(pos[1]+1,pos[2]+1,pos[3]+1)
        # even line
        else
            # diagonal neighbors
            n5 = CI(pos[1]+1,pos[2]+1,pos[3])
            n6 = CI(pos[1]-1,pos[2]+1,pos[3])
            # lower neighbors
            n7 = CI(pos[1],pos[2],pos[3]-1)
            n8 = CI(pos[1],pos[2]+1,pos[3]-1)
            n9 = CI(pos[1]+1,pos[2],pos[3]-1)
            # upper neighbors
            n10 = CI(pos[1],pos[2],pos[3]+1)
            n11 = CI(pos[1],pos[2]+1,pos[3]+1)
            n12 = CI(pos[1]+1,pos[2],pos[3]+1)
        end
    # first layer
    else
        # even line
        if mod(pos[1],2) == 0
            n5 = CI(pos[1]+1,pos[2]+1,pos[3])
            n6 = CI(pos[1]-1,pos[2]+1,pos[3])
            # lower neighbors
            n7 = CI(pos[1],pos[2],pos[3]-1)
            n8 = CI(pos[1],pos[2]-1,pos[3]-1)
            n9 = CI(pos[1]-1,pos[2],pos[3]-1)
            # upper neighbors
            n10 = CI(pos[1],pos[2],pos[3]+1)
            n11 = CI(pos[1],pos[2]-1,pos[3]+1)
            n12 = CI(pos[1]-1,pos[2],pos[3]+1)
        # odd line
        else
            n5 = CI(pos[1]+1,pos[2]-1,pos[3])
            n6 = CI(pos[1]-1,pos[2]-1,pos[3])
            # lower neighbors
            n7 = CI(pos[1],pos[2],pos[3]-1)
            n8 = CI(pos[1],pos[2]-1,pos[3]-1)
            n9 = CI(pos[1]-1,pos[2]-1,pos[3]-1)
            # upper neighbors
            n10 = CI(pos[1],pos[2],pos[3]+1)
            n11 = CI(pos[1],pos[2]-1,pos[3]+1)
            n12 = CI(pos[1]-1,pos[2]-1,pos[3]+1)
        end
    end
    return [n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12]
end

function centeredHcpCoordinates(pos)
    r3 = sqrt(3)
    r2 = sqrt(2)
    unitVectors = [[1/2,r3/2,0],[1,0,0],[1/2,1/(2*r3),r2/r3]]
    s = div(pos[3]-1,2)
    t = div(pos[1]-1,2)
    #posLatt[pos] = (pos[1])*unitVectors[1] + pos[2]*unitVectors[2]            
    return   pos[1]*unitVectors[1] + (pos[2]-t)*unitVectors[2] + pos[3]*unitVectors[3] - s*[1,1/r3,0]
end


# 3D bcc lattice
function bccNeighbors(pos)
    # lower neighbors
    n1 = CI(pos[1],pos[2],pos[3]-1)
    n2 = CI(pos[1]+1,pos[2],pos[3]-1)
    n3 = CI(pos[1],pos[2]+1,pos[3]-1)
    n4 = CI(pos[1]+1,pos[2]+1,pos[3]-1)
    # upper neighbors
    n5 = CI(pos[1],pos[2],pos[3]+1)
    n6 = CI(pos[1]-1,pos[2],pos[3]+1)
    n7 = CI(pos[1],pos[2]-1,pos[3]+1)
    n8 = CI(pos[1]-1,pos[2]-1,pos[3]+1)
    return [n1,n2,n3,n4,n5,n6,n7,n8]
end

function bccCoordinates(pos)
    unitVectors = [[1,0,0],[0,1,0],[1/2,1/2,1/2]]
    return pos[1]*unitVectors[1] + pos[2]*unitVectors[2] + pos[3]*unitVectors[3]
end

# 3D bcc lattice
function centeredBccNeighbors(pos)
    # second layer
    if mod(pos[3],2)==0
        # lower neighbors
        n1 = CI(pos[1],pos[2],pos[3]-1)
        n2 = CI(pos[1]+1,pos[2],pos[3]-1)
        n3 = CI(pos[1],pos[2]+1,pos[3]-1)
        n4 = CI(pos[1]+1,pos[2]+1,pos[3]-1)
        # upper neighbors
        n5 = CI(pos[1],pos[2],pos[3]+1)
        n6 = CI(pos[1]+1,pos[2],pos[3]+1)
        n7 = CI(pos[1],pos[2]+1,pos[3]+1)
        n8 = CI(pos[1]+1,pos[2]+1,pos[3]+1)
    else
    # first layer
        # lower neighbors
        n1 = CI(pos[1],pos[2],pos[3]-1)
        n2 = CI(pos[1]-1,pos[2],pos[3]-1)
        n3 = CI(pos[1],pos[2]-1,pos[3]-1)
        n4 = CI(pos[1]-1,pos[2]-1,pos[3]-1)
        # upper neighbors
        n5 = CI(pos[1],pos[2],pos[3]+1)
        n6 = CI(pos[1]-1,pos[2],pos[3]+1)
        n7 = CI(pos[1],pos[2]-1,pos[3]+1)
        n8 = CI(pos[1]-1,pos[2]-1,pos[3]+1)
    end
    return [n1,n2,n3,n4,n5,n6,n7,n8]
end

function centeredBccCoordinates(pos)
    unitVectors = [[1,0,0],[0,1,0],[1/2,1/2,1/2]]
    return (pos[1]-div(pos[3]-1,2))*unitVectors[1] + (pos[2]- div(pos[3]-1,2))*unitVectors[2] + pos[3]*unitVectors[3]
end


# honeycomb 2D lattice
function honeycombNeighbors(pos)
    n1 = CI(pos[1]-1, pos[2])
    n2 = CI(pos[1]+1, pos[2])
    # even line
    if mod(pos[1],2) == 0
        if mod(pos[2],2) == 0
            n3 = CI(pos[1],pos[2]-1)
        else
            n3 = CI(pos[1],pos[2]+1)
        end
    # odd line
    else
        if mod(pos[2],2) == 0
            n3 = CI(pos[1],pos[2]+1)
        else
            n3 = CI(pos[1],pos[2]-1)
        end
    end
    return [n1,n2,n3]
end

function honeycombCoordinates(pos)
    a = div(pos[2]+1,2)
    b = mod(pos[2],2)
    if mod(pos[1],2)==1    
        return centeredTriangularCoordinates(CI(pos[1],3*a-2*b))
    else
        return centeredTriangularCoordinates(CI(pos[1],3*a-(b+1)))
    end
end

function layeredHoneycombNeighbors(pos)
    n1 = CI(pos[1]-1, pos[2],pos[3])
    n2 = CI(pos[1]+1, pos[2],pos[3])
    # even line
    if mod(pos[1],2) == 0
        if mod(pos[2],2) == 0
            n3 = CI(pos[1],pos[2]-1,pos[3])
        else
            n3 = CI(pos[1],pos[2]+1,pos[3])
        end
    # odd line
    else
        if mod(pos[2],2) == 0
            n3 = CI(pos[1],pos[2]+1,pos[3])
        else
            n3 = CI(pos[1],pos[2]-1,pos[3])
        end
    end
    n4 = CI(pos[1],pos[2],pos[3]+1)
    n5 = CI(pos[1],pos[2],pos[3]-1)
    return [n1,n2,n3,n4,n5]
end

function layeredHoneycombCoordinates(pos)
    if mod(pos[1],2)==1
        a = div(pos[2]+1,2)
        b = mod(pos[2],2)
        return centeredHexagonalCoordinates(CI(pos[1],3*a-2*b,pos[3]))
    else
        a = div(pos[2]+1,2)
        b = mod(pos[2],2)
        return centeredHexagonalCoordinates(CI(pos[1],3*a-(b+1),pos[3]))
    end
end

"""
makeLattice(name::AbstractString,shape;periodic=false)

Using a name-compatible string `name` and a tuple `shape`, 
it creates a multidimensional "adyacency array", the normal adyacency list
and the coordinates of a lattice. 

`periodic` indicates if the boundary conditions are periodic or not. 
valid values for `name`, for example, are `layeredHoneycomb` or `fcc`
"""
function makeLattice(name::AbstractString,shape;periodic=false)
    # search for the neighbors and coordinates functions
    neighFunc = getfield(Main,Symbol(string(name,"Neighbors")))
    coordFunc = getfield(Main,Symbol(string(name,"Coordinates")))
    d = length(shape)
    # linear or natural indices for the points of the lattice
    linIdx = LinearIndices(shape)
    # initialize arrays for values
    neighLatt = Array{Array{CI{d},1},d}(undef,shape)
    linearNeighLatt = Array{Array{Int32,1},1}(undef,prod(shape))
    coordinates = Array{Float64,2}(undef,(prod(shape),d))
    # iterate over cartesian indices in shape range
    for pos in CartesianIndices(shape)
        # create neighbors
        neighs = neighFunc(pos)
        # modify due to boundary conditions
        if periodic
            neighs = [periodicBoundary(shape,neig) for neig in neighs]
        else
            neighs = [neig for neig in neighs if inBounds(shape,neig)]
        end
        neighLatt[pos] = neighs
        # get linear index of position
        l = linIdx[pos]
        # get linear indexes of neighbors
        linearNeighLatt[l] = [linIdx[neig] for neig in neighs]
        coordinates[l,:] = coordFunc(pos)
    end
    return neighLatt, linearNeighLatt, coordinates
end


"""
plotLattice!(linearNeighLatt,coordinates;la=0.3,ma=0.5,lw=1,ms=10)

Plots a lattice given its adyacency list `linearNeighLatt` of length N 
and a 3xN array of coordinates `coordinates`. 
"""
function plotLattice!(linearNeighLatt,coordinates;la=0.3,ma=0.5,lw=1,ms=10)
    d = size(coordinates)[2]
    if d == 3
        scatter!(coordinates[:,1],coordinates[:,2],coordinates[:,3],markeriize=ms,alpha=ma)
        edges = [Set([u,v]) for u in 1:length(linearNeighLatt) for v in linearNeighLatt[u]]
        edges = unique(edges)
        for e in edges
            vs = [i for i in e]
            plot!([coordinates[vs[1],1],coordinates[vs[2],1]],[coordinates[vs[1],2],coordinates[vs[2],2]],[coordinates[vs[1],3],coordinates[vs[2],3]],color="black",alpha=la,lw=lw)
        end
    elseif d == 2
        scatter!(coordinates[:,1],coordinates[:,2],markeriize=ms,alpha=ma,label=false)
        edges = [Set([u,v]) for u in 1:length(linearNeighLatt) for v in linearNeighLatt[u]]
        edges = unique(edges)
        for e in edges
            vs = [i for i in e]
            plot!([coordinates[vs[1],1],coordinates[vs[2],1]],[coordinates[vs[1],2],coordinates[vs[2],2]],color="black",alpha=la,lw=lw)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    makeLattice("square",(10,10),periodic=true)
end
