"""
    EdgesInCycles(edgList::Array{Array{T,1},1})
    
    Function to return number of all edges in cycles of a graph represented by edgelist. 
    
    Uses [DFS-based algorithm](https://www.geeksforgeeks.org/bridge-in-a-graph/)
"""
function EdgesInCycles(edgList::Array{Array{T,1},1}) where {T<:Integer}
    global count = 0
    n = length(edgList)
    visited = falses(n)
    disc = [typemax(Int16) for i in 1:n]
    low = [typemax(Int16) for i in 1:n]
    parent = zeros(Int16,n)
    global bridges = 0
    global bridgeArray = Array{Array{Int16,1},1}()
    for i in 1:n
        if !visited[i]
            DFS(i,visited,parent,low,disc,edgList)
        end
    end
    normal = Int64(sum([length(y) for y in edgList])/2)
    return (normal-bridges)
end


"""
    BridgesArray(edgList)

    function to return an array of all the bridges of a graph with in edge list representation. 
    
    Uses [DFS-based algorithm](https://www.geeksforgeeks.org/bridge-in-a-graph/)
"""
function BridgesArray(edgList)
    global count = 0
    n = length(edgList) 
    visited = falses(n)
    disc = [typemax(Int16) for i in 1:n]
    low = [typemax(Int16) for i in 1:n]
    parent = 0 * ones(Int16,n)
    global bridges = 0
    global bridgeArray = Array{Array{Int16,1},1}()
    for i in 1:n
        if !visited[i]
            DFS(i,visited,parent,low,disc,edgList)
        end
    end
    return bridgeArray
end


"""
DFS(i,visited,parent,low,disc,edgList)

Constructs and transverses the DFS tree of a graph represented by `edgList`
searching for bridge edges. It stores only keeps the count of bridge edges, 
but not the edges.
"""
function DFS(i,visited,parent,low,disc,edgList)
    # globals must exist before the function call
    global count
    global bridges
    visited[i] = true
    disc[i] = count
    low[i] = count
    count = count+1
    for j in edgList[i]
        if !visited[j]
            parent[j] = i
            DFS(j,visited,parent,low,disc,edgList)
            low[i] = minimum([low[j],low[i]])
            if low[j] > disc[i]
                #push!(bridgeArray,[i,j])
                bridges = bridges+1
            end
        elseif j != parent[i]
            low[i] = minimum([low[i],disc[j]])
        end
    end
end

"""
DFS2(i,visited,parent,low,disc,edgList)

Constructs and transverses the DFS tree of a graph represented by `edgList`
searching for bridge edges. It stores both the number of bridge edges and 
the edges themselves.
"""
function DFS2(i,visited,parent,low,disc,edgList)
    # globals must exist before the function call
    global count
    global bridges
    global bridgeArray
    visited[i] = true
    disc[i] = count
    low[i] = count
    count = count+1
    for j in edgList[i]
        if !visited[j]
            parent[j] = i
            DFS2(j,visited,parent,low,disc,edgList)
            low[i] = minimum([low[j],low[i]])
            if low[j] > disc[i]
                push!(bridgeArray,[i,j])
                bridges = bridges+1
            end
        elseif j != parent[i]
            low[i] = minimum([low[i],disc[j]])
        end
    end
end




if abspath(PROGRAM_FILE) == @__FILE__
    #   Testing
    include("systems.jl")
end