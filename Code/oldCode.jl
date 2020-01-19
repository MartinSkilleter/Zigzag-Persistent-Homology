using Eirene, CSV, LinearAlgebra, SparseArrays, GraphPlot, LightGraphs

#= Yossi Bokor was instrumental in developing this code. He gave
an explicit example, which I then abstracted =#
cd("C:\\Users\\marty\\Downloads\\Zigzag Persistent Homology\\example_graphs")

# compute adjacency matrix for graph from .txt file
adjacency = convert(Matrix, CSV.read("graph_1.txt", header=false))

# number of 0-cells = vertices
n = size(adjacency, 1)

#= list of cells: each 0-cell is entered as a list because 1-cells
will be lists consisting of their endpoints =#
cells = [[i] for i = 1:n]
for i = 1:n
    adjacency[i, i] = 0 #= there are no edges from a vertex to itself;
    this shouldn't affect if the vertex is a bridge =#
    for j = i:n
        if adjacency[i, j] == 1
            push!(cells, [i, j]) # add 1-cells as lists of endpoints
        end
    end
end

nc = length(cells) # total number of cells in complex
f_v = zeros(nc, 1) # birth time of simplices

#= Define the distance function f : X → ℝ as follows:
fix a vertex v ∈ X. For any vertex w ∈ X, f(w) is the length
of the shortest path from w to v (computed using Dijkstra's Algorithm
from the Graphs package). For an edge (a,b), f((a,b)) = max{f(a), f(b)}.
We are using max so that we can consider sublevelsets i.e. go from -∞
to some real value. =#

v = 101
g = LightGraphs.Graph(adjacency)
ds = dijkstra_shortest_paths(g, v)
for i = 1:n
    f_v[i] = ds.dists[i]
end

# boundary matrix (block matrix where only top-right corner has nonzero entries)
boundary = zeros(nc, nc)

# set end vertices of a 1-cell to be boundaries
for i = n+1 : nc
    boundary[cells[i][1], i] = 1
    boundary[cells[i][2], i] = 1

    #= since we're in the for loop already, set value of f on an edge
    to be the max of values on the endpoints =#
    f_v[i] = max(f_v[cells[i][1]], f_v[cells[i][2]])
end

# black box of stuff Eirene needs
B = sparse(convert(Matrix, boundary))
r_v = B.rowval
c_p = B.colptr
d_v = [0 for i = 1:nc] # dimension of simplices
for i = (n+1) : nc
    d_v[i] = 1
end

C = eirene(rv = r_v, cp = c_p, dv = d_v, fv = f_v, maxdim = 2)
println(barcode(C, dim=1))

#= plot the graph corresponding to the given adjacency matrix =#
function plotGraph()
    G = LightGraphs.Graph(adjacency)
    nodeLabel = [1 : n;]
    gplot(G, nodelabel = nodeLabel)
end

plotGraph()

#= local neighbourhood test to see if a vertex v is a bridge :
let d be the smallest distance such that the subgraph H of G
consisting of all vertices / edges distance d or less from v
contains a 1-cycle. If removing v from H increases the number
of path-components then v was a bridge =#

function nbhdTest()
    oneBar = barcode(C, dim=1)
    d = oneBar[size(oneBar, 1)]
    newCells = []
    index = 0 # index of v in list of newCells
    for i = 1 : n
        if ds.dists[i] <= d
            push!(newCells, [i])
            for J in newCells
                j = J[1]
                if adjacency[i, j] == 1
                    push!(newCells, [j, i])
                end
            end
        end
        println(newCells)
        if i == v
            index = length(newCells)
        end
    end
end

nbhdTest()
