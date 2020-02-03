using Eirene, CSV, LinearAlgebra, SparseArrays, GraphPlot, LightGraphs, Plots

#= Yossi Bokor was instrumental in developing this code. He gave
an explicit example, which I then abstracted =#
cd("C:\\Users\\marty\\Downloads\\Zigzag Persistent Homology\\example_graphs")

# compute adjacency matrix for graph from .txt file
adjacency = convert(Matrix, CSV.read("graph_1.txt", header=false))

# number of 0-cells = vertices
n = size(adjacency, 1)

# a simple graph has no edges from a vertex to itself
for v = 1 : n
    adjacency[v, v] = 0
end

G = Graph(adjacency)

e = ne(G)

# construct the Laplacian matrix
L = zeros(n, n)
for i = 1 : n
    L[i, i] = degree(G, i)
end

L -= adjacency

eigenDecomp = eigen(L)
spectralGap = eigenDecomp.values[2]
Fiedler = eigenDecomp.vectors[:,2]

avg = sum(abs, Fiedler) / n

function isClusterPoint(v)
    return(abs(Fiedler[v]) > avg)
end

#= a matrix where the i^th column represents the Djikstra distances from vertex
v_i =#
dists = []

#= a list of dictionary. The i^th entry is the dictionary corresponding to the
vertex v_i, which takes a distance d as a key and returns a list of vertices
of distance d from v_i =#
distdict = []

#= Define the distance function f : X → ℝ as follows:
fix a vertex v ∈ X. For any vertex w ∈ X, f(w) is the length
of the shortest path from w to v (computed using Dijkstra's Algorithm
from the Graphs package). For an edge (a,b), f((a,b)) = max{f(a), f(b)}.
We are using max so that we can consider sublevelsets i.e. go from -∞
to some real value. =#

for v = 1 : n
    ds = dijkstra_shortest_paths(G, v)
    push!(dists, ds)
    dict = Dict()
    for w = 1 : n
        d = ds.dists[w]
        if haskey(dict, d)
            push!(dict[d], w)
        else
            dict[d] = [w]
        end
    end
    push!(distdict, dict)
end

# list of the 1-dimensional persistence barcodes
barcodes = []
# list of the cumulative 1-dim. cycles born for each vertex
cumulativeBarcodes = []

# given the vertex v, generate the 1-dimensional barcode
function generateOneBar(v)
    oneBar = [0 0]
    dict = distdict[v]
    d = maximum(keys(dict))
    cumul = zeros(Int64, d)
    previouslyDone = 0
    for ρ = 0 : d
        (H, _) = induced_subgraph(G, neighborhood(G, v, ρ))
        new = 1 - nv(H) + ne(H)
        for row = previouslyDone + 1 : new
            oneBar = vcat([ρ Inf], oneBar)
        end
        if ρ > 0
            cumul[ρ] = new - previouslyDone
        end
        previouslyDone = new
    end
    push!(cumulativeBarcodes, cumul)
    return(oneBar[1 : size(oneBar, 1)-1, :])
end

for v = 1 : n
    oneBar = generateOneBar(v)
    push!(barcodes, oneBar)
end

clusterpoints = ones(Int64, n) # 1 for bridges, 2 for cluster points

#= given a vertex v known to be a cluster point, look at the neighbourhood of v
consisting of all vertices distance <= d from v. For a given vertex w, if the
Wasserstein distance between their 1-dimensional persistence diagrams is <= wd
then call w also a cluster point =#
function cluster(v, d, wd)
    vertices = neighborhood(G, v, d)
    for w in vertices
        if clusterpoints[w] == 1 && wasserstein_distance(barcodes[v], barcodes[w]) <= wd
            clusterpoints[w] = 2
        end
    end
end

function isClusterPointHelper(v)
    cumul = cumulativeBarcodes[v]
    lastBorn = size(cumul, 1)
    currentMax = cumul[1]
    for ρ = 2 : lastBorn
        if cumul[ρ] >= currentMax
            return((ρ-1, lastBorn))
        end
        currentMax = max(currentMax, cumul[ρ])
    end
    return((lastBorn, lastBorn))
end

function isClusterPoint2(v, jumpLimit)
    (jump, lastBorn) = isClusterPointHelper(v)
    return(jump < jumpLimit * lastBorn)
end

function findJumpLimit()
    jumpLimit = 0.05
    found = false
    while !found
        count = 0
        for v = 1 : n
            if isClusterPoint2(v, jumpLimit)
                count += 1
            end
        end
        if count > 0
            found = true
        else
            jumpLimit += 0.05
        end
    end
    return(jumpLimit)
end

#= plot the graph corresponding to the given adjacency matrix =#
function plotGraph()
    nodecolor = [GraphPlot.colorant"red", GraphPlot.colorant"blue"] # red for bridges, blue for cluster points
    nodelabel = [1 : n;]

    jumpLimit = findJumpLimit()
    println(jumpLimit)

    # identify the cluster points
    for v = 1 : n
        if isClusterPoint2(v, jumpLimit)
            cluster(v, 4, 23)
            clusterpoints[v] = 2
        end
    end

    nodefillc = nodecolor[clusterpoints]

    gplot(G, nodelabel = nodelabel, nodefillc = nodefillc)
end

plotGraph()

#= p = histogram(randn(1000))
xlabel!("Time, t")
ylabel!("Number of 1-cycles Born at Time t")
display(p) =#
