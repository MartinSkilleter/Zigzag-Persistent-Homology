using GraphPlot, LightGraphs

# number of vertices
n = 300

P = [   [0.7 0.2];
        [0.2 0.7]]

adjacency = zeros(n, n)

for i = 1 : n
    for j = i : n
        # identify which communities i and j lie in
        c1 = i <= n / 2 ? 1 : 2
        c2 = j <= n / 2 ? 1 : 2
        p = rand(Float64,1)[1]
        if p <= P[c1, c2]
            adjacency[i, j] = adjacency[j, i] = 1
        end
    end
end

G = Graph(adjacency)

gplot(G)
