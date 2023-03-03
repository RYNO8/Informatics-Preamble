#ifndef GRAPH_H
#define GRAPH_H
#include <type_traits>
#include <vector>
#include <unordered_set>
#include <set>
#include <iostream>
#include <cassert>
#include <stack>
#include <queue>
#include <map>
#include <limits>
#include <algorithm>
#include <numeric>

namespace DS {
    template<
        // weight data held for each edge
        // TODO incorporate with isWeighted
        typename T,
        // Whether edges are weighted
        // Unweighted edges are assigned unit weights
        bool isWeighted = false,
        // Whether edges are directed or bidirectional
        bool isDirected = false
    >
    class Graph {
        using node_index = size_t;
        using edge_index = size_t;

        /************************************************
         *                 INITIALISATION               *
         ************************************************/

    private:
        // Number of nodes
        size_t V;
        size_t wantedM = 0;

        // `inEdges[node]` is the multiset of edges leading into `node`
        std::vector<std::multiset<std::pair<node_index, T>>> inEdges;
        // `outEdges[node]` is the multiset of edges leading out of `node`
        std::vector<std::multiset<std::pair<node_index, T>>> outEdges;
        // A multiset of all edges
        std::multiset<std::pair<std::pair<node_index, node_index>, T>> edges;

        inline bool validNode(node_index node) const {
            return 1 <= node && node <= V;
        }

    public:
        // O(V)
        // Initialises an empty graph, assuming 1 indexed
        Graph(size_t _V): V(_V) {
            inEdges.resize(V + 1);
            outEdges.resize(V + 1);
        }

        // O(V + E log E)
        // Initalises an unweighted graph from vector of edges, assuming 1 indexed
        Graph(size_t _V, std::vector<std::pair<node_index, node_index>>& _edges): V(_V) {
            inEdges.resize(V + 1);
            outEdges.resize(V + 1);

            for (std::pair<node_index, node_index> edge: _edges) addEdge(edge.first, edge.second);
        }

        // O(V + E log E)
        // Initalises a weighted graph from vector of edges, assuming 1 indexed
        Graph(size_t _V, std::vector<std::pair<std::pair<node_index, node_index>, T>>& _edges): V(_V) {
            inEdges.resize(V + 1);
            outEdges.resize(V + 1);

            for (std::pair<std::pair<node_index, node_index>, T> edge: _edges) addEdge(edge.first.first, edge.first.second, edge.second);
        }

        /************************************************
         *                    DISPLAY                   *
         ************************************************/

        void reserve(size_t M) {
            wantedM = M;
        }

        friend std::istream& operator>>(std::istream &in, Graph &g) {
            for (size_t u, v, i = 0; i < g.wantedM; ++i) {
                in >> u >> v;
                if (isWeighted) {
                    T w;
                    in >> w;
                    g.addEdge(u, v, w);
                }
                else g.addEdge(u, v);
            }
            return in;
        }

        // O(V + E)
        // Displays the graph, showing the outwards edge connections of each edge in lexographic order
        // @param `out` The string representation of the graph is piped to this output stream
        friend std::ostream& operator<<(std::ostream &out, const Graph<T, isWeighted, isDirected> &graph) {
            for (size_t u = 1; u <= graph.V; ++u) {
                out << u << ':';
                for (std::pair<node_index, T> edge: graph.outEdges[u]) {
                    // on an undirected graph, don't print edges twice
                    if (!isDirected && u > edge.first) continue;

                    out << ' ' << edge.first;
                    if (isWeighted) out << " (w = " << edge.second << ')';
                }
                out << '\n';
            }
            return out;
        }

        /************************************************
         *                EDGE UTILITIES                *
         ************************************************/

        // O(1)
        // @returns `V`, the number of nodes
        inline size_t size() const {
            return V;
        }

        // O(log E)
        // Add the specified edge to the graph - this may create a double edge
        // @note Edge weight is ignored if the graph is unweighted
        void addEdge(node_index u, node_index v, T w = T(1)) {
            if (!isWeighted) assert(w == T(1) && "Unweighted graphs should only have edges of unit weight");
            assert(validNode(u) && validNode(v) && "Node index out of range");

            outEdges[u].insert({ v, w });
            inEdges[v].insert({ u, w });
            edges.insert({ {u, v}, w });
            if (!isDirected && u != v) {
                outEdges[v].insert({ u, w });
                inEdges[u].insert({ v, w });
                edges.insert({ {v, u}, w });
            }
        }

        // O(1)
        // Adds a new node with no edges
        // @returns The index of this node
        inline const node_index addNode() {
            ++V;
            inEdges.push_back({});
            outEdges.push_back({});
            return V;
        }

        // O(log E + answer)
        // @returns the number of edges from `u` to `v`
        // @note doesn't count number of edges in the opposite direction
        inline node_index numEdges(node_index u, node_index v) const {
            assert(validNode(u) && validNode(v) && "Node index out of range");

            auto lowest = edges.lower_bound({ { u, v }, std::numeric_limits<T>::min() });
            auto highest = edges.upper_bound({ { u, v }, std::numeric_limits<T>::max() });
            return distance(lowest, highest);
        }

        // O(log E)
        // @returns Whether the specified edge (with any edge weight) is contained in the graph
        inline bool containsEdge(node_index u, node_index v) const {
            assert(validNode(u) && validNode(v) && "Node index out of range");

            auto iter = edges.lower_bound({ { u, v }, std::numeric_limits<T>::min() });
            return iter != edges.end() && iter->first.first == u && iter->first.second == v;
        }

        // O(log E)
        // @returns Whether the specified edge is contained in the graph
        inline bool containsEdge(node_index u, node_index v, T w) const {
            assert(validNode(u) && validNode(v) && "Node index out of range");

            assert(isWeighted && "Can't find weighted edges in an unweighted graph");
            return edges.find({ {u, v}, w }) != edges.end();
        }

        // O(log E)
        // If the specified edge (with any edge weight) is present, remove it, otherwise silently do nothing
        inline void removeEdge(node_index u, node_index v) {
            assert(validNode(u) && validNode(v) && "Node index out of range");

            if (!containsEdge(u, v)) return;
            for (node_index rep = 0; rep <= node_index(!isDirected && u != v); ++rep) {
                outEdges[u].erase(
                    outEdges[u].lower_bound({ v, std::numeric_limits<T>::min() }),
                    outEdges[u].upper_bound({ v, std::numeric_limits<T>::max() })
                );
                inEdges[v].erase(
                    inEdges[v].lower_bound({ u, std::numeric_limits<T>::min() }),
                    inEdges[v].upper_bound({ u, std::numeric_limits<T>::max() })
                );
                edges.erase(
                    edges.lower_bound({ { u, v }, std::numeric_limits<T>::min() }),
                    edges.upper_bound({ { u, v }, std::numeric_limits<T>::max() })
                );

                std::swap(u, v);
            }
        }

        // O(log E)
        // If the specified edge is present, remove it, otherwise silently do nothing
        inline void removeEdge(node_index u, node_index v, T w) {
            assert(isWeighted && "Can't remove weighted edges in an unweighted graph");
            assert(validNode(u) && validNode(v) && "Node index out of range");

            if (!containsEdge(u, v, w)) return;
            for (node_index rep = 0; rep <= node_index(!isDirected && u != v); ++rep) {
                //auto iter = outEdges[u].lower_bound({ v, w });
                outEdges[u].erase({ v, w });
                inEdges[v].erase({ u, w });
                edges.erase({ { u, v }, w });
                std::swap(u, v);
            }
        }

        // O(E)
        // Finds the edges that travel into `node`
        // @returns A multiset of (neighbour, weight) pairs
        inline const std::multiset<std::pair<node_index, T>> getEdgesIn(node_index node) const {
            assert(validNode(node) && "Node index out of range");
            return inEdges[node];
        }

        // O(E)
        // Finds the edges that travel out of `node`
        // @returns A multiset of (neighbour, weight) pairs
        template<std::enable_if<isWeighted>>
        inline const std::multiset<std::pair<node_index, T>> getEdgesOut(node_index node) const {
            assert(validNode(node) && "Node index out of range");
            return outEdges[node];
        }
        template<std::enable_if<!isWeighted>>
        inline const std::multiset<node_index> getEdgesOut(node_index node) const {
            assert(validNode(node) && "Node index out of range");
            return outEdges[node];
        }

        // O(E log E)
        // Finds the edges that travel into OR out of `node`
        // @returns A multiset of (neighbour, weight) pairs
        inline const std::multiset<std::pair<node_index, T>> getEdges(node_index node) const {
            assert(validNode(node) && "Node index out of range");
            std::multiset<std::pair<node_index, T>> output = inEdges[node];
            if (isDirected) {
                output.insert(outEdges[node].begin(), outEdges[node].end());
            }
            return output;
        }

        // O(E)
        // @returns The multiset of all edges
        inline const std::multiset<std::pair<std::pair<node_index, node_index>, T>> getEdges() const {
            return edges;
        }

        /************************************************
         *                NODE UTILITIES                *
         ************************************************/

        // O(1)
        // @returns The in degree of `node`
        inline const size_t degreeIn(node_index node) const {
            assert(validNode(node) && "Node index out of range");
            return inEdges[node].size();
        }

        // O(1)
        // @returns The out degree of `node`
        inline const size_t degreeOut(node_index node) const {
            assert(validNode(node) && "Node index out of range");
            return outEdges[node].size();
        }

        // O(1)
        // @returns The total degree of `node` (combination in in degree and out degree)
        inline const size_t degree(node_index node) const {
            assert(validNode(node) && "Node index out of range");
            if (isDirected) {
                assert(inEdges[node].size() == outEdges[node].size() && "Sanity check!");
                return inEdges[node].size();
            }
            else {
                return inEdges[node].size() + outEdges[node].size();
            }
        }

        // O(V)
        // @returns All nodes, in ascending order
        inline const std::vector<node_index> getNodes() const {
            std::vector<node_index> output(V);
            std::iota(output.begin(), output.end(), 1);
            return output;
        }

        // O(V_component log V_component + E_component)
        // @returns All nodes in the same component as `root`, in ascending order
        std::vector<node_index> getComponentNodes(node_index root = 0) const {
            if (root == 0) return getNodes();

            assert(validNode(root) && "Node index out of range");
            std::unordered_set<node_index> seen;
            seen.insert(root);
            std::queue<node_index> q;
            q.push(root);

            while (!q.empty()) {
                node_index u = q.front();
                q.pop();

                for (node_index type = 0; type <= 1; ++type) {
                    for (std::pair<node_index, T> edge: (type ? inEdges: outEdges)[u]) {
                        node_index v = edge.first;
                        if (seen.find(v) == seen.end()) {
                            seen.insert(v);
                            q.push(v);
                        }
                    }
                }
            }

            return std::vector<node_index>(seen.begin(), seen.end());
        }

        // O(V + E)
        // @returns All nodes grouped by their component
        std::vector<std::vector<node_index>> getComponentsNodes() const {
            std::vector<std::vector<node_index>> output;
            std::vector<bool> seen = std::vector<bool>(V + 1, false);
            std::queue<node_index> q;

            for (node_index node = 1; node <= V; ++node) {
                if (!seen[node]) {
                    seen[node] = true;
                    q.push(node);
                    output.push_back({ node });

                    while (!q.empty()) {
                        node_index u = q.front();
                        q.pop();

                        for (std::pair<node_index, T> edge: inEdges[u]) {
                            node_index v = edge.first;
                            if (!seen[v]) {
                                seen[v] = true;
                                q.push(v);
                                output.back().push_back(v);
                            }
                        }
                    }
                }
            }
            return output;
        }

        // O(V_component log V_component + E_component)
        // Copy the component at `root` to a new graph, maintaining the same graph properties and the same node numberings
        // @note `V` does not change, so there may be disconnected nodes
        Graph<T> getSubgraph(node_index root) const {
            assert(validNode(root) && "Node index out of range");
            Graph<T> output(V, isWeighted, isDirected);

            for (node_index u: getComponentNodes(root)) {
                for (std::pair<node_index, T> edge: outEdges[u]) {
                    if (!isDirected && edge.first > u) break;
                    output.addEdge(u, edge.first, edge.second);
                }
            }

            return output;
        }

        /************************************************
         *        ALGORITHMS (COMPLEXITLY CLASS P)      *
         ***********************************************/

        // O(V + E log E)
        // @returns The union of 2 graphs
        // @note May cause edge doubling
        Graph<T> operator+(const Graph<T>& o) const {
            assert(V == o.V && isDirected == o.isDirected && isWeighted == o.isWeighted && "Can't find union of graphs with different properties");

            Graph<T>* output = new Graph(*this); // make a copy
            for (std::pair<std::pair<node_index, node_index>, T> edge: o.edges) output->edges.insert({ {edge.first.second, edge.first.first}, edge.second });
            for (node_index node = 1; node <= V; ++node) {
                output->inEdges[node].insert(output->inEdges[node].end(), o.inEdges[node].begin(), o.inEdges[node].end());
                output->outEdges[node].insert(output->outEdges[node].end(), o.outEdges[node].begin(), o.outEdges[node].end());
            }
            return *output;
        }

        // O(V + E log E)
        Graph<T, isWeighted, isDirected> operator~() const {
            Graph<T, isWeighted, isDirected> output(V);
            for (std::pair<std::pair<node_index, node_index>, T> edge: edges) output.edges.insert({ {edge.first.second, edge.first.first}, edge.second });
            output.inEdges = outEdges;
            output.outEdges = inEdges;
            return output;
        }

        inline void clear() {
            edges.clear();
            for (node_index node = 1; node <= V; ++node) {
                inEdges[node].clear();
                outEdges[node].clear();
            }
        }

    private:
        // O(V + E log V), or O(VE) if negative edge weights
        // Standard Dijkstra algorithm - shortest path from `node` to all other nodes for non-negative edge weights
        // @param `backwards` Indicates whether to consider `inEdges` or `outEdges`
        // @note If the graph is undirected, this becomes equivalent to a bfs (but slightly slower)
        // @note If negative edge weights encountered, switches to `BellmanFord()`
        void Dijkstra(node_index node, std::vector<node_index>& prevNode, std::vector<T>& dist, bool backwards = false) const {
            dist[node] = T(0);
            prevNode[node] = node;
            std::priority_queue<std::pair<T, node_index>, std::vector<std::pair<T, node_index>>, std::greater<std::pair<T, node_index>>> pq;
            pq.push({ 0, node });

            while (!pq.empty()) {
                node_index currNode = pq.top().second;
                T cost = pq.top().first;
                pq.pop();
                if (cost > dist[currNode]) continue;

                for (std::pair<node_index, T> child: (backwards ? inEdges : outEdges)[currNode]) {
                    if (child.second < 0) {
                        // negative edge weight, try `BellmanFord()` isntead
                        fill(prevNode.begin(), prevNode.end(), 0);
                        fill(dist.begin(), dist.end(), std::numeric_limits<T>::max());
                        BellmanFord(node, prevNode, dist, backwards);
                        return;
                    }

                    T newCost = cost + child.second;
                    if (newCost < dist[child.first]) {
                        dist[child.first] = newCost;
                        prevNode[child.first] = currNode;
                        pq.push({ newCost, child.first });

                    }
                }
            }
        }

        // O(VE) Standard Bellman-Ford algorithm - shortest path from `node` to all other nodes
        // @param `backwards` Indicates whether to consider `inEdges` or `outEdges`
        void BellmanFord(node_index node, std::vector<node_index>& prevNode, std::vector<node_index>& dist, bool backwards = false) const {
            for (node_index rep = 0; rep < V; ++rep) {
                for (std::pair<std::pair<node_index, node_index>, T> edge: edges) {
                    node_index u = edge.first.first, v = edge.first.second;
                    if (backwards) std::swap(u, v);

                    if (dist[u] == std::numeric_limits<T>::max()) continue;
                    T newCost = dist[u] + edge.second;
                    if (newCost < dist[v]) {
                        // last rep, but the minimum distance can still be reduced
                        assert(rep != V - 1 && "negative weight cycle");

                        dist[v] = newCost;
                        prevNode[v] = u;
                    }
                }
            }
        }

    public:
        // O(V + E log V), but O(VE) if negative edge weights. Finds the shortest path from `u` to `v`
        // @returns Distance of the path
        // @note Component containing `u` and `v` cannot have negative weight cycles
        T shortestDist(node_index u, node_index v) const {
            assert(validNode(u) && validNode(v) && "Node index out of range");
            std::vector<node_index> prevNode = std::vector<node_index>(V + 1, 0);
            std::vector<T> dist = std::vector<T>(V + 1, std::numeric_limits<T>::max());
            Dijkstra(u, prevNode, dist);

            return dist[v];
        }

        // O(V + E log V), but O(VE) if negative edge weights
        // @returns Distance of the shortest path from `u` to `v`, and a vector of nodes representing the path
        // @note Component containing `u` and `v` cannot have negative weight cycles
        std::pair<T, std::vector<node_index>> shortestPath(node_index u, node_index v) const {
            assert(validNode(u) && validNode(v) && "Node index out of range");
            std::vector<node_index> prevNode = std::vector<node_index>(V + 1, 0);
            std::vector<T> dist = std::vector<T>(V + 1, std::numeric_limits<T>::max());
            Dijkstra(u, prevNode, dist);

            std::vector<node_index> path;
            if (dist[v] != std::numeric_limits<T>::max()) {
                for (node_index node = v; node != u; node = prevNode[node]) path.push_back(node);
                path.push_back(u);
            }
            reverse(path.begin(), path.end());
            return { dist[v], path };
        }

        // O(V + E log V), but O(VE) if negative edge weights
        // @returns The eccentricity - the greatest distance between `root` and any other node in the same component
        // @note Component containing `root` cannot have negative weight cycles
        T eccentricity(node_index root) const {
            assert(validNode(root) && "Node index out of range");
            // uh should i return the furthest node too?
            std::vector<node_index> prevNode = std::vector<node_index>(V + 1, 0);
            std::vector<node_index> dist = std::vector<node_index>(V + 1, std::numeric_limits<T>::max());
            Dijkstra(root, prevNode, dist);
            T furthest = 0;
            for (node_index node = 1; node <= V; ++node) {
                if (dist[node] != std::numeric_limits<T>::max()) furthest = std::max(furthest, dist[node]);
            }
            return furthest;
        }

        // O(V + E log V), but O(VE) if negative edge weights. Finds the diameter - the simple path with the maximum distance
        // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
        // @returns Distance of the path, and a vector of nodes representing the path
        // @note Component containing `root` cannot have negative weight cycles
        std::vector<node_index> diameter(node_index root = 0) {
            if (root == 0) {
                std::vector<node_index> output;
                for (std::vector<node_index>& component: getComponentsNodes()) {
                    std::vector<node_index> curr = diameter(component[0]);
                    if (curr.size() > output.size()) output = curr;
                }
                return output;
            }

            assert(validNode(root) && "Node index out of range");
            std::vector<node_index> prevNode = std::vector<node_index>(V + 1, 0);
            std::vector<node_index> dist = std::vector<node_index>(V + 1, std::numeric_limits<T>::max());
            Dijkstra(root, prevNode, dist, true);

            T furthestDist1 = std::numeric_limits<T>::min();
            node_index furthestNode1 = root;
            for (node_index node = 1; node <= V; ++node) {
                if (dist[node] != std::numeric_limits<T>::max() && dist[node] > furthestDist1) {
                    furthestDist1 = dist[node];
                    furthestNode1 = node;
                }
            }

            prevNode = std::vector<node_index>(V + 1, 0);
            dist = std::vector<node_index>(V + 1, std::numeric_limits<T>::max());
            Dijkstra(furthestNode1, prevNode, dist, false);

            T furthestDist2 = std::numeric_limits<T>::min();
            node_index furthestNode2 = root;
            for (node_index node = 1; node <= V; ++node) {
                if (dist[node] != std::numeric_limits<T>::max() && dist[node] > furthestDist2) {
                    furthestDist2 = dist[node];
                    furthestNode2 = node;
                }
            }

            std::vector<node_index> path;
            for (node_index node = furthestNode2; node != prevNode[node]; node = prevNode[node]) {
                path.push_back(node);
            }
            path.push_back(furthestNode1);
            return path;
        }

    private:
        // O(V^3 + E) Standard Floyd-Warshall algorithm - shortest path between every 2 pairs of nodes
        void FloydWarshall(std::vector<std::vector<T>>& dists) const {
            for (node_index node = 1; node <= V; ++node) dists[node][node] = 0;
            for (std::pair<std::pair<node_index, node_index>, T>& edge: edges) {
                node_index u = edge.first.first, v = edge.first.second;
                dists[u][v] = std::min(dists[u][v], edge.second);
            }

            for (node_index mid = 1; mid <= V; ++mid) {
                for (node_index u = 1; u <= V; ++u) {
                    for (node_index v = 1; v <= V; ++v) {
                        if (std::max(dists[u][mid], dists[mid][v]) != std::numeric_limits<T>::max()) {
                            dists[u][v] = std::min(dists[u][v], dists[u][mid] + dists[mid][v]);
                        }
                    }
                }
            }
        }

    public:
        // O(V^3 + E) Finds the minimum distance between every pair of nodes
        // @returns A vector of vectors, with `dist[u][v]` being the minimim distance on the path from `u` to `v`
        std::vector<std::vector<T>> allShortestDist() const {
            std::vector<std::vector<node_index>> dists = std::vector<std::vector<node_index>>(V + 1, std::vector<node_index>(V + 1, std::numeric_limits<T>::max()));
            FloydWarshall();
            return dists;
        }

    private:
        // Amortised O(log V)
        // @returns The parent node of the current disjoint set
        node_index DSUgetParent(node_index node, std::vector<node_index>& parent) const {
            if (parent[node] == node) return node;
            return parent[node] = DSUgetParent(parent[node], parent);
        }

        // Amortised O(E log V)
        // Standard Kruskal's algorithm - finds the cost of the minimum spanning tree, using union find (DSU)
        T KruskalCost() const {
            std::vector<node_index> parent = std::vector<node_index>(V + 1);
            std::iota(parent.begin(), parent.end(), 0);

            std::vector<std::pair<T, std::pair<node_index, node_index>>> sortedEdges;
            for (std::pair<std::pair<node_index, node_index>, T> edge: edges) {
                sortedEdges.push_back({ edge.second, edge.first });
            }
            sort(sortedEdges.begin(), sortedEdges.end());

            T cost = 0;
            for (std::pair<T, std::pair<node_index, node_index>>& edge: sortedEdges) {
                node_index parentU = DSUgetParent(edge.second.first, parent), parentV = DSUgetParent(edge.second.second, parent);
                if (parentU != parentV) {
                    cost += edge.first;
                    parent[parentU] = parentV;
                }
            }
            return cost;
        }

        // Amortised O(E log V)
        // Standard Kruskal's algorithm - finds the minimum spanning tree, using union find (DSU)
        Graph<T> Kruskal() const {
            std::vector<node_index> parent = std::vector<node_index>(V + 1);
            std::iota(parent.begin(), parent.end(), 0);

            std::vector<std::pair<T, std::pair<node_index, node_index>>> sortedEdges;
            for (std::pair<std::pair<node_index, node_index>, T> edge: edges) {
                sortedEdges.push_back({ edge.second, edge.first });
            }
            sort(sortedEdges.begin(), sortedEdges.end());

            Graph<T> output(V, isWeighted, isDirected);
            for (std::pair<T, std::pair<node_index, node_index>>& edge: sortedEdges) {
                node_index parentU = DSUgetParent(edge.second.first, parent), parentV = DSUgetParent(edge.second.second, parent);
                if (parentU != parentV) {
                    output.addEdge(parentU, parentV, edge.first);
                    parent[parentU] = parentV;
                }
            }
            return output;
        }

    public:
        // Amortised O(E log V)
        // @returns The cost of the MST
        T MSTcost() const {
            return KruskalCost();
        }

        // Amortised O(E log V)
        // @returns A Graph object of MST
        Graph<T> MST() const {
            return Kruskal();
        }

    private:
        // O(V + E)
        // Standard Tarjan algorithm - find strongly connected components
        node_index Tarjan(node_index u, node_index index, std::vector<std::vector<node_index>>& components, std::vector<node_index>& s, std::vector<bool>& seen, std::vector<int>& getIndex, std::vector<int>& lowLink) {
            getIndex[u] = lowLink[u] = index;
            ++index;
            s.push_back(u);
            seen[u] = true;

            for (std::pair<node_index, T> edge: outEdges[u]) {
                node_index v = edge.first;
                if (getIndex[v] == -1) {
                    index = Tarjan(v, index, components, s, seen, getIndex, lowLink);
                    lowLink[u] = std::min(lowLink[u], lowLink[v]);
                }
                else if (seen[v]) {
                    lowLink[u] = std::min(lowLink[u], getIndex[v]);
                }
            }

            if (getIndex[u] == lowLink[u]) {
                components.push_back({});

                node_index node;
                do {
                    node = s.back();
                    s.pop_back();
                    seen[node] = false;
                    components.back().push_back(node);
                } while (node != u);
            }

            return index;
        }

    public:
        // O(V + E) 
        // @returns Vector of all strongly connected components, where each component is a vector of nodes
        std::vector<std::vector<node_index>> getSCCnodes() const {
            std::vector<bool> seen = std::vector<bool>(V + 1, false);
            std::vector<int> getIndex = std::vector<int>(V + 1, -1);
            std::vector<int> lowLink = std::vector<int>(V + 1, 0);
            std::vector<std::vector<node_index>> components;
            std::vector<node_index> s;

            for (node_index node = 1; node <= V; ++node) {
                if (getIndex[node] == -1) Tarjan(node, 1, components, s, seen, getIndex, lowLink);
            }
            return components;
        }

        // O(V + E)
        // Finds a strongly connected component
        // @param `root` Consider its connected component
        // @returns The connected components reachable from root, where each component is represented as a vector of nodes
        std::vector<node_index> getSCCnodes(node_index root) const {
            assert(isDirected && "In an undirected graph, all components are strongly connected components");
            assert(validNode(root) && "Node index out of range");

            std::vector<bool> seen = std::vector<bool>(V + 1, false);
            std::vector<int> getIndex = std::vector<int>(V + 1, -1);
            std::vector<node_index> lowLink = std::vector<node_index>(V + 1, 0);
            std::vector<std::vector<node_index>> components;
            std::vector<node_index> s;

            Tarjan(root, 1, components, s, seen, getIndex, lowLink);
            assert(components.size() == 1 && "Sanity check!");
            return components[0];
        }

    private:
        // O(V + E)
        // Standard DFS Bridge finding algorithm
        node_index _getBridges(node_index u, node_index parent, node_index t, std::vector<std::pair<std::pair<node_index, node_index>, T>>& bridges, std::vector<bool>& seen, std::vector<node_index>& minTime, std::vector<node_index>& entryTime) const {
            seen[u] = true;
            minTime[u] = entryTime[u] = ++t;

            for (std::pair<node_index, T> child: outEdges[u]) {
                node_index v = child.first;
                T w = child.second;
                if (v == parent) continue;
                if (seen[v]) {
                    minTime[u] = std::min(minTime[u], entryTime[v]);
                }
                else {
                    t = _getBridges(v, u, t, bridges, seen, minTime, entryTime);
                    minTime[u] = std::min(minTime[u], minTime[v]);

                    // bridge iff minTime[childNode] > entryTime[node] and the edge is not a double edge
                    if (minTime[v] > entryTime[u]) {
                        if (numEdges(u, v) == 1) bridges.push_back({ { u, v }, w });
                    }
                }
            }
            return t;
        }

    public:
        // O(V + E)
        // Find all bridges
        // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
        // @return A vector of edges, where each edge is represented by ((u, v), w)
        std::vector<std::pair<std::pair<node_index, node_index>, T>> getBridges(node_index root = 0) const {
            std::vector<bool> seen = std::vector<bool>(V + 1, false);
            std::vector<node_index> minTime = std::vector<node_index>(V + 1, 0), entryTime = std::vector<node_index>(V + 1, 0);

            std::vector<std::pair<std::pair<node_index, node_index>, T>> bridges;
            if (root == 0) {
                for (node_index node = 1; node <= V; ++node) {
                    if (!seen[node]) _getBridges(node, 0, 0, bridges, seen, minTime, entryTime);
                }
            }
            else {
                assert(validNode(root) && "Node index out of range");
                _getBridges(root, 0, 0, bridges, seen, minTime, entryTime);
            }
            return bridges;
        }

    private:
        // O(V + E log E)
        // Greedy colouring
        void _greedyColouring(node_index node, std::vector<int>& colours) const {
            if (colours[node] != -1) return;

            // calculate mex of neighbours colours
            std::set<int> usedColours;
            for (std::pair<node_index, T> edge: getEdges(node)) usedColours.insert(colours[edge.first]);
            for (colours[node] = 0; usedColours.find(colours[node]) != usedColours.end(); ++colours[node]);

            for (std::pair<node_index, T> edge: getEdges(node)) _greedyColouring(edge.first, colours);
        }

    public:
        // O(V_component log V_component + E log E)
        // @returns A possible greedy colouring of the graph - `0 <= colour[node] < max_colours`
        // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
        // @note Ignore colour[0]
        std::vector<int> greedyColouring(node_index root = 0) const {
            assert((root == 0 || validNode(root)) && "Node index out of range");
            std::vector<int> colours(V + 1, -1);
            for (node_index node: getComponentNodes(root)) _greedyColouring(node, colours);
            return colours;
        }

    private:
        // O(V_component log V_component + E_component)
        // Standard Kahn's algorithm - topological sort
        std::vector<node_index> Kahn(node_index root = 0) const {
            std::vector<node_index> depths = std::vector<node_index>(V + 1, 0);

            std::vector<node_index> q, topSort;
            for (node_index node: getComponentNodes(root)) {
                depths[node] = inEdges[node].size();
                if (inEdges[node].empty()) q.push_back(node);
            }

            while (!q.empty()) {
                node_index u = q.back();
                q.pop_back();
                topSort.push_back(u);

                for (std::pair<node_index, T> edge: outEdges[u]) {
                    node_index v = edge.first;
                    depths[v]--;
                    if (depths[v] == 0) q.push_back(v);
                }
            }
            return topSort;
        }

        // O(V_component)
        // Standard DFS topsort
        void dfsTopSort(node_index node, std::vector<node_index>& topSort, std::vector<node_index>& state) const {
            if (state[node] == 2) return;
            assert(state[node] == 1 && "Graph contains cycle");

            state[node] = 1;
            for (std::pair<node_index, T> edge: outEdges[node]) dfsTopSort(edge.first, topSort, state);
            state[node] = 2;
            topSort.push_back(node);
        }

    public:
        // O(V_component log V_component + E_component)
        // @returns A topological sorted order of nodes
        // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
        // @note Graph cannot be bidirectional or cyclical
        std::vector<node_index> getTopSort(node_index root = 0, bool doKahn = false) const {
            assert(isDirected && "Can't get a topological sort for a bidirectional graph");
            assert((root == 0 || validNode(root)) && "Node index out of range");

            if (doKahn) {
                std::vector<node_index> topSort = Kahn(root);
                assert(topSort.size() <= V && "Sanity check!");
                assert(topSort.size() == V && "Graph contains cycle");
                return topSort;
            }

            else {
                std::vector<node_index> state = std::vector<node_index>(V + 1, 0);

                std::vector<node_index> topSort;
                for (node_index node: getComponentNodes(root)) dfsTopSort(node, topSort, state);
                reverse(topSort.begin(), topSort.end());
                return topSort;
            }
        }

    private:
        T addFlow(node_index source, node_index sink, const std::vector<node_index> &prev, std::map<std::pair<node_index, node_index>, T> &remainingCap) const {
            // find max flow of augmenting path
            T newFlow = std::numeric_limits<T>::max();
            for (node_index u = sink; u != source; u = prev[u]) {
                newFlow = std::min(newFlow, remainingCap[{prev[u], u}]);
            }

            // update path
            for (node_index u = sink; u != source; u = prev[u]) {
                remainingCap[{prev[u], u}] -= newFlow;
                remainingCap[{u, prev[u]}] += newFlow;
            }
            return newFlow;
        }

    public:
        // O(V E^2 log E) or O(V^2 E log E)
        // @returns The value of max flow, using (fast) Dinic's or (slow) Edmonds-Karp
        T maxFlow(node_index source, node_index sink, bool doSlow = false) const {
            assert(validNode(source) && validNode(sink) && "Node index out of range");
            if (source == sink) return std::numeric_limits<T>::max();
            std::map<std::pair<node_index, node_index>, T> remainingCap;
            for (std::pair<std::pair<node_index, node_index>, T> edge: edges) {
                remainingCap[edge.first] += edge.second;
            }

            T totalFlow = T(0);

            while (true) {
                std::vector<int> prev = std::vector<int>(V + 1, -1);
                std::vector<int> levels = std::vector<int>(V + 1, -1);
                std::queue<node_index> q;
                levels[source] = 0;
                q.push(source);

                while (!q.empty()) {
                    node_index u = q.front();
                    q.pop();
                    for (node_index type = 0; type <= 1; ++type) {
                        for (std::pair<node_index, T> edge: (type ? inEdges : outEdges)[u]) {
                            node_index v = edge.first;
                            if (prev[v] == -1 && remainingCap[{u, v}] > 0) {
                                prev[v] = u;
                                levels[v] = levels[u] + 1;
                                q.push(v);
                            }
                        }
                    }
                }

                // no augmenting path, so max flow has been found
                if (prev[sink] == -1) return totalFlow;

                if (doSlow) {
                    // Edmonds-Karp (implmentation of Ford-Fulkerson) - add any augmenting path
                    totalFlow += addFlow(source, sink, prev, remainingCap);
                }
                else {
                    // Dinic's - find best augmenting path in residual "level" graph
                    node_index newFlow;
                    do {
                        std::stack<node_index> s;
                        s.push(source);
                        while (!s.empty() && s.top() != sink) {
                            node_index u = s.top();
                            s.pop();

                            for (node_index type = 0; type <= 1; ++type) {
                                for (std::pair<node_index, T> edge: (type ? inEdges : outEdges)[u]) {
                                    node_index v = edge.first;
                                    if (levels[v] == levels[u] + 1 && remainingCap[{u, v}] > 0) {
                                        prev[v] = u;
                                        q.push(v);
                                    }
                                }
                            }
                        }

                        newFlow = addFlow(source, sink, prev, remainingCap);
                        totalFlow += newFlow;
                    } while (newFlow != 0);
                }
            }
        }

        // TODO: is bipartite graph?
        // TODO: edge matching
        // TODO: see: https://usaco.guide/adv/

        /************************************************
         *                    PROPERTIES                *
         ***********************************************/

        // O(1)
        // @returns Whether the graph is directed
        inline const bool isDirectedGraph() const {
            return isDirected;
        }

        // O(1)
        // @returns Whether the graph is weighted
        inline const bool isWeightedGraph() const {
            return isWeighted;
        }

        // O(V_component log V_component + V_component log E + E_component)
        // @returns Whether there exists a edge which connects a node to itself
        // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
        bool hasSelfEdges(node_index root = 0) const {
            assert((root == 0 || validNode(root)) && "Node index out of range");
            for (node_index node: getComponentNodes(root)) {
                if (containsEdge(node, node)) return true;
            }
            return false;
        }

        // O(V_component log V_component + E_component)
        // @returns Whether there exists "redundant" edges
        // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
        bool hasDoubleEdges(node_index root = 0) const {
            for (node_index node: getComponentNodes(root)) {
                for (node_index type = 0; type <= 1; ++type) {
                    const std::multiset<std::pair<node_index, T>> *currEdges = &(type ? inEdges[node] : outEdges[node]);
                    if (currEdges->empty()) continue;
                    for (auto it1 = currEdges->begin(), it2 = ++currEdges->begin(); it2 != currEdges->end(); ++it1, ++it2) {
                        if (it1->first == it2->first) return true;
                    }
                }
            }
            
            return false;
        }

        // O(V_component log V_component + E_component)
        // @returns Whether the graph has a cycle
        // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
        bool hasCycle(node_index root = 0) const {
            assert((root == 0 || validNode(root)) && "Node index out of range");
            return (node_index)Kahn(root).size() < V;
        }

        // O(V_component log V_component + V_component log E + E_component)
        // @returns Whether the graph is a simple graph
        // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
        bool isSimpleGraph(node_index root = 0) const {
            assert((root == 0 || validNode(root)) && "Node index out of range");
            return !hasSelfEdges(root) && !hasDoubleEdges(root);
        }

        // O(V_component log V_component + E_component)
        // @returns Whether the graph is a directed acrylic graph
        // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
        bool isDAG(node_index root = 0) const {
            assert((root == 0 || validNode(root)) && "Node index out of range");
            return isDirected && !hasCycle(root);
        }

        // O(V_component log V_component + E_component)
        // @returns Whether the graph is a tree
        // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
        bool isTree(node_index root = 0) const {
            assert((root == 0 || validNode(root)) && "Node index out of range");
            if (hasCycle(root) || hasDoubleEdges(root)) return false;

            std::vector<node_index> component = getComponentNodes(root);
            node_index numEdges = 0, numLeaves = 0, numRoots = 0;
            for (node_index node: component) {
                for (std::pair<node_index, T> edge: outEdges[node]) {
                    if (!isDirected && edge.first > node) break;
                    else ++numEdges;
                }
                numLeaves += outEdges[node].empty();
                numRoots += inEdges[node].empty();
            }

            if (isDirected && numLeaves > 1 && numRoots > 1) return false;
            return numEdges == component.size() - 1; // E == V - 1
        }

        // O(V + E) Determines whether the graph is a forest
        bool isForest() const {
            for (std::vector<node_index>& component: getComponentsNodes()) {
                if (!isTree(component[0])) return false;
            }
            return true;
        }

        /************************************************
         *       ALGORITHMS (COMPLEXITLY CLASS NP)      *
         ************************************************/

        // TODO: hamiltonian path / tour
        // TODO: euler path / tour
        // TODO: covering set
        // TODO: optimal colouring
        // TODO: planar embedding?
    };
};

#endif
