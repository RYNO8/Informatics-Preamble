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
    struct UnitWeight {
        int w = 1;
        UnitWeight(): w(1) {}
        UnitWeight(int _w): w(_w) {}
        friend UnitWeight operator+(const UnitWeight &w1, const UnitWeight &w2) {
            return UnitWeight(w1.w + w2.w);
        }
        friend bool operator<(const UnitWeight &a, const UnitWeight &b) {
            return a.w < b.w;
        }
        friend bool operator==(const UnitWeight &a, const UnitWeight &b) {
            return a.w == b.w;
        }
        // lol
        friend std::istream& operator>>(std::istream &in, const UnitWeight &w) {
            return in;
        }
        friend std::ostream& operator<<(std::ostream &out, const UnitWeight &w) {
            out << w.w;
            return out;
        }
    };

    // TODO: define restrictions for weight for each algorithm
    // Weight needs a default constructor `Weight()`, which can initialise to anything
    // Weight needs the following methods, which need to behave properly
    //  - `Weight::operator==`
    //  - `Weight::operator<`
    //  - `Weight::operator+`
    //  - read from istream, display to ostream
    // A insufficiently defined Weight might still produce correct results sometimes
    template<
        size_t maxV,
        // weight data held for each edge
        typename Weight,
        // Whether edges are directed or bidirectional
        bool isDirected = false
    >
    class Graph {
public:
        /************************************************
         *                 INITIALISATION               *
         ************************************************/

        // id of node is an unsigned integer
        using Node = size_t;
        using isWeighted = std::negation<std::is_same<Weight, UnitWeight>>;

        struct Edge {
            using Node = size_t;
            Node u, v;
            Weight w;

            Edge() {
                
            }

            Edge(Node u_, Node v_, Weight w_) : u(u_), v(v_), w(w_) {
                assert(u != v && "Self cycles not supported yet");
                assert(u < maxV && v < maxV && "Node index out of range");
                if (!isDirected && u > v) std::swap(u, v);
            }

            // #if __cplusplus >= 202002L
            // template<std::enable_if_t<isWeighted::value, bool> = false>
            // #endif
            Edge(Node u_, Node v_) : u(u_), v(v_), w(Weight()) {
                assert(u_ != v_ && "Self cycles not supported yet");
                assert(u < maxV && v < maxV && "Node index out of range");
                if (!isDirected && u > v) std::swap(u, v);
            }

            Edge flip() const {
                assert(isDirected && "cannot flip undirected edge");
                return Edge(v, u, w);
            }

            friend bool operator<(const Edge &a, const Edge &b) {
                if (a.u != b.u) return a.u < b.u;
                else if (a.v != b.v) return a.v < b.v;
                else return a.w < b.w;
            }

            friend bool operator==(const Edge &a, const Edge &b) {
                return a.u == b.u && a.v == b.v && a.w == b.w;
            }

            friend std::istream& operator>>(std::istream &in, Edge &e) {
                in >> e.u >> e.v >> e.w;
                return in;
            }
        };

        // comparison ignoring weights - use this in containers
        struct EdgeComp {
            bool operator() (const Edge &a, const Edge &b) const {
                if (a.u != b.u) return a.u < b.u;
                return a.v < b.v;
            }
        };

    private:
        std::vector<Node> nodes; // not garunteed sorted!
        bool validNode[maxV + 1];

        // TODO: is it reasonable to support multiset?
        // TODO: `inEdges` and `outEdges` store copies of each edge - not memory efficient?
        // A set of all edges
        std::set<Edge, EdgeComp> edges;
        // `inEdges[node]` is the set of edges leading into `node`
        std::set<Edge, EdgeComp> inEdges[maxV + 1];
        // `outEdges[node]` is the set of edges leading out of `node`
        std::set<Edge, EdgeComp> outEdges[maxV + 1];

    public:
        // O(V)
        // Initialises an empty graph with no nodes
        Graph() {
            std::fill(std::begin(validNode), std::end(validNode), false);
        }

        // O(V)
        // Initialises an empty graph
        Graph(size_t V) {
            std::fill(std::begin(validNode), std::end(validNode), false);
            for (Node node = 1; node <= V; ++node) {
                nodes.push_back(node);
                validNode[node] = true;
            }
        }

        // O(V + |nodes| log E)
        // Initialise this graph as the graph induced from `g` by `nodes`
        // @note `V` does not change, so there may be disconnected nodes
        Graph(const Graph<maxV, Weight, isDirected> &g, const std::vector<Node> &nodes) {
            for (Node node : nodes) assert(g.containsNode(node) && "Node index out of range");
            for (Node node : nodes) {
                for (Edge incident : g.getEdges(node)) {
                    if (!containsNode(incident.u)) pushNode(incident.u);
                    if (!containsNode(incident.v)) pushNode(incident.v);
                    if (!containsEdge(incident)) insertEdge(incident);
                }
            }
        }

        // O(V + E log E)
        // Initalises an unweighted graph from vector of edges, assuming 1 indexed
        Graph(size_t V, const std::vector<Edge> &_edges) {
            std::fill(std::begin(validNode), std::end(validNode), false);
            for (Node node = 1; node <= V; ++node) {
                nodes.push_back(node);
                validNode[node] = true;
            }

            for (Edge edge: _edges) insertEdge(edge);
        }

        /************************************************
         *                    DISPLAY                   *
         ************************************************/

        // O(V + E)
        // Displays the graph, showing the outwards edge connections of each edge in lexographic order
        // @param `out` The string representation of the graph is piped to this output stream
        friend std::ostream& operator<<(std::ostream &out, const Graph<maxV, Weight, isDirected> &graph) {
            for (Node u : graph.nodes) {
                out << u << ':';
                for (Edge edge: graph.outEdges[u]) {
                    // on an undirected graph, don't print edges twice
                    if (!isDirected && u > edge.u) continue;

                    out << ' ' << edge.v;
                    if (isWeighted::value) out << " (w = " << edge.w << ')';
                }
                out << '\n';
            }
            return out;
        }


        /************************************************
         *                  PROPERTIES                  *
         ************************************************/

        // O(1)
        // @returns `V`, the number of nodes
        inline size_t size() const {
            return nodes.size();
        }
        // @returns `V`, the number of nodes
        inline size_t V() const {
            return nodes.size();
        }
        // @returns `V`, the number of nodes
        inline size_t N() const {
            return nodes.size();
        }

        // O(1)
        // @returns `E`, the number of edges
        inline size_t E() const {
            return edges.size();
        }
        // O(1)
        // @returns `E`, the number of edges
        inline size_t M() const {
            return edges.size();
        }

        // @ returns the imutable set of all ndoes
        inline const std::vector<Node>& getNodes() const {
            return nodes;
        }

        // O(1)
        // @returns the imutable set of all edges
        const std::vector<Edge> getEdges() const {
            return std::vector<Edge>(edges.begin(), edges.end());
        }

        /************************************************
         *                INCIDENT DATA                 *
         ************************************************/

        // O(num_neighbours log num_neighbours) = O(N log N)
        // @returns the imutable set of all edges
        const std::vector<Edge> getEdges(Node node) const {
            assert(containsNode(node) && "Node index out of range");
            std::vector<Edge> out(inEdges[node].size() + outEdges[node].size());
            auto it = std::set_union(inEdges[node].begin(), inEdges[node].end(), outEdges[node].begin(), outEdges[node].end(), out.begin());
            out.resize(it - out.begin());
            return out;
        }

        // O(E)
        // Finds the edges that travel into `node`
        // @returns A multiset of (neighbour, weight) pairs
        inline std::vector<Edge> getEdgesIn(Node node) const {
            assert(containsNode(node) && "Node index out of range");
            return std::vector<Edge>(inEdges[node].begin(), inEdges[node].end());
        }

        // O(E)
        // Finds the edges that travel out of `node`
        // @returns A multiset of (neighbour, weight) pairs
        inline std::vector<Edge> getEdgesOut(Node node) const {
            assert(containsNode(node) && "Node index out of range");
            return std::vector<Edge>(outEdges[node].begin(), outEdges[node].end());
        }

        // O(1)
        // @returns The total degree of `node` (combination in in degree and out degree)
        inline const size_t degree(Node node) const {
            assert(containsNode(node) && "Node index out of range");
            if (isDirected) {
                assert(inEdges[node].size() == outEdges[node].size() && "Sanity check!");
                return inEdges[node].size();
            } else {
                return inEdges[node].size() + outEdges[node].size();
            }
        }

        // O(1)
        // @returns The in degree of `node`
        inline const size_t degreeIn(Node node) const {
            assert(containsNode(node) && "Node index out of range");
            return inEdges[node].size();
        }

        // O(1)
        // @returns The out degree of `node`
        inline const size_t degreeOut(Node node) const {
            assert(containsNode(node) && "Node index out of range");
            return outEdges[node].size();
        }

        /************************************************
         *               GRAPH MODIFICATIONS            *
         ************************************************/

        // O(1)
        // Adds a new node with no edges
        // @returns The index of this node
        // TODO: should I have a variation where silently passes when node exists?
        inline size_t pushNode(Node node) {
            assert(/*0 <= node &&*/ node < maxV + 1 && "Node index out of range");
            assert(!containsNode(node) && "node already exists");
            nodes.push_back(node);
            validNode[node] = true;
            return size();
        }

        // O(log E)
        // Add the specified edge to the graph - this may create a double edge
        // @note Edge weight is ignored if the graph is unweighted
        // TODO: should I have a variation where silently passes when edge exists?
        void insertEdge(Edge e) {
            assert(containsNode(e.u) && containsNode(e.v) && "Node index out of range");
            assert(!edges.count(e) && "Edge already in graph, is this an unintended operation?");

            edges.insert(e);
            // NOTE: this is preventing support of multiedges
            // when duplicaate edge is added, how do you find the most recently added edge
            inEdges[e.v].insert(e);
            outEdges[e.u].insert(e);
            if (!isDirected) {
                inEdges[e.u].insert(e);
                outEdges[e.v].insert(e);
            }
        }

        // O(log E)
        // If the specified edge (with any edge weight) is present, remove it, otherwise silently do nothing
        inline void eraseEdge(Edge e) {
            assert(containsNode(e.u) && containsNode(e.v) && "Node index out of range");
            assert(containsEdge(e) && "Cannot remove edge which isn't in graph");

            edges.erase(e);
            inEdges[e.v].erase(e);
            outEdges[e.u].erase(e);
            if (!isDirected) {
                inEdges[e.u].erase(e);
                outEdges[e.v].erase(e);
            }
        }


        /************************************************
         *                    CONTAINS                  *
         ************************************************/

        inline bool containsNode(Node node) const {
            return /*0 <= node &&*/ node < maxV + 1 && validNode[node];
        }

        // O(log E)
        // @returns Whether the specified edge (with any edge weight) is contained in the graph
        inline bool containsEdge(Edge e) const {
            assert(containsNode(e.u) && containsNode(e.v) && "Node index out of range");
            auto edgeIt = edges.find(e);
            return edgeIt != edges.end() && edgeIt->w == e.w;
        }

        // O(log E)
        // @returns Whether the specified edge is contained in the graph
        inline bool containsEdgeUnweighted(Edge e) const {
            assert(containsNode(e.u) && containsNode(e.v) && "Node index out of range");
            return edges.count(e);
        }

        
        /************************************************
         *                   COMPONENTS                 *
         ************************************************/

        // O(V_component log V_component + E_component log E_component)
        // @returns All nodes in the same component as `root`, in sorted order
        std::vector<Node> getComponent(Node at) const {
            assert(containsNode(at) && "Node index out of range");
            std::unordered_set<Node> seen;
            std::queue<Node> q;
            seen.insert(at);
            q.push(at);

            while (!q.empty()) {
                Node node = q.front();
                q.pop();

                for (Edge incident : getEdges(node)) {
                    if (!seen.count(incident.u)) {
                        seen.insert(incident.u);
                        q.push(incident.u);
                    }
                    if (!seen.count(incident.v)) {
                        seen.insert(incident.v);
                        q.push(incident.v);
                    }
                };
            }
            std::vector<Node> component = std::vector<Node>(seen.begin(), seen.end());
            sort(component.begin(), component.end());
            return component;
        }

        // O(maxV + E)
        // @returns All nodes grouped by their component, in arbitary order
        std::vector<std::vector<Node>> getComponents() const {
            std::vector<std::vector<Node>> output;
            bool seen[maxV];
            std::fill(std::begin(seen), std::end(seen), false);

            for (Node node : getNodes()) {
                if (!seen[node]) {
                    seen[node] = true;
                    output.push_back({ node });
                    for (size_t i = 0; i < output.back().size(); ++i) {
                        for (Edge incident : getEdges(output.back()[i])) {
                            if (!seen[incident.u]) {
                                seen[incident.u] = true;
                                output.back().push_back(incident.u);
                            }
                            if (!seen[incident.v]) {
                                seen[incident.v] = true;
                                output.back().push_back(incident.v);
                            }
                        };
                    }
                    // this adds a log factor
                    // sort(output.begin(), output.end());
                }
            }
            return output;
        }

        /************************************************
         *        ALGORITHMS (COMPLEXITLY CLASS P)      *
         ***********************************************/

        // O(V + E log E)
        // @returns The union of 2 graphs
        // @note Does not cause edge doubling
        friend Graph<maxV, Weight, isDirected> operator+(const Graph<maxV, Weight, isDirected> &a, const Graph<maxV, Weight, isDirected> &b) {
            Graph<maxV, Weight, isDirected> out = a; // make a copy
            for (Edge edge: b.getEdges()) {
                if (!out.containsNode(edge.u)) out.pushNode(edge.u);
                if (!out.containsNode(edge.v)) out.pushNode(edge.v);
                if (!out.containsEdge(edge)) out.insertEdge(edge);
            }
            return out;
        }

        // O(V + E)
        // @returns whether the 2 graphs are identical in edge (including weights!) and vertex sets
        // @note has short circuiting
        friend bool operator==(const Graph<maxV, Weight, isDirected> &a, const Graph<maxV, Weight, isDirected> &b) {
            if (a.E() != b.E()) return false;
            if (a.V() != b.V()) return false;
            for (Node node = 1; node <= maxV; ++node) {
                if (a.containsNode(node) != b.containsNode(node)) return false;
            }
            return true;
        }

        friend bool operator!=(const Graph<maxV, Weight, isDirected> &a, const Graph<maxV, Weight, isDirected> &b) {
            return !(a == b);
        }

        // O(V + E log E)
        Graph<maxV, Weight, isDirected> flip() const {
            assert(isDirected && "cannot flip undirected graph");
            Graph<maxV, Weight, isDirected> out;
            for (Node node: getNodes()) out.pushNode(node);
            for (Edge edge: getEdges()) out.insertEdge(edge.flip());
            return out;
        }

        // O(V + E)
        inline void clear() {
            edges.clear();
            for (Node node : getNodes()) {
                inEdges[node].clear();
                outEdges[node].clear();
            }
            std::fill(std::begin(validNode), std::end(validNode), false);
            nodes.clear();
        }

    // private:
    //     // O(V + E log V), or O(VE) if negative edge weights
    //     // Standard Dijkstra algorithm - shortest path from `node` to all other nodes for non-negative edge weights
    //     // @param `backwards` Indicates whether to consider `inEdges` or `outEdges`
    //     // @note If the graph is undirected, this becomes equivalent to a bfs (but slightly slower)
    //     // @note If negative edge weights encountered, switches to `BellmanFord()`
    //     void Dijkstra(Node node, std::vector<Node>& prevNode, std::vector<T>& dist, bool backwards = false) const {
    //         dist[node] = T(0);
    //         prevNode[node] = node;
    //         std::priority_queue<std::pair<T, Node>, std::vector<std::pair<T, Node>>, std::greater<std::pair<T, Node>>> pq;
    //         pq.push({ 0, node });

    //         while (!pq.empty()) {
    //             Node currNode = pq.top().second;
    //             T cost = pq.top().first;
    //             pq.pop();
    //             if (cost > dist[currNode]) continue;

    //             for (std::pair<Node, T> child: (backwards ? inEdges : outEdges)[currNode]) {
    //                 if (child.second < 0) {
    //                     // negative edge weight, try `BellmanFord()` isntead
    //                     fill(prevNode.begin(), prevNode.end(), 0);
    //                     fill(dist.begin(), dist.end(), std::numeric_limits<T>::max());
    //                     BellmanFord(node, prevNode, dist, backwards);
    //                     return;
    //                 }

    //                 T newCost = cost + child.second;
    //                 if (newCost < dist[child.first]) {
    //                     dist[child.first] = newCost;
    //                     prevNode[child.first] = currNode;
    //                     pq.push({ newCost, child.first });

    //                 }
    //             }
    //         }
    //     }

    //     // O(VE) Standard Bellman-Ford algorithm - shortest path from `node` to all other nodes
    //     // @param `backwards` Indicates whether to consider `inEdges` or `outEdges`
    //     void BellmanFord(Node node, std::vector<Node>& prevNode, std::vector<Node>& dist, bool backwards = false) const {
    //         for (Node rep = 0; rep < V; ++rep) {
    //             for (std::pair<std::pair<Node, Node>, T> edge: edges) {
    //                 Node u = edge.u.first, v = edge.u.second;
    //                 if (backwards) std::swap(u, v);

    //                 if (dist[u] == std::numeric_limits<T>::max()) continue;
    //                 T newCost = dist[u] + edge.v;
    //                 if (newCost < dist[v]) {
    //                     // last rep, but the minimum distance can still be reduced
    //                     assert(rep != V - 1 && "negative weight cycle");

    //                     dist[v] = newCost;
    //                     prevNode[v] = u;
    //                 }
    //             }
    //         }
    //     }

    // public:
    //     // O(V + E log V), but O(VE) if negative edge weights. Finds the shortest path from `u` to `v`
    //     // @returns Distance of the path
    //     // @note Component containing `u` and `v` cannot have negative weight cycles
    //     T shortestDist(Node u, Node v) const {
    //         assert(containsNode(u) && containsNode(v) && "Node index out of range");
    //         std::vector<Node> prevNode = std::vector<Node>(V + 1, 0);
    //         std::vector<T> dist = std::vector<T>(V + 1, std::numeric_limits<T>::max());
    //         Dijkstra(u, prevNode, dist);

    //         return dist[v];
    //     }

    //     // O(V + E log V), but O(VE) if negative edge weights
    //     // @returns Distance of the shortest path from `u` to `v`, and a vector of nodes representing the path
    //     // @note Component containing `u` and `v` cannot have negative weight cycles
    //     std::pair<T, std::vector<Node>> shortestPath(Node u, Node v) const {
    //         assert(containsNode(u) && containsNode(v) && "Node index out of range");
    //         std::vector<Node> prevNode = std::vector<Node>(V + 1, 0);
    //         std::vector<T> dist = std::vector<T>(V + 1, std::numeric_limits<T>::max());
    //         Dijkstra(u, prevNode, dist);

    //         std::vector<Node> path;
    //         if (dist[v] != std::numeric_limits<T>::max()) {
    //             for (Node node = v; node != u; node = prevNode[node]) path.push_back(node);
    //             path.push_back(u);
    //         }
    //         reverse(path.begin(), path.end());
    //         return { dist[v], path };
    //     }

    //     // O(V + E log V), but O(VE) if negative edge weights
    //     // @returns The eccentricity - the greatest distance between `root` and any other node in the same component
    //     // @note Component containing `root` cannot have negative weight cycles
    //     T eccentricity(Node root) const {
    //         assert(containsNode(root) && "Node index out of range");
    //         // uh should i return the furthest node too?
    //         std::vector<Node> prevNode = std::vector<Node>(V + 1, 0);
    //         std::vector<Node> dist = std::vector<Node>(V + 1, std::numeric_limits<T>::max());
    //         Dijkstra(root, prevNode, dist);
    //         T furthest = 0;
    //         for (Node node = 1; node <= V; ++node) {
    //             if (dist[node] != std::numeric_limits<T>::max()) furthest = std::max(furthest, dist[node]);
    //         }
    //         return furthest;
    //     }

    //     // O(V + E log V), but O(VE) if negative edge weights. Finds the diameter - the simple path with the maximum distance
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
    //     // @returns Distance of the path, and a vector of nodes representing the path
    //     // @note Component containing `root` cannot have negative weight cycles
    //     std::vector<Node> diameter(Node root = 0) {
    //         if (root == 0) {
    //             std::vector<Node> output;
    //             for (std::vector<Node>& component: getComponentsNodes()) {
    //                 std::vector<Node> curr = diameter(component[0]);
    //                 if (curr.size() > output.size()) output = curr;
    //             }
    //             return output;
    //         }

    //         assert(containsNode(root) && "Node index out of range");
    //         std::vector<Node> prevNode = std::vector<Node>(V + 1, 0);
    //         std::vector<Node> dist = std::vector<Node>(V + 1, std::numeric_limits<T>::max());
    //         Dijkstra(root, prevNode, dist, true);

    //         T furthestDist1 = std::numeric_limits<T>::min();
    //         Node furthestNode1 = root;
    //         for (Node node = 1; node <= V; ++node) {
    //             if (dist[node] != std::numeric_limits<T>::max() && dist[node] > furthestDist1) {
    //                 furthestDist1 = dist[node];
    //                 furthestNode1 = node;
    //             }
    //         }

    //         prevNode = std::vector<Node>(V + 1, 0);
    //         dist = std::vector<Node>(V + 1, std::numeric_limits<T>::max());
    //         Dijkstra(furthestNode1, prevNode, dist, false);

    //         T furthestDist2 = std::numeric_limits<T>::min();
    //         Node furthestNode2 = root;
    //         for (Node node = 1; node <= V; ++node) {
    //             if (dist[node] != std::numeric_limits<T>::max() && dist[node] > furthestDist2) {
    //                 furthestDist2 = dist[node];
    //                 furthestNode2 = node;
    //             }
    //         }

    //         std::vector<Node> path;
    //         for (Node node = furthestNode2; node != prevNode[node]; node = prevNode[node]) {
    //             path.push_back(node);
    //         }
    //         path.push_back(furthestNode1);
    //         return path;
    //     }

    // private:
    //     // O(V^3 + E) Standard Floyd-Warshall algorithm - shortest path between every 2 pairs of nodes
    //     void FloydWarshall(std::vector<std::vector<T>>& dists) const {
    //         for (Node node = 1; node <= V; ++node) dists[node][node] = 0;
    //         for (std::pair<std::pair<Node, Node>, T>& edge: edges) {
    //             Node u = edge.u.first, v = edge.u.second;
    //             dists[u][v] = std::min(dists[u][v], edge.v);
    //         }

    //         for (Node mid = 1; mid <= V; ++mid) {
    //             for (Node u = 1; u <= V; ++u) {
    //                 for (Node v = 1; v <= V; ++v) {
    //                     if (std::max(dists[u][mid], dists[mid][v]) != std::numeric_limits<T>::max()) {
    //                         dists[u][v] = std::min(dists[u][v], dists[u][mid] + dists[mid][v]);
    //                     }
    //                 }
    //             }
    //         }
    //     }

    // public:
    //     // O(V^3 + E) Finds the minimum distance between every pair of nodes
    //     // @returns A vector of vectors, with `dist[u][v]` being the minimim distance on the path from `u` to `v`
    //     std::vector<std::vector<T>> allShortestDist() const {
    //         std::vector<std::vector<Node>> dists = std::vector<std::vector<Node>>(V + 1, std::vector<Node>(V + 1, std::numeric_limits<T>::max()));
    //         FloydWarshall();
    //         return dists;
    //     }

    // private:
    //     // Amortised O(log V)
    //     // @returns The parent node of the current disjoint set
    //     Node DSUgetParent(Node node, std::vector<Node>& parent) const {
    //         if (parent[node] == node) return node;
    //         return parent[node] = DSUgetParent(parent[node], parent);
    //     }

    //     // Amortised O(E log V)
    //     // Standard Kruskal's algorithm - finds the cost of the minimum spanning tree, using union find (DSU)
    //     T KruskalCost() const {
    //         std::vector<Node> parent = std::vector<Node>(V + 1);
    //         std::iota(parent.begin(), parent.end(), 0);

    //         std::vector<std::pair<T, std::pair<Node, Node>>> sortedEdges;
    //         for (std::pair<std::pair<Node, Node>, T> edge: edges) {
    //             sortedEdges.push_back({ edge.v, edge.u });
    //         }
    //         sort(sortedEdges.begin(), sortedEdges.end());

    //         T cost = 0;
    //         for (std::pair<T, std::pair<Node, Node>>& edge: sortedEdges) {
    //             Node parentU = DSUgetParent(edge.v.first, parent), parentV = DSUgetParent(edge.v.second, parent);
    //             if (parentU != parentV) {
    //                 cost += edge.u;
    //                 parent[parentU] = parentV;
    //             }
    //         }
    //         return cost;
    //     }

    //     // Amortised O(E log V)
    //     // Standard Kruskal's algorithm - finds the minimum spanning tree, using union find (DSU)
    //     Graph<T> Kruskal() const {
    //         std::vector<Node> parent = std::vector<Node>(V + 1);
    //         std::iota(parent.begin(), parent.end(), 0);

    //         std::vector<std::pair<T, std::pair<Node, Node>>> sortedEdges;
    //         for (std::pair<std::pair<Node, Node>, T> edge: edges) {
    //             sortedEdges.push_back({ edge.v, edge.u });
    //         }
    //         sort(sortedEdges.begin(), sortedEdges.end());

    //         Graph<T> output(V, isWeighted, isDirected);
    //         for (std::pair<T, std::pair<Node, Node>>& edge: sortedEdges) {
    //             Node parentU = DSUgetParent(edge.v.first, parent), parentV = DSUgetParent(edge.v.second, parent);
    //             if (parentU != parentV) {
    //                 output.insertEdge(parentU, parentV, edge.u);
    //                 parent[parentU] = parentV;
    //             }
    //         }
    //         return output;
    //     }

    // public:
    //     // Amortised O(E log V)
    //     // @returns The cost of the MST
    //     T MSTcost() const {
    //         return KruskalCost();
    //     }

    //     // Amortised O(E log V)
    //     // @returns A Graph object of MST
    //     Graph<T> MST() const {
    //         return Kruskal();
    //     }

    // private:
    //     // O(V + E)
    //     // Standard Tarjan algorithm - find strongly connected components
    //     Node Tarjan(Node u, Node index, std::vector<std::vector<Node>>& components, std::vector<Node>& s, std::vector<bool>& seen, std::vector<int>& getIndex, std::vector<int>& lowLink) {
    //         getIndex[u] = lowLink[u] = index;
    //         ++index;
    //         s.push_back(u);
    //         seen[u] = true;

    //         for (std::pair<Node, T> edge: outEdges[u]) {
    //             Node v = edge.u;
    //             if (getIndex[v] == -1) {
    //                 index = Tarjan(v, index, components, s, seen, getIndex, lowLink);
    //                 lowLink[u] = std::min(lowLink[u], lowLink[v]);
    //             }
    //             else if (seen[v]) {
    //                 lowLink[u] = std::min(lowLink[u], getIndex[v]);
    //             }
    //         }

    //         if (getIndex[u] == lowLink[u]) {
    //             components.push_back({});

    //             Node node;
    //             do {
    //                 node = s.back();
    //                 s.pop_back();
    //                 seen[node] = false;
    //                 components.back().push_back(node);
    //             } while (node != u);
    //         }

    //         return index;
    //     }

    // public:
    //     // O(V + E) 
    //     // @returns Vector of all strongly connected components, where each component is a vector of nodes
    //     std::vector<std::vector<Node>> getSCCnodes() const {
    //         std::vector<bool> seen = std::vector<bool>(V + 1, false);
    //         std::vector<int> getIndex = std::vector<int>(V + 1, -1);
    //         std::vector<int> lowLink = std::vector<int>(V + 1, 0);
    //         std::vector<std::vector<Node>> components;
    //         std::vector<Node> s;

    //         for (Node node = 1; node <= V; ++node) {
    //             if (getIndex[node] == -1) Tarjan(node, 1, components, s, seen, getIndex, lowLink);
    //         }
    //         return components;
    //     }

    //     // O(V + E)
    //     // Finds a strongly connected component
    //     // @param `root` Consider its connected component
    //     // @returns The connected components reachable from root, where each component is represented as a vector of nodes
    //     std::vector<Node> getSCCnodes(Node root) const {
    //         assert(isDirected && "In an undirected graph, all components are strongly connected components");
    //         assert(containsNode(root) && "Node index out of range");

    //         std::vector<bool> seen = std::vector<bool>(V + 1, false);
    //         std::vector<int> getIndex = std::vector<int>(V + 1, -1);
    //         std::vector<Node> lowLink = std::vector<Node>(V + 1, 0);
    //         std::vector<std::vector<Node>> components;
    //         std::vector<Node> s;

    //         Tarjan(root, 1, components, s, seen, getIndex, lowLink);
    //         assert(components.size() == 1 && "Sanity check!");
    //         return components[0];
    //     }

    // private:
    //     // O(V + E)
    //     // Standard DFS Bridge finding algorithm
    //     Node _getBridges(Node u, Node parent, Node t, std::vector<std::pair<std::pair<Node, Node>, T>>& bridges, std::vector<bool>& seen, std::vector<Node>& minTime, std::vector<Node>& entryTime) const {
    //         seen[u] = true;
    //         minTime[u] = entryTime[u] = ++t;

    //         for (std::pair<Node, T> child: outEdges[u]) {
    //             Node v = child.first;
    //             T w = child.second;
    //             if (v == parent) continue;
    //             if (seen[v]) {
    //                 minTime[u] = std::min(minTime[u], entryTime[v]);
    //             }
    //             else {
    //                 t = _getBridges(v, u, t, bridges, seen, minTime, entryTime);
    //                 minTime[u] = std::min(minTime[u], minTime[v]);

    //                 // bridge iff minTime[childNode] > entryTime[node] and the edge is not a double edge
    //                 if (minTime[v] > entryTime[u]) {
    //                     if (numEdges(u, v) == 1) bridges.push_back({ { u, v }, w });
    //                 }
    //             }
    //         }
    //         return t;
    //     }

    // public:
    //     // O(V + E)
    //     // Find all bridges
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
    //     // @return A vector of edges, where each edge is represented by ((u, v), w)
    //     std::vector<std::pair<std::pair<Node, Node>, T>> getBridges(Node root = 0) const {
    //         std::vector<bool> seen = std::vector<bool>(V + 1, false);
    //         std::vector<Node> minTime = std::vector<Node>(V + 1, 0), entryTime = std::vector<Node>(V + 1, 0);

    //         std::vector<std::pair<std::pair<Node, Node>, T>> bridges;
    //         if (root == 0) {
    //             for (Node node = 1; node <= V; ++node) {
    //                 if (!seen[node]) _getBridges(node, 0, 0, bridges, seen, minTime, entryTime);
    //             }
    //         }
    //         else {
    //             assert(containsNode(root) && "Node index out of range");
    //             _getBridges(root, 0, 0, bridges, seen, minTime, entryTime);
    //         }
    //         return bridges;
    //     }

    // private:
    //     // O(V + E log E)
    //     // Greedy colouring
    //     void _greedyColouring(Node node, std::vector<int>& colours) const {
    //         if (colours[node] != -1) return;

    //         // calculate mex of neighbours colours
    //         std::set<int> usedColours;
    //         for (std::pair<Node, T> edge: getEdges(node)) usedColours.insert(colours[edge.u]);
    //         for (colours[node] = 0; usedColours.find(colours[node]) != usedColours.end(); ++colours[node]);

    //         for (std::pair<Node, T> edge: getEdges(node)) _greedyColouring(edge.u, colours);
    //     }

    // public:
    //     // O(V_component log V_component + E log E)
    //     // @returns A possible greedy colouring of the graph - `0 <= colour[node] < max_colours`
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
    //     // @note Ignore colour[0]
    //     std::vector<int> greedyColouring(Node root = 0) const {
    //         assert((root == 0 || containsNode(root)) && "Node index out of range");
    //         std::vector<int> colours(V + 1, -1);
    //         for (Node node: getComponentNodes(root)) _greedyColouring(node, colours);
    //         return colours;
    //     }

    // private:
    //     // O(V_component log V_component + E_component)
    //     // Standard Kahn's algorithm - topological sort
    //     std::vector<Node> Kahn(Node root = 0) const {
    //         std::vector<Node> depths = std::vector<Node>(V + 1, 0);

    //         std::vector<Node> q, topSort;
    //         for (Node node: getComponentNodes(root)) {
    //             depths[node] = inEdges[node].size();
    //             if (inEdges[node].empty()) q.push_back(node);
    //         }

    //         while (!q.empty()) {
    //             Node u = q.back();
    //             q.pop_back();
    //             topSort.push_back(u);

    //             for (std::pair<Node, T> edge: outEdges[u]) {
    //                 Node v = edge.u;
    //                 depths[v]--;
    //                 if (depths[v] == 0) q.push_back(v);
    //             }
    //         }
    //         return topSort;
    //     }

    //     // O(V_component)
    //     // Standard DFS topsort
    //     void dfsTopSort(Node node, std::vector<Node>& topSort, std::vector<Node>& state) const {
    //         if (state[node] == 2) return;
    //         assert(state[node] == 1 && "Graph contains cycle");

    //         state[node] = 1;
    //         for (std::pair<Node, T> edge: outEdges[node]) dfsTopSort(edge.u, topSort, state);
    //         state[node] = 2;
    //         topSort.push_back(node);
    //     }

    // public:
    //     // O(V_component log V_component + E_component)
    //     // @returns A topological sorted order of nodes
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
    //     // @note Graph cannot be bidirectional or cyclical
    //     std::vector<Node> getTopSort(Node root = 0, bool doKahn = false) const {
    //         assert(isDirected && "Can't get a topological sort for a bidirectional graph");
    //         assert((root == 0 || containsNode(root)) && "Node index out of range");

    //         if (doKahn) {
    //             std::vector<Node> topSort = Kahn(root);
    //             assert(topSort.size() <= V && "Sanity check!");
    //             assert(topSort.size() == V && "Graph contains cycle");
    //             return topSort;
    //         }

    //         else {
    //             std::vector<Node> state = std::vector<Node>(V + 1, 0);

    //             std::vector<Node> topSort;
    //             for (Node node: getComponentNodes(root)) dfsTopSort(node, topSort, state);
    //             reverse(topSort.begin(), topSort.end());
    //             return topSort;
    //         }
    //     }

    // private:
    //     T addFlow(Node source, Node sink, const std::vector<Node> &prev, std::map<std::pair<Node, Node>, T> &remainingCap) const {
    //         // find max flow of augmenting path
    //         T newFlow = std::numeric_limits<T>::max();
    //         for (Node u = sink; u != source; u = prev[u]) {
    //             newFlow = std::min(newFlow, remainingCap[{prev[u], u}]);
    //         }

    //         // update path
    //         for (Node u = sink; u != source; u = prev[u]) {
    //             remainingCap[{prev[u], u}] -= newFlow;
    //             remainingCap[{u, prev[u]}] += newFlow;
    //         }
    //         return newFlow;
    //     }

    // public:
    //     // O(V E^2 log E) or O(V^2 E log E)
    //     // @returns The value of max flow, using (fast) Dinic's or (slow) Edmonds-Karp
    //     T maxFlow(Node source, Node sink, bool doSlow = false) const {
    //         assert(containsNode(source) && containsNode(sink) && "Node index out of range");
    //         if (source == sink) return std::numeric_limits<T>::max();
    //         std::map<std::pair<Node, Node>, T> remainingCap;
    //         for (std::pair<std::pair<Node, Node>, T> edge: edges) {
    //             remainingCap[edge.u] += edge.v;
    //         }

    //         T totalFlow = T(0);

    //         while (true) {
    //             std::vector<int> prev = std::vector<int>(V + 1, -1);
    //             std::vector<int> levels = std::vector<int>(V + 1, -1);
    //             std::queue<Node> q;
    //             levels[source] = 0;
    //             q.push(source);

    //             while (!q.empty()) {
    //                 Node u = q.front();
    //                 q.pop();
    //                 for (Node type = 0; type <= 1; ++type) {
    //                     for (std::pair<Node, T> edge: (type ? inEdges : outEdges)[u]) {
    //                         Node v = edge.u;
    //                         if (prev[v] == -1 && remainingCap[{u, v}] > 0) {
    //                             prev[v] = u;
    //                             levels[v] = levels[u] + 1;
    //                             q.push(v);
    //                         }
    //                     }
    //                 }
    //             }

    //             // no augmenting path, so max flow has been found
    //             if (prev[sink] == -1) return totalFlow;

    //             if (doSlow) {
    //                 // Edmonds-Karp (implmentation of Ford-Fulkerson) - add any augmenting path
    //                 totalFlow += addFlow(source, sink, prev, remainingCap);
    //             }
    //             else {
    //                 // Dinic's - find best augmenting path in residual "level" graph
    //                 Node newFlow;
    //                 do {
    //                     std::stack<Node> s;
    //                     s.push(source);
    //                     while (!s.empty() && s.top() != sink) {
    //                         Node u = s.top();
    //                         s.pop();

    //                         for (Node type = 0; type <= 1; ++type) {
    //                             for (std::pair<Node, T> edge: (type ? inEdges : outEdges)[u]) {
    //                                 Node v = edge.u;
    //                                 if (levels[v] == levels[u] + 1 && remainingCap[{u, v}] > 0) {
    //                                     prev[v] = u;
    //                                     q.push(v);
    //                                 }
    //                             }
    //                         }
    //                     }

    //                     newFlow = addFlow(source, sink, prev, remainingCap);
    //                     totalFlow += newFlow;
    //                 } while (newFlow != 0);
    //             }
    //         }
    //     }

    //     // TODO: is bipartite graph?
    //     // TODO: edge matching
    //     // TODO: see: https://usaco.guide/adv/

    //     /************************************************
    //      *                    PROPERTIES                *
    //      ***********************************************/

    //     // O(1)
    //     // @returns Whether the graph is directed
    //     inline const bool isDirectedGraph() const {
    //         return isDirected;
    //     }

    //     // O(1)
    //     // @returns Whether the graph is weighted
    //     inline const bool isWeightedGraph() const {
    //         return isWeighted;
    //     }

    //     // O(V_component log V_component + V_component log E + E_component)
    //     // @returns Whether there exists a edge which connects a node to itself
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
    //     bool hasSelfEdges(Node root = 0) const {
    //         assert((root == 0 || containsNode(root)) && "Node index out of range");
    //         for (Node node: getComponentNodes(root)) {
    //             if (containsEdge(node, node)) return true;
    //         }
    //         return false;
    //     }

    //     // O(V_component log V_component + E_component)
    //     // @returns Whether there exists "redundant" edges
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
    //     bool hasDoubleEdges(Node root = 0) const {
    //         for (Node node: getComponentNodes(root)) {
    //             for (Node type = 0; type <= 1; ++type) {
    //                 const std::multiset<std::pair<Node, T>> *currEdges = &(type ? inEdges[node] : outEdges[node]);
    //                 if (currEdges->empty()) continue;
    //                 for (auto it1 = currEdges->begin(), it2 = ++currEdges->begin(); it2 != currEdges->end(); ++it1, ++it2) {
    //                     if (it1->first == it2->first) return true;
    //                 }
    //             }
    //         }
            
    //         return false;
    //     }

    //     // O(V_component log V_component + E_component)
    //     // @returns Whether the graph has a cycle
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
    //     bool hasCycle(Node root = 0) const {
    //         assert((root == 0 || containsNode(root)) && "Node index out of range");
    //         return (Node)Kahn(root).size() < V;
    //     }

    //     // O(V_component log V_component + V_component log E + E_component)
    //     // @returns Whether the graph is a simple graph
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
    //     bool isSimpleGraph(Node root = 0) const {
    //         assert((root == 0 || containsNode(root)) && "Node index out of range");
    //         return !hasSelfEdges(root) && !hasDoubleEdges(root);
    //     }

    //     // O(V_component log V_component + E_component)
    //     // @returns Whether the graph is a directed acrylic graph
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
    //     bool isDAG(Node root = 0) const {
    //         assert((root == 0 || containsNode(root)) && "Node index out of range");
    //         return isDirected && !hasCycle(root);
    //     }

    //     // O(V_component log V_component + E_component)
    //     // @returns Whether the graph is a tree
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole graph
    //     bool isTree(Node root = 0) const {
    //         assert((root == 0 || containsNode(root)) && "Node index out of range");
    //         if (hasCycle(root) || hasDoubleEdges(root)) return false;

    //         std::vector<Node> component = getComponentNodes(root);
    //         Node numEdges = 0, numLeaves = 0, numRoots = 0;
    //         for (Node node: component) {
    //             for (std::pair<Node, T> edge: outEdges[node]) {
    //                 if (!isDirected && edge.u > node) break;
    //                 else ++numEdges;
    //             }
    //             numLeaves += outEdges[node].empty();
    //             numRoots += inEdges[node].empty();
    //         }

    //         if (isDirected && numLeaves > 1 && numRoots > 1) return false;
    //         return numEdges == component.size() - 1; // E == V - 1
    //     }

    //     // O(V + E) Determines whether the graph is a forest
    //     bool isForest() const {
    //         for (std::vector<Node>& component: getComponentsNodes()) {
    //             if (!isTree(component[0])) return false;
    //         }
    //         return true;
    //     }

    //     /************************************************
    //      *       ALGORITHMS (COMPLEXITLY CLASS NP)      *
    //      ************************************************/

    //     // TODO: hamiltonian path / tour
    //     // TODO: euler path / tour
    //     // TODO: covering set
    //     // TODO: optimal colouring
    //     // TODO: planar embedding?
    };
};

#endif
