#ifndef GRAPH_H
#define GRAPH_H
#include <algorithm>
#include <cassert>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <stack>
#include <type_traits>
#include <unordered_set>
#include <vector>

namespace DS {
// I made the decision to use a class template argument for Graph, rather than
// a default value
// this makes reading input more convinient, and also saves some memory (this struct takes 1 byte)
struct UnitEdgeWeight {
    UnitEdgeWeight() {
    }

    template<typename T>
    UnitEdgeWeight(T _) {
    }

    template<
        typename PathWeight,
        std::enable_if_t<
            std::is_integral_v<PathWeight> || std::is_floating_point_v<PathWeight>, bool> = true>
    friend PathWeight operator+(const PathWeight &w1, const UnitEdgeWeight &w2) {
        return w1 + PathWeight(1);
    }

    template<
        typename PathWeight,
        std::enable_if_t<
            std::is_integral_v<PathWeight> || std::is_floating_point_v<PathWeight>, bool> = true>
    friend PathWeight operator+=(PathWeight &w1, const UnitEdgeWeight &w2) {
        return ++w1;
    }

    // lol
    friend bool operator==(const UnitEdgeWeight &a, const UnitEdgeWeight &b) {
        return true;
    }

    // required for Kruskal's MST sort edges by weights
    friend bool operator<(const UnitEdgeWeight &a, const UnitEdgeWeight &b) {
        return false;
    }

    // theres a way to determine if weights are negative
    friend bool operator<(const UnitEdgeWeight &a, int b) {
        assert(b == 0);
        return false;
    }

    // lol
    friend std::istream &operator>>(std::istream &in, const UnitEdgeWeight &w) {
        return in;
    }

    friend std::ostream &operator<<(std::ostream &out, const UnitEdgeWeight &w) {
        out << 1;
        return out;
    }
};

// @TODO support conversion from directed to undirected and vice versa
// @TODO support mapping of edges e.g. (u, v, w) -> (u, v, 1), so conversion from weighted to
// unweighted
// @TODO breaks for maxV = 1e5 (too slow? too much memory?)

// Graph
// @note A insufficiently defined EdgeWeight might still produce correct results sometimes
template<
    size_t maxV,
    // weight data held for each edge
    typename EdgeWeight_,
    // accumulation of weight data over many edges
    typename PathWeight_,
    // Whether edges are directed or bidirectional
    bool isDirectedEdges                                                   = false,
    std::enable_if_t<std::is_arithmetic_v<PathWeight_> && maxV != 0, bool> = true>
class Graph {
public:
    using EdgeWeight = EdgeWeight_;
    using PathWeight = PathWeight_;
    // id of node is an unsigned integer
    using Node                             = size_t;
    using MyGraph                          = Graph<maxV, EdgeWeight, PathWeight, isDirectedEdges>;
    constexpr static PathWeight MAX_WEIGHT = std::numeric_limits<PathWeight>::max();

    // Edges are represented internaly as (u, v) with weight w
    // Undirected edges are either (u, v) or (v, u), which are interfaced with as such
    // even though they compare equal
    struct Edge {
        Node u = 0, v = 0;
        EdgeWeight w;

        Edge() {
        }

        Edge(Node u_, Node v_, EdgeWeight w_) : u(u_), v(v_), w(w_) {
            assert(u != v && "Self cycles not supported yet");
            assert(u < maxV && v < maxV && "Node index out of range");
        }

        // #if __cplusplus >= 202002L
        // template<std::enable_if_t<isWeighted::value, bool> = false>
        // #endif
        Edge(Node u_, Node v_) : u(u_), v(v_), w(EdgeWeight()) {
            assert(u_ != v_ && "Self cycles not supported yet");
            assert(u < maxV && v < maxV && "Node index out of range");
        }

        // If this edge is (u, v), returns edge (v, u)
        // @note Undirected edges can be flipped
        Edge flip() const {
            return Edge(v, u, w);
        }

        Edge directFrom(Node node) {
            assert((node == v || node == u) && "Must provide an endpoint of this edge");
            if (node == u)
                return Edge(*this);
            else
                return flip();
        }

        Edge directTo(Node node) {
            assert((node == v || node == u) && "Must provide an endpoint of this edge");
            if (node == v)
                return Edge(*this);
            else
                return flip();
        }

        // Given u, return v
        // Given v, return u
        Node otherSide(Node node) const {
            assert((node == v || node == u) && "Must provide an endpoint of this edge");
            return u ^ v ^ node;
            // if (node == u) return v;
            // else if (node == v) return u;
            // else assert(false && "Must provide an endpoint of this edge");
        }

        // friend bool operator<(const Edge &a, const Edge &b) {
        //     if (a.u != b.u) return a.u < b.u;
        //     else if (a.v != b.v) return a.v < b.v;
        //     else return a.w < b.w;
        // }

        friend bool operator==(const Edge &a, const Edge &b) {
            return a.w == b.w &&
                   ((a.u == b.u && a.v == b.v) || (isUndirected && a.u == b.v && a.v == b.u));
        }

        friend bool operator!=(const Edge &a, const Edge &b) {
            return !(a == b);
        }

        friend std::istream &operator>>(std::istream &in, Edge &e) {
            in >> e.u >> e.v >> e.w;
            return in;
        }

        friend std::ostream &operator<<(std::ostream &out, Edge &e) {
            if constexpr (isDirectedEdges) {
                out << e.u << "->" << e.v << " (w = " << e.w << ')';
            } else {
                out << e.u << "--" << e.v << " (w = " << e.w << ')';
            }
            return out;
        }
    };

    using Path = std::vector<Edge>;

private:
    // comparison ignoring weights - use this in containers
    struct EdgeComp {
        bool operator()(const Edge &a, const Edge &b) const {
            Node au = a.u, av = a.v;
            if (isUndirected() && au > av) std::swap(au, av);
            Node bu = b.u, bv = b.v;
            if (isUndirected() && bu > bv) std::swap(bu, bv);

            if (au != bu) return au < bu;
            return av < bv;
        }
    };

private:
    std::vector<Node> nodes; // not garunteed sorted!
    std::array<bool, maxV> validNode;

    // @TODO is it reasonable to support multiset?
    // @TODO `edgesIn` and `edgesOut` store copies of each edge - not memory efficient?
    // A set of all edges
    std::set<Edge, EdgeComp> edges;
    // `edgesIn[node]` is the set of edges leading into `node`
    std::array<std::set<Edge, EdgeComp>, maxV> edgesIn;
    // `edgesOut[node]` is the set of edges leading out of `node`
    std::array<std::set<Edge, EdgeComp>, maxV> edgesOut;

public:
    /************************************************
     *                 INITIALISATION               *
     ************************************************/

    // O(V)
    // Initialises an empty graph with no nodes
    Graph() {
        validNode.fill(false);
    }

    // O(V + E log E)
    // Initalises an graph from vector of edges, assuming 1 indexed
    Graph(size_t V, const std::vector<Edge> &_edges = {}) {
        assert(V < maxV && "Not enough capacity");
        validNode.fill(false);
        for (Node node = 1; node <= V; ++node)
            pushNode(node);
        for (Edge edge : _edges)
            insertEdge(edge);
    }

    // O(V + E log E)
    // Initialses a G with node set `_nodes` and edge set `_edges`
    Graph(const std::vector<Node> &_nodes, const std::vector<Edge> &_edges = {}) {
        validNode.fill(false);
        for (Node node : _nodes)
            pushNode(node);
        for (Edge edge : _edges)
            insertEdge(edge);
    }

    // O(V^2)
    // Initialises a weighted graph from adjacency matrix of edge weights
    // @note Assuming nodes 1 indexed

    Graph(const std::vector<std::vector<EdgeWeight>> &_adj) {
        size_t V = _adj.size();
        assert(V < maxV && "Not enough capacity");
        validNode.fill(false);
        for (Node node = 1; node <= V; ++node)
            pushNode(node);
        for (Node u = 1; u <= V; ++u) {
            assert(_adj[u - 1].size() == V && "Not square matrix");
            for (Node v = 1; v <= V; ++v) {
                if (_adj[u - 1][v - 1] != MAX_WEIGHT) {
                    insertEdge(Edge(u, v, _adj[u - 1][v - 1]));
                }
            }
        }
    }

    // Graph(const Matrix &_adj) {
    // @TODO
    // }

    // template<typename Weight = PathWeight>
    // Graph(const Tree<maxV, Weight> &_tree) {
    //     // @TODO
    // }

    // O(|nodes| + num_induced_edges log E) = O(|nodes| + E log E)
    // Initialise this graph as the graph induced from `g` by `nodes`
    // @note `V` does not change, so there may be disconnected nodes
    Graph(const MyGraph &_g, const std::vector<Node> &_nodes) {
        validNode.fill(false);
        for (Node node : _nodes)
            pushNode(node);
        for (Node node : _nodes) {
            for (Edge incident : _g.getEdges(node)) {
                if (containsNode(incident.otherSide(node)) && !containsEdge(incident))
                    insertEdge(incident);
            }
        }
    }

    /************************************************
     *                    DISPLAY                   *
     ************************************************/

    // O(V + E)
    // Displays the graph, showing the outwards edge connections of each edge in lexographic order
    // @param `out` The string representation of the graph is piped to this output stream
    friend std::ostream &operator<<(std::ostream &out, const MyGraph &graph) {
        for (Node u : graph.nodes) {
            out << u << ':';
            for (Edge edge : graph.edgesOut[u]) {
                // on an undirected graph, don't print edges twice
                if (isUndirected() && u > edge.u) continue;

                out << ' ' << edge.v;
                if (isWeighted()) out << " (w = " << edge.w << ')';
            }
            out << '\n';
        }
        return out;
    }

    // O(maxV^2)
    // @TODO alignment ahh
    void displayAdj(std::ostream &out) {
        for (Node u = 0; u < maxV; ++u) {
            for (Node v = 0; v < maxV; ++v) {
                if (u == v || !containsEdgeUnweighted(Edge(u, v))) {
                    out << "- ";
                } else {
                    out << getEdge(u, v, EdgeWeight()) << ' ';
                }
            }
            out << '\n';
        }
    }

    /************************************************
     *                  PROPERTIES                  *
     ************************************************/

    // @note No `size()` because ambiguous, use the following functions instead

    // @returns `V`, the number of nodes
    size_t V() const {
        return nodes.size();
    }
    // @returns `V`, the number of nodes
    size_t N() const {
        return nodes.size();
    }

    // O(1)
    // @returns `E`, the number of edges
    size_t E() const {
        return edges.size();
    }
    // O(1)
    // @returns `E`, the number of edges
    size_t M() const {
        return edges.size();
    }

    // @returns the imutable set of all ndoes
    const std::vector<Node> &getNodes() const {
        return nodes;
    }

    // O(E)
    // @returns the imutable set of all edges
    // @note Edges sorted by endpoints
    const std::vector<Edge> getEdges() const {
        return std::vector<Edge>(edges.begin(), edges.end());
    }

    /************************************************
     *                INCIDENT DATA                 *
     ************************************************/

    // O(E)
    // @returns The imutable collection of unique edges incident to this node
    const std::vector<Edge> getEdges(Node node) const {
        assert(containsNode(node) && "Node index out of range");
        std::vector<Edge> out(edgesIn[node].size() + edgesOut[node].size());
        auto it = std::set_union(
            edgesIn[node].begin(), edgesIn[node].end(), edgesOut[node].begin(),
            edgesOut[node].end(), out.begin(), EdgeComp()
        );
        out.resize(it - out.begin());
        return out;
    }

    // O(E)
    // Finds the edges that travel into `node`
    // @returns The imutable collection of incoming edges incident to this node
    const std::vector<Edge> getEdgesIn(Node node) const {
        assert(containsNode(node) && "Node index out of range");
        return std::vector<Edge>(edgesIn[node].begin(), edgesIn[node].end());
    }

    // O(E)
    // Finds the edges that travel out of `node`
    // @returns The imutable collection of outgoing edges incident to this node
    const std::vector<Edge> getEdgesOut(Node node) const {
        assert(containsNode(node) && "Node index out of range");
        return std::vector<Edge>(edgesOut[node].begin(), edgesOut[node].end());
    }

    // O(E)
    // @returns The imutable colelction of unique neighbours, in sorted order
    const std::vector<Node> getNeighbours(Node node) const {
        assert(containsNode(node) && "Node index out of range");
        std::vector<Node> out;
        for (const Edge e : getEdges(node)) {
            out.push_back(e.otherSide(node));
        }
        return out;
    }

    // O(E)
    // @returns The imutable collection of neighbours reachable from `node` via 1 edge
    const std::vector<Node> getNeighboursIn(Node node) const {
        assert(containsNode(node) && "Node index out of range");
        std::vector<Node> out;
        for (const Edge e : edgesIn[node]) {
            out.push_back(e.otherSide(node));
        }
        return out;
    }

    // O(E)
    // @returns The imutable collection of neighbours which reach `node` via 1 edge
    const std::vector<Node> getNeighboursOut(Node node) const {
        assert(containsNode(node) && "Node index out of range");
        std::vector<Node> out;
        for (const Edge e : getEdgesOut(node)) {
            out.push_back(e.otherSide(node));
        }
        return out;
    }

    // O(E)
    // @returns The imutable collection of nodes which can be the start of a topsort
    const std::vector<Node> getSources() const {
        assert(isDirected() && "cannot find sources of undirected graph");
        std::vector<Node> out;
        for (const Node node : nodes) {
            if (degreeIn(node) == 0) {
                out.push_back(node);
            }
        }
        return out;
    }

    // O(E)
    // @returns The imutable collection of nodes which can be the end of a topsort
    const std::vector<Node> getSinks() const {
        assert(isDirected() && "cannot find sources of undirected graph");
        std::vector<Node> out;
        for (const Node node : nodes) {
            if (degreeOut(node) == 0) {
                out.push_back(node);
            }
        }
        return out;
    }

    // O(1)
    // @returns The total degree of `node` (combination in in degree and out degree)
    // @note In an undirected graph, edge uv and vu will add 2 to the degree count
    const size_t degree(Node node) const {
        assert(containsNode(node) && "Node index out of range");
        if (isDirected()) {
            assert(edgesIn[node].size() == edgesOut[node].size() && "Sanity check!");
            return edgesIn[node].size();
        } else {
            return edgesIn[node].size() + edgesOut[node].size();
        }
    }

    // O(1)
    // @returns The in degree of `node`
    const size_t degreeIn(Node node) const {
        assert(containsNode(node) && "Node index out of range");
        return edgesIn[node].size();
    }

    // O(1)
    // @returns The out degree of `node`
    const size_t degreeOut(Node node) const {
        assert(containsNode(node) && "Node index out of range");
        return edgesOut[node].size();
    }

    /************************************************
     *               GRAPH MODIFICATIONS            *
     ************************************************/

    // O(1)
    // Adds a new node with no edges
    // @returns The index of this node
    // @TODO should I have a variation where silently passes when node exists?
    size_t pushNode(Node node) {
        assert(node < maxV + 1 && "Node index out of range");
        assert(!containsNode(node) && "node already exists");
        nodes.push_back(node);
        validNode[node] = true;
        return V();
    }

    // O(log E)
    // Add the specified edge to the graph, or error if this will become a double edge
    // @note Edge weight is ignored if the graph is unweighted
    // @TODO should I have a variation where silently passes when edge exists?
    void insertEdge(Edge e) {
        // @note this check is duplicated, but included for clear error messages
        assert(containsNode(e.u) && containsNode(e.v) && "Node index out of range");
        assert(
            !containsEdgeUnweighted(e) && "Edge already in graph, is this an unintended operation?"
        );

        edges.insert(e);
        // @note this is preventing support of multiedges
        // when duplicaate edge is added, how do you find the most recently added edge
        edgesIn[e.v].insert(e);
        edgesOut[e.u].insert(e);
        if (isUndirected()) {
            edgesIn[e.u].insert(e);
            edgesOut[e.v].insert(e);
        }
    }

    // O(log E)
    // If the specified edge (with any edge weight) is present, remove it, otherwise silently do
    // nothing
    void eraseEdge(Edge e) {
        if (!containsEdgeUnweighted(e)) return;

        edges.erase(e); // @TODO if edge weight is different is this logic error?
        edgesIn[e.v].erase(e);
        edgesOut[e.u].erase(e);
        if (isUndirected()) {
            edgesIn[e.u].erase(e);
            edgesOut[e.v].erase(e);
        }
    }

    /************************************************
     *                    CONTAINS                  *
     ************************************************/

    // O(1)
    // @returns whether `node` is within the acceptable bounds
    bool isNode(Node node) const {
        // @note dont need to check 0 <= node, because Node=size_t
        return node < maxV;
    }

    // O(1)
    // @returns whether `node` is in the vertex set of this graph
    bool containsNode(Node node) const {
        // @note dont need to check 0 <= node, because Node=size_t
        return node < maxV && validNode[node];
    }

    // O(log E)
    // @returns Whether the specified edge (with any edge weight) is contained in the graph
    bool containsEdge(Edge e) const {
        if (!containsNode(e.u) || !containsNode(e.v)) return false;
        auto edgeIt = edges.find(e);
        return edgeIt != edges.end() && edgeIt->w == e.w;
    }

    // O(log E)
    // @returns Whether the specified edge is contained in the graph
    bool containsEdgeUnweighted(Edge e) const {
        if (!containsNode(e.u) || !containsNode(e.v)) return false;
        return edges.count(e);
    }

    // O(log E)
    // @returns the weight of the edge between nodes u and v, if it exists
    // @returns default otherwise
    EdgeWeight getEdge(Node u, Node v, EdgeWeight defaultWeight) {
        if (!containsNode(u) || !containsNode(v)) return defaultWeight;
        auto edgeIt = edges.find(Edge(u, v));
        if (edgeIt == edges.end())
            return defaultWeight;
        else
            return edgeIt->w;
    }

    /************************************************
     *                   COMPONENTS                 *
     ************************************************/

    // O(V_component log V_component + E_component log E_component)
    // @returns All nodes in the same component as `root`, in sorted order
    // @note Big constant factor due to `std::unordered_set` to store seen nodes
    MyGraph getComponent(Node at) const {
        assert(containsNode(at) && "Node index out of range");
        std::unordered_set<Node> seen;
        std::queue<Node> q;
        seen.insert(at);
        q.push(at);

        while (!q.empty()) {
            Node node = q.front();
            q.pop();

            for (Node child : getNeighbours(node)) {
                if (!seen.count(child)) {
                    seen.insert(child);
                    q.push(child);
                }
            };
        }
        return MyGraph(*this, std::vector<Node>(seen.begin(), seen.end()));
    }

    // O(maxV + E log E)
    // @returns All nodes grouped by their component, in arbitary order
    std::vector<MyGraph> getComponents() const {
        std::vector<MyGraph> output;
        std::array<bool, maxV> seen;
        seen.fill(false);

        for (Node node : nodes) {
            if (!seen[node]) {
                seen[node]                  = true;
                std::vector<Node> component = {node};
                for (size_t i = 0; i < component.size(); ++i) {
                    for (Node child : getNeighbours(component[i])) {
                        if (!seen[child]) {
                            seen[child] = true;
                            component.push_back(child);
                        }
                    };
                }
                output.push_back(MyGraph(*this, component));
                // this adds a log factor
                // sort(output.begin(), output.end());
            }
        }
        return output;
    }

    /************************************************
     *   TRIVIAL ALGORITHMS (COMPLEXITLY CLASS P)   *
     ***********************************************/

    // O(V + E log E)
    // @returns The union of 2 graphs
    // @note Does not cause edge doubling
    friend MyGraph operator+(const MyGraph &a, const MyGraph &b) {
        MyGraph out = a; // make a copy
        for (Edge edge : b.getEdges()) {
            if (!out.containsNode(edge.u)) out.pushNode(edge.u);
            if (!out.containsNode(edge.v)) out.pushNode(edge.v);
            if (!out.containsEdge(edge)) out.insertEdge(edge);
        }
        return out;
    }

    // O(V + E)
    // @returns whether the 2 graphs are identical in edge (including weights!) and vertex sets
    // @note has short circuiting
    friend bool operator==(const MyGraph &a, const MyGraph &b) {
        return (
            a.E() == b.E() && a.V() == b.V() &&
            std::equal(std::begin(a.validNode), std::end(a.validNode), std::begin(b.validNode)) &&
            std::equal(a.edges.begin(), a.edges.end(), b.edges.begin())
        );
    }

    friend bool operator!=(const MyGraph &a, const MyGraph &b) {
        return !(a == b);
    }

    // O(V + E log E)
    MyGraph flip() const {
        assert(isDirected() && "cannot flip undirected graph");
        MyGraph out;
        for (Node node : getNodes())
            out.pushNode(node);
        for (Edge edge : getEdges())
            out.insertEdge(edge.flip());
        return out;
    }

    // O(V + E)
    void clear() {
        edges.clear();
        for (Node node : nodes) {
            edgesIn[node].clear();
            edgesOut[node].clear();
        }
        validNode.fill(false);
        nodes.clear();
    }

    /************************************************
     * NON TRIVIAL ALGORITHMS (COMPLEXITLY CLASS P) *
     ***********************************************/

    // leave as placeholder, not worthwhile implementing this
    // https://en.wikipedia.org/wiki/Shortest_path_problem#Directed_graphs_with_arbitrary_weights_with_negative_cycles
    std::pair<std::array<PathWeight, maxV>, std::array<Edge, maxV>>
    NegativeCycleShortestPath(Node node) const {
        assert("negative weight cycle, cannot SSSP");
        throw std::exception();
    }

    // O(maxV + E log V), or O(maxV E) if negative edge weights
    // Standard Dijkstra algorithm
    // shortest path from `node` to all other nodes for non-negative weighted graph (either directed
    // or undirected)
    // @returns { cost, pred }
    // @note If the graph uses `UnitWeight`, this becomes equivalent to a bfs (but slightly slower)
    // @note If negative edge weights encountered, throws
    // @note Will check and error if negative weight cycle
    std::pair<std::array<PathWeight, maxV>, std::array<Edge, maxV>> Dijkstra(Node node) const {
        std::array<PathWeight, maxV> cost;
        cost.fill(MAX_WEIGHT);
        cost[node] = PathWeight(0);
        std::array<Edge, maxV> pred; // predecessor

        if (containsNode(node)) {
            std::priority_queue<
                std::pair<PathWeight, Node>, std::vector<std::pair<PathWeight, Node>>,
                std::greater<std::pair<PathWeight, Node>>>
                pq;
            pq.push({PathWeight(0), node});

            while (!pq.empty()) {
                Node currNode       = pq.top().second;
                PathWeight currCost = pq.top().first;
                pq.pop();
                if (cost[currNode] < currCost) continue;

                for (Edge incident : edgesOut[currNode]) {
                    assert(!(incident.w < 0) && "negative weights, use BellmanFord algo instead");

                    PathWeight newCost = currCost + incident.w;
                    Node newNode       = incident.otherSide(currNode);
                    if (newCost < cost[newNode]) {
                        cost[newNode] = newCost;
                        pred[newNode] = incident;
                        pq.push({newCost, newNode});
                    }
                }
            }
        }

        return {cost, pred};
    }

    // O(maxV E)
    // Standard Bellman-Ford algorithm
    // shortest path from `node` to all other nodes for any weighted directed graph
    // @returns { cost, pred }
    // @note Will check and error if negative weight cycle
    std::pair<std::array<PathWeight, maxV>, std::array<Edge, maxV>> BellmanFord(Node node) const {
        std::array<PathWeight, maxV> cost;
        cost.fill(MAX_WEIGHT);
        cost[node] = PathWeight(0);
        std::array<Edge, maxV> pred; // predecessor

        if (containsNode(node)) {
            for (size_t rep = 0; rep < V() - 1; ++rep) {
                for (Edge edge : getEdges()) {
                    if (cost[edge.u] != MAX_WEIGHT && cost[edge.u] + edge.w < cost[edge.v] &&
                        pred[edge.u] != edge) {
                        cost[edge.v] = cost[edge.u] + edge.w;
                        pred[edge.v] = edge;
                    }
                    if (isUndirected() && cost[edge.v] != MAX_WEIGHT &&
                        cost[edge.v] + edge.w < cost[edge.u] && pred[edge.v] != edge) {
                        cost[edge.u] = cost[edge.v] + edge.w;
                        pred[edge.u] = edge;
                    }
                }
            }

            // last rep has finished, but the minimum distance can still be reduced
            // indicating negative cycles
            for (Edge edge : getEdges()) {
                if (cost[edge.u] != MAX_WEIGHT && cost[edge.v] != MAX_WEIGHT) {
                    // any more simplification of this expression assume properties of PathWeight
                    if (isDirected()) {
                        if (cost[edge.u] + edge.w < cost[edge.v])
                            return NegativeCycleShortestPath(node);
                    } else if (edge != pred[edge.u]) {
                        if (cost[edge.u] + edge.w < cost[edge.v])
                            return NegativeCycleShortestPath(node);
                    } else if (edge != pred[edge.v]) {
                        if (cost[edge.v] + edge.w < cost[edge.u])
                            return NegativeCycleShortestPath(node);
                    }
                }
            }
        }

        return {cost, pred};
    }

    // O(maxV + E log V), or O(maxV E) if negative edge weights
    // shortest path from `node` to all other nodes (using either Dijkstra of BellmanFord)
    // @returns { cost, pred }
    // @note If the graph uses `UnitWeight`, this becomes equivalent to a bfs (but slightly slower)
    // @note Will check and error if negative weight cycle
    std::pair<std::array<PathWeight, maxV>, std::array<Edge, maxV>> SSSP(Node node) const {
        if (!containsNode(node)) {
            std::array<PathWeight, maxV> cost;
            cost.fill(MAX_WEIGHT);
            cost[node] = PathWeight(0);
            std::array<Edge, maxV> pred; // predecessor

            return {cost, pred};
        }

        for (const Edge edge : edges) {
            if (edge.w < 0) return BellmanFord(node);
        }
        return Dijkstra(node);
    }

    // O(maxV + E log V), or O(maxV E) if negative edge weights
    // @returns Whether there exists a path from u to v
    // @note Component containing `u` and `v` cannot have negative weight cycles
    bool isReachable(Node u, Node v) {
        return containsNode(u) && containsNode(v) && SSSP(u).first[v] != MAX_WEIGHT;
    }

    // O(maxV + E log V), or O(maxV E) if negative edge weights
    // @returns Distance of the path
    // @note Component containing `u` and `v` cannot have negative weight cycles
    // @note return MAX_WEIGHT when nodes are unreachable/out of range
    PathWeight shortestDist(Node u, Node v) const {
        if (isNode(u) && u == v)
            return PathWeight(0);
        else if (!containsNode(u) || !containsNode(v))
            return MAX_WEIGHT;
        else
            return SSSP(u).first[v];
    }

    // O(maxV + E log V), or O(maxV E) if negative edge weights
    // @returns Distance of the shortest path from `u` to `v`, and a vector of nodes representing
    // the path
    // @note Component containing `u` and `v` cannot have negative weight cycles
    Path shortestPath(Node u, Node v) const {
        assert(containsNode(u) && containsNode(v) && "Node index out of range");

        // @note Important to start from `u` and reverse path later, because digraphs
        std::pair<std::array<PathWeight, maxV>, std::array<Edge, maxV>> res = SSSP(u);
        std::array<PathWeight, maxV> &cost                                  = res.first;
        std::array<Edge, maxV> &pred                                        = res.second;

        Path path;
        if (cost[v] != MAX_WEIGHT) {
            for (Node node = v; node != u; node = pred[node].otherSide(node)) {
                path.push_back(pred[node].directTo(node));
            }
        }
        reverse(path.begin(), path.end());
        return path;
    }

    // O(maxV + E log V), but O(maxV + VE) if negative edge weights
    // @returns The eccentricity - the greatest distance between `root` and any other node in the
    // same component
    // @note Component containing `root` cannot have negative weight cycles
    PathWeight eccentricity(Node root) const {
        assert(containsNode(root) && "Node index out of range");
        std::array<PathWeight, maxV> cost = SSSP(root).first;
        PathWeight furthest               = 0;
        for (Node node = 0; node < maxV; ++node) {
            if (cost[node] != MAX_WEIGHT) furthest = std::max(furthest, cost[node]);
        }
        return furthest;
    }

    // O(V + E log V), but O(VE) if negative edge weights
    // Diameter is defined as the shortest simple path with maximal weight
    // @returns Distance of the path, and a vector of nodes representing the path
    // @note Component containing `root` cannot have negative weight cycles
    Path diameter(Node root) {
        assert(containsNode(root) && "Node index out of range");
        std::array<PathWeight, maxV> cost1 = SSSP(root).first;

        Node furthestNode1                 = root;
        for (const Node &node : nodes) {
            if (cost1[node] != MAX_WEIGHT && cost1[node] > cost1[furthestNode1]) {
                furthestNode1 = node;
            }
        }

        std::pair<std::array<PathWeight, maxV>, std::array<Edge, maxV>> res2 = SSSP(furthestNode1);
        std::array<PathWeight, maxV> &cost2                                  = res2.first;
        std::array<Edge, maxV> &pred2                                        = res2.second;

        Node furthestNode2                                                   = root;
        for (const Node &node : nodes) {
            if (cost2[node] != MAX_WEIGHT && cost2[node] > cost2[furthestNode2]) {
                furthestNode2 = node;
            }
        }

        Path path;
        for (Node node = furthestNode2; node != furthestNode1; node = pred2[node].otherSide(node)) {
            path.push_back(pred2[node].directFrom(node));
        }

        return path;
    }

    // leave as placeholder, not worthwhile implementing this
    // https://en.wikipedia.org/wiki/Shortest_path_problem#Directed_graph
    std::function<PathWeight(Node, Node)> NegativeCyclesAllPairsShortestPaths() const {
        assert(false && "not worth implementing");
        throw std::exception();
    }

    // O(maxV^3 + E)
    // Standard Floyd-Warshall algorithm - shortest path between every 2 pairs of nodes
    // @note all edge weights must be non negative
    // @note No path reconstruction
    std::function<PathWeight(Node, Node)> FloydWarshall() const {
        // assert()
        std::array<std::array<PathWeight, maxV>, maxV> dists;
        for (Node node = 0; node < maxV; ++node) {
            dists[node].fill(MAX_WEIGHT);
            dists[node][node] = PathWeight(0);
        }

        for (const Edge edge : edges) {
            dists[edge.u][edge.v] = std::min(dists[edge.u][edge.v], PathWeight(0) + edge.w);
            if (isUndirected()) {
                dists[edge.v][edge.u] = std::min(dists[edge.v][edge.u], PathWeight(0) + edge.w);
            }
        }

        for (Node k = 0; k < maxV; ++k) {
            for (Node u = 0; u < maxV; ++u) {
                for (Node v = 0; v < maxV; ++v) {
                    if (dists[u][k] != MAX_WEIGHT && dists[k][v] != MAX_WEIGHT) {
                        dists[u][v] = std::min(dists[u][v], dists[u][k] + dists[k][v]);
                    }
                }
            }
        }
        return [this, dists](Node u, Node v) {
            if (isNode(u) && u == v) return PathWeight(0);
            if (!containsNode(u) || !containsNode(v)) return MAX_WEIGHT;
            return dists[u][v];
        };
    }

    // O(maxV^3 + E) Finds the minimum distance between every pair of nodes
    // @returns A vector of vectors, with `dist[u][v]` being the minimim distance on the path from
    // `u` to `v`
    std::function<PathWeight(Node, Node)> allShortestDist() const {
        for (const Edge edge : edges) {
            if (edge.w < 0 && isUndirected()) return NegativeCyclesAllPairsShortestPaths();
        }

        return FloydWarshall();
    }

    // Amortised O(E log V)
    // Standard Kruskal's algorithm - finds the cost of the minimum spanning tree, using union find
    // (DSU)
    // @note If k components, return lowest cost to make k disjoint spanning trees
    // @TODO return a tree
    PathWeight KruskalsCost() const {
        std::array<Node, maxV> parent;
        std::iota(parent.begin(), parent.end(), 0);

        // Amortised O(log V)
        // @returns The parent node of the current disjoint set
        std::function<Node(Node)> getParent;
        getParent = [&parent, &getParent](Node node) {
            if (parent[node] == node) return node;
            return parent[node] = getParent(parent[node]);
        };

        std::vector<Edge> sortedEdges(edges.begin(), edges.end());
        sort(sortedEdges.begin(), sortedEdges.end(), [](Edge a, Edge b) { return a.w < b.w; });

        PathWeight cost = PathWeight(0);
        for (const Edge &edge : sortedEdges) {
            Node parentU = getParent(edge.u), parentV = getParent(edge.v);
            if (parentU != parentV) {
                cost += edge.w;
                parent[parentU] = parentV;
            }
        }
        return cost;
    }

    // O(V + E)
    // Standard Kahn's algorithm - topological sort
    std::vector<Node> Kahns() const {
        std::array<int, maxV> depths;
        depths.fill(0);

        std::vector<Node> q, topsort;
        for (const Node &node : nodes) {
            depths[node] = edgesIn[node].size();
            if (edgesIn[node].empty()) q.push_back(node);
        }

        while (!q.empty()) {
            Node u = q.back();
            q.pop_back();
            topsort.push_back(u);

            for (Node v : getNeighboursOut(u)) {
                if (--depths[v] == 0) q.push_back(v);
            }
        }
        assert(topsort.size() == V() && "Graph contains cycle");
        return topsort;
    }

    // O(V + E)
    // Standard DFS topsort
    std::vector<Node> dfsTopsort() const {
        std::vector<Node> topsort;
        std::array<int, maxV> state; // @TODO actually only needs to hold 3 states
        state.fill(0);

        std::function<void(int)> dfs;
        dfs = [&](Node node) {
            if (state[node] == 2) return;
            assert(state[node] != 1 && "Graph contains cycle");

            state[node] = 1;
            for (Node v : getNeighboursOut(node))
                dfs(v);
            state[node] = 2;
            topsort.push_back(node);
        };
        for (const Node &node : nodes) {
            if (state[node] == 0) dfs(node);
        }
        reverse(topsort.begin(), topsort.end());
        return topsort;
    }

    // O(V + E)
    // Find strongly connected components
    // https://cp-algorithms.com/graph/strongly-connected-components.html#implementation
    std::vector<std::vector<Node>> Kosaraju() const {
        std::vector<Node> order;
        std::array<bool, maxV> seen;
        seen.fill(false);
        std::function<void(int)> dfs1;
        dfs1 = [&](Node node) {
            if (seen[node]) return;
            seen[node] = true;
            for (Node v : getNeighboursOut(node))
                dfs1(v);
            order.push_back(node);
        };
        for (const Node &node : nodes)
            dfs1(node);
        reverse(order.begin(), order.end());

        seen.fill(false);
        std::vector<std::vector<Node>> components;
        std::function<void(int)> dfs2;
        dfs2 = [&](Node node) {
            seen[node] = true;
            components.back().push_back(node);
            for (Node v : getNeighboursIn(node)) {
                if (!seen[v]) dfs2(v);
            }
        };
        for (const Node &node : nodes) {
            if (!seen[node]) {
                components.push_back({});
                dfs2(node);
            }
        }

        return components;
    }

    // O(V + E)
    // Standard Tarjan algorithm - find strongly connected components
    // @returns Vector of all strongly connected components, where each component is a vector of
    // nodes
    std::vector<std::vector<Node>> Tarjan() const {
        assert(isDirected() && "SCC only applies to directed graphs");
        std::array<bool, maxV> seen;
        seen.fill(false);
        std::array<int, maxV> getIndex;
        getIndex.fill(-1);
        std::array<int, maxV> lowLink;
        lowLink.fill(0);
        std::vector<std::vector<Node>> components;
        std::vector<Node> s;
        int index = 0;

        std::function<Node(Node)> dfs;
        dfs = [&](Node u) {
            getIndex[u] = lowLink[u] = index++;
            s.push_back(u);
            seen[u] = true;

            for (Edge edge : edgesOut[u]) {
                Node v = edge.otherSide(u);
                if (getIndex[v] == -1) {
                    index      = dfs(v);
                    lowLink[u] = std::min(lowLink[u], lowLink[v]);
                } else if (seen[v]) {
                    lowLink[u] = std::min(lowLink[u], getIndex[v]);
                }
            }

            if (getIndex[u] == lowLink[u]) {
                components.push_back({});

                Node node;
                do {
                    node = s.back();
                    s.pop_back();
                    seen[node] = false;
                    components.back().push_back(node);
                } while (node != u);
            }

            return index;
        };
        for (const Node &node : nodes) {
            if (getIndex[node] == -1) dfs(node);
        }
        return components;
    }

    // O(V + E)
    // Find strongly connecsted components
    // @returns a vector of components, where a component is a vector of nodes
    // @TODO test whether Tarjans or dfs is faster
    // should I instead be returning `std::vector<MyGraph>` ?
    std::vector<std::vector<Node>> SSC() const {
        return Tarjan();
    }

    // O(maxV + E)
    // Standard DFS Bridge finding algorithm
    // Find all bridges
    // @note Breaks when we introduce multiple edges
    // @param `root` If specified, consider its connected component. Otherwise consider the whole
    // graph
    // @returns A vector of edges, where each edge is represented by ((u, v), w)
    std::vector<Edge> bridgesDFS() const {
        assert(isUndirected() && "Bridges are a property of only undirected graphs");
        std::array<bool, maxV> seen;
        seen.fill(false);
        std::array<int, maxV> minTime;
        minTime.fill(0);
        std::array<int, maxV> entryTime;
        entryTime.fill(0);
        std::vector<Edge> bridges;
        int t = 0;

        // huh this looks very similar to Tarjan's
        std::function<void(Node, Node)> dfs;
        dfs = [&](Node u, Node parent) {
            seen[u]    = true;
            minTime[u] = entryTime[u] = ++t;

            for (Edge edge : edgesOut[u]) {
                Node v = edge.otherSide(u);
                if (v == parent)
                    continue;
                else if (seen[v]) {
                    minTime[u] = std::min(minTime[u], entryTime[v]);
                } else {
                    dfs(v, u);
                    minTime[u] = std::min(minTime[u], minTime[v]);

                    // bridge iff minTime[childNode] > entryTime[node]
                    // @note When graph can support multiedges, check that the edge is not a double
                    // edge
                    if (minTime[v] > entryTime[u]) {
                        bridges.push_back(edge);
                    }
                }
            }
        };

        for (const Node &node : nodes) {
            if (!seen[node]) dfs(node, maxV);
        }
        return bridges;
    }

    // private:
    //     // O(V + E log E)
    //     // Greedy colouring
    //     void _greedyColouring(Node node, std::vector<int>& colours) const {
    //         if (colours[node] != -1) return;

    //         // calculate mex of neighbours colours
    //         std::set<int> usedColours;
    //         for (std::pair<Node, T> edge: getEdges(node)) usedColours.insert(colours[edge.u]);
    //         for (colours[node] = 0; usedColours.find(colours[node]) != usedColours.end();
    //         ++colours[node]);

    //         for (std::pair<Node, T> edge: getEdges(node)) _greedyColouring(edge.u, colours);
    //     }

    // public:
    //     // O(V_component log V_component + E log E)
    //     // @returns A possible greedy colouring of the graph - `0 <= colour[node] < max_colours`
    //     // @param `root` If specified, consider its connected component. Otherwise consider the
    //     whole graph
    //     // @note Ignore colour[0]
    //     std::vector<int> greedyColouring(Node root = 0) const {
    //         assert((root == 0 || containsNode(root)) && "Node index out of range");
    //         std::vector<int> colours(V + 1, -1);
    //         for (Node node: getComponentNodes(root)) _greedyColouring(node, colours);
    //         return colours;
    //     }

    // O(V E^2 log E) or O(V^2 E log E)
    // Edmonds-Karp max flow algorithm
    // @note Using `PathWeight` to represent flow, will this cause issues
    // @returns The value of max flow
    PathWeight EdmondsKarp(Node source, Node sink, bool doSlow = false) const {
        assert(containsNode(source) && containsNode(sink) && "Node index out of range");

        if (source == sink) return MAX_WEIGHT;

        std::map<std::pair<Node, Node>, PathWeight> remainingCap;
        std::function<PathWeight(std::array<Node, maxV>)> addFlow;
        addFlow = [&](std::array<Node, maxV> prev) {
            // find flow of augmenting path
            PathWeight newFlow = MAX_WEIGHT;
            for (Node u = sink; u != source; u = prev[u]) {
                newFlow = std::min(newFlow, remainingCap[{prev[u], u}]);
            }

            // update path
            for (Node u = sink; u != source; u = prev[u]) {
                remainingCap[{prev[u], u}] -= newFlow;
                remainingCap[{u, prev[u]}] += newFlow;
            }
            return newFlow;
        };

        for (Edge edge : edges) {
            remainingCap[{edge.u, edge.v}] += edge.w;
            if (isUndirected()) {
                remainingCap[{edge.v, edge.u}] += edge.w;
            }
        }

        PathWeight totalFlow = PathWeight(0);

        std::array<Node, maxV> prev;
        while (true) {
            prev.fill(maxV);
            prev[source] = source;
            std::queue<Node> q;
            q.push(source);

            while (!q.empty()) {
                Node u = q.front();
                q.pop();
                for (Edge edge : getEdgesOut(u)) {
                    Node v = edge.otherSide(u);
                    if (prev[v] == maxV && remainingCap[{u, v}] > 0) {
                        prev[v] = u;
                        q.push(v);
                    }
                }
            }

            // no augmenting path, so max flow has been found
            if (prev[sink] == maxV) return totalFlow;

            // add any augmenting path
            totalFlow += addFlow(prev);
        }
    }

    // O(V^2 E), or O(min(V^(2/3), sqrt E) E) if `!isWeightedGraph()`
    // Dinics max flow algorithm
    // @note Using `PathWeight` to represent flow, will this cause issues
    // @returns The value of max flow
    PathWeight Dinics(Node source, Node sink) {
        assert(containsNode(source) && containsNode(sink) && "Node index out of range");

        if (source == sink) return MAX_WEIGHT;

        std::array<PathWeight, maxV> level;
        std::map<std::pair<Node, Node>, PathWeight> remainingCap;
        std::array<typename std::set<Edge, EdgeComp>::iterator, maxV> lastSeen;

        std::function<PathWeight(Node, PathWeight)> dfs;
        dfs = [&](Node u, PathWeight pushed) {
            if (pushed == 0) return 0;
            if (u == sink) return pushed;

            for (; lastSeen[u] != edgesOut[u].end(); lastSeen[u]++) {
                Node v = lastSeen[u]->otherSide(u);
                if (level[u] + 1 == level[v] && remainingCap[{u, v}] > 0) {
                    PathWeight tr = dfs(v, std::min(pushed, remainingCap[{u, v}]));
                    if (tr != 0) {
                        remainingCap[{u, v}] -= tr;
                        remainingCap[{v, u}] += tr;
                        return tr;
                    }
                }
            }
            return 0;
        };

        std::function<void(void)> bfs = [&]() {
            level.fill(PathWeight(0));
            level[source] = PathWeight(1);
            std::queue<Node> q;
            q.push(source);

            while (!q.empty()) {
                Node u = q.front();
                q.pop();
                for (Edge edge : getEdgesOut(u)) {
                    Node v = edge.otherSide(u);
                    if (remainingCap[{u, v}] > 0 && level[v] == PathWeight(0)) {
                        level[v] = level[u] + 1;
                        q.push(v);
                    }
                }
            }
        };

        for (Edge edge : edges) {
            remainingCap[{edge.u, edge.v}] += edge.w;
            if (isUndirected()) {
                remainingCap[{edge.v, edge.u}] += edge.w;
            }
        }

        PathWeight f = 0;
        while (true) {
            bfs();
            if (level[sink] == PathWeight(0)) break;
            for (const Node &node : nodes)
                lastSeen[node] = edgesOut[node].begin();
            for (PathWeight pushed; (pushed = dfs(source, MAX_WEIGHT)) != PathWeight(0);
                 f += pushed)
                ;
        }
        return f;
    }

    //     // @TODO edge matching
    //     // @TODO see: https://usaco.guide/adv/

    /************************************************
     *                    PROPERTIES                *
     ***********************************************/

    // O(1)
    // @returns Whether the graph is directed
    static constexpr bool isDirected() {
        return isDirectedEdges;
    }

    // O(1)
    // @returns Whether the graph is undirected
    static constexpr bool isUndirected() {
        return !isDirectedEdges;
    }

    // O(1)
    // @returns Whether the graph is weighted
    static constexpr bool isWeighted() {
        return !std::is_same_v<EdgeWeight, UnitEdgeWeight>;
    }

    // // @TODO is bipartite graph? (try greedy colouring)
    // bool isBipartite() const {
    //     return greedyColouring() <= 2;
    // }

    //     // O(V_component log V_component + V_component log E + E_component)
    //     // @returns Whether there exists a edge which connects a node to itself
    //     // @param `root` If specified, consider its connected component. Otherwise consider the
    //     whole graph bool hasSelfEdges(Node root = 0) const {
    //         assert((root == 0 || containsNode(root)) && "Node index out of range");
    //         for (Node node: getComponentNodes(root)) {
    //             if (containsEdge(node, node)) return true;
    //         }
    //         return false;
    //     }

    //     // O(V_component log V_component + E_component)
    //     // @returns Whether there exists "redundant" edges
    //     // @param `root` If specified, consider its connected component. Otherwise consider the
    //     whole graph bool hasDoubleEdges(Node root = 0) const {
    //         for (Node node: getComponentNodes(root)) {
    //             for (Node type = 0; type <= 1; ++type) {
    //                 const std::multiset<std::pair<Node, T>> *currEdges = &(type ? edgesIn[node] :
    //                 edgesOut[node]); if (currEdges->empty()) continue; for (auto it1 =
    //                 currEdges->begin(), it2 = ++currEdges->begin(); it2 != currEdges->end();
    //                 ++it1, ++it2) {
    //                     if (it1->first == it2->first) return true;
    //                 }
    //             }
    //         }

    //         return false;
    //     }

    // O(V + E)
    // @returns Whether the graph has a cycle
    // @param `root` If specified, consider its connected component. Otherwise consider the whole
    // graph
    bool isCyclic() const {
        std::array<int, maxV> state; // @TODO actually only needs to hold 3 states
        state.fill(0);

        std::function<bool(Node)> dfs;
        dfs = [&](Node u) {
            if (state[u] == 2)
                return false;
            else if (state[u] == 1)
                return true;

            state[u] = 1;
            for (Node v : getNeighboursOut(u)) {
                if (dfs(v)) return true;
            }
            state[u] = 2;
            return false;
        };
        for (const Node &node : nodes) {
            if (dfs(node)) return true;
        }
        return false;
    }

    bool isAcyclic() const {
        return !isCyclic();
    }

    //     // O(V_component log V_component + V_component log E + E_component)
    //     // @returns Whether the graph is a simple graph
    //     // @param `root` If specified, consider its connected component. Otherwise consider the
    //     whole graph bool isSimpleGraph(Node root = 0) const {
    //         assert((root == 0 || containsNode(root)) && "Node index out of range");
    //         return !hasSelfEdges(root) && !hasDoubleEdges(root);
    //     }

    // O(V + E)
    // @returns Whether the graph is a directed acrylic graph
    // @param `root` If specified, consider its connected component. Otherwise consider the whole
    // graph
    bool isDAG() const {
        return isDirected() && isAcyclic();
    }

    // O(V + E) Determines whether the graph is a forest
    bool isForest() const {
        return isAcyclic();
    }

    // O(V + E)
    // @returns Whether the graph is a tree
    // @param `root` If specified, consider its connected component. Otherwise consider the whole
    // graph
    bool isTree() const {
        return isForest() && E() == V() - 1;
    }

    /************************************************
     *       ALGORITHMS (COMPLEXITLY CLASS NP)      *
     ************************************************/

    // @TODO hamiltonian path / tour
    // @TODO euler path / tour
    // @TODO covering set
    // @TODO optimal colouring
    // @TODO planar embedding?
};

template<size_t maxV>
using UnweightedDiGraph = Graph<maxV, UnitEdgeWeight, int, true>;
template<size_t maxV>
using DirectedUnweightedGraph = Graph<maxV, UnitEdgeWeight, int, true>;

template<size_t maxV>
using WeightedDiGraph = Graph<maxV, int, int, true>;
template<size_t maxV>
using DirectedWeightedGraph = Graph<maxV, int, int, true>;

template<size_t maxV>
using WeightedGraph = Graph<maxV, int, int, false>;
template<size_t maxV>
using UndirectedWeightedGraph = Graph<maxV, int, int, false>;

template<size_t maxV>
using UnWeightedGraph = Graph<maxV, UnitEdgeWeight, int, false>;
template<size_t maxV>
using UndirectedUnweightedGraph = Graph<maxV, UnitEdgeWeight, int, false>;
}; // namespace DS

#endif
