#pragma once
#include "Util.h"
#include "Constants.h"
#include "Graph.h"
#include <limits>
#include <type_traits>

constexpr size_t numberOfBits(size_t x) {
    return x < 2 ? x : 1 + numberOfBits(x >> 1);
}

namespace DS {
    constexpr size_t TREE_PRINT_SPACING = 2;

    // Tree
    // @note: does not preserve/store order of children nodes
    // @note: there is no distinction between weighted and unweighted trees
    // because all algorithms can be applied regardless
    // @note: asserts that Weight() is the 0 (and minimum) weight
    // @note: if one child in a binary tree, it is considered left child
    template<
        // Max number of nodes allowed
        size_t maxV,
        // accumulation of weight data over many edges
        typename Weight_,
        std::enable_if_t<
            std::is_arithmetic_v<Weight_> &&
            Weight_() == std::numeric_limits<Weight_>::min() &&
            maxV != 0
        , bool> = true
    >
    class Tree {
        /************************************************
         *                 INITIALISATION               *
         ************************************************/

public:
        using Weight = Weight_;
        // id of node is an unsigned integer
        using Node = size_t;
        using Depth = size_t;

        constexpr static Weight INVALID_WEIGHT = 0;
        constexpr static Depth logMaxV = numberOfBits(maxV);
        constexpr static Depth INVALID_DEPTH = std::numeric_limits<Depth>::max();

private:
        using MyTree = Tree<maxV, Weight>;

public:
        // Edges are represented internaly as (u, v) with weight w
        // u is the jmp, v is the child
        // NOTE: allowing self cycles, because the root has an edge to itself in lca
        struct UndirectedPath {
            Node u = 0, v = 0;
            Weight w = Weight(1);

            UndirectedPath() {}

            UndirectedPath(Node u_, Node v_) : u(u_), v(v_) {
                assert(u < maxV && v < maxV && "Node index out of range");
            }

            template<typename ForeignWeight>
            UndirectedPath(Node u_, Node v_, ForeignWeight w_) : u(u_), v(v_), w(Weight() + w_) {
                assert(u < maxV && v < maxV && "Node index out of range");
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

            // friend bool operator<(const UndirectedPath &a, const UndirectedPath &b) {
            //     if (a.u != b.u) return a.u < b.u;
            //     else if (a.v != b.v) return a.v < b.v;
            //     else return a.w < b.w;
            // }

            friend bool operator==(const UndirectedPath &a, const UndirectedPath &b) {
                return a.w == b.w && (
                    (a.u == b.u && a.v == b.v) ||
                    (a.u == b.v && a.v == b.u)
                );
            }

            friend bool operator!=(const UndirectedPath &a, const UndirectedPath &b) {
                return !(a == b);
            }

            friend std::istream& operator>>(std::istream &in, UndirectedPath &e) {
                in >> e.u >> e.v >> e.w;
                return in;
            }

            friend std::ostream& operator<<(std::ostream &out, UndirectedPath &e) {
                out << e.u << "->" << e.v << " (w = " << e.w << ')';
                return out;
            }
        };

        struct Edge: UndirectedPath {
            using UndirectedPath::u;
            using UndirectedPath::v;
            using UndirectedPath::w;
            using UndirectedPath::UndirectedPath;

            // If this edge is (u, v), returns edge (v, u)
            // @note Undirected edges can be flipped
            Edge flip() const {
                return Edge(v, u, w);
            }

            friend bool operator==(const Edge &a, const Edge &b) {
                return a.u == b.u && a.v == b.v && a.w == b.w;
            }

            friend bool operator!=(const Edge &a, const Edge &b) {
                return !(a == b);
            }
        };

        // arbitary total order of edges
        // struct UnitEdge: Edge {
        //     friend std::istream& operator>>(std::istream &in, Edge &e) {
        //         in >> e.u >> e.v;
        //         return in;
        //     }
        // };

        // Directed edges can represent simple (possible empty) directed paths
        // because every simple path is defined uniquely by its endpoints
        using DirectedPath = Edge;

private:
        struct EdgeComp {
            bool operator() (const Edge &a, const Edge &b) const {
                if (a.u != b.u) return a.u < b.u;
                return a.v < b.v;
            }
        };

        // root node
        // must be a valid node
        Node root;

        // garunteed sorted!
        std::set<Node> nodes;

        // `depth[node]` is the number of edges between `node` and `root`
        // depth of invalid nodes is `INVALID_DEPTH`
        std::array<Depth, maxV> depth;

        // height of root is 0, height of invalid node should not be accessed, +w for each child
        // height of invalid nodes should not be accessed
        std::array<Weight, maxV> height;

        // `children[node]` is a vector of children, where each child is a (node, weight) pair
        // children of invalid nodes should not be accessed
        std::array<std::set<Edge, EdgeComp>, maxV> children;

        // `jmp[d][node]` is the jmp which is 2^d steps higher than `node`, so `jmp[0][node]` is the direct jmp
        // jump of invalid nodes should not be accessed
        Edge jmp[logMaxV + 1][maxV];

        // O(N log N) Builds `jmps` from [1..LOG_MAXN], assuming `jmps[0]` has been built
        void buildJumpPtrs(const std::set<Node> &todo) {
            for (Depth d = 0; d < logMaxV; ++d) {
                for (Node node : todo) {
                    jmp[d + 1][node] = Edge(
                        jmp[d][jmp[d][node].u].u,
                        node,
                        jmp[d][node].w + jmp[d][jmp[d][node].u].w
                    );
                }
            }
        }

public:
        // O(V)
        // Initialises a singleton tree with the node labeled `root`
        Tree(Node root_): root(root_) {
            assert(isNode(root) && "Invalid root");

            depth.fill(INVALID_DEPTH);
            height.fill(INVALID_WEIGHT);

            nodes = { root };
            depth[root] = 0;
            height[root] = Weight();
            children[root] = {};
            jmp[0][root] = Edge(root, root, Weight());
            buildJumpPtrs(nodes);
        }

        // O(V log V)
        // Initialises tree from a graph, given a root
        // @note G is treated as undirected
        template<typename EdgeWeight, typename PathWeight = Weight, bool isDirected>
        Tree(const Graph<maxV, EdgeWeight, PathWeight, isDirected> &G, Node root_): root(root_) {
            assert(G.containsNode(root) && "Wanted root is not a node in input graph");
            assert(G.E() == G.V() - 1 && "Definitely not a tree"); // redundant check, remove?
            
            depth.fill(INVALID_DEPTH);
            height.fill(INVALID_WEIGHT);

            std::function<void(Node, Depth, Weight)> buildDfs;
            buildDfs = [&](Node node, Depth d, Weight h) {
                nodes.insert(node);
                depth[node] = d;
                height[node] = h;
                for (auto edge : G.getEdges(node)) {
                    Node child = edge.otherSide(node);
                    if (depth[child] == INVALID_DEPTH) {
                        Edge treeEdge = Edge(node, child, Weight() + edge.w);
                        children[node].insert(treeEdge);
                        jmp[0][child] = treeEdge;
                        buildDfs(child, d + 1, h + edge.w);
                    }
                }
            };
            jmp[0][root] = Edge(root, root, Weight());
            buildDfs(root, 0, Weight());
            
            buildJumpPtrs(nodes);
        }

        // O(V log V)
        // Initialise path graph with nodes 1..V, rooted at 1
        static MyTree PathGraph(size_t V, Weight w) {
            assert(V < maxV && "Insufficient capacity");
            MyTree T(1);
            for (Node node = 2; node <= V; ++node) {
                T.insert(Edge(node - 1, node, w));
            }
            return T;
        }

        // O(V log V)
        // Initialise star graph with nodes 1..V, rooted at 1
        static MyTree StarGraph(size_t V, Weight w) {
            assert(V < maxV && "Insufficient capacity");
            MyTree T(1);
            for (Node node = 2; node <= V; ++node) {
                T.insert(Edge(1, node, w));
            }
            return T;
        }

        /************************************************
         *                    DISPLAY                   *
         ************************************************/

    public:
        // O(V) Displays the graph to `out` with fancy indenting for depth, showing the jmp of each node
        friend std::ostream& operator<<(std::ostream &out, const MyTree &tree) {
            std::function<void(Node, std::string)> print;
            print = [&](Node node, std::string prefix) {
                size_t i = 0;
                for (Edge edge : tree.children[node]) {
                    bool isLast = i++ == tree.children[node].size() - 1;
                    Node child = edge.otherSide(node);
                    out << prefix << (isLast ? "└" : "├") << std::string("─") * edge.w << child << '\n';
                    print(child, prefix + (isLast ? " " : "│") + std::string(" ") * edge.w);
                }
            };
            out << tree.root << '\n';
            print(tree.root, "");
            return out;
        }

        // O(maxV)
        // trying to do it with minimal additional memory is too painful
        // friend std::ostream& operator<<(std::ostream &out, const MyTree *tree) {
        //     size_t spacing = 1;
        //     std::vector<size_t> offset(maxV, 0);
        //     std::function<size_t(Node, size_t)> findSize;
        //     findSize = [&](Node node, size_t currOffset) {
        //         if (tree->numChildren(node) == 0) {
        //             offset[node] = currOffset;
        //             return repr(node).size();
        //         } else if (tree->numChildren(node) == 1) {
        //             Node child = tree->getChildren(node)[0];
        //             size_t childSize = findSize(child, currOffset);
        //             offset[node] = currOffset;
        //             return std::max(childSize, repr(node).size());
        //         } else if (tree->numChildren(node) == 2) {
        //             Node leftChild = tree->getChildren(node)[0];
        //             Node rightChild = tree->getChildren(node)[1];
        //             size_t leftSize = findSize(leftChild, currOffset);
        //             offset[node] = currOffset + leftSize + spacing;
        //             size_t rightSize = findSize(rightChild, currOffset + leftSize + spacing + repr(node).size() + spacing);
        //             return leftSize + spacing + repr(node).size() + spacing + rightSize;
        //         } else {
        //             assert(false && "Cannot print non binary tree in binary format, use `out << tree` instead");
        //         }
        //     };
        //     findSize(tree->root, 0);

        //     std::vector<Node> bfsOrder = tree->getBfsOrder();
        //     size_t i = 0;
        //     for (Depth d = 0; i < bfsOrder.size(); ++d) {
        //         size_t currOffset = 0;
        //         for ( ; i < bfsOrder.size() && tree->depth[bfsOrder[i]] == d; ++i) {
        //             while (currOffset < offset[bfsOrder[i]]) {
        //                 out << ' ';
        //                 currOffset++;
        //             }
        //             out << bfsOrder[i];
        //             currOffset += repr(bfsOrder[i]).size();
        //         }
        //         out << '\n';
        //     }

        //     return out;
        // }

private:
        // additional requirement that all lines have equal length
        struct Display {
            size_t width;
            std::deque<std::string> display;
            size_t offset = 0;

            void operator+=(Display o) {
                while (display.size() < o.display.size()) {
                    display.push_back(std::string(" ") * width);
                }
                for (size_t i = 0; i < display.size(); ++i) {
                    display[i] += i < o.display.size() ? o.display[i] : std::string(" ") * o.width;
                }
                width += o.width;
            }

            friend std::ostream& operator<<(std::ostream& out, Display d) {
                for (std::string line: d.display) {
                    out << line.substr(0, line.size() - TREE_PRINT_SPACING) << '\n';
                }
                return out;
            }
        };


        Display print(Node node) const {
            std::string nodeData = repr(node);
            if (children[node].empty()) {
                return Display{
                    nodeData.size() + TREE_PRINT_SPACING, 
                    { nodeData + std::string(" ") * TREE_PRINT_SPACING },
                    nodeData.size() / 2
                };
            }

            std::vector<size_t> offsets;
            Display d = Display{0, {}};
            for (Node child: getChildren(node)) {
                Display childDisplay = print(child);
                offsets.push_back(d.width + childDisplay.offset);
                d += childDisplay;
            }

            d.offset = (offsets.front() + offsets.back()) / 2;

            std::string header;
            size_t i = 0;
            for (; i < offsets[0]; ++i) header += " ";
            for (size_t offsetI = 0; offsetI < offsets.size(); i++, offsetI++) {
                for (; i < offsets[offsetI]; ++i) header += i == d.offset ? "┴" : "─";
                if (offsets.size() == 1) {
                    header += i == d.offset ? "│" : "│";
                } else if (offsetI == 0) {
                    header += i == d.offset ? "├" : "┌";
                } else if (offsetI == offsets.size() - 1) {
                    header += i == d.offset ? "┤" : "┐";
                } else {
                    header += i == d.offset ? "┼" : "┬";
                }
            }
            for (; i < d.width; ++i) header += " ";

            d.display.push_front(header);
            d.display.push_front(
                std::string(" ") * (d.offset - nodeData.size() / 2) + 
                nodeData +
                std::string(" ") * (d.width - d.offset + nodeData.size() / 2 - nodeData.size())
            );
            return d;
        };

public:
        // O(maxV ^ 2) probably
        friend std::ostream& operator<<(std::ostream &out, const MyTree *tree) {
            out << tree->print(tree->root);
            return out;
        }


        /************************************************
         *                     NODES                    *
         ************************************************/

        // O(1)
        // @returns `V`, the number of nodes
        size_t size() const {
            return nodes.size();
        }

        // O(1)
        // @returns `V`, the number of nodes
        size_t V() const {
            return nodes.size();
        }

        // O(1)
        // @returns `V`, the number of nodes
        size_t N() const {
            return nodes.size();
        }

        // O(1)
        // @returns the root of the tree
        Node getRoot() const {
            return root;
        }

        // O(1)
        // @returns whether `node` is within the acceptable bounds
        bool isNode(Node node) const {
            // NOTE: dont need to check 0 <= node, because Node=size_t
            return node < maxV;
        }

        // O(1)
        // @returns whether `node` is in the vertex set of this graph
        bool containsNode(Node node) const {
            // NOTE: dont need to check 0 <= node, because Node=size_t
            return node < maxV && depth[node] != INVALID_DEPTH;
        }

        // O(1)
        // @returns the imutable set of all ndoes
        const std::vector<Node> getNodes() const {
            return std::vector<Node>(nodes.begin(), nodes.end());
        }

        
        // O(1)
        // Determines if `node` is a leaf node (has no children)
        bool isLeaf(Node node) const {
            return containsNode(node) && children[node].empty();
        }

        // O(V)
        // Finds the leaves
        // @returns vector of leaves, in increasing order
        std::vector<Node> getLeaves() const {
            std::vector<Node> leaves;
            for (Node node: nodes) {
                if (isLeaf(node)) leaves.push_back(node);
            }
            return leaves;
        }

        /************************************************
         *                    EDGES                     *
         ************************************************/

        // O(1)
        // @returns `E`, the number of edges
        size_t E() const {
            return nodes.size() - 1;
        }
        // O(1)
        // @returns `E`, the number of edges
        size_t M() const {
            return nodes.size() - 1;
        }

        // O(V)
        // @returns the imutable set of all edges
        // @note Edges sorted by endpoints
        std::vector<Edge> getEdges() const {
            std::vector<Edge> out;
            out.reserve(E()); // optimisation
            for (Node node: nodes) {
                out.insert(out.end(), children[node].begin(), children[node].end());
            }
            return out;
        }

        // O(num children) = O(V)
        // @returns The imutable collection of edges (to children) incident to this node
        const std::vector<Edge> getEdges(Node node) const {
            assert(containsNode(node) && "Node index out of range");
            return std::vector<Edge>(children[node].begin(), children[node].end());
        }

        // O(log E)
        // @returns the weight of the edge between nodes u and v, if it exists
        // @returns default otherwise
        Weight getEdge(Node u, Node v, Weight defaultWeight) {
            if (!containsNode(u) || !containsNode(v)) return defaultWeight;
            auto edgeIt = children[u].find(Edge(u, v));
            return edgeIt == children[u].end() ? defaultWeight: edgeIt->w;
        }

        // O(log E)
        // @returns Whether the specified edge (with any weight) is contained in the graph
        bool containsEdge(Edge e) const {
            if (!containsNode(e.u) || !containsNode(e.v)) return false;
            auto edgeIt = children[e.u].find(e);
            return edgeIt != children[e.u].end() && edgeIt->w == e.w;
        }

        // O(log E)
        // @returns Whether the specified edge (with the specific weight) is contained in the graph
        bool containsEdgeUnweighted(Edge e) const {
            if (!containsNode(e.u) || !containsNode(e.v)) return false;
            return children[e.u].count(e);
        }


        // void initEdge(Edge e) const {
        //     assert(containsNode(e.u) && containsNode(e.v) && "Endpoint out of range");
        //     if (depth[e.u] > depth[e.v]) std::swap(e.u, e.v);
        //     e.w = dist(e.u, e.v);
        // }

        // O(num children) = O(V)
        // @returns The imutable colelction of children, in sorted order
        const std::vector<Node> getChildren(Node node) const {
            assert(containsNode(node) && "Node index out of range");
            std::vector<Node> out;
            for (const Edge e : getEdges(node)) {
                out.push_back(e.otherSide(node));
            }
            return out;
        }

        // O(1)
        // @returns The total degree of `node`, including the parent (if there is one)
        size_t degree(Node node) const {
            assert(containsNode(node) && "Node index out of range");
            return children[node].size() + (node != root);
        }

        // O(1)
        // @returns The number of children of `node`
        size_t numChildren(Node node) const {
            assert(containsNode(node) && "Node index out of range");
            return children[node].size();
        }

        /************************************************
         *               TREE MODIFICATIONS             *
         ************************************************/

        // O(log V)
        // Add a edge from an existing node (parent) to a new node, or error if this isn't possible
        void insert(Edge e) {
            // NOTE: this check is duplicated, but included for clear error messages
            assert(containsNode(e.u) && "Parent node of edge already in graph");
            assert(!containsNode(e.v) && "Child node of edge already in graph");

            nodes.insert(e.v);
            depth[e.v] = depth[e.u] + 1;
            height[e.v] = height[e.u] + e.w;
            children[e.u].insert(e);
            jmp[0][e.v] = e;
            buildJumpPtrs({ e.v });
        }

        // O(log V)
        // If the specified edge (with any edge weight) is present, remove it, otherwise silently do nothing
        // @note Can only remove edges incident to non-root leaves
        // TODO: support remove proper subtrees
        void erase(Node node) {
            assert(isLeaf(node) && node != root && "Only supports removing edges incident to leaves");

            nodes.erase(node);
            Node par = jmp[0][node].u;
            depth[node] = INVALID_DEPTH;
            height[node] = INVALID_WEIGHT;
            children[par].erase(children[par].lower_bound(Edge(par, node, Weight())));
            // NOTE: leave jump pointers as is
            // remember to not touch jump pointers if node is invalid
        }


        /************************************************
         *               ANCESTOR RELATED               *
         ************************************************/

        // O(1)
        // @returns The edge to the immediate parent of node, where root is its own parent
        Edge getParent(Node node) const {
            assert(containsNode(node) && "Node index out of range");
            return jmp[0][node];
        }

        // O(log V)
        // @returns The `n`th parent of node, where root is the parent of itself
        // @note A node's 0th parent is itself
        Edge getParent(Node node, size_t n) const {
            assert(containsNode(node) && "Node index out of range");
            size_t bits = std::min(n, E());

            Node par = node;
            Weight w = Weight();
            for (Depth d = 0; bits; bits /= 2, d++) {
                if (bits & 1) {
                    w += jmp[d][par].w;
                    par = jmp[d][par].u;
                }
            }
            return Edge(par, node, w);
        }
 
        // O(log V)
        // Finds whether `ancestor` is actually an ancestor of `node`
        bool isAncestor(Node ancestor, Node node) const {
            assert(containsNode(ancestor) && "Node index out of range"); // or should return false?
            assert(containsNode(node) && "Node index out of range"); // or should return false?

            // TODO: test which implementation is faster

            // if (depth[ancestor] > depth[node]) return false;
            // return getParent(node, depth[node] - depth[ancestor]) == node;
            
            for (Depth d = logMaxV; d-- > 0; ) {
                if (depth[jmp[d][node].u] >= depth[ancestor]) {
                    node = jmp[d][node].u;
                }
            }
            return node == ancestor;
        }

        // O(log V)
        // @returns The lowest common ancestor of `u` and `v`
        // uses jump pointers with `depth` and `jmp`
        Node lca(Node u, Node v) const {
            assert(containsNode(u) && containsNode(v) && "Node index out of range");
            for (Depth d = logMaxV; d-- > 0; ) {
                if (depth[jmp[d][u].u] >= depth[v]) u = jmp[d][u].u;
                if (depth[jmp[d][v].u] >= depth[u]) v = jmp[d][v].u;
            }

            for (Depth d = logMaxV; d-- > 0; ) {
                if (jmp[d][u].u != jmp[d][v].u) {
                    u = jmp[d][u].u;
                    v = jmp[d][v].u;
                }
            }

            assert(jmp[0][u].u == jmp[0][v].u);
            return u == v ? u : jmp[0][u].u;
        }

        // O(log N)
        // Finds the distance between nodes `u` and `v`, by traversing `u` -> `lca(u, v)` -> `v`
        // @returns non negative weight
        Weight dist(Node u, Node v) const {
            assert(containsNode(u) && containsNode(v) && "Node index out of range");
            Depth targetDepth = depth[lca(u, v)];
            Weight total = Weight();
            for (Depth d = logMaxV; d-- > 0; ) {
                if (depth[jmp[d][u].u] >= targetDepth) {
                    total += jmp[d][u].w;
                    u = jmp[d][u].u;
                }
                if (depth[jmp[d][v].u] >= targetDepth) {
                    total += jmp[d][v].w;
                    v = jmp[d][v].u;
                }
            }

            return total;
        }

        // O(N)
        // Find the path with the shortest distance, between nodes `u` and `v`
        // @returns u vector of nodes on the path, starting from `u` and ending with `v`
        std::vector<Node> getPathNodes(Edge e) const {
            std::vector<Node> path;
            Node mid = lca(e.u, e.v);
            for (Node node = e.u; node != mid; node = jmp[0][node].u) path.push_back(node);
            path.push_back(mid);

            Node midLen = path.size();
            for (Node node = e.v; node != mid; node = jmp[0][node].u) path.push_back(node);
            reverse(path.begin() + midLen, path.end());

            return path;
        }

        // O(N)
        // Find the path with the shortest distance, between nodes `u` and `v`
        // @returns u vector of nodes on the path, starting from `u` and ending with `v`
        std::vector<Edge> getPathEdges(Edge e) const {
            std::vector<Edge> path;
            Node mid = lca(e.u, e.v);
            for (Node node = e.u; node != mid; node = jmp[0][node].u) path.push_back(jmp[0][node].flip());

            Node midLen = path.size();
            for (Node node = e.v; node != mid; node = jmp[0][node].u) path.push_back(jmp[0][node]);
            reverse(path.begin() + midLen, path.end());

            return path;
        }

        /************************************************
         *                   ALGORITHMS                 *
         ************************************************/

        Weight getHeight() const {
            Weight maxHeight = Weight();
            for (Node node: nodes) {
                maxHeight = std::max(maxHeight, height[node]);
            }
            return maxHeight;
        }
        
        // O(maxV)
        // @returns array, where the ith weight is the height of node i from the root, or INVALID_WEIGHT if node i is invalid
        const std::array<Weight, maxV>& getHeightFromRoot() const {
            return height;
        }

        // O(maxV)
        // @returns array, where the ith weight is the height of node i from the furthest leaf, or INVALID_WEIGHT if node i is invalid
        std::array<Weight, maxV> getHeightFromLeaves() const {
            std::array<Weight, maxV> leafHeight;
            leafHeight.fill(INVALID_WEIGHT);
            std::function<void(Node)> dfs;
            dfs = [&](Node node) {
                leafHeight[node] = Weight();
                for (Edge edge: children[node]) {
                    dfs(edge.v);
                    leafHeight[node] = std::max(leafHeight[node], leafHeight[edge.v] + edge.w);
                }
            };
            dfs(root);
            return leafHeight;
        }

        // O(N)
        // Finds the diameter - the simple path with the maximum distance
        // @returns Distance of the path, and a vector of nodes representing the path
        Weight getDiameter() const {
            std::array<Weight, maxV> leafHeight = getHeightFromLeaves();
            Weight best = Weight();
            
            for (Node node: nodes) {
                Weight furthest1 = Weight(), furthest2 = Weight();
                for (Edge edge : children[node]) {
                    Weight curr = leafHeight[edge.v] + edge.w;
                    if (curr >= furthest1) {
                        furthest2 = furthest1;
                        furthest1 = curr;
                    } else if (curr > furthest2) {
                        furthest2 = curr;
                    }
                }
                if (furthest1 + furthest2 > best) {
                    best = furthest1 + furthest2;
                }
            }
            return best;
        }

        /************************************************
         *                  ORDERINGS                   *
         ************************************************/

        // O(V)
        // @returns the bfs order of the tree
        // @note children in increasing order
        std::vector<Node> getBfsOrder() const {
            std::vector<Node> traversal = { root };
            traversal.reserve(V());
            for (size_t i = 0; i < traversal.size(); ++i) {
                for (Node child: getChildren(traversal[i])) {
                    traversal.push_back(child);
                }
            }
            assert(traversal.size() == V());
            return traversal;
        }

        // O(V)
        // @returns the in-order traversal of all nodes of a binary tree, or throws
        std::vector<Node> getInOrder() const {
            std::vector<Node> traversal;
            traversal.reserve(V());
            std::function<void(Node)> dfs;
            dfs = [&](Node node) {
                assert(numChildren(node) <= 2 && "Can only get in order traversal of binary tree");
                if (numChildren(node) >= 1) dfs(children[node].begin()->v);
                traversal.push_back(node);
                if (numChildren(node) >= 2) dfs((++children[node].begin())->v);
            };
            dfs(root);
            assert(traversal.size() == V());
            return traversal;
        }

        // O(V)
        // @returns the in-order traversal of all nodes
        std::vector<Node> getPreOrder() const {
            std::vector<Node> traversal;
            traversal.reserve(V());
            std::function<void(Node)> dfs;
            dfs = [&](Node node) {
                traversal.push_back(node);
                for (Node child: getChildren(node)) {
                    dfs(child);
                }
            };
            dfs(root);
            assert(traversal.size() == V());
            return traversal;
        }

        // O(V)
        // @returns the in-order traversal of all nodes
        std::vector<Node> getPostOrder() const {
            std::vector<Node> traversal;
            traversal.reserve(V());
            std::function<void(Node)> dfs;
            dfs = [&](Node node) {
                for (Node child: getChildren(node)) {
                    dfs(child);
                }
                traversal.push_back(node);
            };
            dfs(root);
            assert(traversal.size() == V());
            return traversal;
        }

        // O(V)
        // @returns the in-order traversal of all nodes
        std::vector<Node> getEulerTour() const {
            std::vector<Node> traversal;
            traversal.reserve(2 * V());
            std::function<void(Node)> dfs;
            dfs = [&](Node node) {
                assert(numChildren(node) <= 2 && "In-order traversals are only valid in binary trees");
                traversal.push_back(node);
                for (Node child: getChildren(node)) {
                    dfs(child);
                }
                traversal.push_back(node);
            };
            dfs(root);
            assert(traversal.size() == 2 * V());
            return traversal;
        }

        /************************************************
         *                  TREE TYPES                  *
         ************************************************/

        // O(V)
        // Finds whether the tree is a height balanced tree - depth of every left and right subtree differ by not more than 1
        bool isHeightBalanced() const {
            bool balanced = true;
            std::function<Depth(Node)> dfs;
            dfs = [&](Node node) {
                Depth minDepth = std::numeric_limits<Depth>::max();
                Depth maxDepth = 0;
                for (Node child: getChildren(node)) {
                    Depth currDepth = dfs(child);
                    minDepth = std::min(minDepth, currDepth);
                    maxDepth = std::max(maxDepth, currDepth);
                }
                balanced &= maxDepth - minDepth <= 1;
                return maxDepth;
            };
            return balanced;
        }

        // O(V)
        // Finds whether this this is a node balanced tree - whether V is in the range [2^ddepth .. 2*2^depth)
        // Should this still apply to non binary trees?
        bool isNodeBalanced() const {
            Depth depth = getHeight();
            return (1U << depth) <= V() && V() < 2 * (1U << depth);
        }

        // O(V)
        // @returns Whether this tree is a binary tree - there are at most a left and right child for each node
        bool isBinary() const {
            for (Node node : nodes) {
                if (numChildren(node) > 2) return false;
            }
            return true;
        }

        // O(V)
        // @returns Whether this tree is a path graph
        bool isPath() const {
            for (Node node : nodes) {
                if (!isLeaf(node) && numChildren(node) != 1) return false;
            }
            return true;
        }

        // O(1)
        bool isStar() const {
            Node node = *nodes.begin();
            return degree(node) == V() - 1 || (degree(node) == 1 && degree(getChildren(node)[0]) == V() - 1);
        }

        /************************************************
         *                     MISC                     *
         ************************************************/


        // void rootAt() {
        //     // TODO
        // }

        // MyTree getSubtree(Node node) {
        //     // TODO
        // }
    };
};