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
    // @note: My definition: root is never a leaf
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
        
        constexpr static Weight MAX_WEIGHT = std::numeric_limits<Weight>::max(); 
        constexpr static Depth logMaxV = numberOfBits(maxV);
        constexpr static Depth INVALID_DEPTH = std::numeric_limits<Depth>::max();

private:
        using MyTree = Tree<maxV, Weight>;

public:
        // Edges are represented internaly as (u, v) with weight w
        // u is the jmp, v is the child
        // NOTE: allowing self cycles, because the root has an edge to itself in lca
        struct Edge {
            Node u = 0, v = 0;
            Weight w = Weight(1);

            Edge() {}

            Edge(Node u_, Node v_) : u(u_), v(v_) {
                assert(u < maxV && v < maxV && "Node index out of range");
            }

            template<typename ForeignWeight>
            Edge(Node u_, Node v_, ForeignWeight w_) : u(u_), v(v_), w(Weight() + w_) {
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

            // friend bool operator<(const Edge &a, const Edge &b) {
            //     if (a.u != b.u) return a.u < b.u;
            //     else if (a.v != b.v) return a.v < b.v;
            //     else return a.w < b.w;
            // }

            friend bool operator==(const Edge &a, const Edge &b) {
                return a.u == b.u && a.v == b.v && a.w == b.w;
            }

            friend bool operator!=(const Edge &a, const Edge &b) {
                return !(a == b);
            }

            friend std::istream& operator>>(std::istream &in, Edge &e) {
                in >> e.u >> e.v >> e.w;
                return in;
            }

            friend std::ostream& operator<<(std::ostream &out, Edge &e) {
                out << e.u << "->" << e.v << " (w = " << e.w << ')';
                return out;
            }
        };
        
        // struct UnitEdge: Edge {
        //     friend std::istream& operator>>(std::istream &in, Edge &e) {
        //         in >> e.u >> e.v;
        //         return in;
        //     }
        // };

        // arbitary total order of edges

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

        // not garunteed sorted!
        std::vector<Node> nodes;

        // `depth[node]` is the number of edges between `node` and `root`
        // depth of invalid nodes is `INVALID_DEPTH`
        Depth depth[maxV];

        // height of root is 0, height of invalid node should not be accessed, +w for each child
        // height of invalid nodes should not be accessed
        Weight height[maxV];

        // `children[node]` is a vector of children, where each child is a (node, weight) pair
        // children of invalid nodes should not be accessed
        std::set<Edge, EdgeComp> children[maxV];

        // `jmp[d][node]` is the jmp which is 2^d steps higher than `node`, so `jmp[0][node]` is the direct jmp
        // jump of invalid nodes should not be accessed
        Edge jmp[logMaxV + 1][maxV];

        // O(N log N) Builds `jmps` from [1..LOG_MAXN], assuming `jmps[0]` has been built
        void buildJumpPtrs(const std::vector<Node> &todo) {
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
            std::fill(std::begin(depth), std::end(depth), INVALID_DEPTH);
            nodes = { root };
            depth[root] = 0;
            height[root] = Weight();
            children[root] = {};
            jmp[0][root] = 0;
            buildJumpPtrs(nodes);
        }

        // O(V log V)
        // Initialises tree from a graph, given a root
        // @note G is treated as undirected
        template<typename EdgeWeight, typename PathWeight = Weight, bool isDirected>
        Tree(const Graph<maxV, EdgeWeight, PathWeight, isDirected> &G, Node root_): root(root_) {
            assert(G.containsNode(root) && "Wanted root is not a node in input graph");
            assert(G.E() == G.V() - 1 && "Definitely not a tree"); // redundant check, remove?
            
            std::fill(std::begin(depth), std::end(depth), -1);
            
            std::function<void(Node, Depth, Weight)> buildDfs;
            buildDfs = [&](Node node, Depth d, Weight h) {
                nodes.push_back(node);
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
        Tree(size_t V, Weight w) {
            assert(V < maxV && "Insufficient capacity");

            std::fill(std::begin(depth), std::end(depth), -1);

            nodes.push_back(1);
            depth[root] = 0;
            height[root] = Weight();
            jmp[0][root] = 0;
            for (Node node = 2; node <= V; ++node) {
                nodes.push_back(node);
                depth[node] = depth[node - 1] + 1;
                height[node] = height[node - 1] + w;
                children[node - 1].push(Edge(node - 1, node, w));
                jmp[0][node] = Edge(node - 1, node, w);
            }

            buildJumpPtrs(nodes);
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
                    out << line << '\n';
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
         *                  PROPERTIES                  *
         ************************************************/

        size_t size() const {
            return nodes.size();
        }

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
            return nodes.size() - 1;
        }
        // O(1)
        // @returns `E`, the number of edges
        size_t M() const {
            return nodes.size() - 1;
        }

        // O(1)
        // @returns the root of the tree
        Node getRoot() const {
            return root;
        }

        // O(1)
        // @returns the imutable set of all ndoes
        const std::vector<Node> getNodes() const {
            return nodes;
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

        /************************************************
         *                INCIDENT DATA                 *
         ************************************************/

        // O(num children) = O(V)
        // @returns The imutable collection of edges (to children) incident to this node
        const std::vector<Edge> getEdges(Node node) const {
            assert(containsNode(node) && "Node index out of range");
            return std::vector<Edge>(children[node].begin(), children[node].end());
        }

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

            depth[e.v] = depth[e.u] + 1;
            height[e.v] = height[e.u] + e.w;
            children[e.u].insert(e);
            jmp[0][e.v] = e;
            buildJumpPtrs({ e.v });
        }

        // O(log V)
        // If the specified edge (with any edge weight) is present, remove it, otherwise silently do nothing
        // @note Can only remove edges incident to non-root leaves
        void erase(Node node) {
            assert(isLeaf(node) && "Only supports removing edges incident to leaves");

            Node par = children[jmp[0][node]];
            depth[node] = INVALID_DEPTH;
            children[par].erase(children[par].lower_bound(Edge(par, node, Weight())));
            // NOTE: left jump pointers as is
            // remember to not touch jump pointers if node is invalid
        }

        /************************************************
         *                    CONTAINS                  *
         ************************************************/

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

        // O(log E)
        // @returns the weight of the edge between nodes u and v, if it exists
        // @return default otherwise
        Weight getEdge(Node u, Node v, Weight defaultWeight) {
            if (!containsNode(u) || !containsNode(v)) return defaultWeight;
            auto edgeIt = children[u].find(Edge(u, v));
            return edgeIt == children[u].end() ? defaultWeight: edgeIt->w;
        }








        // TODO: stub pls fix
        bool isBinary() const {
            for (Node node : nodes) {
                if (children[node].size() > 2) return false;
            }
            return true;
        }

        // TODO: move to the right place, write documentation
        std::vector<Node> getBfsOrder() const {
            std::vector<Node> out = { root };
            for (size_t i = 0; i < out.size(); ++i) {
                for (Node child: getChildren(out[i])) {
                    out.push_back(child);
                }
            }
            return out;
        }

    //     // O(N) Finds all descendants of `node`
    //     // @returns vector of nodes
    //     std::vector<int> getDescendants(int node = 1) {
    //         return preOrder(node); // lmao
    //     }

        // O(1) Determines if `node` is a leaf node (has no children)
        // @note From my own definition, root is never a leaf
        bool isLeaf(Node node) {
            return children[node].empty() && node != root;
        }

    //     // /O(N) Finds the leaves (nodes that have no children) in ascending order
    //     std::vector<int> getLeaves() {
    //         std::vector<int> leaves;
    //         for (int node = 1; node <= N; ++node) {
    //             if (isLeaf(node)) leaves.push_back(node);
    //         }
    //         return leaves;
    //     }

    //     // O(N) Finds the breath of the tree - the number of leaves
    //     int getBreadth() {
    //         int total = 0;
    //         for (int node = 1; node <= N; ++node) total += children[node].size() == 0;
    //         return total;
    //     }


    //     /************************************************
    //      *                   PROPERTIES                 *
    //      ************************************************/

    //     // [ O(N log N) ] Returns whether the tree is a binary tree - there are at most a left and right child for each node
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole tree
    //     bool isBinary(int root = 1) {
    //         for (int node = 1; node <= N; ++node) {
    //             if (children[node].size() > 2 && isAncestor(root, node)) return false;
    //         }
    //         return true;
    //     }

    //     // O(N log N) Finds whether the tree is a balanced tree - height of the left and right subtree differ by not more than 1
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole tree
    //     bool isBalanced(int root = 1) {
    //         int subtreeSize = getDescendants(root).size();
    //         int maxDepth = 0;
    //         for (int node = 1; node <= N; ++node) {
    //             if (isAncestor(root, node)) maxDepth = std::max(maxDepth, depth[node] - depth[root]);
    //         }
    //         return (1 << maxDepth) <= subtreeSize && subtreeSize < 2 * (1 << maxDepth);
    //     }

    //     /************************************************
    //      *                   ALGORITHMS                 *
    //      ************************************************/

    //     // O(log N) Finds the lowest common ancestor of `u` and `v`, using jump pointers with `depth` and `jmp`
    //     int lca(int u, int v) {
    //         for (int d = LOG_MAXN; d >= 0; --d) {
    //             if (depth[jmp[d][u].first] >= depth[v]) u = jmp[d][u].first;
    //             if (depth[jmp[d][v].first] >= depth[u]) v = jmp[d][v].first;
    //         }

    //         for (int d = LOG_MAXN; d >= 0; --d) {
    //             if (jmp[d][u].first != jmp[d][v].first) {
    //                 u = jmp[d][u].first;
    //                 v = jmp[d][v].first;
    //             }
    //         }

    //         if (u == v) return u;
    //         assert(jmp[0][u].first == jmp[0][v].first);
    //         return jmp[0][u].first;
    //     }

    //     // O(log N) Finds whether `ancestor` is actually an ancestor of `node`. Note that u node is its own ancestor
    //     bool isAncestor(int ancestor, int node) {
    //         return lca(ancestor, node) == ancestor;
    //     }

    //     // O(log n) Finds the `nth` jmp of `node`
    //     // @note The `0th` jmp of `node` is `node`, and `1` is the jmp of itself
    //     int getParent(int node, int n = 1) {
    //         assert(n >= 0 && "Can't find a negative jmp");
    //         n = std::min(n, N);
    //         int par = node;
    //         for (int d = 0; n; ++d, n /= 2) {
    //             if (n & 1) par = jmp[d][par].first;
    //         }
    //         return par;
    //     }

    //     // O(log N) Finds the shortest distance between nodes `u` and `v`, by traversing `u` -> `lca(u, v)` -> `v`
    //     T dist(int u, int v) {
    //         int targetDepth = depth[lca(u, v)];
    //         T total = T(0);
    //         for (int d = LOG_MAXN; d >= 0; --d) {
    //             if (depth[jmp[d][u].first] >= targetDepth) {
    //                 total += jmp[d][u].second;
    //                 u = jmp[d][u].first;
    //             }
    //             if (depth[jmp[d][v].first] >= targetDepth) {
    //                 total += jmp[d][v].second;
    //                 v = jmp[d][v].first;
    //             }
    //         }

    //         return total;
    //     }

    //     // O(N) Find the path with the shortest distance, between nodes `u` and `v`
    //     // @returns u vector of nodes on the path, starting from `u` and ending with `v`
    //     std::vector<int> getPath(int u, int v) {
    //         std::vector<int> path;
    //         int mid = lca(u, v);
    //         for (int node = u; node != mid; node = jmp[0][node].first) path.push_back(node);
    //         path.push_back(mid);

    //         int midLen = path.size();
    //         for (int node = v; node != mid; node = jmp[0][node].first) path.push_back(node);
    //         reverse(path.begin() + midLen, path.end());

    //         return path;
    //     }

    // private:
    //     // O(N) Finds the diameter - the simple path with the maximum distance
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole tree
    //     // @returns Distance of the path, and a vector of nodes representing the path
    //     // @note Cannot have negative weight cycles
    //     std::pair<int, T> _diameterPath(std::pair<std::pair<int, int>, T>& best, int root = 1) {
    //         T furthest1 = T(0), furthest2 = T(0);
    //         int node1 = root, node2 = -1;

    //         for (std::pair<int, T> child : children[root]) {
    //             std::pair<int, T> path = _diameterPath(best, child.first);
    //             path.second += child.second;
    //             if (path.second > furthest1) {
    //                 furthest2 = furthest1;
    //                 node2 = node1;

    //                 furthest1 = path.second;
    //                 node1 = path.first;
    //             }
    //             else if (path.second > furthest2) {
    //                 furthest2 = path.second;
    //                 node2 = path.first;
    //             }
    //         }

    //         if (node1 != -1 && node2 != -1) {
    //             if (furthest1 + furthest2 > best.second) {
    //                 best = { {node1, node2}, furthest1 + furthest2 };
    //             }
    //         }

    //         return { node1, furthest1 };
    //     }

    // public:
    //     // O(N) Returns the end points of the diameter (the longest simple path through the tree)
    //     std::pair<int, int> diameterPath() {
    //         std::pair<std::pair<int, int>, T> best = { {-1, -1}, T(0) };
    //         _diameterPath(best);
    //         return best.first;
    //     }

    //     // O(N) Returns the length of the diameter (the longest simple path through the tree) 
    //     int diameterDist() {
    //         std::pair<std::pair<int, int>, T> best = { {-1, -1}, T(0) };
    //         _diameterPath(best);
    //         return best.second;
    //     }

    // private:
    //     // O(N) DFS to find in-order traversal, assuming tree is a binary tree
    //     std::vector<int> _inOrder(std::vector<int>& traversal, int node = 1) {
    //         assert(children[node].size() <= 2 && "In-order traversals are only valid in binary trees");
    //         if (children[node].size() > 0) _inOrder(traversal, children[node][0].first);
    //         traversal.push_back(node);
    //         if (children[node].size() > 1) _inOrder(traversal, children[node][1].first);
    //         return traversal;
    //     }

    //     // O(N) DFS to find a tree traversal
    //     // @param `pushPre` Whether to add the node to the traversal before exploring children
    //     // @param `pushPost` Whether to add the node to the traversal after exploring children
    //     std::vector<int> _multiOrder(std::vector<int>& traversal, bool pushPre, bool pushPost, int node = 1) {
    //         if (pushPre) traversal.push_back(node);
    //         for (std::pair<int, T> child : children[node]) _multiOrder(traversal, pushPre, pushPost, child.first);
    //         if (pushPost) traversal.push_back(node);
    //         return traversal;
    //     }

    // public:
    //     // O(N) Finds in-order traversal, assuming tree is a binary tree
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole tree
    //     std::vector<int> inOrder(int root = 1) {
    //         std::vector<int> traversal;
    //         return _inOrder(traversal, root);
    //     }

    //     // O(N) Finds pre-order traversal
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole tree
    //     std::vector<int> preOrder(int root = 1) {
    //         std::vector<int> traversal;
    //         return _multiOrder(traversal, true, false, root);
    //     }

    //     // O(N) Finds post-order traversal
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole tree
    //     std::vector<int> postOrder(int root = 1) {
    //         std::vector<int> traversal;
    //         return _multiOrder(traversal, false, true, root);
    //     }

    //     // O(N) Finds euler tour traversal
    //     // @param `root` If specified, consider its connected component. Otherwise consider the whole tree
    //     std::vector<int> eulerTour(int root = 1) {
    //         std::vector<int> traversal;
    //         return _multiOrder(traversal, true, true, root);
    //     }
    // };
    };
};