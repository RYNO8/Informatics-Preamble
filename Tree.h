template<typename T> class Tree {
	/************************************************
	 *                INITIALISATION                *
	 ************************************************/

private:
	// Number of nodes
	int N;
	// Whether edges are weighted
	// Unweighted edges are assigned unit weights
	bool isWeighted;
	// `parent[d][node]` is the parent which is 2^d steps higher than `node`, so `parent[0][node]` is the direct parent
	vector<pair<int, T>> parent[LOG_MAXN + 5];
	// `children[node]` is a vector of children, where each child is a (node, weight) pair
	vector<vector<pair<int, T>>> children;
	// `depth[node]` is the distance from node 1, with `depth[1] = 1`
	vector<int> depth;

	// O(N log N) Initialises all variables
	void init(int _N, bool _isWeighted) {
		assert(1 <= _N && _N <= MAXN);
		N = _N;
		isWeighted = _isWeighted;
		for (int d = 0; d <= LOG_MAXN; ++d) parent[d].resize(N + 1);
		parent[0][1] = { 1, T(0) };
		children.resize(N + 1);
		depth.resize(N + 1);
	}
	
	// O(N) Find direction of edge each such that the tree can be rooted at 1
	void makeDirected(Graph<T>& g, int node = 1, int par = -1) {
		for (pair<int, T> child : g.getOutEdges(node)) {
			if (child.first == par) continue;
			parent[0][child.first] = { node, child.second };
			children[node].push_back(child);
			makeDirected(g, child.first, node);
		}

		if (g.isDirectedGraph()) {
			for (pair<int, T> child : g.getInEdges(node)) {
				if (child.first == par) continue;
				parent[0][child.first] = { node, child.second };
				children[node].push_back(child);
				makeDirected(g, child.first, node);
			}
		}
	}

public:
	// O(N log N) Initialises a tree from a graph, assuming rooted at 1
	Tree(Graph<T>& g) {
		assert(g.isTree() && "Initialiser graph must be a tree");
		init(g.size(), g.isWeightedGraph());

		makeDirected(g);

		sortChildren();
		buildJumpPtrs();
		buildDepths();
	}

	// O(N log N) Initialises a tree from cin, assuming rooted at 1*/
	Tree(int _N, istream& in = cin, bool _isWeighted = false) {
		init(_N, _isWeighted);

		for (int i = 2; i <= N; ++i) {
			in >> parent[0][i].first;
			if (isWeighted) in >> parent[0][i].second;
			else parent[0][i].second = T(1);
			children[parent[0][i].first].push_back({ i, parent[0][i].second });
		}

		sortChildren();
		buildJumpPtrs();
		buildDepths();
	}

	// O(1) Initialises a tree with a single node (node 1)
	Tree(bool _isWeighted = false) {
		init(1, _isWeighted);

		sortChildren();
		buildJumpPtrs();
		buildDepths();
	}

private:
	// O(N log N) Ensure that the children are always in sorted order
	void sortChildren() {
		for (int node = 1; node <= N; ++node) sort(children[node].begin(), children[node].end());
	}

	// O(N log N) Builds `parents` from [1..LOG_MAXN], assuming `parents[0]` has been built
	void buildJumpPtrs() {
		for (int d = 0; d < LOG_MAXN; ++d) {
			for (int i = 1; i <= N; ++i) {
				parent[d + 1][i] = {
					parent[d][parent[d][i].first].first,
					parent[d][i].second + parent[d][parent[d][i].first].second
				};
			}
		}
	}

	// O(N) Builds `depth` using a dfs
	void buildDepths(int node = 1) {
		depth[node] = depth[parent[0][node].first] + 1;
		for (pair<int, T> child : children[node]) buildDepths(child.first);
	}

	/************************************************
	 *                    DISPLAY                   *
	 ************************************************/

public:
	// O(N) Displays the graph, showing the parent of each node
	// @param `out` The string representation of the graph is piped to this output stream
	// @param `newLine` Indicates whether to end with a trailing `\\n`
	void print(ostream& out = cout, bool newLine = true) const {
		for (int i = 1; i <= N; ++i) out << parent[0][i].first << ' ';
		out << '\n';
		if (newLine) out << '\n';
	}

	// O(N) Displays the graph to `out` with fancy indenting for depth, showing the parent of each node
	// @param `out` The string representation of the graph is piped to this output stream
	// @param `newLine` Indicates whether to end with a trailing `\\n`
	// @param `showWeight` Whether to indent each node according to their edge weights
	void pprint(ostream& out = cout, bool newLine = true, bool showWeight = true, int node = 1, string prefix = "", bool isLast = true) const {
		out << prefix << "|\n" << prefix << (isLast ? "`" : "|");
		if (showWeight) {
			out << string(parent[0][node].second, '-');
		}
		else {
			out << '-';
		}
		out << node << '\n';

		prefix += (isLast ? " " : "|");
		if (showWeight) prefix += string(parent[0][node].second, ' ');
		for (pair<int, T> child : children[node]) {
			pprint(out, false, showWeight, child.first, prefix, child == children[node].back());
		}

		if (newLine) out << '\n';
	}

	/************************************************
	 *                    UTILITY                   *
	 ************************************************/

	// O(1) Gets the size of the tree - the number of nodes
	size_t size() {
		return N;
	}

	// O(1) Gets the depth of `node`
	int getDepth(int node) {
		return depth[node];
	}
	// O(N) Gets the children of `node`
	vector<pair<int, T>> getChildren(int node) {
		return children[node];
	}

	// O(N) Gets the neighbours (children and parent) of `node`
	vector<pair<int, T>> getNeigbours(int node) {
		vector<pair<int, T>> neighbours = children[node];
		neighbours.push_back(parent[node]);
		return neighbours;
	}

	// O(N) Finds all descendants of `node`
	// @returns vector of nodes
	vector<int> getDescendants(int node = 1) {
		return preOrder(node); // lmao
	}

	// O(1) Determines if `node` is a leaf node (has no children)
	bool isLeaf(int node) {
		return children[node].empty();
	}

	// /O(N) Finds the leaves (nodes that have no children) in ascending order
	vector<int> getLeaves() {
		vector<int> leaves;
		for (int node = 1; node <= N; ++node) {
			if (isLeaf(node)) leaves.push_back(node);
		}
		return leaves;
	}

	// O(N) Finds the breath of the tree - the number of leaves
	int getBreadth() {
		int total = 0;
		for (int node = 1; node <= N; ++node) total += children[node].size() == 0;
		return total;
	}

	// O(log N) Add a node to the tree, and return the index of the new node */
	// @note Edge weight is ignored if the graph is unweighted
	int addNode(int par, T w = T(1)) {
		if (!isWeighted) w = T(1);

		++N;
		parent[0].push_back({ par, w });
		for (int d = 0; d < LOG_MAXN; ++d) {
			parent[d + 1].push_back({
				parent[d][parent[d][N].first].first,
				parent[d][N].second + parent[d][parent[d][N].first].second
				});
		}
		children[par].push_back({ N, w }); // make sure children[par] are still sorted
		children.push_back({});
		depth.push_back(depth[par] + 1);
		return N;
	}

	/************************************************
	 *                   PROPERTIES                 *
	 ************************************************/

	// [ O(N log N) ] Returns whether the tree is a binary tree - there are at most a left and right child for each node
	// @param `root` If specified, consider its connected component. Otherwise consider the whole tree
	bool isBinary(int root = 1) {
		for (int node = 1; node <= N; ++node) {
			if (children[node].size() > 2 && isAncestor(root, node)) return false;
		}
		return true;
	}

	// O(N log N) Finds whether the tree is a balanced tree - height of the left and right subtree differ by not more than 1
	// @param `root` If specified, consider its connected component. Otherwise consider the whole tree
	bool isBalanced(int root = 1) {
		int subtreeSize = getDescendants(root).size();
		int maxDepth = 0;
		for (int node = 1; node <= N; ++node) {
			if (isAncestor(root, node)) maxDepth = max(maxDepth, depth[node] - depth[root]);
		}
		return (1 << maxDepth) <= subtreeSize && subtreeSize < 2 * (1 << maxDepth);
	}

	/************************************************
	 *                   ALGORITHMS                 *
	 ************************************************/

	// O(log N) Finds the lowest common ancestor of `u` and `v`, using jump pointers with `depth` and `parent`
	int lca(int u, int v) {
		for (int d = LOG_MAXN; d >= 0; --d) {
			if (depth[parent[d][u].first] >= depth[v]) u = parent[d][u].first;
			if (depth[parent[d][v].first] >= depth[u]) v = parent[d][v].first;
		}

		for (int d = LOG_MAXN; d >= 0; --d) {
			if (parent[d][u].first != parent[d][v].first) {
				u = parent[d][u].first;
				v = parent[d][v].first;
			}
		}

		if (u == v) return u;
		assert(parent[0][u].first == parent[0][v].first);
		return parent[0][u].first;
	}

	// O(log N) Finds whether `ancestor` is actually an ancestor of `node`. Note that u node is its own ancestor
	bool isAncestor(int ancestor, int node) {
		return lca(ancestor, node) == ancestor;
	}

	// O(log n) Finds the `nth` parent of `node`
	// @note The `0th` parent of `node` is `node`, and `1` is the parent of itself
	int getParent(int node, int n = 1) {
		assert(n >= 0 && "Can't find a negative parent");
		n = min(n, N);
		int par = node;
		for (int d = 0; n; ++d, n /= 2) {
			if (n & 1) par = parent[d][par].first;
		}
		return par;
	}

	// O(log N) Finds the shortest distance between nodes `u` and `v`, by traversing `u` -> `lca(u, v)` -> `v`
	T dist(int u, int v) {
		int targetDepth = depth[lca(u, v)];
		T total = T(0);
		for (int d = LOG_MAXN; d >= 0; --d) {
			if (depth[parent[d][u].first] >= targetDepth) {
				total += parent[d][u].second;
				u = parent[d][u].first;
			}
			if (depth[parent[d][v].first] >= targetDepth) {
				total += parent[d][v].second;
				v = parent[d][v].first;
			}
		}

		return total;
	}

	// O(N) Find the path with the shortest distance, between nodes `u` and `v`
	// @returns u vector of nodes on the path, starting from `u` and ending with `v`
	vector<int> getPath(int u, int v) {
		vector<int> path;
		int mid = lca(u, v);
		for (int node = u; node != mid; node = parent[0][node].first) path.push_back(node);
		path.push_back(mid);

		int midLen = path.size();
		for (int node = v; node != mid; node = parent[0][node].first) path.push_back(node);
		reverse(path.begin() + midLen, path.end());

		return path;
	}

private:
	// O(N) Finds the diameter - the simple path with the maximum distance
	// @param `root` If specified, consider its connected component. Otherwise consider the whole tree
	// @returns Distance of the path, and a vector of nodes representing the path
	// @note Cannot have negative weight cycles
	pair<int, T> _diameterPath(pair<pair<int, int>, T>& best, int root = 1) {
		T furthest1 = T(0), furthest2 = T(0);
		int node1 = root, node2 = -1;

		for (pair<int, T> child : children[root]) {
			pair<int, T> path = _diameterPath(best, child.first);
			path.second += child.second;
			if (path.second > furthest1) {
				furthest2 = furthest1;
				node2 = node1;

				furthest1 = path.second;
				node1 = path.first;
			}
			else if (path.second > furthest2) {
				furthest2 = path.second;
				node2 = path.first;
			}
		}

		if (node1 != -1 && node2 != -1) {
			if (furthest1 + furthest2 > best.second) {
				best = { {node1, node2}, furthest1 + furthest2 };
			}
		}

		return { node1, furthest1 };
	}

public:
	// O(N) Returns the end points of the diameter (the longest simple path through the tree)
	pair<int, int> diameterPath() {
		pair<pair<int, int>, T> best = { {-1, -1}, T(0) };
		_diameterPath(best);
		return best.first;
	}

	// O(N) Returns the length of the diameter (the longest simple path through the tree) 
	int diameterDist() {
		pair<pair<int, int>, T> best = { {-1, -1}, T(0) };
		_diameterPath(best);
		return best.second;
	}

private:
	// O(N) DFS to find in-order traversal, assuming tree is a binary tree
	vector<int> _inOrder(vector<int>& traversal, int node = 1) {
		assert(children[node].size() <= 2 && "In-order traversals are only valid in binary trees");
		if (children[node].size() > 0) _inOrder(traversal, children[node][0].first);
		traversal.push_back(node);
		if (children[node].size() > 1) _inOrder(traversal, children[node][1].first);
		return traversal;
	}

	// O(N) DFS to find a tree traversal
	// @param `pushPre` Whether to add the node to the traversal before exploring children
	// @param `pushPost` Whether to add the node to the traversal after exploring children
	vector<int> _multiOrder(vector<int>& traversal, bool pushPre, bool pushPost, int node = 1) {
		if (pushPre) traversal.push_back(node);
		for (pair<int, T> child : children[node]) _multiOrder(traversal, pushPre, pushPost, child.first);
		if (pushPost) traversal.push_back(node);
		return traversal;
	}

public:
	// O(N) Finds in-order traversal, assuming tree is a binary tree
	// @param `root` If specified, consider its connected component. Otherwise consider the whole tree
	vector<int> inOrder(int root = 1) {
		vector<int> traversal;
		return _inOrder(traversal, root);
	}

	// O(N) Finds pre-order traversal
	// @param `root` If specified, consider its connected component. Otherwise consider the whole tree
	vector<int> preOrder(int root = 1) {
		vector<int> traversal;
		return _multiOrder(traversal, true, false, root);
	}

	// O(N) Finds post-order traversal
	// @param `root` If specified, consider its connected component. Otherwise consider the whole tree
	vector<int> postOrder(int root = 1) {
		vector<int> traversal;
		return _multiOrder(traversal, false, true, root);
	}

	// O(N) Finds euler tour traversal
	// @param `root` If specified, consider its connected component. Otherwise consider the whole tree
	vector<int> eulerTour(int root = 1) {
		vector<int> traversal;
		return _multiOrder(traversal, true, true, root);
	}
};

/************************************************
 *                    DISPLAY                   *
 ************************************************/
template<typename T> ostream& operator<<(ostream& out, Tree<T> tree) {
	tree.print(out);
	return out;
}