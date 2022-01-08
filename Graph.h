template<typename T> class Graph {
	/************************************************
	 *                 INITIALISATION               *
	 ************************************************/

private:
	// Number of nodes
	int V;
	// Whether edges are weighted
	// Unweighted edges are assigned unit weights
	bool isWeighted;
	// Whether edges are directed or bidirectional
	bool isDirected;
	// `inEdges[node]` is the multiset of edges leading into `node`
	vector<multiset<pair<int, T>>> inEdges;
	// `outEdges[node]` is the multiset of edges leading out of `node`
	vector<multiset<pair<int, T>>> outEdges;
	// A multiset of all edges
	multiset<pair<pair<int, int>, T>> edges;

	// O(V) Initialises all non-temporary variables.
	void init(int _V, bool _isWeighted = false, bool _isDirected = false) {
		V = _V;
		isDirected = _isDirected;
		isWeighted = _isWeighted;
		inEdges.resize(V + 1);
		outEdges.resize(V + 1);
	}

public:
	// O(V) Initialises an empty graph, assuming 1 indexed
	Graph(int _V, bool _isWeighted = false, bool _isDirected = false) {
		init(_V, _isWeighted, _isDirected);
	}

	// O(V) Initalises a graph from cin, assuming 1 indexed
	Graph(int _V, int M, istream& in = cin, bool _isWeighted = false, bool _isDirected = false) {
		init(_V, _isWeighted, _isDirected);

		for (int u, v, i = 0; i < M; ++i) {
			in >> u >> v;

			if (isWeighted) {
				T w;
				in >> w;
				addEdge(u, v, w);
			}
			else addEdge(u, v);
		}
	}

	// O(V + E) Initalises an unweighted graph from vector of edges, assuming 1 indexed
	Graph(int _V, vector<pair<int, int>> _edges, bool _isDirected = false) {
		init(_V, false, _isDirected);

		for (pair<int, int> edge : _edges) addEdge(edge.first, edge.second);
	}

	// O(V + E) Initalises a weighted graph from vector of edges, assuming 1 indexed
	Graph(int _V, vector<pair<pair<int, int>, T>> _edges, bool _isDirected = false) {
		init(_V, true, _isDirected);

		for (pair<pair<int, int>, T> edge : _edges) addEdge(edge.first.first, edge.first.second, edge.second);
	}


	/************************************************
	 *                    DISPLAY                   *
	 ************************************************/

	 // O(V + E) Displays the graph, showing the outwards edge connections of each edge in lexographic order
	 // @param `out` The string representation of the graph is piped to this output stream
	 // @param `newLine` Indicates whether to end with a trailing `\\n`
	void print(ostream& out = cout, bool newLine = true) {
		for (int u = 1; u <= V; ++u) {
			out << u << ":\n";
			for (pair<int, T> edge : outEdges[u]) {
				if (!isDirected && u > edge.first) continue;

				out << ' ' << edge.first;
				if (isWeighted) out << " (w = " << edge.second << ')';
				out << '\n';
			}
		}
		if (newLine) out << '\n';
	}

	/************************************************
	 *                   OPERATIONS                 *
	 ************************************************/

	 // O(V + E log E) Finds the union of 2 graphs
	 // @note May cause edge doubling
	Graph operator+(Graph<T> o) {
		assert(V == o.V && isDirected == o.isDirected && isWeighted == o.isWeighted && "Can't find union of graphs with different properties");

		Graph* output = new Graph(*this); // make a copy
		for (pair<pair<int, int>, T> edge : o.edges) output->addEdge(edge.first.first, edge.first.second, edge.second);
		return *output;
	}

	// O(1) Finds `V`, the number of nodes
	size_t size() {
		return V;
	}

	// O(log E + log V) Add the specified edge to the graph - this may create a double edge
	// @note Edge weight is ignored if the graph is unweighted
	void addEdge(int u, int v, T w = T(1)) {
		if (!isWeighted) w = 1;

		for (int rep = 0; rep <= int(!isDirected && u != v); ++rep) {
			outEdges[u].insert({ v, w });
			inEdges[v].insert({ u, w });
			edges.insert({ {u, v}, w });

			swap(u, v);
		}
	}

	// O(1) Adds a new node with no edges
	// @returns The index of this node
	int addNode() {
		++V;
		inEdges.push_back({});
		outEdges.push_back({});
		return V;
	}

	// O(log E) Finds the number of edges from `u` to `v`
	int numEdges(int u, int v) {
		int total = 0;
		for (auto iter = edges.lower_bound({ { u, v }, numeric_limits<T>::min() }); iter != edges.end(); ++iter) {
			if (iter->first.first != u || iter->first.second != v) break;
			++total;
		}

		return total;
	}

	// O(log E) Finds whether the specified edge (with any edge weight) is contained in the graph
	bool containsEdge(int u, int v) {
		auto iter = edges.lower_bound({ { u, v }, numeric_limits<T>::min() });
		if (iter == edges.end()) return false;
		return iter->first.first == u && iter->first.second == v;
	}

	// O(log E) Finds whether the specified edge is contained in the graph
	bool containsEdge(int u, int v, T w) {
		assert(isWeighted && "Can't find weighted edges in an unweighted graph");
		return edges.find({ {u, v}, w }) != edges.end();
	}

	// O(log E) If the specified edge (with any edge weight) is present, remove it, otherwise silently do nothing
	void removeEdge(int u, int v) {
		if (!containsEdge(u, v)) return;
		for (int rep = 0; rep <= int(!isDirected && u != v); ++rep) {
			outEdges[u].erase(
				outEdges[u].lower_bound({ v, numeric_limits<T>::min() }),
				outEdges[u].upper_bound({ v, numeric_limits<T>::max() })
			);
			inEdges[v].erase(
				inEdges[v].lower_bound({ u, numeric_limits<T>::min() }),
				inEdges[v].upper_bound({ u, numeric_limits<T>::max() })
			);
			edges.erase(
				edges.lower_bound({ { u, v }, numeric_limits<T>::min() }),
				edges.upper_bound({ { u, v }, numeric_limits<T>::max() })
			);

			swap(u, v);
		}
	}

	// O(log E) If the specified edge is present, remove it, otherwise silently do nothing
	void removeEdge(int u, int v, T w) {
		if (!containsEdge(u, v, w)) return;
		for (int rep = 0; rep <= int(!isDirected && u != v); ++rep) {
			//auto iter = outEdges[u].lower_bound({ v, w });
			outEdges[u].erase({ v, w });
			inEdges[v].erase({ u, w });
			edges.erase({ { u, v }, w });
			swap(u, v);
		}
	}

	// O(E) Finds the edges that travel into `node`
	// @returns A multiset of (neighbour, weight) pairs
	multiset<pair<int, T>> getInEdges(int node) {
		return inEdges[node];
	}

	// O(E) Finds the edges that travel out of `node`
	// @returns A multiset of (neighbour, weight) pairs
	multiset<pair<int, T>> getOutEdges(int node) {
		return outEdges[node];
	}

	// O(E log E) Finds the edges that travel into OR out of `node`
	// @returns A multiset of (neighbour, weight) pairs
	multiset<pair<int, T>> getEdges(int node) {
		multiset<pair<int, T>> output = inEdges[node];
		if (isDirected) {
			for (pair<int, T> edge : outEdges[node]) output.insert(edge);
		}
		return output;
	}

	// O(V + E) Finds the multiset of all edges
	// @param `root` If specified, consider its connected component. Otherwise consider the whole graph
	// @returns A multiset of ((u, v), weight) pairs
	multiset<pair<pair<int, int>, T>> getAllEdges(int root = -1) {
		if (root == -1) return edges;

		multiset<pair<pair<int, int>, T>> output;
		for (int u : getComponent(root)) {
			for (pair<int, T> edge : outEdges[u]) {
				if (!isDirected && edge.first > u) break;
				output.insert({ { u, edge.first }, edge.second });
			}
		}
		return edges;
	}

	// O(1) Gets the in degree of `node`
	size_t inDegree(int node) {
		return inEdges[node].size();
	}

	// O(1) Gets the out degree of `node`
	size_t outDegree(int node) {
		return outEdges[node].size();
	}

	// O(1) Gets the total degree of `node`
	size_t degree(int node) {
		if (isDirected) return (inEdges[node].size() + outEdges[node].size()) / 2;
		return inEdges[node].size() + outEdges[node].size();
	}

	/************************************************
	 *                    PROPERTIES
	 ***********************************************/

	 // O(1) Gets whether the graph has directed edges
	bool isDirectedGraph() {
		return isDirected;
	}

	// O(1) Gets whether the graph has weighted edges
	bool isWeightedGraph() {
		return isWeighted;
	}

	// O(V + E) Determines whether there exists edges to itself i.e. `u -> u`
	// @param `root` If specified, consider its connected component. Otherwise consider the whole graph
	bool hasSelfEdges(int root = -1) {
		for (int node : getComponent(root)) {
			if (containsEdge(node, node)) return true;
		}
		return false;
	}

	// O(V + E) Determines whether there exists multiple edges between any 2 nodes
	// @param `root` If specified, consider its connected component. Otherwise consider the whole graph
	bool hasDoubleEdges(int root = -1) {
		multiset<pair<pair<int, int>, T>> currEdges = getAllEdges(root);
		if (edges.empty()) return false;
		for (auto it1 = edges.begin(), it2 = ++edges.begin(); it2 != edges.end(); ++it1, ++it2) {
			if (it1->first == it2->first) return true;
		}
		return false;
	}

	// O(V + E) Determines whether the graph is a simple graph
	bool isSimpleGraph() {
		return !hasSelfEdges() && !hasDoubleEdges();
	}

	// O(V) Determines whether the graph has a cycle
	bool hasCycle() {
		return Kahn().size() != V;
	}

	// O(V) Determines whether the graph is a directed acrylic graph
	bool isDAG() {
		return isDirected && !hasCycle();
	}

	// O(V + E) Determines hether the graph is a tree
	// @param `root` If specified, consider its connected component. Otherwise consider the whole graph
	bool isTree(int root = -1) {
		if (hasSelfEdges(root)) return false;
		if (hasDoubleEdges(root)) return false;

		vector<int> component = getComponent(root);
		int numEdges = 0, numLeaves = 0, numRoots = 0;
		for (int node : component) {
			for (pair<int, T> edge : outEdges[node]) {
				if (!isDirected && edge.first > node) break;
				else ++numEdges;
			}
			numLeaves += outEdges[node].empty();
			numRoots += inEdges[node].empty();
		}

		if (isDirected && numLeaves > 1 && numRoots > 1) return false;
		return numEdges == component.size() - 1; // E == V - 1
	}

	// O(V + E) Determines hether the graph is a forest
	bool isForest() {
		for (vector<int>& component : getComponents()) {
			if (!isTree(component[0])) return false;
		}
		return true;
	}

	/************************************************
	 *        ALGORITHMS (COMPLEXITLY CLASS P)
	 ***********************************************/

private:
	// O(V + E log V) Standard Dijkstra algorithm - shortest path from `node` to all other nodes for non-negative edge weights
	// @param `backwards` Indicates whether to consider `inEdges` or `outEdges`
	// @note If the graph is undirected, this becomes equivalent to a bfs (but slightly slower)
	// @note If negative edge weights encountered, switches to `BellmanFord()`
	void Dijkstra(int node, vector<int>& prevNode, vector<T>& dist, bool backwards = false) {
		dist[node] = T(0);
		prevNode[node] = node;
		priority_queue<pair<T, int>, vector<pair<T, int>>, greater<pair<T, int>>> pq;
		pq.push({ 0, node });

		while (!pq.empty()) {
			int currNode = pq.top().second;
			T cost = pq.top().first;
			pq.pop();
			if (cost > dist[currNode]) continue;

			for (pair<int, T> child : (backwards ? inEdges : outEdges)[currNode]) {
				if (child.second < 0) {
					// negative edge weight, try `BellmanFord()` isntead
					fill(prevNode.begin(), prevNode.end(), 0);
					fill(dist.begin(), dist.end(), numeric_limits<T>::max());
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
	void BellmanFord(int node, vector<int>& prevNode, vector<int>& dist, bool backwards = false) {
		for (int rep = 0; rep < V; ++rep) {
			for (pair<pair<int, int>, T> edge : edges) {
				if (backwards) swap(edge.first.first, edge.first.second);

				if (dist[edge.first.first] == numeric_limits<T>::max()) continue;
				T newCost = dist[edge.first.first] + edge.second;
				if (newCost < dist[edge.first.second]) {
					// last rep and the minimum distance can still be reduced
					assert(rep == V - 1 && "negative weight cycle");

					dist[edge.first.second] = newCost;
					prevNode[edge.first.second] = edge.first.first;
				}
			}
		}
	}

	// O(V^3 + E) Standard Floyd-Warshall algorithm - shortest path between every 2 pairs of nodes
	void FloydWarshall(vector<vector<T>>& dists) {
		for (int node = 1; node <= V; ++node) dists[node][node] = 0;
		for (pair<pair<int, int>, T>& edge : edges) {
			int u = edge.first.first, v = edge.first.second;
			dists[u][v] = min(dists[u][v], edge.second);
		}

		for (int mid = 1; mid <= V; ++mid) {
			for (int u = 1; u <= V; ++u) {
				for (int v = 1; v <= V; ++v) {
					if (max(dists[u][mid], dists[mid][v]) != numeric_limits<T>::max()) dists[u][v] = min(dists[u][v], dists[u][mid] + dists[mid][v]);
				}
			}
		}
	}

public:
	// O(E log V), but O(VE) if negative edge weights. Finds the shortest path from `u` to `v`
	// @returns Distance of the path, and a vector of nodes representing the path
	// @note Component containing `u` and `v` cannot have negative weight cycles
	pair<T, vector<int>> shortestPath(int u, int v) {
		vector<int> prevNode = vector<int>(V + 1, 0);
		vector<T> dist = vector<T>(V + 1, numeric_limits<T>::max());
		Dijkstra(u, prevNode, dist);

		vector<int> path;
		if (dist[v] != numeric_limits<T>::max()) {
			for (int node = v; prevNode[node] != node; node = prevNode[node]) path.push_back(node);
			path.push_back(u);

			reverse(path.begin(), path.end());
		}

		return { dist[v], path };
	}

	// O(E log V), but O(VE) if negative edge weights. Finds the shortest path from `u` to `v`
	// @returns Distance of the path
	// @note Component containing `u` and `v` cannot have negative weight cycles
	T shortestDist(int u, int v) {
		vector<int> prevNode = vector<int>(V + 1, 0);
		vector<T> dist = vector<T>(V + 1, numeric_limits<T>::max());
		Dijkstra(u, prevNode, dist);

		return dist[v];
	}

	// O(V^3 + E) Finds the minimum distance between every pair of nodes
	// @returns A vector of vectors, with `dist[u][v]` being the minimim distance on the path from `u` to `v`
	vector<vector<T>> allShortestDist() {
		vector<vector<int>> dists = vector<vector<int>>(V + 1, vector<int>(V + 1, numeric_limits<T>::max()));
		FloydWarshall();
		return dists;
	}

	// O(E log V), but O(VE) if negative edge weights. Finds the eccentricity - the greatest distance between `root` and any other node in the same component
	// @note Component containing `root` cannot have negative weight cycles
	T eccentricity(int root) {
		// uh should i return the furthest node too?
		vector<int> prevNode = vector<int>(V + 1, 0);
		vector<int> dist = vector<int>(V + 1, numeric_limits<T>::max());
		Dijkstra(root, prevNode, dist);
		T furthest = 0;
		for (int node = 1; node <= V; ++node) {
			if (dist[node] != numeric_limits<T>::max()) furthest = max(furthest, dist[node]);
		}
		return furthest;
	}

	// O(E log V), but O(VE) if negative edge weights. Finds the diameter - the simple path with the maximum distance
	// @param `root` If specified, consider its connected component. Otherwise consider the whole graph
	// @returns Distance of the path, and a vector of nodes representing the path
	// @note Cannot have negative weight cycles
	vector<int> diameter(int root = -1) {
		if (root == -1) {
			vector<int> output;
			for (vector<int> component : getComponents()) {
				vector<int> curr = diameter(component[0]);
				if (curr.size() > output.size()) output = curr;
			}
			return output;
		}

		vector<int> prevNode = vector<int>(V + 1, 0);
		vector<int> dist = vector<int>(V + 1, numeric_limits<T>::max());
		Dijkstra(root, prevNode, dist, true);

		int furthestDist1 = numeric_limits<T>::min();
		int furthestNode1 = root;
		for (int node = 1; node <= V; ++node) {
			if (dist[node] != numeric_limits<T>::max() && dist[node] > furthestDist1) {
				furthestDist1 = dist[node];
				furthestNode1 = node;
			}
		}

		prevNode = vector<int>(V + 1, 0);
		dist = vector<int>(V + 1, numeric_limits<T>::max());
		Dijkstra(furthestNode1, prevNode, dist, false);

		int furthestDist2 = numeric_limits<T>::min();
		int furthestNode2 = root;
		for (int node = 1; node <= V; ++node) {
			if (dist[node] != numeric_limits<T>::max() && dist[node] > furthestDist2) {
				furthestDist2 = dist[node];
				furthestNode2 = node;
			}
		}

		vector<int> path;
		for (int node = furthestNode2; node != prevNode[node]; node = prevNode[node]) {
			path.push_back(node);
		}
		path.push_back(furthestNode1);
		return path;
	}

private:
	// Amortised O(log V)
	// @returns The parent node of the current disjoint set
	int _getParent(int node, vector<int>& parent) {
		if (parent[node] == node) return node;
		return parent[node] = _getParent(parent[node], parent);
	}

	// Amortised O(log V) Merges the disjoint sets of `a` and `b`
	void _merge(int a, int b, vector<int>& parent) {
		parent[_getParent(a, parent)] = _getParent(b, parent);
	}

public:
	// Amortised O(E log V) Standard Kruskal's algorithm - finds the minimum spanning tree, using union find (DSU)
	// @returns Graph object of MST
	Graph Kruskal() {
		vector<int> parent = vector<int>(V + 1);
		iota(parent.begin(), parent.end(), 0);

		vector<pair<T, pair<int, int>>> sortedEdges;
		for (pair<pair<int, int>, T> edge : edges) {
			sortedEdges.push_back({ edge.second, edge.first });
			_merge(edge.first.first, edge.first.second, parent);
		}
		sort(sortedEdges.begin(), sortedEdges.end());

		Graph output(V);
		for (pair<T, pair<int, int>>& edge : sortedEdges) {
			int parentU = _getParent(edge.second.first, parent), parentV = _getParent(edge.second.second, parent);
			if (parentU != parentV) {
				output.addEdge(parentU, parentV, edge.first);
				_merge(parentU, parentV, parent);
			}
		}
		return output;
	}
	Graph MST() {
		return Kruskal();
	}

private:
	// O(V + E) Standard Tarjan algorithm - find strongly connected components
	// @note If the graph is bidirectional, this finds all components
	int Tarjan(int u, int index, vector<vector<int>>& components, vector<int>& s, vector<bool>& seen, vector<int>& getIndex, vector<int>& lowLink) {
		getIndex[u] = lowLink[u] = index;
		++index;
		s.push_back(u);
		seen[u] = true;

		for (pair<int, T> edge : outEdges[u]) {
			int v = edge.first;
			if (getIndex[v] == -1) {
				index = Tarjan(v, index, components, s, seen, getIndex, lowLink);
				lowLink[u] = min(lowLink[u], lowLink[v]);
			}
			else if (seen[v]) {
				lowLink[u] = min(lowLink[u], getIndex[v]);
			}
		}

		if (getIndex[u] == lowLink[u]) {
			components.push_back({});

			int node;
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
	// O(V + E) Finds strongly connected components. In a bidirectional graph, this returns all connected components
	// @returns Vector of components, where each component is a vector of nodes
	vector<vector<int>> getComponents() {
		vector<bool> seen = vector<bool>(V + 1, false);
		vector<int> getIndex = vector<int>(V + 1, -1);
		vector<int> lowLink = vector<int>(V + 1, 0);
		vector<vector<int>> components;
		vector<int> s;

		for (int node = 1; node <= V; ++node) {
			if (getIndex[node] == -1) Tarjan(node, 1, components, s, seen, getIndex, lowLink);
		}
		return components;
	}

	// O(V + E) Finds a strongly connected component. In a bidirectional graph, this returns all connected components
	// @param `root` If specified, consider its connected component. Otherwise consider the whole graph
	// @returns The connected components reachable from root, where each component is represented as a vector of nodes
	vector<int> getComponent(int root = -1) {
		if (root == -1) {
			// return all nodes
			vector<int> output(V);
			iota(output.begin(), output.end(), 1);
			return output;
		}

		vector<bool> seen = vector<bool>(V + 1, false);
		vector<int> getIndex = vector<int>(V + 1, -1);
		vector<int> lowLink = vector<int>(V + 1, 0);
		vector<vector<int>> components;
		vector<int> s;

		Tarjan(root, 1, components, s, seen, getIndex, lowLink);
		assert(components.size() == 1 && "Sanity check!");
		return components[0];
	}

	// O(V + E) Copy the component at root to a new graph, maintaining the same graph properties and the same node numberings
	Graph<T> getSubgraph(int root) {
		Graph output<T>(V, isWeighted, isDirected);
		for (pair<pair<int, int>, T> edge : getAllEdges(root)) output.addEdge(edge.first.first, edge.first.second, edge.second);
		return output;
	}

private:
	// O(V + E) Standard DFS Bridge finding algorithm
	int _getBridges(int u, int parent, int t, vector<pair<pair<int, int>, T>>& bridges, vector<bool>& seen, vector<int>& minTime, vector<int>& entryTime) {
		seen[u] = true;
		minTime[u] = entryTime[u] = ++t;

		for (pair<int, T> child : outEdges[u]) {
			int v = child.first;
			T w = child.second;
			if (v == parent) continue;
			if (seen[v]) {
				minTime[u] = min(minTime[u], entryTime[v]);
			}
			else {
				t = _getBridges(v, u, t, bridges, seen, minTime, entryTime);
				minTime[u] = min(minTime[u], minTime[v]);

				// bridge iff minTime[childNode] > entryTime[node] and the edge is not a double edge
				if (minTime[v] > entryTime[u]) {
					if (numEdges(u, v) == 1) bridges.push_back({ { u, v }, w });
				}
			}
		}
		return t;
	}

public:
	// O(V + E) Find all bridges
	// @param `root` If specified, consider its connected component. Otherwise consider the whole graph
	// @return A vector of edges, where each edge is represented by ((u, v), w)
	vector<pair<pair<int, int>, T>> getBridges(int root = -1) {
		vector<bool> seen = vector<bool>(V + 1, false);
		vector<int> minTime = vector<int>(V + 1, 0), entryTime = vector<int>(V + 1, 0);

		vector<pair<pair<int, int>, T>> bridges;
		if (root == -1) {
			for (int node = 1; node <= V; ++node) {
				if (!seen[node]) _getBridges(node, -1, 0, bridges, seen, minTime, entryTime);
			}
		}
		else {
			_getBridges(root, -1, 0, bridges, seen, minTime, entryTime);
		}
		return bridges;
	}

private:
	// O(V + E log E) Greedy colouring
	void _greedyColouring(int node, vector<int>& colours) {
		if (colours[node] != -1) return;

		// calculate mex of neighbours colours
		unordered_set<int> usedColours; // BEWARE! ORAC IS STUPID!
		for (pair<int, T> edge : getEdges(node)) usedColours.insert(colours[edge.first]);
		for (colours[node] = 0; usedColours.find(colours[node]) != usedColours.end(); ++colours[node]);

		for (pair<int, T> edge : getEdges(node)) _greedyColouring(edge.first, colours);
	}

public:
	// O(V + E log E) Performs a greedy colouring of nodes - assigning the minimum number of colours such that adjacent nodes have different colours
	// @returns A vector where `0 <= colour[node] < max_colours`
	// @note Ignore colour[0]
	vector<int> greedyColouring() {
		vector<int> colours(V + 1, -1);
		for (int node = 1; node <= V; ++node) _greedyColouring(node, colours);
		return colours;
	}

private:
	// O(V) Standard Kahn's algorithm - topological sort
	vector<int> Kahn() {
		vector<int> depths = vector<int>(V + 1, 0);

		vector<int> q, topSort;
		for (int node = 1; node <= V; ++node) {
			depths[node] = inEdges[node].size();
			if (depths[node] == 0) q.push_back(node);
		}
		while (!q.empty()) {
			int node = q.back();
			q.pop_back();
			topSort.push_back(node);

			for (pair<int, T> edge : outEdges[node]) {
				depths[edge.first]--;
				if (depths[edge.first] == 0) q.push_back(edge.first);
			}
		}
		return topSort;
	}

	// O(V) Stanard DFS topsort algorithm
	void _dfsTopSort(int node, vector<int>& topSort, vector<int>& state) {
		if (state[node] == 2) return;
		assert(state[node] == 1 && "Graph contains cycle");

		state[node] = 1;
		for (pair<int, T> edge : outEdges[node]) _dfsTopSort(edge.first, topSort, state);
		state[node] = 2;
		topSort.push_back(node);
	}

public:
	// O(V) Find a topological sort - either Kahn's algorithm of DFS topsort
	// @note Graph cannot be bidirectional or cyclical
	// @returns vector of nodes in topological sort order
	vector<int> getTopSort(bool doKahn = false) {
		assert(isDirected && "Can't get a topological sort for a bidirectional graph");

		if (doKahn) {
			vector<int> topSort = Kahn();
			assert(topSort.size() != V && "Graph contains cycle");
			return topSort;
		}
		else {
			vector<int> state = vector<int>(V + 1, 0);

			vector<int> topSort;
			for (int node = 1; node <= V; ++node) _dfsTopSort(node, topSort, state);
			reverse(topSort.begin(), topSort.end());
			return topSort;
		}
	}

	// TODO: is bipartite graph?
	// TODO: find cycles
	// TODO: bipartite matching
	// TODO: edge matching
	// TODO: max flow
	// TODO: see: https://usaco.guide/adv/

	/************************************************
	 *       ALGORITHMS (COMPLEXITLY CLASS NP)      *
	 ************************************************/

	 // TODO: hamiltonian path / tour
	 // TODO: euler path / tour
	 // TODO: covering set
	 // TODO: optimal colouring
	 // TODO: planar embedding?
};

/************************************************
 *                    DISPLAY                   *
 ************************************************/
template<typename T> ostream& operator<<(ostream& out, Graph<T> graph) {
	graph.print(out);
	return out;
}