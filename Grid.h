// @returns The number of characters of the printed representation of `obj` */
template<typename T> int reprLen(T obj) {
	stringstream s;
	s << obj;
	return s.str().size();
}

template<typename T> class Grid {
	/************************************************
	 *                INITIALISATION                *
	 ************************************************/

private:
	// Number of rows
	int R = 0;
	// Number of columns;
	int C = 0;
	// Grid
	vector<vector<T>> grid;
	// Prefix sum of grid
	vector<vector<T>> prefix;
	// Whether the graph hasn't been modified after the last precomp
	bool hasPrecomp = false;

	// O(RC) Initialises all variables
	void init(int _R, int _C) {
		R = _R;
		C = _C;
		grid = vector<vector<T>>(R);
		for (int r = 0; r < R; ++r) grid[r] = vector<T>(C);
		prefix = vector<vector<T>>(R + 1);
		for (int r = 0; r < R + 1; ++r) prefix[r] = vector<T>(C + 1);
	}

public:
	// O(1) initialises an empty grid
	Grid() {

	}

	// O(RC) Initialises a grid filled with `val`
	Grid(int _R, int _C, T val = T(0)) {
		init(_R, _C);
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) {
				grid[r][c] = val;
			}
		}
	}

	// O(RC) Initialises grid from cin, assuming 0 indexed
	Grid(int _R, int _C, istream& in = cin) {
		init(_R, _C);
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) {
				in >> grid[r][c];
			}
		}
		buildPrefix();
	}

	// O(RC) Initialises grid from a 2d array, assuming 0 indexed
	template<size_t rows, size_t cols> Grid(T(&_grid)[rows][cols]) {
		init(sizeof(_grid) / sizeof(_grid[0]), sizeof(_grid[0]) / sizeof(T));

		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) {
				grid[r][c] = _grid[r][c];
			}
		}
		buildPrefix();
	}

	// O(RC) Initialises grid from a 2d vector, assuming 0 indexed
	Grid(vector<vector<T>> _grid) {
		init((int)_grid.size(), (int)_grid[0].size());
		for (int r = 0; r < R; ++r) {
			assert(_grid[r].size() == grid[0].size() && "Invalid grid shape");
			grid[r] = _grid[r];
		}
		buildPrefix();
	}

	// Precomp prefix sum
	// @note automatically called during `init()`
	void buildPrefix() {
		hasPrecomp = true;
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) {
				prefix[r + 1][c + 1] = grid[r][c] + (prefix[r + 1][c] + prefix[r][c + 1] - prefix[r][c]);
			}
		}
	}

	/************************************************
	 *                    DISPLAY                   *
	 ************************************************/

	// O(RC) Displays the grid
	// @param `out` The string representation of the graph is piped to this output stream
	// @param `newLine` Indicates whether to end with a trailing `\\n`
	void print(ostream& out = cout, bool newLine = true) {
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) {
				out << grid[r][c] << ' ';
			}
			out << '\n';
		}
		if (newLine) out << '\n';
	}

	// O(RC) Displays the grid
	// @param `out` The string representation of the graph is piped to this output stream
	// @param `newLine` Indicates whether to end with a trailing `\\n`
	// @param `spacing` Indicates how large each cell is (use spacing = -1 for automatic spacing)
	void pprint(ostream& out = cout, bool newLine = true, int spacing = -1) {
		if (spacing == -1) {
			// find max repr length of grid values
			for (int r = 0; r < R; ++r) {
				for (int c = 0; c < C; ++c) spacing = max(spacing, (reprLen(grid[r][c])+3) / 4);
			}
		}

		string sep = "+", space = "|";
		for (int r = 0; r < R; ++r) {
			sep += string(4 * spacing + 1, '-') + '+';
			space += string(4 * spacing + 1, ' ') + '|';
		}
		sep += '\n';
		space += '\n';
		
		for (int r = 0; r < R; ++r) {
			out << sep;

			for (int rep = 0; rep < spacing; ++rep) out << space;

			out << '|';
			for (int c = 0; c < C; ++c) {
				int paddingL = max(0, 2 * spacing - reprLen(grid[r][c]) / 2);
				int paddingR = max(0, 2 * spacing - (reprLen(grid[r][c]) - 1) / 2);

				out << string(paddingL, ' ') << grid[r][c] << string(paddingR, ' ') << '|';
			}
			out << '\n';

			for (int rep = 0; rep < spacing; ++rep) out << space;
		}
		out << sep;

		if (newLine) out << '\n';
	}

	/************************************************
	 *             COMPARISON FUNCTIONS             *
	 ************************************************/

	// O(1) Filter generator, to determine whether cell values are equal to `val`
	// @returns A function which acts to "filter" cells 
	function<bool(int, int)> compEQ(T val) {
		return [&, val](int r, int c) {
			if (validCoord(r, c)) return getVal(r, c) == val;
			return false;
		};
	}

	// O(1) Filter generator, to determine whether cell values are not equal to `val`
	// @returns A function which acts to "filter" cells 
	function<bool(int, int)> compNE(T val) {
		return [&, val](int r, int c) {
			if (validCoord(r, c)) return getVal(r, c) != val;
			return false;
		};
	}

	// O(1) Filter generator, to determine whether cell values are greater than `val`
	// @returns A function which acts to "filter" cells 
	function<bool(int, int)> compGT(T val) {
		return [&, val](int r, int c) {
			if (validCoord(r, c)) return getVal(r, c) > val;
			return false;
		};
	}

	// O(1) Filter generator, to determine whether cell values are greater or equal to `val`
	// @returns A function which acts to "filter" cells 
	function<bool(int, int)> compGE(T val) {
		return [&, val](int r, int c) {
			if (validCoord(r, c)) return getVal(r, c) >= val;
			return false;
		};
	}

	// O(1) Filter generator, to determine whether cell values are less than `val`
	// @returns A function which acts to "filter" cells 
	function<bool(int, int)> compLT(T val) {
		return [&, val](int r, int c) {
			if (validCoord(r, c)) return getVal(r, c) < val;
			return false;
		};
	}

	// O(1) Filter generator, to determine whether cell values are less or equal to `val`
	// @returns A function which acts to "filter" cells 
	function<bool(int, int)> compLE(T val) {
		return [&, val](int r, int c) {
			if (validCoord(r, c)) return getVal(r, c) <= val;
			return false;
		};
	}

	// O(1) Filter generator, which returns a truthy function - to consider any cell
	function<bool(int, int)> compAny() {
		return [](int r, int c) {
			return true;
		};
	}

	/************************************************
	 *             COORDINATES & VALUES             *
	 ************************************************/

	// O(1) Gets R - the number of rows
	int getR() {
		return R;
	}

	// O(1) Gets C - the number of columns
	int getC() {
		return C;
	}

	// O(1) Gets the size of the grid
	size_t size() {
		return R * C;
	}

	// O(1)
	// @returns whether a coordinate is within the grid
	bool validCoord(int r, int c) {
		return 0 <= r && r < R && 0 <= c && c < C;
	}
	bool validCoord(pair<int, int> coord) {
		return validCoord(coord.first, coord.second);
	}

	// O(RC)
	// @returns Whether there exists a coordinate that satisfies `isValid`
	template<typename Condition> bool containsVal(Condition isValid) {
		return findVal(isValid).first != -1;
	}

	// O(RC)
	// @returns The first (lexigraphically smallest) coordinate which satisfies `isValid`, otherwise returns {-1, -1}
	template<typename Condition> pair<int, int> findVal(Condition isValid) {
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) {
				if (isValid(r, c)) return { r, c };
			}
		}
		return { -1, -1 };
	}

	// O(RC)
	// @returns A vector of coordinates which satisfies `isValid`
	template<typename Condition> vector<pair<int, int>> findAllVal(Condition isValid) {
		vector<pair<int, int>> output;
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) {
				if (isValid(r, c)) output.push_back({ r, c });
			}
		}
		return output;
	}

	// O(RC)
	// @returns the number of cells that satisfy `isValid`
	template<typename Condition> int count(Condition isValid) {
		int total = 0;
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) total += (bool)isValid(r, c);
		}
		return total;
	}

	// O(1)
	// @returns If the provided coordinate is valid, set its value to `val`, otherwise silently do nothing
	void setVal(int r, int c, T val) {
		if (validCoord(r, c)) {
			hasPrecomp = false;
			grid[r][c] = val;
		}
	}
	void setVal(pair<int, int> coord, T val) {
		setVal(coord.first, coord.second, val);
	}

	// O(area of rectangular region) <= O(RC) Sets the value of all cells within a rectangular region (inclusive)
	// @note The rectangular region will be capped so it stays within the grid
	void setVals(int r1, int c1, int r2, int c2, T val) {
		hasPrecomp = false;
		if (r1 > r2) swap(r1, r2);
		if (c1 > c2) swap(c1, c2);
		r1 = max(0, r1);
		c1 = max(0, c1);
		r2 = min(R - 1, r2);
		c2 = max(C - 1, c2);

		for (int r = r1; r <= r2; ++r) {
			for (int c = c1; c <= c2; ++c) grid[r][c] = val;
		}
	}
	void setVals(pair<int, int> coord1, pair<int, int> coord2, T val) {
		setVals(coord1.first, coord1.second, coord2.first, coord2.second, val);
	}
	
	// O(1) If the provided coordinate is valid, increment its value by `val`
	void incrVal(int r, int c, T val) {
		if (validCoord(r, c)) {
			hasPrecomp = false;
			grid[r][c] += val;
		}
	}
	void incrVal(pair<int, int> coord, T val) {
		incrVal(coord.first, coord.second, val);
	}

	// O(area of rectangular region) <= O(RC) Increments the value of all cells within a rectangular region (inclusive)
	// @note The rectangular region will be capped so it stays within the grid
	void incrVals(int r1, int c1, int r2, int c2, T val) {
		hasPrecomp = false;
		if (r1 > r2) swap(r1, r2);
		if (c1 > c2) swap(c1, c2);
		r1 = max(0, r1);
		c1 = max(0, c1);
		r2 = min(R - 1, r2);
		c2 = max(C - 1, c2);

		for (int r = r1; r <= r2; ++r) {
			for (int c = c1; c <= c2; ++c) grid[r][c] += val;
		}
	}
	void incrVals(pair<int, int> coord1, pair<int, int> coord2, T val) {
		incrVals(coord1.first, coord1.second, coord2.first, coord2.second, val);
	}

	// O(1)
	// @returns The value at the coordinate (r, c) if it is valid, otherwise return `defaultVal`
	T getVal(int r, int c, T defaultVal = T()) {
		if (validCoord(r, c)) return grid[r][c];
		return defaultVal;
	}
	T getVal(pair<int, int> coord, T defaultVal = T()) {
		return getVal(coord.first, coord.second, defaultVal);
	}

	// O(1) if precomp, otherwise O(RC)
	// @returns The sum of all cells within a rectangular region (inclusive)
	// @note This is only garunteed to be correct after initPrefix() is called, and the grid is not modified after that
	// @note The rectangular region will be capped so it stays within the grid
	T sumVals(int r1, int c1, int r2, int c2) {
		if (r1 > r2) swap(r1, r2);
		if (c1 > c2) swap(c1, c2);
		r1 = max(0, r1);
		c1 = max(0, c1);
		r2 = min(R - 1, r2);
		c2 = max(C - 1, c2);

		if (hasPrecomp) {
			return prefix[r2 + 1][c2 + 1] + prefix[r1][c1] - prefix[r2 + 1][c1] - prefix[r1][c2 + 1];
		}
		else {
			T total = 0;
			for (int r = r1; r <= r2; ++r) {
				for (int c = c1; c <= c2; ++c) total += grid[r][c];
			}
			return total;
		}
	}
	T sumVals(pair<int, int> coord1, pair<int, int> coord2) {
		return sumVals(coord1.first, coord1.second, coord2.first, coord2.second);
	}

	/************************************************
	 *                 OPERATIONS                   *
	 ************************************************/

	// O(RC)
	bool operator==(Grid o) {
		if (R != o.getR() || C != o.getC()) return false;
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) {
				if (grid[r][c] != o.getVal(r, c)) return false;
			}
		}
		return true;
	}

	/************************************************
	 *               TRANSFORMATIONS                *
	 ************************************************/

	// O(RC) Flips the grid over its primary diagonal (top left to bottom right)
	// @returns grid object of the modified grid
	Grid<T> flipPrimaryDiag() {
		Grid<T> output(C, R, 0);
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) output.setVal(c, r, grid[r][c]);
		}
		return output;
	}

	// O(RC) Flips the grid over its secondary diagonal (top right to bottom left)
	// @returns grid object of the modified grid
	Grid<T> flipSecondaryDiag() {
		return rot180().flipPrimaryDiag();
	}

	// O(RC) Flips the order of rows - first row becomes last row and vice versa
	// @returns grid object of the modified grid
	Grid<T> flipRows() {
		Grid<T> output(R, C, 0);
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) output.setVal(R - r - 1, c, grid[r][c]);
		}
		return output;
	}

	// O(RC) Flips the order of columns - first column becomes last row and vice versa
	// @returns grid object of the modified grid
	Grid<T> flipCols() {
		Grid<T> output(R, C, 0);
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) output.setVal(r, C - c - 1, grid[r][c]);
		}
		return output;
	}

	// O(RC) Rotates the grid clockwise/anticlockwise by 180 degrees, about its center
	// @returns grid object of the modified grid
	Grid<T> rot180() {
		return flipRows().flipCols();
	}

	// O(RC) Rotates the grid clockwise by 90 degrees, about its center
	// @returns grid object of the modified grid
	Grid<T> rotClockwise() {
		return flipPrimaryDiag().flipCols();
	}

	// O(RC) Rotates the grid anticlockwise by 90 degrees, about its center
	// @returns grid object of the modified grid
	Grid<T> rotAnticlockwise() {
		return flipPrimaryDiag().flipRows();
	}

	// O(RC) Replaces every instance of `prevVal` with `newVal`
	// @returns grid object of the modified grid`
	Grid<T> replace(T prevVal, T newVal) {
		Grid<T> *output = new Grid<T>(*this); // make a copy
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) {
				if (output->getVal(r, c) == prevVal) output->setVal(r, c, newVal);
			}
		}
		return *output;
	}

	/************************************************
	 *           GROUPINGS OF COORDINATES           *
	 ************************************************/

	// O(C)
	// @returns All coordinates in a particular row, sorted in increasing order
	vector<pair<int, int>> getRowCoords(int r) {
		vector<pair<int, int>> output;
		for (int c = 0; c < C; ++c) output.push_back({ r, c });
		return output;
	}

	// O(R)
	// @returns All coordinates in a particular column, sorted in increasing order
	vector<pair<int, int>> getColCoords(int c) {
		vector<pair<int, int>> output;
		for (int r = 0; r < R; ++r) output.push_back({ r, c });
		return output;
	}

	// O(RC)
	// @returns All coordinates in a vector, sorted by lexicographic order
	vector<pair<int, int>> getAllCoords() {
		//return getSomeCoords(compAny);
		vector<pair<int, int>> output;
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) {
				output.push_back({ r, c });
			}
		}
		return output;
	}

	// O(RC) Finds all coordinates that satisfies `isValid`
	// @param `dirs` the valid directions to move to adjacent nodes
	// @returns aAl coordinates in a vector, sorted by lexicographic order
	vector<pair<int, int>> getSomeCoords(function<bool(int, int)> isValid) {
		vector<pair<int, int>> output;
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) {
				if (isValid(r, c)) output.push_back({ r, c });
			}
		}
		return output;
	}

	// O(1)
	// @param `dirs` the valid directions to move to adjacent nodes
	// @returns coordinates which are adjacent to `(r, c)`
	vector<pair<int, int>> getAdjCoords(int r, int c, vector<pair<int, int>> dirs) {
		vector<pair<int, int>> output;
		for (pair<int, int> dir : dirs) {
			int newR = r + dir.first, newC = c + dir.second;
			if (validCoord(newR, newC)) output.push_back({ newR, newC });
		}
		return output;
	}
	vector<pair<int, int>> getAdjCoords(pair<int, int> coord, vector<pair<int, int>> dirs) {
		return getAdjCoords(coord.first, coord.second, dirs);
	}

	/************************************************
	 *               BFS & VARIATIONS               *
	 ************************************************/

private:
	// O(size of region) <= O(RC)
	// Standard Breath first search algorithm, starting from (r, c)
	// @param `isValid` Specifies which coordinates to explore that satisfy in the directions of `dirs`, taking a maximum of `maxD` steps
	// @param `dirs` the valid directions to move to adjacent nodes
	// @param `depth` Stores the shortest distance to the starting node, and it is also used as the seen array
	// @param `maxD` Indicates the number of BFS steps taken
	void bfs(int r, int c, function<bool(int, int)> isValid, vector<pair<int, int>> dirs, vector<vector<int>> &depth, int maxD = MAXN) {
		if (!validCoord(r, c)) return;
		if (!isValid(r, c)) return;

		queue<pair<int, int>> q;
		q.push({ r, c });
		depth[r][c] = 0;

		for (int d = 1; !q.empty() && d < maxD; ++d) {
			int size = q.size();
			for (int i = 0; i < size; ++i) {
				auto curr = q.front();
				q.pop();

				for (pair<int, int> adj : getAdjCoords(curr.first, curr.second, dirs)) {
					if (depth[adj.first][adj.second] == -1 && isValid(adj.first, adj.second)) {
						q.push(adj);
						depth[adj.first][adj.second] = d;
					}
				}
			}
		}
	}
	void bfs(pair<int, int> coord, function<bool(int, int)> isValid, vector<pair<int, int>> dirs, vector<vector<int>>& depth, int maxD = MAXN) {
		return bfs(coord, isValid, dirs, depth);
	}

public:
	// O(RC) Finds coordinates that satisfy `isValid`, and are connected to (startR, startC)
	// @param `dirs` the valid directions to move to adjacent nodes
	// @returns a vector of coordinates, where each coordinate is a pair of (r, c)
	vector<pair<int, int>> getRegion(int startR, int startC, function<bool(int, int)> isValid, vector<pair<int, int>> dirs) {
		vector<vector<int>> depth = vector<vector<int>>(R, vector<int>(C, -1));
		bfs(startR, startC, isValid, dirs, depth);

		vector<pair<int, int>> output;
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) {
				if (depth[r][c] >= 0) output.push_back({ r, c });
			}
		}
		return output;
	}
	vector<pair<int, int>> getRegion(pair<int, int> coord, function<bool(int, int)> isValid, vector<pair<int, int>> dirs) {
		return getRegion(coord.first, coord.second, isValid, dirs);
	}

	// O(RC) Finds regions whcih satisfy `isValid`
	// @param `dirs` the valid directions to move to adjacent nodes
	// @returns a vector of regions, where each region is a vector of coordinates
	vector<vector<pair<int, int>>> getRegions(function<bool(int, int)> isValid, vector<pair<int, int>> dirs) {
		vector<vector<int>> depth = vector<vector<int>>(R, vector<int>(C, -1));

		vector<vector<pair<int, int>>> output;
		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) {
				if (depth[r][c] == -1 && isValid(r, c)) {
					bfs(r, c, isValid, dirs, depth);
					output.push_back(getRegion(r, c, isValid, dirs));
				}
			}
		}
		return output;
	}

	// O(RC) Determines whether a region of coordinates which satisfy `isValid` are commpletely surrounded coordinates which satisfy `isBoundary`
	bool isSurroundedBy(int startR, int startC, function<bool(int, int)> isValid, function<bool(int, int)> isBoundary, vector<pair<int, int>> dirs) {
		vector<vector<int>> depth = vector<vector<int>>(R, vector<int>(C, -1));
		bfs(startR, startC, isValid, dirs, depth);

		for (int r = 0; r < R; ++r) {
			for (int c = 0; c < C; ++c) {
				if (isValid(r, c)) {
					for (pair<int, int> adj : getAdjCoords(r, c, dirs)) {
						if (!isValid(adj.first, adj.second) && !isBoundary(adj.first, adj.second)) return false;
					}
				}
			}
		}
		return true;
	}
	bool isSurroundedBy(pair<int, int> coord, T val, function<bool(int, int)> isValid, function<bool(int, int)> isBoundary, vector<pair<int, int>> dirs) {
		return isSurroundedBy(coord.first, coord.second, val, isValid, isBoundary, dirs);
	}

	// Finds the shortest path from (r1, c1) to (r2, c2) inclusive, where each coordinate on the path satisfies `isValid`
	// @param `dirs` the valid directions to move to adjacent nodes
	vector<pair<int, int>> shortestPath(int r1, int c1, int r2, int c2, function<bool(int, int)> isValid, vector<pair<int, int>> dirs) {
		vector<vector<int>> depth = vector<vector<int>>(R, vector<int>(C, -1));
		bfs(r2, c2, isValid, dirs, depth);
		assert(depth[r1][c1] != -1 && "No path found");

		vector<pair<int, int>> output;
		for (int r = r1, c = c1; !(r == r2 && c == c2); ) {
			output.push_back({ r, c });
			for (pair<int, int> adj : getAdjCoords(r, c, dirs)) {
				int newDepth = depth[adj.first][adj.second];
				if (newDepth != -1 && newDepth < depth[r][c]) {
					r = adj.first;
					c = adj.second;
				}
			}
		}
		output.push_back({ r2, c2 });

		return output;
	}
	vector<pair<int, int>> shortestPath(pair<int, int> coord1, pair<int, int> coord2, function<bool(int, int)> isValid, vector<pair<int, int>> dirs) {
		return shortestPath(coord1.first, coord1.second, coord2.first, coord2.second, isValid, dirs);
	}
};

/************************************************
 *                    DISPLAY                   *
 ************************************************/
template<typename T> ostream& operator<<(ostream& out, Grid<T> grid) {
	grid.print(out);
	return out;
}
