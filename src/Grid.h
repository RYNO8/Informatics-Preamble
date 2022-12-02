#pragma once
#include "Constants.h"

namespace DS {
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
		std::vector<std::vector<T>> grid;
		// Prefix sum of grid
		std::vector<std::vector<T>> prefix;
		// Whether the graph hasn't been modified after the last precomp
		bool hasPrecomp = false;

		// O(RC)
		// Initialises all variables
		void init(int _R, int _C) {
			R = _R;
			C = _C;
			grid = std::vector<std::vector<T>>(R);
			for (int r = 0; r < R; ++r) grid[r] = std::vector<T>(C);
			prefix = std::vector<std::vector<T>>(R + 1);
			for (int r = 0; r < R + 1; ++r) prefix[r] = std::vector<T>(C + 1);
		}

	public:
		// O(1)
		// Initialises an empty grid
		Grid() {

		}

		// O(RC)
		// Initialises a grid filled with `val`
		Grid(int _R, int _C, T val = T(0)) {
			init(_R, _C);
			for (int r = 0; r < R; ++r) {
				for (int c = 0; c < C; ++c) {
					grid[r][c] = val;
				}
			}
		}

		// O(RC)
		// Initialises grid from cin, assuming 0 indexed
		Grid(int _R, int _C, std::istream& in) {
			init(_R, _C);
			for (int r = 0; r < R; ++r) {
				for (int c = 0; c < C; ++c) {
					in >> grid[r][c];
				}
			}
			buildPrefix();
		}

		// O(RC)
		// Initialises grid from a 2d array, assuming 0 indexed
		template<size_t rows, size_t cols> Grid(T(&_grid)[rows][cols]) {
			init(sizeof(_grid) / sizeof(_grid[0]), sizeof(_grid[0]) / sizeof(T));

			for (int r = 0; r < R; ++r) {
				for (int c = 0; c < C; ++c) {
					grid[r][c] = _grid[r][c];
				}
			}
			buildPrefix();
		}

		// O(RC)
		//Initialises grid from a 2d vector, assuming 0 indexed
		Grid(std::vector<std::vector<T>> _grid) {
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

		// O(RC)
		// Displays the grid
		// @param `out` The string representation of the graph is piped to this output stream
		friend std::ostream& operator<<(std::ostream& out, const Grid<T> grid) {
			for (int r = 0; r < grid.R; ++r) {
				for (int c = 0; c < grid.C; ++c) {
					out << grid.grid[r][c] << ' ';
				}
				out << '\n';
			}
			return out;
		}

		// O(RC)
		// Displays the grid with pretty borders
		// @param `out` The string representation of the graph is piped to this output stream
		// @param `spacing` Indicates how large each cell is (use spacing = -1 for automatic spacing)
		void pprint(std::ostream& out = std::cout, int spacing = -1) {
			if (spacing == -1) {
				// find max repr length of grid values
				for (int r = 0; r < R; ++r) {
					for (int c = 0; c < C; ++c) spacing = std::max(spacing, (reprLen(grid[r][c])+3) / 4);
				}
			}

			std::string sep = "+", space = "|";
			for (int r = 0; r < R; ++r) {
				sep += std::string(4 * spacing + 1, '-') + '+';
				space += std::string(4 * spacing + 1, ' ') + '|';
			}
			sep += '\n';
			space += '\n';
			
			for (int r = 0; r < R; ++r) {
				out << sep;

				for (int rep = 0; rep < spacing; ++rep) out << space;

				out << '|';
				for (int c = 0; c < C; ++c) {
					int paddingL = std::max(0, 2 * spacing - reprLen(grid[r][c]) / 2);
					int paddingR = std::max(0, 2 * spacing - (reprLen(grid[r][c]) - 1) / 2);

					out << std::string(paddingL, ' ') << grid[r][c] << std::string(paddingR, ' ') << '|';
				}
				out << '\n';

				for (int rep = 0; rep < spacing; ++rep) out << space;
			}
			out << sep;
		}

		/************************************************
		 *             COMPARISON FUNCTIONS             *
		 ************************************************/

		// O(1)
		// Filter generator, to determine whether cell values are equal to `val`
		// @returns A std::function which acts to "filter" cells 
		std::function<bool(int, int)> compEQ(T val) {
			return [&, val](int r, int c) {
				if (validCoord(r, c)) return getVal(r, c) == val;
				return false;
			};
		}

		// O(1)
		// Filter generator, to determine whether cell values are not equal to `val`
		// @returns A std::function which acts to "filter" cells 
		std::function<bool(int, int)> compNE(T val) {
			return [&, val](int r, int c) {
				if (validCoord(r, c)) return getVal(r, c) != val;
				return false;
			};
		}

		// O(1)
		// Filter generator, to determine whether cell values are greater than `val`
		// @returns A std::function which acts to "filter" cells 
		std::function<bool(int, int)> compGT(T val) {
			return [&, val](int r, int c) {
				if (validCoord(r, c)) return getVal(r, c) > val;
				return false;
			};
		}

		// O(1)
		// Filter generator, to determine whether cell values are greater or equal to `val`
		// @returns A std::function which acts to "filter" cells 
		std::function<bool(int, int)> compGE(T val) {
			return [&, val](int r, int c) {
				if (validCoord(r, c)) return getVal(r, c) >= val;
				return false;
			};
		}
		
		// O(1)
		// Filter generator, to determine whether cell values are less than `val`
		// @returns A std::function which acts to "filter" cells 
		std::function<bool(int, int)> compLT(T val) {
			return [&, val](int r, int c) {
				if (validCoord(r, c)) return getVal(r, c) < val;
				return false;
			};
		}
		
		// O(1)
		// Filter generator, to determine whether cell values are less or equal to `val`
		// @returns A std::function which acts to "filter" cells 
		std::function<bool(int, int)> compLE(T val) {
			return [&, val](int r, int c) {
				if (validCoord(r, c)) return getVal(r, c) <= val;
				return false;
			};
		}

		// O(1)
		// Filter generator, which returns a truthy std::function - to consider any cell
		std::function<bool(int, int)> compAny() {
			return [](int r, int c) {
				return true;
			};
		}

		/************************************************
		 *             COORDINATES & VALUES             *
		 ************************************************/

		// O(1)
		// @returns R - the number of rows
		inline const int getR() const {
			return R;
		}

		// O(1)
		// @returns C - the number of columns
		inline const int getC() const {
			return C;
		}

		// O(1)
		// @returns The size of the grid
		inline const int size() const {
			return R * C;
		}

		// O(1)
		// @returns Whether a coordinate is within the grid
		inline bool validCoord(int r, int c) const {
			return 0 <= r && r < R && 0 <= c && c < C;
		}
		inline bool validCoord(std::pair<int, int> coord) const {
			return validCoord(coord.first, coord.second);
		}

		// O(RC)
		// @returns The first (lexigraphically smallest) coordinate which satisfies `isValid`, otherwise returns {-1, -1}
		std::pair<int, int> findVal(std::function<bool(int, int)> isValid) {
			for (int r = 0; r < R; ++r) {
				for (int c = 0; c < C; ++c) {
					if (isValid(r, c)) return { r, c };
				}
			}
			return { -1, -1 };
		}

		// O(RC)
		// @returns Whether there exists a coordinate that satisfies `isValid`
		inline bool containsVal(std::function<bool(int, int)> isValid) const {
			return findVal(isValid).first != -1;
		}

		// O(RC)
		// @returns A vector of coordinates which satisfies `isValid`
		std::vector<std::pair<int, int>> findAllVal(std::function<bool(int, int)> isValid) const {
			std::vector<std::pair<int, int>> output;
			for (int r = 0; r < R; ++r) {
				for (int c = 0; c < C; ++c) {
					if (isValid(r, c)) output.push_back({ r, c });
				}
			}
			return output;
		}

		// O(RC)
		// @returns the number of cells that satisfy `isValid`
		int count(std::function<bool(int, int)> isValid) const {
			int total = 0;
			for (int r = 0; r < R; ++r) {
				for (int c = 0; c < C; ++c) total += (bool)isValid(r, c);
			}
			return total;
		}

		// O(1)
		// @returns If the provided coordinate is valid, set its value to `val`, otherwise silently do nothing
		inline void setVal(int r, int c, T val) {
			if (validCoord(r, c)) {
				hasPrecomp = false;
				grid[r][c] = val;
			}
		}
		void setVal(std::pair<int, int> coord, T val) {
			setVal(coord.first, coord.second, val);
		}

		// O(area of rectangular region) <= O(RC) Sets the value of all cells within a rectangular region (inclusive)
		// @note The rectangular region will be capped so it stays within the grid
		void setVals(int r1, int c1, int r2, int c2, T val) {
			hasPrecomp = false;
			if (r1 > r2) std::swap(r1, r2);
			if (c1 > c2) std::swap(c1, c2);
			r1 = std::max(0, r1);
			c1 = std::max(0, c1);
			r2 = std::min(R - 1, r2);
			c2 = std::max(C - 1, c2);

			for (int r = r1; r <= r2; ++r) {
				for (int c = c1; c <= c2; ++c) grid[r][c] = val;
			}
		}
		void setVals(std::pair<int, int> coord1, std::pair<int, int> coord2, T val) {
			setVals(coord1.first, coord1.second, coord2.first, coord2.second, val);
		}
		
		// O(1) If the provided coordinate is valid, increment its value by `val`
		void incrVal(int r, int c, T val) {
			if (validCoord(r, c)) {
				hasPrecomp = false;
				grid[r][c] += val;
			}
		}
		void incrVal(std::pair<int, int> coord, T val) {
			incrVal(coord.first, coord.second, val);
		}

		// O(area of rectangular region) <= O(RC) Increments the value of all cells within a rectangular region (inclusive)
		// @note The rectangular region will be capped so it stays within the grid
		void incrVals(int r1, int c1, int r2, int c2, T val) {
			hasPrecomp = false;
			if (r1 > r2) std::swap(r1, r2);
			if (c1 > c2) std::swap(c1, c2);
			r1 = std::max(0, r1);
			c1 = std::max(0, c1);
			r2 = std::min(R - 1, r2);
			c2 = std::max(C - 1, c2);

			for (int r = r1; r <= r2; ++r) {
				for (int c = c1; c <= c2; ++c) grid[r][c] += val;
			}
		}
		void incrVals(std::pair<int, int> coord1, std::pair<int, int> coord2, T val) {
			incrVals(coord1.first, coord1.second, coord2.first, coord2.second, val);
		}

		// O(1)
		// @returns The value at the coordinate (r, c) if it is valid, otherwise return `defaultVal`
		T getVal(int r, int c, T defaultVal = T()) const {
			if (validCoord(r, c)) return grid[r][c];
			return defaultVal;
		}
		T getVal(std::pair<int, int> coord, T defaultVal = T()) const {
			return getVal(coord.first, coord.second, defaultVal);
		}

		// O(1) if precomp, otherwise O(RC)
		// @returns The sum of all cells within a rectangular region (inclusive)
		// @note This is only garunteed to be correct after initPrefix() is called, and the grid is not modified after that
		// @note The rectangular region will be capped so it stays within the grid
		T sumVals(int r1, int c1, int r2, int c2) const {
			if (r1 > r2) std::swap(r1, r2);
			if (c1 > c2) std::swap(c1, c2);
			r1 = std::max(0, r1);
			c1 = std::max(0, c1);
			r2 = std::min(R - 1, r2);
			c2 = std::max(C - 1, c2);

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
		T sumVals(std::pair<int, int> coord1, std::pair<int, int> coord2) const {
			return sumVals(coord1.first, coord1.second, coord2.first, coord2.second);
		}

		/************************************************
		 *                 OPERATIONS                   *
		 ************************************************/

		// O(RC)
		// @returns Whether 2 grids are identical
		bool operator==(Grid o) const {
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

		// O(RC)
		// Flips the order of rows - first row becomes last row and vice versa
		// @returns Grid object of the modified grid
		Grid<T> flipRows() const {
			Grid<T> output(R, C, 0);
			for (int r = 0; r < R; ++r) {
				for (int c = 0; c < C; ++c) output.setVal(R - r - 1, c, grid[r][c]);
			}
			return output;
		}

		// O(RC)
		// Flips the order of columns - first column becomes last row and vice versa
		// @returns grid object of the modified grid
		Grid<T> flipCols() const {
			Grid<T> output(R, C, 0);
			for (int r = 0; r < R; ++r) {
				for (int c = 0; c < C; ++c) output.setVal(r, C - c - 1, grid[r][c]);
			}
			return output;
		}

		// O(RC) Rotates the grid clockwise/anticlockwise by 180 degrees, about its center
		// @returns grid object of the modified grid
		Grid<T> rot180() const {
			return flipRows().flipCols();
		}

		// O(RC)
		// Flips the grid over its primary diagonal (top left to bottom right)
		// @returns Grid object of the modified grid
		Grid<T> flipPrimaryDiag() const {
			Grid<T> output(C, R, 0);
			for (int r = 0; r < R; ++r) {
				for (int c = 0; c < C; ++c) output.setVal(c, r, grid[r][c]);
			}
			return output;
		}

		// O(RC)
		// Flips the grid over its secondary diagonal (top right to bottom left)
		// @returns Grid object of the modified grid
		Grid<T> flipSecondaryDiag() const {
			return rot180().flipPrimaryDiag();
		}

		// O(RC)
		// Rotates the grid clockwise by 90 degrees, about its center
		// @returns grid object of the modified grid
		Grid<T> rotClockwise() {
			return flipPrimaryDiag().flipCols();
		}

		// O(RC)
		// Rotates the grid anticlockwise by 90 degrees, about its center
		// @returns grid object of the modified grid
		Grid<T> rotAnticlockwise() const {
			return flipPrimaryDiag().flipRows();
		}

		// O(RC)
		// Replaces every instance of `prevVal` with `newVal`
		// @returns Grid object of the modified grid`
		Grid<T> replace(T prevVal, T newVal) const {
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
		// @returns All coordinates in a particular row, sorted by increasing order
		std::vector<std::pair<int, int>> getRowCoords(int r) const {
			std::vector<std::pair<int, int>> output;
			for (int c = 0; c < C; ++c) output.push_back({ r, c });
			return output;
		}

		// O(R)
		// @returns All coordinates in a particular column, sorted by increasing order
		std::vector<std::pair<int, int>> getColCoords(int c) const {
			std::vector<std::pair<int, int>> output;
			for (int r = 0; r < R; ++r) output.push_back({ r, c });
			return output;
		}

		// O(RC)
		// @returns All coordinates in a vector, sorted by lexicographic order
		std::vector<std::pair<int, int>> getAllCoords() const {
			//return getSomeCoords(compAny);
			std::vector<std::pair<int, int>> output;
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
		std::vector<std::pair<int, int>> getSomeCoords(std::function<bool(int, int)> isValid) const {
			std::vector<std::pair<int, int>> output;
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
		std::vector<std::pair<int, int>> getAdjCoords(int r, int c, std::vector<std::pair<int, int>> &dirs) const {
			std::vector<std::pair<int, int>> output;
			for (std::pair<int, int> dir : dirs) {
				int newR = r + dir.first, newC = c + dir.second;
				if (validCoord(r, c) && validCoord(newR, newC)) output.push_back({ newR, newC });
			}
			return output;
		}
		std::vector<std::pair<int, int>> getAdjCoords(std::pair<int, int> coord, std::vector<std::pair<int, int>> &dirs) const {
			return getAdjCoords(coord.first, coord.second, dirs);
		}

		/************************************************
		 *               BFS & VARIATIONS               *
		 ************************************************/

	private:
		// O(RC)
		// Standard Breadth First Search, starting from (r, c)
		// @param `isValid` Specifies which coordinates to explore that satisfy in the directions of `dirs`, taking a maximum of `maxD` steps
		// @param `dirs` the valid directions to move to adjacent nodes
		// @param `depth` Stores the shortest distance to the starting node, and it is also used as the seen array
		// @param `maxD` Indicates the number of BFS steps taken
		void bfs(int r, int c, std::function<bool(int, int)> isValid, std::vector<std::pair<int, int>> dirs, std::vector<std::vector<int>> &depth, int maxD = INT_MAX) const {
			if (!validCoord(r, c)) return;
			if (!isValid(r, c)) return;

			std::queue<std::pair<int, int>> q;
			q.push({ r, c });
			depth[r][c] = 0;

			for (int d = 1; !q.empty() && d < maxD; ++d) {
				int size = q.size();
				for (int i = 0; i < size; ++i) {
					auto curr = q.front();
					q.pop();

					for (std::pair<int, int> adj : getAdjCoords(curr.first, curr.second, dirs)) {
						if (depth[adj.first][adj.second] == -1 && isValid(adj.first, adj.second)) {
							q.push(adj);
							depth[adj.first][adj.second] = d;
						}
					}
				}
			}
		}
		void bfs(std::pair<int, int> coord, std::function<bool(int, int)> isValid, std::vector<std::pair<int, int>> &dirs, std::vector<std::vector<int>>& depth, int maxD = INT_MAX) const {
			return bfs(coord, isValid, dirs, depth);
		}

	public:
		// O(RC)
		// Finds coordinates that satisfy `isValid`, and are connected to (startR, startC)
		// @param `dirs` the valid directions to move to adjacent nodes
		// @returns a vector of coordinates, where each coordinate is a pair of (r, c)
		std::vector<std::pair<int, int>> getRegion(int startR, int startC, std::function<bool(int, int)> isValid, std::vector<std::pair<int, int>>& dirs) const {
			std::vector<std::vector<int>> depth = std::vector<std::vector<int>>(R, std::vector<int>(C, -1));
			bfs(startR, startC, isValid, dirs, depth);

			std::vector<std::pair<int, int>> output;
			for (int r = 0; r < R; ++r) {
				for (int c = 0; c < C; ++c) {
					if (depth[r][c] >= 0) output.push_back({ r, c });
				}
			}
			return output;
		}
		std::vector<std::pair<int, int>> getRegion(std::pair<int, int> coord, std::function<bool(int, int)> isValid, std::vector<std::pair<int, int>> &dirs) const {
			return getRegion(coord.first, coord.second, isValid, dirs);
		}

		// O(RC)
		// Finds regions whcih satisfy `isValid`
		// @param `dirs` the valid directions to move to adjacent nodes
		// @returns a vector of regions, where each region is a vector of coordinates
		std::vector<std::vector<std::pair<int, int>>> getRegions(std::function<bool(int, int)> isValid, std::vector<std::pair<int, int>> &dirs) const {
			std::vector<std::vector<int>> depth = std::vector<std::vector<int>>(R, std::vector<int>(C, -1));

			std::vector<std::vector<std::pair<int, int>>> output;
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

		// O(RC)
		// Determines whether a region of coordinates which satisfy `isValid` are commpletely surrounded coordinates which satisfy `isBoundary`
		bool isSurroundedBy(int startR, int startC, std::function<bool(int, int)> isValid, std::function<bool(int, int)> isBoundary, std::vector<std::pair<int, int>> &dirs) const {
			std::vector<std::vector<int>> depth = std::vector<std::vector<int>>(R, std::vector<int>(C, -1));
			bfs(startR, startC, isValid, dirs, depth);

			for (int r = 0; r < R; ++r) {
				for (int c = 0; c < C; ++c) {
					if (isValid(r, c)) {
						for (std::pair<int, int> adj : getAdjCoords(r, c, dirs)) {
							if (!isValid(adj.first, adj.second) && !isBoundary(adj.first, adj.second)) return false;
						}
					}
				}
			}
			return true;
		}
		bool isSurroundedBy(std::pair<int, int> coord, T val, std::function<bool(int, int)> isValid, std::function<bool(int, int)> isBoundary, std::vector<std::pair<int, int>> &dirs) const {
			return isSurroundedBy(coord.first, coord.second, val, isValid, isBoundary, dirs);
		}

		// O(RC)
		// Finds the shortest path from (r1, c1) to (r2, c2) inclusive, where each coordinate on the path satisfies `isValid`
		// @param `dirs` the valid directions to move to adjacent nodes
		std::vector<std::pair<int, int>> shortestPath(int r1, int c1, int r2, int c2, std::function<bool(int, int)> isValid, std::vector<std::pair<int, int>> &dirs) const {
			std::vector<std::vector<int>> depth = std::vector<std::vector<int>>(R, std::vector<int>(C, -1));
			bfs(r2, c2, isValid, dirs, depth);
			assert(depth[r1][c1] != -1 && "No path found");

			std::vector<std::pair<int, int>> output;
			for (int r = r1, c = c1; !(r == r2 && c == c2); ) {
				output.push_back({ r, c });
				for (std::pair<int, int> adj : getAdjCoords(r, c, dirs)) {
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
		std::vector<std::pair<int, int>> shortestPath(std::pair<int, int> coord1, std::pair<int, int> coord2, std::function<bool(int, int)> isValid, std::vector<std::pair<int, int>> &dirs) const {
			return shortestPath(coord1.first, coord1.second, coord2.first, coord2.second, isValid, dirs);
		}
	};
};