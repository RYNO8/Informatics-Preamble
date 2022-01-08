// TODO: test all
template<typename T> class Point {
	/************************************************
	 *                INITIALISATION                *
	 ************************************************/

private:
    T x, y;

public:
	// O(1) Initialises a Point
    Point(T x_, T y_) {
		assert(max(abs(x_), abs(y_)) <= sqrt(numeric_limits<T>::max()) && "Operations of points will cause overflow");
		x = x_;
		y = y_;
    }

	/************************************************
	 *                    DISPLAY                   *
	 ************************************************/

	// O(1)
	// @param `out` The string representation of the graph is piped to this output stream
	// @param `newLine` Indicates whether to end with a trailing `\\n`
	void print(ostream& out = cout, bool newLine = false) {
		out << '(' << x << ", " << y << ')';
		if (newLine) out << '\n';
	}

	/************************************************
	 *                   PROPERTIES                 *
	 ************************************************/

	// O(1) Gets x coordinate
	const T getX() {
		return x;
	}

	// O(1) Gets y coordinate
	const T getY() {
		return y;
	}

	/************************************************
	 *               BASCIC OPERATIONS              *
	 ************************************************/
	bool operator==(Point o) const {
		// TODO: if `T` is a floating point type, account for floating point imprecisions
		return x == o.x && y == o.y;
	}

	bool operator!=(Point o) const {
		// TODO: if `T` is a floating point type, account for floating point imprecisions
		return x != o.x || y != o.y;
	}

	// O(1) Standard modular addition
	Point operator+(Point o) const {
		return Point(x + o.x, y + o.y);
	}

	// O(1) Standard modular subtraction
	Point operator-(Point o) const {
		return Point(x - o.x, y - o.y);
	}

	// O(1) Addition assignment
	Point operator+=(Point o) {
		*this = *this + o;
		return *this;
	}

	// O(1) Subtraction assignment
	Point operator-=(Point o) {
		*this = *this - o;
		return *this;
	}
	
	/************************************************
	 *               VECTOR OPERATIONS              *
	 ************************************************/

	// O(1) Scale by `k`
	Point operator*(T k) {
		return Point(x * k, y * k);
	}

	// O(1) Multiplication assignment
	Point operator*=(T k) {
		*this = *this * k;
		return *this;
	}

	// O(1) Dot product
	T dot(Point o) {
		return x * o.x + y * o.y;
	}

	// O(1) Skew product
	T skew(Point o) {
		return x * o.x - y * o.y;
	}

	// O(1) Pythagorean distance from origin
	ld length() {
		return sqrt(x * x + y * y);
	}
	// O(1)
	// @returns anticlockwise order from `*this` to `q` to `r`
	T ccw(Point<T> &q, Point<T> &r) {
		return (q.y - y) * (r.x - q.x) - (q.x - x) * (r.y - q.y);
	}

	// O(1)
	// @returns orientation from `*this` to `q` to `r`
	// 1: clockwise, 0: collinear, -1: anticlockwise
	T orientation(Point<T>& q, Point<T>& r) {
		return sgn(ccw(q, r));
	}

	// O(1)
	// @returns whether `*this`, `q` and `r` are collinear
	bool isCollinear(Point<T>& q, Point<T>& r) {
		return orientation(q, r) == 0;
	}

	// O(1)
	// whether `*this` lies in the bounding box of `pq`
	bool onSegment(Point p, Point q) {
		return min(p.x, q.x) <= x && x <= max(p.x, q.x) && min(p.y, q.y) <= y && y <= max(p.y, q.y);
	}
};

/************************************************
 *                    DISPLAY                   *
 ************************************************/
template<typename T> ostream& operator<<(ostream& out, Point<T> val) {
	val.print(out);
	return out;
}

template<typename T> class Line {
	/************************************************
	 *                INITIALISATION                *
	 ************************************************/
private:
	Point<T> a, b;

public:
	// O(1) Initialises a Line
	Line(Point<T> a_, Point<T> b_) : a(a_), b(b_) {
		assert(a != b && "This is a degenerate line");
	}

	/************************************************
	 *                    DISPLAY                   *
	 ************************************************/
	// O(1)
	// @param `out` The string representation of the graph is piped to this output stream
	// @param `newLine` Indicates whether to end with a trailing `\\n`
	void print(ostream& out = cout, bool newLine = false) {
		out << a << " -- " << b;
		if (newLine) out << '\n';
	}

	/************************************************
	 *                   PROPERTIES                 *
	 ************************************************/
	// O(1) Gets both endpoints of the line, in no particular order
	const pair<Point<T>, Point<T>> getPoints() {
		return { a, b };
	}
	
	/************************************************
	 *               BASCIC OPERATIONS              *
	 ************************************************/
	// O(1)
	// @returns Whether if the line segment `q` was extended infinitely, it would intesect `this`
	bool intersects(Line<T> &l) {
		T denom = (b - a).skew(l.b - l.a);

		T numerator1 = (b - a).skew(a - l.a);
		T numerator2 = (l.a - a).skew(l.b - l.a);
		if (denom < 0) {
			numerator1 = -numerator1;
			numerator2 = -numerator2;
			denom = -denom;
		}

		if (denom == 0) {
			// TODO: collinear, case bash to check for overlap
			return false;
		}
		return 0 < min(numerator1, numerator2) && max(numerator1, numerator2) < denom;
	}

	// TODO: ccw algorithm to determine line line intersection
	bool intersects2(Line<T>& l) {

	}

	// O(1)
	// @returns Whether there is point `p` lies on this line segment
	bool intersects(Point<T>& p) {
		return p.collinear(a, b) && p.orientation(a, b) == 1;
	}

	// O(1)
	bool touches(Point<T>& p) {
		return p == a || p == b;
	}
};

/************************************************
 *                    DISPLAY                   *
 ************************************************/
template<typename T> ostream& operator<<(ostream& out, Line<T> val) {
	val.print(out);
	return out;
}

template<typename T> class Polygon {
	/************************************************
	 *                INITIALISATION                *
	 ************************************************/
private:
	vector<Point<T>> points;

public:
	// O(n^2) Initialises a Polygon
	Polygon(vector<Point<T>> points_, bool makeConvex = true) : points(points_) {
		assert(points.size() >= 3 && "This is a degenerate polygon");
		// also assuming that all points are unique
		if (makeConvex) makeConvexHull();
	}

	/************************************************
	 *                    DISPLAY                   *
	 ************************************************/
	 // O(1)
	 // @param `out` The string representation of the graph is piped to this output stream
	 // @param `newLine` Indicates whether to end with a trailing `\\n`
	void print(ostream& out = cout, bool newLine = false) {
		out << points;
		if (newLine) out << '\n';
	}

	/************************************************
	 *                   OPERATIONS                 *
	 ************************************************/
	 // O(n^2) Gift wrapping algorithm to find convex hull
	 // @note updates in place
	void makeConvexHull() {
		vector<Point<T>> hull;
		int p = min_element(points.begin(), points.end()) - points.begin();
		do {
			hull.push_back(points[p]);

			int best = (p + 1) % points.size();
			for (int i = 0; i < points.size(); i++) {
				if (orientation(points[p], points[i], points[best]) == 2) best = i;

				// handling collinear points: if point q lies in the middle, then also update q
				if (p != i && points[p].collinear(points[i], points[best]) && points[best].onSegment(points[p], points[i])) {
					best = i;
				}
			}
			p = best;
		} while (points[p] != hull[0]);

		points = hull;
	}

	// O(n)
	// @returns whether point is on the boundary of the polygon
	// @note Assuming polygon is a convex hull
	bool isOnBoundary(Point<T> &point) {
		for (int i0 = 0; i0 < points.size(); ++i0) {
			int i1 = (i0 + 1) % points.size();
			if (Line<T>(points[i0], points[i1]).intersects(point)) return true;
		}
		return false;
	}

	// O(n)
	// @returns whether point is stictly within the boundary of the polygon
	// @note Assuming polygon is a convex hull
	bool isContained(Point<T> &point) {
		for (int i0 = 0; i0 < points.size(); ++i0) {
			int i1 = (i0 + 1) % points.size();
			T curr = orientation(points[0], points[1], point);
			if (curr == 0 || curr != orientation(points[i0], points[i1], point)) return false;
		}
		return true;
	}

	// O(n)
	// @returns whether point is stictly within the boundary of the polygon
	bool isContained2(Point<T> point) {
		// TODO: use ray tracing
	}
};

/************************************************
 *                    DISPLAY                   *
 ************************************************/
template<typename T> ostream& operator<<(ostream& out, Polygon<T> val) {
	val.print(out);
	return out;
}
