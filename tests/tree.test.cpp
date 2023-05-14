#include <iostream>
#include <fstream>
#include <sstream>
#include "../src/Tree.h"
#include "../src/Graph.h"
#include "../src/Util.h"
using namespace std;
using namespace DS;

/*
 * 1<----->2
 * ^      /|
 * |     / |
 * |    /  |
 * v   |   v
 * 3<-`    4
 * 
 */
void testBinaryTree() {
    using Graph1 = Graph<10, UnitEdgeWeight, int, true>;
    using Tree1 = Tree<10, UnitEdgeWeight, int>;
    stringstream s1;
    s1 << R"(
7
1 2
1 3
2 4
2 5
3 6
3 7
)";
    size_t N1;
    s1 >> N1;
    vector<Graph1::Edge> e(N1 - 1);
    s1 >> e;
    Graph1 G1(N1, e);
    Tree1 T1(G1, 1);
    assert(repr(T1) == 
"1\n"
"├─2\n"
"│ ├─4\n"
"│ └─5\n"
"└─3\n"
"  ├─6\n"
"  └─7\n"
    );
}

void testWeightedTree() {
    using Graph2 = WeightedDiGraph<12>;
    using Tree2 = Tree<12, int, int>;
    stringstream s2;
    s2 << R"""(
11
1 2 3
1 3 4
3 4 5
4 5 4
5 6 6
1 7 3
7 8 2
7 9 5
9 10 6
9 11 7
)""";
    size_t N2;
    s2 >> N2;
    vector<Graph2::Edge> e(N2 - 1);
    s2 >> e;
    Graph2 G2(N2, e);
    Tree2 T2(G2, 1);
    assert(repr(T2) == 
"1\n"
"├───2\n"
"├────3\n"
"│    └─────4\n"
"│          └────5\n"
"│               └──────6\n"
"└───7\n"
"    ├──8\n"
"    └─────9\n"
"          ├──────10\n"
"          └───────11\n"
    );
}

int main() {
    testBinaryTree();
    testWeightedTree();
}