#include <iostream>
#include <fstream>
#include <sstream>
#include "../src/Tree.h"
#include "../src/Graph.h"
#include "../src/Util.h"
#include "helpers.h"
using namespace std;
using namespace DS;

/*
 *     1     
 *  ┌──┴──┐  
 *  2     3  
 * ┌┴─┐  ┌┴─┐
 * 4  5  6  7
 *          │
 *          8
 *
 *
 * 1
 * ├─2
 * │ ├─4
 * │ └─5
 * │   ├─9
 * │   └─10
 * └─3
 *   ├─6
 *   └─7
 *     └─8
 *
 **/
void testBinaryTree() {
    using Graph1 = UnweightedDiGraph<12>;
    using Tree1 = Tree<12, unsigned int>;
    stringstream s1;
    s1 << R"(
8
1 2
1 3
2 4
2 5
3 6
3 7
7 8
)";
    size_t N1;
    s1 >> N1;
    vector<Graph1::Edge> e(N1 - 1);
    s1 >> e;
    Graph1 G1(N1, e);
    Tree1 T1(G1, 1);

    assert(T1.N() == 8);
    assert(T1.V() == 8);
    assert(T1.M() == 7);
    assert(T1.E() == 7);
    assert(T1.getRoot() == 1);
    assert(unordered_eq(T1.getNodes(), vector<Tree1::Node>({ 1, 2, 3, 4, 5, 6, 7, 8 })));
    assert(unordered_eq(T1.getEdges(), vector<Tree1::Edge>({
        Tree1::Edge(1, 2, 1),
        Tree1::Edge(1, 3, 1),
        Tree1::Edge(2, 4, 1),
        Tree1::Edge(2, 5, 1),
        Tree1::Edge(3, 6, 1),
        Tree1::Edge(3, 7, 1),
        Tree1::Edge(7, 8, 1),
    })));

    assert(unordered_eq(T1.getEdges(2), vector<Tree1::Edge>({
        Tree1::Edge(2, 4, 1),
        Tree1::Edge(2, 5, 1),
    })));
    assert(T1.getEdges(8) == vector<Tree1::Edge>());
    assert(T1.getChildren(2) == vector<Tree1::Node>({ 4, 5 }));
    assert(T1.getChildren(8) == vector<Tree1::Node>({}));
    assert(T1.degree(1) == 2);
    assert(T1.degree(2) == 3);
    assert(T1.degree(8) == 1);
    assert(T1.numChildren(1) == 2);
    assert(T1.numChildren(2) == 2);
    assert(T1.numChildren(8) == 0);


    assert(T1.isNode(0) == true);
    assert(T1.isNode(1) == true);
    assert(T1.isNode(11) == true);
    assert(T1.isNode(12) == false);
    assert(T1.containsNode(0) == false);
    assert(T1.containsNode(1) == true);
    assert(T1.containsNode(11) == false);
    assert(T1.containsNode(12) == false);
    assert(T1.containsEdge(Tree1::Edge(1, 2, 1)) == true);
    assert(T1.containsEdge(Tree1::Edge(1, 2, 3175368)) == false);
    assert(T1.containsEdge(Tree1::Edge(1, 5, 1)) == false);
    assert(T1.containsEdgeUnweighted(Tree1::Edge(1, 2, 1)) == true);
    assert(T1.containsEdgeUnweighted(Tree1::Edge(1, 2, 3175368)) == true);
    assert(T1.containsEdgeUnweighted(Tree1::Edge(1, 5, 1)) == false);
    assert(T1.getEdge(2, 4, 69) == 1);
    assert(T1.getEdge(2, 3, 69) == 69);


    assert(repr(T1) == 
"1\n"
"├─2\n"
"│ ├─4\n"
"│ └─5\n"
"└─3\n"
"  ├─6\n"
"  └─7\n"
"    └─8\n"
    );
    assert(TREE_PRINT_SPACING == 2);
    assert(repr(&T1) ==
"    1     \n"
" ┌──┴──┐  \n"
" 2     3  \n"
"┌┴─┐  ┌┴─┐\n"
"4  5  6  7\n"
"         │\n"
"         8\n"
    );

    T1.insert(Tree1::Edge(5, 9));
    T1.insert(Tree1::Edge(5, 10));
    assert(repr(T1) == 
"1\n"
"├─2\n"
"│ ├─4\n"
"│ └─5\n"
"│   ├─9\n"
"│   └─10\n"
"└─3\n"
"  ├─6\n"
"  └─7\n"
"    └─8\n"
    );
    assert(repr(&T1) ==
"      1       \n"
"  ┌───┴────┐  \n"
"  2        3  \n"
"┌─┴──┐    ┌┴─┐\n"
"4    5    6  7\n"
"   ┌─┴─┐     │\n"
"   9  10     8\n"
    );

    T1.erase(9);
    T1.erase(10);

    assert(T1.isBinary());
    T1.insert(Tree1::Edge(5, 9));
    assert(T1.isBinary());
    T1.insert(Tree1::Edge(5, 10));
    assert(T1.isBinary());
    T1.insert(Tree1::Edge(5, 11));
    assert(!T1.isBinary());
    T1.erase(9);
    T1.erase(10);
    T1.erase(11);

    assert(!T1.isPath());
    // assert(Tree1(10, 1).isPath());
}

/*
 * 1
 * ├───2
 * ├────3
 * │    └─────4
 * │          └────5
 * │               └──────6
 * └───7
 *     ├──8
 *     └─────9
 *         ├──────10
 *         └───────11
 *
 *
 *     1          
 * ┌──┬┴────┐     
 * 2  3     7     
 *    │  ┌──┴──┐  
 *    4  8     9  
 *    │      ┌─┴─┐
 *    5     10  11
 *    │           
 *    6           
 *
 **/
void testWeightedTree() {
    using Graph2 = WeightedDiGraph<12>;
    using Tree2 = Tree<12, unsigned int>;
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
    assert(repr(&T2) == 
"    1          \n"
"┌──┬┴────┐     \n"
"2  3     7     \n"
"   │  ┌──┴──┐  \n"
"   4  8     9  \n"
"   │      ┌─┴─┐\n"
"   5     10  11\n"
"   │           \n"
"   6           \n"
    );
}

int main() {
    testBinaryTree();
    testWeightedTree();
}