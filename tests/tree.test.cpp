#include "../src/Tree.h"
#include "../src/Graph.h"
#include "../src/Util.h"
#include "helpers.h"
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;
using namespace DS;

using MyTree = Tree<12, unsigned int>;

/* Mr B
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
 *
 *       1
 *   ┌───┴────┐
 *   2        3
 * ┌─┴──┐    ┌┴─┐
 * 4    5    6  7
 *    ┌─┴─┐     │
 *    9  10     8
 *
 **/
MyTree createBinary() {
    using MyGraph = UnweightedDiGraph<12>;
    stringstream input;
    input << R"(
    10
    1 2
    1 3
    2 4
    2 5
    3 6
    3 7
    7 8
    5 9
    5 10
    )";
    size_t N;
    input >> N;
    vector<MyGraph::Edge> e(N - 1);
    input >> e;
    MyGraph G(N, e);
    MyTree T(G, 1);
    return T;
}

MyTree B = createBinary();

/* Mr W
 *
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
MyTree createWeighted() {
    using MyGraph = WeightedDiGraph<12>;
    stringstream input;
    input << R"""(
    11
    1 2 0
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
    size_t N;
    input >> N;
    vector<MyGraph::Edge> e(N - 1);
    input >> e;
    MyGraph G(N, e);
    MyTree T(G, 1);
    return T;
}

MyTree W = createWeighted();

// Singleton (unit graph)
MyTree U = MyTree(7);

// Path Graph
MyTree P = MyTree::PathGraph(5, 2);

// Star Graph
MyTree S = MyTree::StarGraph(5, 2);

void testNodeFunc() {
    assert(B.size() == 10);
    assert(B.N() == 10);
    assert(B.V() == 10);
    assert(B.getRoot() == 1);
    assert(unordered_eq(B.getNodes(), vector<MyTree::Node>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10})));

    assert(B.isNode(0) == true);
    assert(B.isNode(1) == true);
    assert(B.isNode(11) == true);
    assert(B.isNode(12) == false);
    assert(B.containsNode(0) == false);
    assert(B.containsNode(1) == true);
    assert(B.containsNode(11) == false);
    assert(B.containsNode(12) == false);

    assert(B.isLeaf(4));
    assert(!B.isLeaf(0));
    assert(!B.isLeaf(2));
    assert(
        B.getLeaves() == vector<MyTree::Node>({
                             4,
                             6,
                             8,
                             9,
                             10,
                         })
    );
}

void testEdgeFunc() {
    assert(B.M() == 9);
    assert(B.E() == 9);
    assert(unordered_eq(
        B.getEdges(), vector<MyTree::Edge>({
                          MyTree::Edge(1, 2, 1),
                          MyTree::Edge(1, 3, 1),
                          MyTree::Edge(2, 4, 1),
                          MyTree::Edge(2, 5, 1),
                          MyTree::Edge(3, 6, 1),
                          MyTree::Edge(3, 7, 1),
                          MyTree::Edge(7, 8, 1),
                          MyTree::Edge(5, 9, 1),
                          MyTree::Edge(5, 10, 1),
                      })
    ));
    assert(unordered_eq(
        B.getEdges(2), vector<MyTree::Edge>({
                           MyTree::Edge(2, 4, 1),
                           MyTree::Edge(2, 5, 1),
                       })
    ));
    assert(B.getEdges(8) == vector<MyTree::Edge>());
    assert(B.getChildren(2) == vector<MyTree::Node>({4, 5}));
    assert(B.getChildren(8) == vector<MyTree::Node>({}));
    assert(B.degree(1) == 2);
    assert(B.degree(2) == 3);
    assert(B.degree(8) == 1);
    assert(B.numChildren(1) == 2);
    assert(B.numChildren(2) == 2);
    assert(B.numChildren(8) == 0);

    assert(B.containsEdge(MyTree::Edge(1, 2, 1)) == true);
    assert(B.containsEdge(MyTree::Edge(1, 2, 3175368)) == false);
    assert(B.containsEdge(MyTree::Edge(1, 5, 1)) == false);
    assert(B.containsEdgeUnweighted(MyTree::Edge(1, 2, 1)) == true);
    assert(B.containsEdgeUnweighted(MyTree::Edge(1, 2, 3175368)) == true);
    assert(B.containsEdgeUnweighted(MyTree::Edge(1, 5, 1)) == false);
    assert(B.getEdge(2, 4, 69) == 1);
    assert(B.getEdge(2, 3, 69) == 69);
}

void testAncestorFunc() {
    assert(B.getParent(1) == MyTree::Edge(1, 1, 0));
    assert(B.getParent(9) == MyTree::Edge(5, 9, 1));
    assert(B.getParent(1, 0) == MyTree::Edge(1, 1, 0));
    assert(B.getParent(5, 0) == MyTree::Edge(5, 5, 0));
    assert(B.getParent(1, 100000) == MyTree::Edge(1, 1, 0));
    assert(B.getParent(8, 3) == MyTree::Edge(1, 8, 3));
    assert(B.getParent(7, 10000) == MyTree::Edge(1, 7, 2));
    assert(B.getParent(9, 1) == MyTree::Edge(5, 9, 1));
    assert(B.isAncestor(1, 2));
    assert(B.isAncestor(1, 4));
    assert(B.isAncestor(2, 10));
    assert(!B.isAncestor(2, 3));
    assert(!B.isAncestor(8, 2));
    assert(B.lca(4, 2) == 2);
    assert(B.lca(1, 1) == 1);
    assert(B.lca(9, 8) == 1);
    assert(B.lca(4, 10) == 2);
    assert(B.dist(4, 2) == 1);
    assert(B.dist(1, 1) == 0);
    assert(B.dist(9, 8) == 6);
    assert(B.dist(4, 10) == 3);
    assert(B.getPathNodes(MyTree::Edge(3, 3)) == vector<MyTree::Node>({3}));
    assert(B.getPathEdges(MyTree::Edge(3, 3)) == vector<MyTree::Edge>({}));
    assert(B.getPathNodes(MyTree::Edge(4, 6)) == vector<MyTree::Node>({4, 2, 1, 3, 6}));
    assert(
        B.getPathEdges(MyTree::Edge(4, 6)) == vector<MyTree::Edge>(
                                                  {MyTree::Edge(4, 2, 1), MyTree::Edge(2, 1, 1),
                                                   MyTree::Edge(1, 3, 1), MyTree::Edge(3, 6, 1)}
                                              )
    );
}

void testAlgos() {
    array<MyTree::Weight, 12> fromRoot{0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 0};
    assert(B.getHeightFromRoot() == fromRoot);
    array<MyTree::Weight, 12> fromLeaves{0, 3, 2, 2, 0, 1, 0, 1, 0, 0, 0, 0};
    assert(B.getHeightFromLeaves() == fromLeaves);
    assert(B.getDiameter() == 6);
}

void testOrderings() {
    assert(B.getBfsOrder() == vector<MyTree::Node>({1, 2, 3, 4, 5, 6, 7, 9, 10, 8}));
    assert(B.getInOrder() == vector<MyTree::Node>({4, 2, 9, 5, 10, 1, 6, 3, 8, 7}));
    assert(B.getPreOrder() == vector<MyTree::Node>({1, 2, 4, 5, 9, 10, 3, 6, 7, 8}));
    assert(B.getPostOrder() == vector<MyTree::Node>({4, 9, 10, 5, 2, 6, 8, 7, 3, 1}));
    assert(B.getEulerTour() == vector<MyTree::Node>({1, 2, 4, 4, 5, 9, 9, 10, 10, 5,
                                                     2, 3, 6, 6, 7, 8, 8, 7,  3,  1}));
}

void testDisplayAndModifications() {
    assert(
        repr(B) == "1\n"
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
    assert(
        repr(&B) == "      1       \n"
                    "  ┌───┴────┐  \n"
                    "  2        3  \n"
                    "┌─┴──┐    ┌┴─┐\n"
                    "4    5    6  7\n"
                    "   ┌─┴─┐     │\n"
                    "   9  10     8\n"
    );
    B.erase(9);
    B.erase(10);
    assert(
        repr(B) == "1\n"
                   "├─2\n"
                   "│ ├─4\n"
                   "│ └─5\n"
                   "└─3\n"
                   "  ├─6\n"
                   "  └─7\n"
                   "    └─8\n"
    );
    assert(TREE_PRINT_SPACING == 2);
    assert(
        repr(&B) == "    1     \n"
                    " ┌──┴──┐  \n"
                    " 2     3  \n"
                    "┌┴─┐  ┌┴─┐\n"
                    "4  5  6  7\n"
                    "         │\n"
                    "         8\n"
    );
    B.insert(MyTree::Edge(5, 9));
    B.insert(MyTree::Edge(5, 10));

    assert(
        repr(W) == "1\n"
                   "├2\n"
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
    assert(
        repr(&W) == "    1          \n"
                    "┌──┬┴────┐     \n"
                    "2  3     7     \n"
                    "   │  ┌──┴──┐  \n"
                    "   4  8     9  \n"
                    "   │      ┌─┴─┐\n"
                    "   5     10  11\n"
                    "   │           \n"
                    "   6           \n"
    );

    assert(repr(U) == "7\n");
    assert(repr(&U) == "7\n");

    assert(
        repr(P) == "1\n"
                   "└──2\n"
                   "   └──3\n"
                   "      └──4\n"
                   "         └──5\n"
    );
    assert(
        repr(&P) == "1\n"
                    "│\n"
                    "2\n"
                    "│\n"
                    "3\n"
                    "│\n"
                    "4\n"
                    "│\n"
                    "5\n"
    );

    assert(
        repr(S) == "1\n"
                   "├──2\n"
                   "├──3\n"
                   "├──4\n"
                   "└──5\n"
    );
    assert(
        repr(&S) == "    1     \n"
                    "┌──┬┴─┬──┐\n"
                    "2  3  4  5\n"
    );
}

void testTreeTypes() {
    assert(B.isHeightBalanced());
    assert(B.isNodeBalanced());
    assert(B.isBinary());
    B.insert(MyTree::Edge(5, 11));
    assert(!B.isBinary());
    B.erase(11);
    assert(!B.isPath());

    assert(U.isHeightBalanced());
    assert(U.isNodeBalanced());
    assert(U.isBinary());
    assert(U.isPath());
    assert(U.isStar());

    assert(P.isHeightBalanced());
    assert(!P.isNodeBalanced());
    assert(P.isBinary());
    assert(P.isPath());
    assert(!P.isStar());

    assert(S.isHeightBalanced());
    assert(S.isNodeBalanced());
    assert(!S.isBinary());
    assert(!S.isPath());
    assert(S.isStar());
}

int main() {
    testNodeFunc();
    testEdgeFunc();
    testAncestorFunc();
    testAlgos();
    testOrderings();
    testDisplayAndModifications();
    testTreeTypes();
}
