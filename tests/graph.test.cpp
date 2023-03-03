#include <iostream>
#include <fstream>
#include <sstream>
#include "../src/Graph.h"
#include "../src/Util.h"
using namespace std;
using namespace DS;

using DiGraph = Graph<20, UnitWeight, true>;
using MyGraph = Graph<20, UnitWeight, false>;
using WeightedGraph = Graph<20, int, false>;

void testDigraph() {
    stringstream s1;
    s1 << R""""(
4 6
1 2
2 1
1 3
3 1
2 3
2 4
)"""";
    size_t N1, M1;
    s1 >> N1 >> M1;
    vector<DiGraph::Edge> e(M1);
    for (size_t i = 0; i < M1; ++i) s1 >> e[i];
    DiGraph G1(N1, e);


    assert(repr(G1) == "1: 2 3\n2: 1 3 4\n3: 1\n4:\n");
    assert(G1.size() == 4 && G1.V() == 4 && G1.N() == 4);
    assert(G1.E() == 6 && G1.M() == 6);
    assert(G1.getNodes() == vector<DiGraph::Node>({1, 2, 3, 4}));
    assert(G1.getEdges() == vector<DiGraph::Edge>({
        DiGraph::Edge(1, 2),
        DiGraph::Edge(1, 3),
        DiGraph::Edge(2, 1),
        DiGraph::Edge(2, 3),
        DiGraph::Edge(2, 4),
        DiGraph::Edge(3, 1),
    }));
    assert(G1.getEdgesOut(3) == vector<DiGraph::Edge>({ DiGraph::Edge(3, 1) }));
    assert(G1.containsEdge(DiGraph::Edge(1, 2)));
    assert(G1.containsEdge(DiGraph::Edge(2, 1)));
    assert(!G1.containsEdge(DiGraph::Edge(3, 4)));
    assert(G1.containsEdge(DiGraph::Edge(2, 1, UnitWeight())));
    assert(G1.degreeIn(4) == 1);

    G1.eraseEdge(DiGraph::Edge(1, 2));
    assert(G1.E() == 5);
    assert(G1.degreeOut(1) == 1);
    G1.insertEdge(DiGraph::Edge(1, 2));

    G1.eraseEdge(DiGraph::Edge(2, 4));
    assert(G1.getComponent(4) == vector<DiGraph::Node>({ 4 }));
    vector<DiGraph::Node> comp = G1.getComponent(1);
    sort(comp.begin(), comp.end());
    assert(comp == vector<DiGraph::Node>({ 1, 2, 3 }));
    assert(G1.getComponents() == vector<vector<DiGraph::Node>>({
        vector<DiGraph::Node>({ 1, 2, 3 }),
        vector<DiGraph::Node>({ 4 }),
    })); // this could actually be any order by cbbs testing that
    G1.insertEdge(DiGraph::Edge(2, 4));
    assert(G1.getComponent(2) == G1.getNodes());

    DiGraph induced(G1, vector<DiGraph::Node>({ 3 }));
    assert(induced.getEdges() == vector<DiGraph::Edge>({
        DiGraph::Edge(1, 3),
        DiGraph::Edge(2, 3),
        DiGraph::Edge(3, 1),
    }));
    
    assert(G1 + induced == G1);
    assert(G1 + DiGraph(2, vector<DiGraph::Edge>({ DiGraph::Edge(1, 2) })) == G1);
    assert(G1 + DiGraph(2, vector<DiGraph::Edge>({ DiGraph::Edge(1, 4) })) != G1);
    assert(repr(G1.flip()) == "1: 2 3\n2:\n3: 1 2\n4: 2\n");

    // G1.addEdge(1, 2);
    // assert(G1.numEdges(1, 2) == 2);
    // assert(G1.numEdges(2, 1) == 1);
    // assert(G1.numEdges(1, 4) == 0);
    // G1.removeEdge(1, 2);
    // assert(G1.getEdgesOut(1) == multiset<size_t>({2, 3}));
    // vector<int> a = G1.greedyColouring();
    // assert(a == vector<int>({-1, 0, 1, 2, 0}));
    G1.clear();
    assert(G1.V() == 0 && G1.E() == 0);
}

void testGraph() {
    stringstream s2;
    s2 << R""""(
10 11
1 7
1 8
1 6
2 8
6 7
5 8
2 5
2 3
2 4
3 4
10 9
)"""";
    size_t N2, M2;
    s2 >> N2 >> M2;
    vector<MyGraph::Edge> e(M2);
    for (size_t i = 0; i < M2; ++i) s2 >> e[i];
    MyGraph G2(N2, e);

    assert(repr(G2) == "1: 6 7 8\n2: 3 4 5 8\n3: 4\n4:\n5: 8\n6: 7\n7:\n8:\n9: 9\n10: 9\n");
    assert(G2.containsEdge(MyGraph::Edge(1, 7)));
    assert(G2.containsEdge(MyGraph::Edge(7, 1)));
}

void testWeightedGraph() {
    stringstream s3;
    s3 << R""""(
3 2
1 2 10
2 3 0
)"""";
    size_t N3, M3;
    s3 >> N3 >> M3;
    vector<WeightedGraph::Edge> e(M3);
    for (size_t i = 0; i < M3; ++i) s3 >> e[i];
    WeightedGraph G3(N3, e);

    assert(repr(G3) == "1: 2 (w = 10)\n2: 3 (w = 0)\n3:\n");
    assert(G3.containsEdge(WeightedGraph::Edge(1, 2, 10)));
    assert(!G3.containsEdge(WeightedGraph::Edge(1, 2, 11)));
    assert(G3.containsEdgeUnweighted(WeightedGraph::Edge(1, 2, 10)));
    assert(G3.containsEdgeUnweighted(WeightedGraph::Edge(1, 2, 11)));
}

int main() {
    // testDigraph();
    testGraph();
    testWeightedGraph();
}