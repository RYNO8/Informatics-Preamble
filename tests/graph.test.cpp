#include <iostream>
#include <fstream>
#include <sstream>
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
void testDigraph() {
    using DiGraph = Graph<10, UnitEdgeWeight, int, true>;
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
    s1 >> e;
    DiGraph G1(N1, e);

    assert(!G1.isWeightedGraph());
    assert(repr(G1) == "1: 2 3\n2: 1 3 4\n3: 1\n4:\n");
    assert(G1.V() == 4 && G1.N() == 4);
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
    assert(G1.containsEdge(DiGraph::Edge(2, 1, UnitEdgeWeight())));
    assert(G1.degreeIn(4) == 1);

    G1.eraseEdge(DiGraph::Edge(1, 2));
    assert(G1.E() == 5);
    assert(G1.degreeOut(1) == 1);
    G1.insertEdge(DiGraph::Edge(1, 2));

    G1.eraseEdge(DiGraph::Edge(2, 4));
    assert(G1.getComponent(4) == DiGraph(G1, vector<DiGraph::Node>({ 4 })));
    assert(G1.getComponent(1) == DiGraph(G1, vector<DiGraph::Node>({ 1, 2, 3 })));
    assert(G1.getComponents() == vector<DiGraph>({
        DiGraph(G1, vector<DiGraph::Node>({ 1, 2, 3 })),
        DiGraph(G1, vector<DiGraph::Node>({ 4 })),
    })); // this could actually be any order but cbbs testing that
    G1.insertEdge(DiGraph::Edge(2, 4));
    assert(G1.getComponent(2) == G1);

    DiGraph induced(G1, vector<DiGraph::Node>({ 3 }));
    assert(induced.getEdges().empty() && induced.getNodes() == vector<DiGraph::Node>({ 3 }));
    assert(induced == DiGraph(vector<DiGraph::Node>({ 3 })));
    
    assert(G1 + induced == G1);
    assert(G1 + DiGraph(3, vector<DiGraph::Edge>({ DiGraph::Edge(1, 2) })) == G1);
    assert(G1 + DiGraph(5, vector<DiGraph::Edge>({ DiGraph::Edge(1, 4) })) != G1);
    assert(repr(G1.flip()) == "1: 2 3\n2: 1\n3: 1 2\n4: 2\n");

    auto a = std::array<int, 10>({
        2147483647, 
        0, 
        1, 
        1, 
        2, 
        2147483647, 
        2147483647, 
        2147483647, 
        2147483647, 
        2147483647, 
    });
    assert(G1.SSSP(1).first == a);
    
    assert(G1.isReachable(2, 4));
    G1.eraseEdge(DiGraph::Edge(2, 4));
    assert(!G1.isReachable(2, 4));
    assert(G1.shortestPath(2, 4) == vector<DiGraph::Edge>());
    G1.insertEdge(DiGraph::Edge(2, 4));

    G1.eraseEdge(DiGraph::Edge(2, 1));
    G1.eraseEdge(DiGraph::Edge(3, 1));
    assert(G1.dfsTopsort() == vector<DiGraph::Node>({ 1, 2, 4, 3 }));
    assert(G1.Kahns() == vector<DiGraph::Node>({ 1, 2, 4, 3 }));
    assert(!G1.hasCycle());
    G1.insertEdge(DiGraph::Edge(2, 1));
    G1.insertEdge(DiGraph::Edge(3, 1));

    assert(G1.SCCdfs() == vector<vector<DiGraph::Node>>({
        vector<DiGraph::Node>({ 1, 2, 3 }),
        vector<DiGraph::Node>({ 4 })
    }));
    assert(G1.Tarjan() == vector<vector<DiGraph::Node>>({
        vector<DiGraph::Node>({ 4 }),
        vector<DiGraph::Node>({ 3, 2, 1 })
    }));

    assert(G1.hasCycle());

    G1.clear();
    assert(G1.V() == 0 && G1.E() == 0);
}

/*
 * 3-----4
 *  \   /
 *   \ /
 *    2
 *   / \
 *  /   \
 * 5-----8
 *       |
 *       |
 * 7-----1
 *  \   /
 *   \ /
 *    6
 *
 * 9------10
 */
void testGraph() {
    using MyGraph = Graph<20, UnitEdgeWeight, int, false>;
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
    s2 >> e;
    MyGraph G2(N2, e);

    assert(repr(G2) == "1: 6 7 8\n2: 3 4 5 8\n3: 4\n4:\n5: 8\n6: 7\n7:\n8:\n9: 9\n10: 9\n");
    assert(G2.containsEdge(MyGraph::Edge(1, 7)));
    assert(G2.containsEdge(MyGraph::Edge(7, 1)));

    assert(G2.shortestPath(7, 5) == vector<MyGraph::Edge>({
        MyGraph::Edge(7, 1),
        MyGraph::Edge(1, 8),
        MyGraph::Edge(8, 5)
    }));

    assert(G2.eccentricity(8) == 2);
    assert(G2.eccentricity(6) == 4);
    assert(G2.diameter(1) == vector<MyGraph::Edge>({
        MyGraph::Edge(6, 1),
        MyGraph::Edge(1, 8),
        MyGraph::Edge(8, 2),
        MyGraph::Edge(2, 3)
    }));; // many valid diameters but again cbbs making generic test

    auto f = G2.FloydWarshall();
    assert(f(3, 4) == 1 && f(5, 7) == 3 && f(7, 5) == 3 && f(6, 3) == 4);
    assert(f(0, 3) == INT_MAX && f(-5, 3) == INT_MAX && f(20, 3) == INT_MAX);
    assert(f(1, 8) == 1);
    G2.eraseEdge(MyGraph::Edge(1, 8));
    assert(f(1, 8) == 1);
    G2.insertEdge(MyGraph::Edge(1, 8));

    assert(G2.bridgesDFS() == vector<MyGraph::Edge>({ 
        MyGraph::Edge(1, 8),
        MyGraph::Edge(10, 9),
    }));

    assert(G2.EdmondsKarp(3, 8) == 2);
    assert(G2.EdmondsKarp(10, 9) == 1);
    assert(G2.EdmondsKarp(9, 10) == 1);
    assert(G2.EdmondsKarp(5, 1) == 1);
    assert(G2.EdmondsKarp(2, 2) == INT_MAX);
    assert(G2.Dinics(3, 8) == 2);
    assert(G2.Dinics(10, 9) == 1);
    assert(G2.Dinics(9, 10) == 1);
    assert(G2.Dinics(5, 1) == 1);
    assert(G2.Dinics(2, 2) == INT_MAX);

    assert(G2.hasCycle());

    assert(G2.KruskalsCost() == 8);
}

void testWeightedGraph() {
    using WeightedGraph = Graph<20, int, int, false>;
    stringstream s3;
    s3 << R""""(
3 2
1 2 10
2 3 -1
)"""";
    size_t N3, M3;
    s3 >> N3 >> M3;
    vector<WeightedGraph::Edge> e(M3);
    s3 >> e;
    WeightedGraph G3(N3, e);

    assert(repr(G3) == "1: 2 (w = 10)\n2: 3 (w = -1)\n3:\n");
    assert(G3.containsEdge(WeightedGraph::Edge(1, 2, 10)));
    assert(!G3.containsEdge(WeightedGraph::Edge(1, 2, 11)));
    assert(G3.containsEdgeUnweighted(WeightedGraph::Edge(1, 2, 10)));
    assert(G3.containsEdgeUnweighted(WeightedGraph::Edge(1, 2, 11)));
}

template<size_t N, typename EdgeWeight, typename PathWeight, bool isWeighted>
void testShortestPath(Graph<N, EdgeWeight, PathWeight, isWeighted> G, array<array<int, N>, N> adj, bool testF = true) {
    for (size_t u = 0; u < N; ++u) {
        auto SSSP = G.SSSP(u).first;
        for (size_t v = 0; v < N; ++v) {
            assert(SSSP[v] == adj[u][v]);
            assert(G.shortestDist(u, v) == adj[u][v]);
        }
    }
    if (!testF) return;
    auto f = G.allShortestDist();
    for (size_t u = 0; u < N; ++u) {
        for (size_t v = 0; v < N; ++v) {
            assert(f(u, v) == adj[u][v]);
        }
    }
}

void testShortestPaths() {
    using DirectedWeightedGraph = Graph<4, int, int, true>;

    testShortestPath(
        DirectedWeightedGraph(4, vector<DirectedWeightedGraph::Edge>({
            DirectedWeightedGraph::Edge(1, 2, 10),
            DirectedWeightedGraph::Edge(2, 3, 1),
        })),
        array<array<int, 4>, 4>{
            array<int, 4>{ 0,       INT_MAX, INT_MAX, INT_MAX },
            array<int, 4>{ INT_MAX, 0,       10,      11,     },
            array<int, 4>{ INT_MAX, INT_MAX, 0,       1,      },
            array<int, 4>{ INT_MAX, INT_MAX, INT_MAX, 0,      },
        }
    );

    testShortestPath(
        DirectedWeightedGraph(4, vector<DirectedWeightedGraph::Edge>({
            DirectedWeightedGraph::Edge(1, 2, 10),
            DirectedWeightedGraph::Edge(2, 3, -11),
        })),
        array<array<int, 4>, 4>{
            array<int, 4>{ 0,       INT_MAX, INT_MAX, INT_MAX },
            array<int, 4>{ INT_MAX, 0,       10,      -1,     },
            array<int, 4>{ INT_MAX, INT_MAX, 0,       -11,    },
            array<int, 4>{ INT_MAX, INT_MAX, INT_MAX, 0,      },
        }
    );

    using UndirectedWeightedGraph = Graph<4, int, int, false>;

    testShortestPath(
        UndirectedWeightedGraph(4, vector<UndirectedWeightedGraph::Edge>({
            UndirectedWeightedGraph::Edge(1, 2, 10),
            UndirectedWeightedGraph::Edge(2, 3, 1),
        })),
        array<array<int, 4>, 4>{
            array<int, 4>{ 0,       INT_MAX, INT_MAX, INT_MAX },
            array<int, 4>{ INT_MAX, 0,       10,      11,     },
            array<int, 4>{ INT_MAX, 10,      0,       1,      },
            array<int, 4>{ INT_MAX, 11,      1,       0,      },
        }
    );

    testShortestPath(
        UndirectedWeightedGraph(4, vector<UndirectedWeightedGraph::Edge>({
        UndirectedWeightedGraph::Edge(1, 2, 10),
            UndirectedWeightedGraph::Edge(2, 3, -11),
        })),
        array<array<int, 4>, 4>{
            array<int, 4>{ 0,       INT_MAX, INT_MAX, INT_MAX },
            array<int, 4>{ INT_MAX, 0,       10,      -1,     },
            array<int, 4>{ INT_MAX, 10,      0,       -11,    },
            array<int, 4>{ INT_MAX, -1,      -11,     0,      },
        },
        false
    );

    using DirectedUnweightedGraph = Graph<4, UnitEdgeWeight, int, true>;

    testShortestPath(
        DirectedUnweightedGraph(4, vector<DirectedUnweightedGraph::Edge>({
            DirectedUnweightedGraph::Edge(1, 2),
            DirectedUnweightedGraph::Edge(2, 3),
        })),
        array<array<int, 4>, 4>{
            array<int, 4>{ 0,       INT_MAX, INT_MAX, INT_MAX },
            array<int, 4>{ INT_MAX, 0,       1,       2,     },
            array<int, 4>{ INT_MAX, INT_MAX, 0,       1,      },
            array<int, 4>{ INT_MAX, INT_MAX, INT_MAX, 0,      },
        }
    );

    using UndirectedUnweightedGraph = Graph<4, UnitEdgeWeight, int, false>;
    
    testShortestPath(
        UndirectedUnweightedGraph(4, vector<UndirectedUnweightedGraph::Edge>({
            UndirectedUnweightedGraph::Edge(1, 2),
            UndirectedUnweightedGraph::Edge(2, 3),
        })),
        array<array<int, 4>, 4>{
            array<int, 4>{ 0,       INT_MAX, INT_MAX, INT_MAX },
            array<int, 4>{ INT_MAX, 0,       1,       2,     },
            array<int, 4>{ INT_MAX, 1,       0,       1,      },
            array<int, 4>{ INT_MAX, 2,       1,       0,      },
        }
    );
}

int main() {
    testDigraph();
    testGraph();
    testWeightedGraph();
    testShortestPaths();

    // TODO: https://muscat2023b.contest.codeforces.com/group/cuULm9FF5q/contest/432252/problem/P
}