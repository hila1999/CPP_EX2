
#include "doctest.h"
#include "Algorithms.hpp"
#include "Graph.hpp"
#include <sstream>
using namespace std;
TEST_CASE("Test graph addition")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    vector<vector<int>> weightedGraph = {
        {0, 1, 1},
        {1, 0, 2},
        {1, 2, 0}};
    g2.loadGraph(weightedGraph);
    ariel::Graph g3 = g1 + g2;
    vector<vector<int>> expectedGraph = {
        {0, 2, 1},
        {2, 0, 3},
        {1, 3, 0}};
    CHECK(g3.getAdjacencyMatrix() == expectedGraph);
}

// TEST_CASE("Test graph addition with different dimensions")
// {
//     ariel::Graph g1;
//     vector<vector<int>> graph = {
//         {0, 1, 0},
//         {1, 0, 1},
//         {0, 1, 0}};
//     g1.loadGraph(graph);
//     ariel::Graph g2;
//     vector<vector<int>> weightedGraph = {
//         {0, 1, 1, 1},
//         {1, 0, 2, 1},
//         {1, 2, 0, 1}};
//     g2.loadGraph(weightedGraph);
//     CHECK_THROWS_AS(g1 + g2, std::runtime_error);
// }
TEST_CASE("Test graph addition with empty graph")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    CHECK_THROWS(g1 + g2);
    CHECK_THROWS(g2 + g1);}

TEST_CASE("Test graph addition with zero matrix")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    vector<vector<int>> weightedGraph = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    g2.loadGraph(weightedGraph);
    CHECK((g1 + g2).getAdjacencyMatrix() == graph);
    CHECK((g2 + g1).getAdjacencyMatrix() == graph);
}

TEST_CASE("Test graph addition with negative values")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, -1, 0},
        {-1, 0, -1},
        {0, -1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    vector<vector<int>> weightedGraph = {
        {0, -1, -1},
        {-1, 0, -2},
        {-1, -2, 0}};
    g2.loadGraph(weightedGraph);
    vector<vector<int>> expectedGraph = {
        {0, -2, -1},
        {-2, 0, -3},
        {-1, -3, 0}};
    CHECK((g1 + g2).getAdjacencyMatrix() == expectedGraph);
    CHECK((g2 + g1).getAdjacencyMatrix() == expectedGraph);
}
TEST_CASE("Test graph addition assignment operator")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, 1, 1},
        {1, 0, 2},
        {1, 2, 0}};
    g2.loadGraph(graph2);
    g1 += g2;
    vector<vector<int>> expectedGraph = {
        {0, 2, 1},
        {2, 0, 3},
        {1, 3, 0}};
    CHECK(g1.getAdjacencyMatrix() == expectedGraph);
}

TEST_CASE("Test graph addition assignment operator with empty graph")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    
    CHECK_THROWS( g1 += g2);
}

TEST_CASE("Test graph addition assignment operator with zero matrix")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    CHECK_THROWS( g1 += g2);
}

TEST_CASE("Test graph addition assignment operator with negative values")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, -1, 0},
        {-1, 0, -1},
        {0, -1, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, -1, -1},
        {-1, 0, -2},
        {-1, -2, 0}};
    g2.loadGraph(graph2);
    g1 += g2;
    vector<vector<int>> expectedGraph = {
        {0, -2, -1},
        {-2, 0, -3},
        {-1, -3, 0}};
    CHECK(g1.getAdjacencyMatrix() == expectedGraph);
}
TEST_CASE("Test graph operator+")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2 = +g1;
    CHECK(g2.getAdjacencyMatrix() == graph);
}


TEST_CASE("Test graph operator+ with zero matrix")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2 = +g1;
    CHECK(g2.getAdjacencyMatrix() == graph);
}

TEST_CASE("Test graph operator+ with negative values")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, -1, 0},
        {-1, 0, -1},
        {0, -1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2 = +g1;
    CHECK(g2.getAdjacencyMatrix() == graph);
}
TEST_CASE("Test graph subtraction")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, 1, 1},
        {1, 0, 2},
        {1, 2, 0}};
    g2.loadGraph(graph2);
    ariel::Graph g3 = g1 - g2;
    vector<vector<int>> expectedGraph = {
        {0, 0, -1},
        {0, 0, -1},
        {-1, -1, 0}};
    CHECK(g3.getAdjacencyMatrix() == expectedGraph);
}

TEST_CASE("Test graph subtraction with empty graph")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    CHECK_THROWS(g1 - g2);
    CHECK_THROWS(g2 - g1);
}

TEST_CASE("Test graph subtraction with zero matrix")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    vector<vector<int>> weightedGraph = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    g2.loadGraph(weightedGraph);
    CHECK((g1 - g2).getAdjacencyMatrix() == graph);
    CHECK((g2 - g1).getAdjacencyMatrix() == graph);
}

TEST_CASE("Test graph subtraction with negative values")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, -1, 0},
        {-1, 0, -1},
        {0, -1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    vector<vector<int>> weightedGraph = {
        {0, -1, -1},
        {-1, 0, -2},
        {-1, -2, 0}};
    g2.loadGraph(weightedGraph);
    vector<vector<int>> expectedGraph = {
        {0, 0, 1},
        {0, 0, 1},
        {1, 1, 0}};
    CHECK((g1 - g2).getAdjacencyMatrix() == expectedGraph);
    
}
TEST_CASE("Test graph negation")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    ariel::Graph negatedGraph = -g1;
    vector<vector<int>> expectedGraph = {
        {0, -1, 0},
        {-1, 0, -1},
        {0, -1, 0}};
    CHECK(negatedGraph.getAdjacencyMatrix() == expectedGraph);
}

TEST_CASE("Test graph negation with zero matrix")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    g1.loadGraph(graph);
    ariel::Graph negatedGraph = -g1;
    vector<vector<int>> expectedGraph = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    CHECK(negatedGraph.getAdjacencyMatrix() == expectedGraph);
}

TEST_CASE("Test graph negation with negative values")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, -1, 0},
        {-1, 0, -1},
        {0, -1, 0}};
    g1.loadGraph(graph);
    ariel::Graph negatedGraph = -g1;
    vector<vector<int>> expectedGraph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    CHECK(negatedGraph.getAdjacencyMatrix() == expectedGraph);
}
TEST_CASE("Test graph equality with same graph")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    vector<vector<int>> weightedGraph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g2.loadGraph(weightedGraph);
    CHECK(g1 == g2);
    CHECK(g2 == g1);
}

TEST_CASE("Test graph equality with different graphs")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, 1, 0},
        {1, 0, 0},
        {0, 0, 0}};
    g2.loadGraph(graph2);
    CHECK_FALSE(g1 == g2);
    CHECK_FALSE(g2 == g1);
}

TEST_CASE("Test graph equality with empty graphs")
{
    ariel::Graph g1;
    ariel::Graph g2;
    CHECK(g1 == g2);
    CHECK(g2 == g1);
}
TEST_CASE("Test graph inequality with different graphs")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph1);

    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 1}};
    g2.loadGraph(graph2);

    CHECK(g1 != g2);
    CHECK(g2 != g1);
}

TEST_CASE("Test graph inequality with same graph")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);

    ariel::Graph g2;
    g2.loadGraph(graph);

    CHECK_FALSE(g1 != g2);
    CHECK_FALSE(g2 != g1);
}

TEST_CASE("Test graph inequality with empty graph")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    g1.loadGraph(graph);

    ariel::Graph g2;

    CHECK(g1 != g2);
    CHECK(g2 != g1);
}
TEST_CASE("Test operator< with contained graph")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, 1, 1},
        {1, 0, 1},
        {1, 1, 0}};
    g2.loadGraph(graph2);
    CHECK(g1 < g2);
    CHECK_FALSE(g2 < g1);
}

TEST_CASE("Test operator< with different number of edges")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, 1, 1,0},
        {1, 0, 1,0},
        {1, 1, 0,0},
        {0, 0, 0,0}};
    g2.loadGraph(graph2);
    CHECK(g1 < g2);
    CHECK_FALSE(g2 < g1);
}

TEST_CASE("Test operator< with same number of edges but different order of adjacency matrices")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, 1, 1},
        {1, 0, 1},
        {0, 1, 0}};
    g2.loadGraph(graph2);
    CHECK(g1 < g2);
    CHECK_FALSE(g2 < g1);
}

// TEST_CASE("Test operator< with same number of edges and same order of adjacency matrices")
// {
//     ariel::Graph g1;
//     vector<vector<int>> graph1 = {
//         {0, 1, 0},
//         {1, 0, 1},
//         {0, 1, 0}};
//     g1.loadGraph(graph1);
//     ariel::Graph g2;
//     vector<vector<int>> graph2 = {
//         {0, 1, 0},
//         {1, 0, 1},
//         {0, 1, 0}};
//     g2.loadGraph(graph2);
//     CHECK_FALSE(g1 < g2);
//     CHECK_FALSE(g2 < g1);
// }

// TEST_CASE("Test graph comparison using operator<=")
// {
//     ariel::Graph g1;
//     vector<vector<int>> graph1 = {
//         {0, 1, 0},
//         {1, 0, 1},
//         {0, 1, 0}};
//     g1.loadGraph(graph1);

//     ariel::Graph g2;
//     vector<vector<int>> graph2 = {
//         {0, 1, 1},
//         {1, 0, 2},
//         {1, 2, 0}};
//     g2.loadGraph(graph2);

//     ariel::Graph g3;
//     vector<vector<int>> graph3 = {
//         {0, 1, 0},
//         {1, 0, 1},
//         {0, 1, 0}};
//     g3.loadGraph(graph3);

//     ariel::Graph g4;
//     vector<vector<int>> graph4 = {
//         {0, 0, 0},
//         {0, 0, 0},
//         {0, 0, 0}};
//     g4.loadGraph(graph4);

//     ariel::Graph g5;
//     vector<vector<int>> graph5 = {
//         {0, -1, 0},
//         {-1, 0, -1},
//         {0, -1, 0}};
//     g5.loadGraph(graph5);

//     CHECK(g1 <= g2); // g1 is less than g2
//     // CHECK_FALSE(g2 <= g3); // g2 is equal to g3 //check if it is correct
//     CHECK(g3 <= g1); // g3 is greater than g1
//     CHECK(g4 <= g4); // g4 is equal to itself
//     CHECK(g5 <= g1); // g5 is less than g1
//     CHECK(g5 <= g2); // g5 is less than g2
//     CHECK(g5 <= g3); // g5 is less than g3
//     CHECK(g5 <= g4); // g5 is less than g4
// }
TEST_CASE("Test graph greater than operator")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, 1, 1},
        {1, 0, 2},
        {1, 2, 0}};
    g2.loadGraph(graph2);
    CHECK(g2 > g1);
    CHECK_FALSE(g1 > g2);
}

TEST_CASE("Test graph greater than operator with zero matrix")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    CHECK_THROWS(g1 > g2);
}

TEST_CASE("Test graph greater than operator with negative values")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, -1, 0},
        {-1, 0, -1},
        {0, -1, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, -1, -1},
        {-1, 0, -2},
        {-1, -2, 0}};
    g2.loadGraph(graph2);
    CHECK(g2 > g1);
    CHECK_FALSE(g1 > g2);
}TEST_CASE("Test graph comparison (>=) with equal graphs")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    vector<vector<int>> weightedGraph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g2.loadGraph(weightedGraph);
    CHECK((g1 >= g2) == true);
    CHECK((g2 >= g1) == true);
}

TEST_CASE("Test graph comparison (>=) with different graphs")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    vector<vector<int>> weightedGraph = {
        {0, 1, 1},
        {1, 0, 2},
        {1, 2, 0}};
    g2.loadGraph(weightedGraph);
    CHECK((g1 >= g2) == false);
    CHECK((g2 >= g1) == true);
}

TEST_CASE("Test graph comparison (>=) with empty graphs")
{
    ariel::Graph g1;
    ariel::Graph g2;
    CHECK((g1 >= g2) == true);
    CHECK((g2 >= g1) == true);
}

TEST_CASE("Test graph comparison (>=) with zero matrix")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    vector<vector<int>> weightedGraph = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    g2.loadGraph(weightedGraph);
    CHECK((g1 >= g2) == true);
    CHECK((g2 >= g1) == true);
}

TEST_CASE("Test graph comparison (>=) with negative values")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, -1, 0},
        {-1, 0, -1},
        {0, -1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    vector<vector<int>> weightedGraph = {
        {0, -1, -1},
        {-1, 0, -2},
        {-1, -2, 0}};
    g2.loadGraph(weightedGraph);
    CHECK((g1 >= g2) == false);
    CHECK((g2 >= g1) == true);
}TEST_CASE("Test graph increment operator")
{
    ariel::Graph g;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g.loadGraph(graph);
    
    // Apply the increment operator
    ++g;
    
    vector<vector<int>> expectedGraph = {
        {1, 2, 1},
        {2, 1, 2},
        {1, 2, 1}};
    
    CHECK(g.getAdjacencyMatrix() == expectedGraph);
}TEST_CASE("Test decrement operator on graph")
{
    ariel::Graph g;
    vector<vector<int>> graph = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}};
    g.loadGraph(graph);
    
    --g;
    
    vector<vector<int>> expectedGraph = {
        {0, 1, 2},
        {3, 4, 5},
        {6, 7, 8}};
    
    CHECK(g.getAdjacencyMatrix() == expectedGraph);
}

TEST_CASE("Test decrement operator on empty graph")
{
    ariel::Graph g;
    
    --g;
    
    vector<vector<int>> expectedGraph;
    
    CHECK_THROWS(g.getAdjacencyMatrix());
}

TEST_CASE("Test decrement operator on graph with negative values")
{
    ariel::Graph g;
    vector<vector<int>> graph = {
        {-1, -2, -3},
        {-4, -5, -6},
        {-7, -8, -9}};
    g.loadGraph(graph);
    
    --g;
    
    vector<vector<int>> expectedGraph = {
        {-2, -3, -4},
        {-5, -6, -7},
        {-8, -9, -10}};
    
    CHECK(g.getAdjacencyMatrix() == expectedGraph);
}TEST_CASE("Test graph division by non-zero scalar")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    int scalar = 2;
    g1 /= scalar;
    vector<vector<int>> expectedGraph = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    CHECK(g1.getAdjacencyMatrix() == expectedGraph);
}

TEST_CASE("Test graph division by zero scalar")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    int scalar = 0;
    g1 /= scalar;
    // Division by zero should not modify the graph
    CHECK(g1.getAdjacencyMatrix() == graph);
}TEST_CASE("Test graph scalar multiplication")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    int scalar = 2;
    g1 *= scalar;
    vector<vector<int>> expectedGraph = {
        {0, 2, 0},
        {2, 0, 2},
        {0, 2, 0}};
    CHECK(g1.getAdjacencyMatrix() == expectedGraph);
}

TEST_CASE("Test graph scalar multiplication with zero matrix")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    g1.loadGraph(graph);
    int scalar = 5;
    g1 *= scalar;
    CHECK(g1.getAdjacencyMatrix() == graph);
}

TEST_CASE("Test graph scalar multiplication with negative values")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, -1, 0},
        {-1, 0, -1},
        {0, -1, 0}};
    g1.loadGraph(graph);
    int scalar = -2;
    g1 *= scalar;
    vector<vector<int>> expectedGraph = {
        {0, 2, 0},
        {2, 0, 2},
        {0, 2, 0}};
    CHECK(g1.getAdjacencyMatrix() == expectedGraph);
}TEST_CASE("Test graph multiplication")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, 1, 1},
        {1, 0, 2},
        {1, 2, 0}};
    g2.loadGraph(graph2);
    ariel::Graph g3 = g1 * g2;
    vector<vector<int>> expectedGraph = {
        {0, 0, 2},
        {1, 0, 1},
        {1, 0, 0}};
    CHECK(g3.getAdjacencyMatrix() == expectedGraph);
}

// TEST_CASE("Test graph multiplication with empty graph")
// {
//     ariel::Graph g1;
//     vector<vector<int>> graph1 = {
//         {0, 1, 0},
//         {1, 0, 1},
//         {0, 1, 0}};
//     g1.loadGraph(graph1);
//     ariel::Graph g2;
//     CHECK((g1 * g2).getAdjacencyMatrix().empty());
//     CHECK((g2 * g1).getAdjacencyMatrix().empty());
// }

TEST_CASE("Test graph multiplication with zero matrix")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    g2.loadGraph(graph2);
    vector<vector<int>> expectedGraph = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}};
    CHECK((g1 * g2).getAdjacencyMatrix() == expectedGraph);
    CHECK((g2 * g1).getAdjacencyMatrix() == expectedGraph);
}


TEST_CASE("Test operator<< with non-empty graph")
{
    ariel::Graph g;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g.loadGraph(graph);
    stringstream ss;
    ss << g;
    string expectedOutput = "[0, 1, 0]\n[1, 0, 1]\n[0, 1, 0]\n";
    CHECK(ss.str() == expectedOutput);
}

TEST_CASE("Test operator<< with negative values")
{
    ariel::Graph g;
    vector<vector<int>> graph = {
        {0, -1, 0},
        {-1, 0, -1},
        {0, -1, 0}};
    g.loadGraph(graph);
    stringstream ss;
    ss << g;
    string expectedOutput = "[0, -1, 0]\n[-1, 0, -1]\n[0, -1, 0]\n";
    CHECK(ss.str() == expectedOutput);
}

TEST_CASE("== operator")
{
    ariel::Graph g1;
    vector<vector<int>> graph = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph);
    ariel::Graph g2;
    vector<vector<int>> weightedGraph = {
        {0, 1, 1},
        {1, 0, 2},
        {1, 2, 0}};
    g2.loadGraph(weightedGraph);
    CHECK(g1 == g1);
    CHECK(g1 != g2);

    ariel::Graph g3;
    vector<vector<int>> graph1 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g3.loadGraph(graph);
    CHECK(g1 == g3);
}

TEST_CASE("Test < operator")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g2.loadGraph(graph2);
    CHECK_FALSE(g1 < g2);

    ariel::Graph g3;
    vector<vector<int>> graph3 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g3.loadGraph(graph3);
    ariel::Graph g4;
    vector<vector<int>> graph4 = {
        {0, 1, 0, 0, 1},
        {1, 0, 1, 0, 0},
        {0, 1, 0, 1, 0},
        {0, 0, 1, 0, 1},
        {1, 0, 0, 1, 0}};
    g4.loadGraph(graph4);
    CHECK(g3 < g4);

    // Test if one graph has more edges than the other
    ariel::Graph g5;
    vector<vector<int>> graph5 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g5.loadGraph(graph5);
    ariel::Graph g6;
    vector<vector<int>> graph6 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 0, 0}};
    g6.loadGraph(graph6);
    CHECK_FALSE(g5 < g6);

    // Test if one graph has more vertices than the other
    ariel::Graph g7;
    vector<vector<int>> graph7 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g7.loadGraph(graph7);
    ariel::Graph g8;
    vector<vector<int>> graph8 = {
        {0, 1, 0, 0},
        {1, 0, 0, 0},
        {0, 0, 0, 1},
        {0, 0, 1, 0}};
    g8.loadGraph(graph8);
    CHECK(g7 < g8);
}

TEST_CASE("Test != operator")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g2.loadGraph(graph2);
    CHECK_FALSE(g1 != g2);

    ariel::Graph g3;
    vector<vector<int>> graph3 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g3.loadGraph(graph3);
    ariel::Graph g4;
    vector<vector<int>> graph4 = {
        {0, 1, 0, 0, 1},
        {1, 0, 1, 0, 0},
        {0, 1, 0, 1, 0},
        {0, 0, 1, 0, 1},
        {1, 0, 0, 1, 0}};
    g4.loadGraph(graph4);
    CHECK(g3 != g4);

    ariel::Graph g5;
    vector<vector<int>> graph5 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g5.loadGraph(graph5);
    ariel::Graph g6;
    vector<vector<int>> graph6 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 0, 0}};
    g6.loadGraph(graph6);
    CHECK(g5 != g6);
}

TEST_CASE("Test <= operator")
{
    ariel::Graph g1;
    vector<vector<int>> graph1 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g1.loadGraph(graph1);
    ariel::Graph g2;
    vector<vector<int>> graph2 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g2.loadGraph(graph2);
    CHECK(g1 <= g2);

    ariel::Graph g3;
    vector<vector<int>> graph3 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g3.loadGraph(graph3);
    ariel::Graph g4;
    vector<vector<int>> graph4 = {
        {0, 1, 0, 0, 1},
        {1, 0, 1, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 1},
        {1, 0, 0, 1, 0}};
    g4.loadGraph(graph4);
    CHECK(g3 <= g4);

    ariel::Graph g5;
    vector<vector<int>> graph5 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 1, 0}};
    g5.loadGraph(graph5);
    ariel::Graph g6;
    vector<vector<int>> graph6 = {
        {0, 1, 0},
        {1, 0, 1},
        {0, 0, 0}};
    g6.loadGraph(graph6);
    CHECK_FALSE(g5 <= g6);
}