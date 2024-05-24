//hila.shamir99@gmail.com  314906983
#include <vector>
#include <stack>
#include <queue>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include "Graph.hpp"
#include "Algorithms.hpp"
#include <climits> // for INT_MAX
#include <cstddef> // for SIZE_MAX
using namespace std;
using namespace ariel;
  void Algorithms::dfs(const vector<std::vector<int>>& adjacencyMatrix, vector<bool>& visited, int vertex) {
    // Mark the current vertex as visited
    visited[static_cast<std::vector<bool>::size_type>(vertex)] = true;

    // Iterate through all possible vertices to check adjacency
    for (size_t i = 0; i < adjacencyMatrix[static_cast<std::vector<std::vector<int>>::size_type>(vertex)].size(); ++i) {
        // If there is an edge to vertex 'i' and it has not been visited yet
        if (adjacencyMatrix[static_cast<vector<std::vector<int>>::size_type>(vertex)][i] != 0 && !visited[i]) {
            // Recursively visit the adjacent vertex
            dfs(adjacencyMatrix, visited, static_cast<int>(i));
        }
    }
}

bool Algorithms::isConnected(const Graph& graph) {
    // Retrieve the adjacency matrix of the graph
    const std::vector<std::vector<int>>& adjacencyMatrix = graph.getAdjacencyMatrix();

    // If the graph is empty (no vertices), it is not connected
    if (adjacencyMatrix.empty()) {
        return false;
    }

    // Create a vector to track visited vertices, initialized to false
    std::vector<bool> visited(adjacencyMatrix.size(), false);

    // Start a depth-first search from the first vertex (vertex 0)
    dfs(adjacencyMatrix, visited, 0);

    // Check if all vertices have been visited
    if (std::all_of(visited.begin(), visited.end(), [](bool isVisited) { return isVisited; })) {
        // If all vertices are visited, the graph is connected
        return true;
    }

    // If any vertex is not visited, the graph is not connected
    return false;
}
string Algorithms::shortestPath(Graph& graph, int start, int end) {
    // If the graph is empty, return "-1" indicating no path
    if(graph.size() == 0) {
        return "-1";
    }
    
    size_t numOfVer = graph.size(); // Number of vertices in the graph
    vector<int> dist(numOfVer, INT_MAX); // Initialize distances to all vertices as infinite
    vector<int> predecessor(numOfVer, -1); // Initialize predecessors for all vertices
    dist[(size_t)start] = 0; // Distance to the start vertex is 0

    vector<vector<int>> adjacencyMatrix = graph.getAdjMatrix(); // Get the adjacency matrix of the graph

    // Relax all edges numOfVer-1 times
    for (size_t uar = 0; uar < numOfVer - 1; ++uar) {
        for (size_t jar = 0; jar < numOfVer; jar++) {
            for (size_t kar = 0; kar < numOfVer; kar++) {
                // If there is an edge and the new path is shorter, update the distance and predecessor
                if (adjacencyMatrix[jar][kar] != 0 && dist[jar] != INT_MAX && dist[jar] + adjacencyMatrix[jar][kar] < dist[kar]) {
                    dist[kar] = dist[jar] + adjacencyMatrix[jar][kar];
                    predecessor[kar] = static_cast<int>(jar);
                }
            }
        }
    }

    // Check for negative weight cycles
    for (size_t jar = 0; jar < numOfVer; jar++) {
        for (size_t k = 0; k < numOfVer; k++) {
            if (adjacencyMatrix[jar][k] != 0 && dist[jar] != INT_MAX && dist[jar] + adjacencyMatrix[jar][k] < dist[k]) {
                return "Negative cycle detected"; // If a cycle is detected, return
            }
        }
    }

    // If there is no path from start to end
    if (dist[(size_t)end] == INT_MAX) {
        return "-1";
    }

    // Reconstruct the path from end to start using predecessors
    string path = to_string(end);
    for (size_t var = (size_t)end; var != start; var = (size_t)predecessor[var]) {
        path = to_string((int)predecessor[var]) + "->" + path;
    }

    return path;
}

std::string Algorithms::isBipartite(const Graph &g) {
    size_t numNodes = g.size(); // Number of nodes in the graph
    vector<vector<int>> adjMatrix = g.getAdjacencyMatrix(); // Get adjacency matrix

    // Make the graph undirected if it is directed
    if (g.isDirected()) {
        for (size_t i = 0; i < numNodes; ++i) {
            for (size_t jar = 0; jar < numNodes; ++jar) {
                if (adjMatrix[i][jar] != 0 && adjMatrix[jar][i] == 0) {
                    adjMatrix[jar][i] = adjMatrix[i][jar];
                } else if (adjMatrix[i][jar] == 0 && adjMatrix[jar][i] != 0) {
                    adjMatrix[i][jar] = adjMatrix[jar][i];
                }
            }
        }
    }

    vector<int> colorArr(numNodes, -1); // Initialize color array for nodes
    std::vector<size_t> setA, setB; // Two sets for bipartite graph

    // BFS to check if graph is bipartite
    for (size_t i = 0; i < numNodes; i++) {
        if (colorArr[i] == -1) {
            std::queue<size_t> q;
            q.push(i);
            colorArr[i] = 1; // Assign first color

            while (!q.empty()) {
                size_t u = q.front();
                q.pop();
                std::vector<size_t> neighbors = g.getNeighbors(u); // Get neighbors of u

                for (size_t v : neighbors) {
                    if (colorArr[v] == -1) {
                        colorArr[v] = 1 - colorArr[u]; // Assign alternate color
                        q.push(v);
                    } else if (colorArr[v] == colorArr[u]) {
                        return "0"; // If same color, graph is not bipartite
                    }
                }
            }
        }
    }

    // Separate nodes into two sets
    for (size_t i = 0; i < numNodes; i++) {
        if (colorArr[i] == 1) {
            setA.push_back(i);
        } else {
            setB.push_back(i);
        }
    }

    // Format the result string
    std::string result = "The graph is bipartite: A={";
    for (size_t i = 0; i < setA.size(); i++) {
        result += std::to_string(setA[i]);
        if (i != setA.size() - 1) {
            result += ", ";
        }
    }
    result += "}, B={";
    for (size_t i = 0; i < setB.size(); i++) {
        result += std::to_string(setB[i]);
        if (i != setB.size() - 1) {
            result += ", ";
        }
    }
    result += "}";
    return result;
}

bool Algorithms::isCyclic(vector<bool>& visited, size_t v, std::vector<bool>& stack, std::vector<size_t>& parent, const Graph& g, std::vector<size_t>& cyclePath) {
    visited[v] = true; // Mark the current node as visited
    stack[v] = true; // Add the current node to the recursion stack
    bool isDirected = g.isDirected(); // Check if the graph is directed
    size_t numNodes = g.size(); // Number of nodes in the graph

    // Check all adjacent vertices
    for (size_t i = 0; i < numNodes; i++) {
        if (g.haveEdge(v, i)) { // If there's an edge from v to i
            if (!visited[i]) { // If the adjacent vertex is not visited
                parent[i] = v; // Set the parent of the adjacent vertex
                if (isCyclic(visited, i, stack, parent, g, cyclePath)) {
                    return true; // If a cycle is detected, return true
                }
            } else if ((isDirected && stack[i]) || (!isDirected && stack[i] && parent[v] != i)) {
                // If a cycle is detected, trace back to print the cycle
                cyclePath.push_back(i);
                for (size_t p = v; p != i; p = parent[p]) {
                    cyclePath.push_back(p);
                }
                cyclePath.push_back(i); // Complete the cycle by adding the start node again
                std::reverse(cyclePath.begin(), cyclePath.end());

                // Print the cycle path
                std::cout << "Cycle detected: ";
                for (size_t jar = 0; jar < cyclePath.size(); jar++) {
                    std::cout << cyclePath[jar];
                    if (jar < cyclePath.size() - 1) std::cout << " -> ";
                }
                std::cout << std::endl;
                return true;
            }
        }
    }

    stack[v] = false; // Remove the vertex from the recursion stack
    return false;
}

bool Algorithms::isContainsCycle(const Graph& g) {
    size_t numNodes = g.size(); // Number of nodes in the graph
    vector<size_t> parent(numNodes, SIZE_MAX); // Parent array to store the parent of each node
    vector<size_t> cyclePath; // To store the path of the cycle
    vector<bool> visited(numNodes, false); // Visited array to mark visited nodes
    vector<bool> stack(numNodes, false); // Stack to keep track of the recursion stack

    // Check for cycles in different DFS trees
    for (size_t i = 0; i < numNodes; i++) {
        if (!visited[i]) {
            if (isCyclic(visited, i, stack, parent, g, cyclePath)) {
                return true; // If a cycle is detected, return true
            }
        }
    }
    return false; // If no cycle is detected, return false
}

std::string Algorithms::negativeCycle(const Graph& g) {
    size_t n = g.size(); // Number of nodes in the graph
    std::vector<int> distance(n, INT_MAX); // Distance vector initialized to infinite
    std::vector<int> pred(n, -1); // Predecessor array initialized to -1
    int cycle_start = -1; // Variable to store the start of the cycle

    // Assume the source vertex is 0 for checking the negative cycle
    distance[0] = 0;

    // Bellman-Ford algorithm to find shortest paths
    for (size_t i = 0; i < n; ++i) {
        cycle_start = -1;
        for (size_t u = 0; u < n; ++u) {
            for (size_t v = 0; v < n; ++v) {
                if (g.getAdjacencyMatrix()[u][v] != 0) { // Check if there is an edge from u to v
                    // Relax the edge
                    int new_distance = distance[u] + g.getAdjacencyMatrix()[u][v];
                    if (new_distance < distance[v]) { // If new distance is shorter, update distance and predecessor
                        distance[v] = new_distance;
                        pred[v] = u;
                        cycle_start = v; // Update cycle start
                    }
                }
            }
        }
    }

    vector<int> cycle; // Vector to store the cycle
    if (cycle_start != -1) { // If a cycle is detected
        // We found a negative cycle
        // Go n steps back to make sure we are in the cycle
        int v = cycle_start;
        for (size_t i = 0; i < n; ++i) {
            v = pred[(size_t)v];
        }

        // Add vertices to the cycle
        for (int u = v;; u = pred[(size_t)u]) {
            cycle.push_back(u);
            if (u == v && cycle.size() > 1) { // Break when we complete the cycle
                break;
            }
        }
        reverse(cycle.begin(), cycle.end()); // Reverse to get the correct order

        // Construct the cycle path string
        std::string cyclePath = "Cycle detected: ";
        for (int i = cycle.size() - 1; i >= 0; --i) {
            cyclePath += std::to_string(cycle[(size_t)i]);
            if (i > 0) {
                cyclePath += " -> ";
            }
        }
        return cyclePath; // Return the cycle path
    }

    return "No negative cycle found."; // Return if no negative cycle is found
}








