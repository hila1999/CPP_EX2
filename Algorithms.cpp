//hila.shamir99@gmail.com  314906983
#include "Graph.hpp"
#include "Algorithms.hpp"
#include <vector>
#include <stack>
#include <queue>
#include <unordered_set>
#include <algorithm>
#include <iostream>

#include <climits> // for INT_MAX
#include <cstddef> // for SIZE_MAX
using namespace std;
using namespace ariel;
    void Algorithms::dfs(const vector<std::vector<int>>& adjacencyMatrix, vector<bool>& visited, int vertex) {
        visited[static_cast<std::vector<bool>::size_type>(vertex)] = true; // Mark the current vertex as visited
        for (size_t i = 0; i < adjacencyMatrix[static_cast<std::vector<std::vector<int>>::size_type>(vertex)].size(); ++i) {
            if (adjacencyMatrix[static_cast<vector<std::vector<int>>::size_type>(vertex)][i]!=0 && !visited[i]) {
                dfs(adjacencyMatrix, visited, static_cast<int>(i)); // Recursively visit adjacent vertices
            }
        }
    }
    
bool Algorithms::isConnected(const Graph& graph) {
    const std::vector<std::vector<int>>& adjacencyMatrix = graph.getAdjacencyMatrix();
    if (adjacencyMatrix.empty()) {
        return false; // An empty graph is not connected
    }
    std::vector<bool> visited(adjacencyMatrix.size(), false);
    dfs(adjacencyMatrix, visited, 0);

    // Check if all vertices are visited using std::all_of()
    if (std::all_of(visited.begin(), visited.end(), [](bool isVisited) { return isVisited; })) {
        return true; // If all vertices are visited, the graph is connected
    } 
        return false; // If any vertex is not visited, the graph is not connected
}
string Algorithms::shortestPath(Graph& graph, int start, int end) {
    if(graph.size()==0){
        return "-1";
    }
    size_t numOfVer = graph.size();
    vector<int> dist(numOfVer, INT_MAX);
     vector<int> predecessor(numOfVer, -1);
    dist[(size_t)start] = 0;
     vector<vector<int>> adjacencyMatrix = graph.getAdjMatrix();
    for (size_t uar = 0; uar < numOfVer-1; ++uar) {
     
       for (size_t jar = 0; jar < numOfVer; jar++)
        {
            for (size_t kar = 0; kar < numOfVer; kar++)
            {
                if (adjacencyMatrix[jar][kar]!=0 && dist[jar] != INT_MAX && dist[jar] + adjacencyMatrix[jar][kar] < dist[kar])
                {
                    dist[kar] = dist[jar] + adjacencyMatrix[jar][kar];
                    predecessor[kar] = static_cast<int>(jar);
                }
            }
        }
    }
   
      for (size_t jar = 0; jar < numOfVer; jar++)
    {
        for (size_t k = 0; k < numOfVer; k++)
        {
            if (adjacencyMatrix[jar][k]!=0 && dist[jar] != INT_MAX && dist[jar] + adjacencyMatrix[jar][k] < dist[k])
            {
                return "Negative cycle detected";
            }
        }
    }
       // If there is no path from start to end
    if (dist[(size_t)end] == INT_MAX)
    {
        return "-1";
    }
                
     string path = to_string(end);
    for (size_t var = (size_t)end; var != start; var = (size_t)predecessor[var])
    {
        path = to_string((int)predecessor[var]) + "->" + path;
    }

    return path;
            }




std::string Algorithms::isBipartite(const Graph &g) {
    size_t numNodes = g.size();
    vector<vector<int>> adjMatrix = g.getAdjacencyMatrix(); // supposed to point on the same matrix

      if (g.isDirected())
    {
        for (size_t i = 0; i < numNodes; ++i)
        {
            for (size_t jar = 0; jar < numNodes; ++jar)
            {
                if (adjMatrix[i][jar] != 0 && adjMatrix[jar][i] == 0)
                {
                    adjMatrix[jar][i] = adjMatrix[i][jar];
                }
                else if (adjMatrix[i][jar] == 0 && adjMatrix[jar][i] != 0)
                {
                    adjMatrix[i][jar] = adjMatrix[jar][i];
                }
            }
        }
        
    }
    vector<int> colorArr(numNodes, -1);
    std::vector<size_t> setA, setB;

    for (size_t i = 0; i < numNodes; i++) {
        if (colorArr[i] == -1) {
            std::queue<size_t> q;
            q.push(i);
            colorArr[i] = 1;
            while (!q.empty()) {
                size_t u = q.front();
                q.pop();
                std::vector<size_t> neighbors = g.getNeighbors(u);
                for (size_t v : neighbors) {
                    if (colorArr[v] == -1) {
                        colorArr[v] = 1 - colorArr[u];
                        q.push(v);
                    } else if (colorArr[v] == colorArr[u]) {
                        return "0";
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < numNodes; i++) {
        if (colorArr[i] == 1) {
            setA.push_back(i);
        } else {
            setB.push_back(i);
        }
    }

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

bool Algorithms::isCyclic(vector<bool>& visited,size_t v,  std::vector<bool>& stack, std::vector<size_t >& parent, const Graph& g, std::vector<size_t >& cyclePath) {
    visited[v] = true;
    stack[v] = true;
    bool isDirected = g.isDirected();
    size_t numNodes = g.size();
    for (size_t i = 0; i < numNodes; i++) {
        if (g.haveEdge(v, i)) {  // Check if there is an edge from v to i
            if (!visited[i]) {
                parent[i] = v;
                if (isCyclic(visited,i, stack, parent, g, cyclePath)) {
                    return true;
                }
            } else if ((isDirected && stack[i]) || (!isDirected && stack[i] && parent[v] != i)) {
                // If a cycle is detected, trace back to print the cycle
                cyclePath.push_back(i);
                for (size_t p = v; p != i; p = parent[p]) {
                    cyclePath.push_back(p);
                }
                cyclePath.push_back(i);  // Complete the cycle by adding the start node again
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

    stack[v] = false;
    return false;
}
bool Algorithms::isContainsCycle(const Graph& g) {
    size_t numNodes = g.size();
    vector<size_t > parent(numNodes, SIZE_MAX);
    vector<size_t> cyclePath;  // To store the path of the cycle
    vector<bool> visited(numNodes, false);
    vector<bool> stack(numNodes, false);
    

    for (size_t i = 0; i < numNodes; i++) {
        if (!visited[i]) {
            if (isCyclic(visited,i,  stack, parent, g, cyclePath)) {
                return true;
            }
        }
    }
    return false;
}
std::string Algorithms::negativeCycle(const Graph& g) {
    size_t n = g.size();
    std::vector<int> distance(n, INT_MAX);
    std::vector<int> pred(n, -1);
    int cycle_start = -1;

    // Assume the source vertex is 0, for checking the negative cycle we'll run bellman-ford from the first vertex
    distance[0] = 0;

    for (size_t i = 0; i < n; ++i) {
        cycle_start = -1;
        for (size_t u = 0; u < n; ++u) {
            for (size_t v = 0; v < n; ++v) {
                if (g.getAdjacencyMatrix()[u][v] != 0) {
                    // relax the edge
                    int new_distance = distance[u] + g.getAdjacencyMatrix()[u][v];
                    if (new_distance < distance[v]) {
                        distance[v] = new_distance;
                        pred[v] = u;
                        cycle_start = v;
                    }
                }
            }
        }
    }

       vector<int> cycle;
    if (cycle_start != -1)
    {
        // We found a negative cycle
        // Go n steps back to make sure we are in the cycle
        int v = cycle_start;
        for (size_t i = 0; i < n; ++i)
        {
            v = pred[(size_t)v];
        }

        // Add vertices to the cycle
        for (int u = v;; u = pred[(size_t)u])
        {
            cycle.push_back(u);
            if (u == v && cycle.size() > 1)
            {
                break;
            }
        }
        reverse(cycle.begin(), cycle.end());
    

        // Print the cycle path
        std::string cyclePath = "Cycle detected: ";
        for (int i = cycle.size() - 1; i >= 0; --i)
{
            cyclePath += std::to_string(cycle[(size_t)i]);
            if (i > 0) {
                cyclePath += " -> ";
            }
        }
        return cyclePath;
    }

    return "No negative cycle found.";
}








