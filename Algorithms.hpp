//hila.shamir99@gmail.com  314906983
#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP
#include <string> 
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include "Graph.hpp"
namespace ariel {
    // Declare the Algorithms namespace or class
    class Algorithms {
    public:
    static bool isConnected(const Graph& graph);
    static std::string shortestPath(Graph& graph, int start, int end) ;
    static std::string isBipartite(const Graph& graph) ;
    static bool isContainsCycle(const Graph& graph);
    static std::string negativeCycle(const Graph& graph);

    private:
    static void dfs(const std::vector<std::vector<int> >&, std::vector<bool>&, int);
    static bool dss( Graph& g, int u, std::vector<int>& parent, std::unordered_set<int>& visited, std::vector<int>& cycle);
    static bool isCyclic(std::vector<bool>& visited,size_t v,  std::vector<bool>& recStack, std::vector<size_t >& parent, const Graph& g, std::vector<size_t>& cycle);  

    };
}

 

#endif /* ALGORITHMS_HPP */
