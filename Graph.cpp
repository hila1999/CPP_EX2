// hila.shamir99@gmail.com  314906983
#include "Graph.hpp"
#include <iostream>

#include <stdexcept>
using namespace std;
// using namespace ariel;
namespace ariel
{
    void Graph::loadGraph(const std::vector<std::vector<int>> &matrix)
    {
        if (matrix.empty())
        {
            throw std::invalid_argument("The graph cannot be empty");
        }
        // Check if the matrix is square
        size_t size = matrix.size();
        for (const auto &row : matrix)
        {
            if (row.size() != size)
            {
                throw std::invalid_argument("Invalid graph: The graph is not a square matrix.");
            }
        }

        // If the matrix is square, load it into the adjacencyMatrix
        this->adjacencyMatrix = matrix;
    }
    void Graph::printGraph() const
    {
        unsigned int numVertices = static_cast<unsigned int>(adjacencyMatrix.size());
        int numEdges = 0;
        for (size_t i = 0; i < adjacencyMatrix.size(); ++i)
        {
            for (size_t j = 0; j < adjacencyMatrix.size(); ++j)
            {
                if (adjacencyMatrix[i][j] != 0)
                {
                    numEdges += 1;
                }
            }
        }

        if (!isDirected())
        {
            numEdges /= 2;
        }
        std::cout << "Graph with " << numVertices << " vertices and " << numEdges << " edges." << std::endl;
    }
    bool Graph::isDirected() const
    {
        for (size_t i = 0; i < adjacencyMatrix.size(); ++i)
        {
            for (size_t j = 0; j < adjacencyMatrix[i].size(); ++j)
            {
                if (adjacencyMatrix[i][j] != adjacencyMatrix[j][i])
                {
                    return true;
                }
            }
        }
        return false;
    }
    std::vector<size_t> Graph::getNeighbors(size_t node) const
    {
        std::vector<size_t> neighbors;
        if (node >= adjacencyMatrix.size())
        {
            throw std::out_of_range("Node index out of range");
        }
        for (size_t i = 0; i < adjacencyMatrix[node].size(); ++i)
        {
            if (adjacencyMatrix[node][i] != 0)
            {
                neighbors.push_back(i);
            }
        }
        return neighbors;
    }
    size_t Graph::size() const
    {
        return adjacencyMatrix.size();
    }
    std::vector<std::vector<int>> Graph::getAdjMatrix()
    {
        return this->adjacencyMatrix;
    }

    void Graph::check_same_size(const Graph &other) const
    {
        if (adjacencyMatrix.size() != other.adjacencyMatrix.size())
        {
            throw std::invalid_argument("Graphs must be of the same size");
        }
    }

    Graph Graph::operator+(const Graph &other) const
    {
        if(other.adjacencyMatrix.empty())
        {
            throw std::invalid_argument("The graph cannot be empty");
        }
        check_same_size(other);
        size_t n = adjacencyMatrix.size();
        std::vector<std::vector<int>> result(n, std::vector<int>(n, 0));
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                result[i][j] = adjacencyMatrix[i][j] + other.adjacencyMatrix[i][j];
            }
        }
        return Graph(result);
    }

    Graph &Graph::operator+=(const Graph &other)
    {
        if(other.adjacencyMatrix.empty())
        {
            throw std::invalid_argument("The graph cannot be empty");
        }
        check_same_size(other);
        size_t n = adjacencyMatrix.size();
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                adjacencyMatrix[i][j] += other.adjacencyMatrix[i][j];
            }
        }
        return *this;
    }

    Graph Graph::operator+() const
    {
        return *this;
    }

    Graph Graph::operator-(const Graph &other) const
    {
        check_same_size(other);
        size_t n = adjacencyMatrix.size();
        std::vector<std::vector<int>> result(n, std::vector<int>(n, 0));
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                result[i][j] = adjacencyMatrix[i][j] - other.adjacencyMatrix[i][j];
            }
        }
        return Graph(result);
    }

    Graph &Graph::operator-=(const Graph &other)
    {
        check_same_size(other);
        size_t n = adjacencyMatrix.size();
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                adjacencyMatrix[i][j] -= other.adjacencyMatrix[i][j];
            }
        }
        return *this;
    }

    Graph Graph::operator-() const
    {
        size_t n = adjacencyMatrix.size();
        std::vector<std::vector<int>> result(n, std::vector<int>(n, 0));
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                result[i][j] = -adjacencyMatrix[i][j];
            }
        }
        return Graph(result);
    }

  bool Graph::operator==(const Graph &other) const {
    // Check if the adjacency matrices are equal
    if (adjacencyMatrix == other.adjacencyMatrix) {
        return true;
    }

    // Check if neither graph is less than the other
    if (!(adjacencyMatrix < other.adjacencyMatrix) && !(other.adjacencyMatrix < adjacencyMatrix)) {
        return true;
    }

    return false;
}


    bool Graph::operator!=(const Graph &other) const
    {
         return !(*this == other);
    }
    // bool Graph::isContained(const Graph &other) const
    // {
    //    bool isContained = true;
    //     for (size_t i = 0; i < adjacencyMatrix.size(); ++i) {
    //         for (size_t j = 0; j < adjacencyMatrix[i].size(); ++j) {
    //             if (adjacencyMatrix[i][j] != 0 && other.adjacencyMatrix[i][j] != 1) {
    //                 isContained = false;
    //                 break;
    //             }
    //         }
    //         if (!isContained) break;
    //     }

    //     return isContained;
    // }
    // bool Graph::edgesNumber(const Graph &other) const
    // {
    //     int thisEdges = 0;
    //     int otherEdges = 0;
    //     for (const auto& row : adjacencyMatrix) {
    //         for (int val : row) {
    //             if (val != 0) thisEdges++;
    //         }
    //     }
    //     for (const auto& row : other.adjacencyMatrix) {
    //         for (int val : row) {
    //             if (val != 0) otherEdges++;
    //         }
    //     }

    //     if (thisEdges != otherEdges) return thisEdges < otherEdges;

    //     // If number of edges is the same, compare the order of the adjacency matrices
    //     return adjacencyMatrix.size() < other.adjacencyMatrix.size();
    // }
//      bool isContained(const Graph &g1, const Graph &g2)
//     {
//         // Checking if g1 is contained in g2, we'll go over g1 and check if all the values are in g2
//         for (size_t i = 0; i < g1.size(); i++)
//         {
//             for (size_t j = 0; j < g1.size(); j++)
//             {
//                 if (g1.adjacencyMatrix[i][j] != 0 && g2.adjacencyMatrix[i][j] == 0)
//                 {
//                     return false;
//                 }
//             }
//         }
//         return true;
//     }
//     bool Graph::operator<(const Graph &other) const
//          { 
//             bool result ;
            
//              if (this->adjacencyMatrix.size() > other.adjacencyMatrix.size())
//         {
//             return false;
//         }
//          if (this->adjacencyMatrix.size() < other.adjacencyMatrix.size()){
//        bool result = isContained(*this, other);
//     }
//      if (!result){
//     // If not contained, compare the number of edges
//     int thisEdges = 0;
//     int otherEdges = 0;
//     for (const auto& row : adjacencyMatrix) {
//         for (int val : row) {
//             if (val != 0) thisEdges++;
//         }
//     }
//     for (const auto& row : other.adjacencyMatrix) {
//         for (int val : row) {
//             if (val != 0) otherEdges++;
//         }
//     }

//     if (thisEdges < otherEdges) return true;
//     }

//     // If number of edges is the same, compare the order of the adjacency matrices
//     return false;
// }
size_t Graph::countEdges() const
    {
        size_t count = 0;
        for (const auto &row : this->adjacencyMatrix)
        {
            for (int val : row)
            {
                if (val != 0)
                {
                    count++;
                }
            }
        }
        return count / 2; // Assuming undirected graph, each edge is counted twice
    }
  bool Graph::operator<(const Graph &other) const
    {
        // Compare number of edges
        size_t thisEdges = this->countEdges();
        size_t otherEdges = other.countEdges();
        if (thisEdges != otherEdges)
        {
            return thisEdges < otherEdges;
        }

        // If number of edges is the same, compare number of vertices
        if (this->size() != other.size())
        {
            return this->size() < other.size();
        }

        // If number of edges and vertices are the same, perform element-wise comparison
        for (size_t i = 0; i < this->size(); ++i)
        {
            for (size_t j = 0; j < this->size(); ++j)
            {
                if (this->adjacencyMatrix[i][j] != other.adjacencyMatrix[i][j])
                {
                    return this->adjacencyMatrix[i][j] < other.adjacencyMatrix[i][j];
                }
            }
        }

        // If all else is equal, the graphs are not less than each other
        return false;
    }

        

    bool Graph::operator<=(const Graph &other) const
    {
        return *this < other || *this == other ;
    }

    bool Graph::operator>(const Graph &other) const //the upside of the < operator
    {
         if (other.adjacencyMatrix.empty()) {
        throw std::invalid_argument("The graph to compare with does not exist.");
    }
        return !(*this <= other);
    }

    bool Graph::operator>=(const Graph &other) const
    {
        return !(*this < other);
    }

    Graph &Graph::operator++()
    {
        size_t n = adjacencyMatrix.size();
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    adjacencyMatrix[i][j]++;
                }
            }
            return *this;
    }

    Graph Graph::operator++(int)
    {
          Graph temp = *this;
            ++(*this);
            return temp;
    }

    Graph &Graph::operator--()
    {
       size_t n = adjacencyMatrix.size();
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    adjacencyMatrix[i][j]--;
                }
            }
            return *this;
    }

    Graph Graph::operator--(int)
    {
        if(adjacencyMatrix.empty()) throw std::invalid_argument("The graph cannot be empty");
         Graph temp = *this;
            --(*this);
            return temp;
    }
      Graph& Graph::operator/=(int scalar) {
        if (scalar != 0) {
            for (size_t i = 0; i < adjacencyMatrix.size(); ++i) {
                for (size_t j = 0; j < adjacencyMatrix[i].size(); ++j) {
                    adjacencyMatrix[i][j] /= scalar;
                }
            }
        }
        return *this;
    }
        // Multiply a graph by a scalar
   Graph& Graph::operator*=(int scalar) {
        for (size_t i = 0; i < adjacencyMatrix.size(); ++i) {
            for (size_t j = 0; j < adjacencyMatrix[i].size(); ++j) {
                adjacencyMatrix[i][j] *= scalar;
            }
        }
        return *this;
    }

    Graph Graph::operator*(const Graph &other) const
    {
           check_same_size(other);
        size_t n = adjacencyMatrix.size();
        std::vector<std::vector<int>> result(n, std::vector<int>(n, 0));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if(i!=j){
                for (size_t k = 0; k < n; ++k) {
                    result[i][j] += adjacencyMatrix[i][k] * other.adjacencyMatrix[k][j];}
                }
            }
        }
        return Graph(result);
    }

    ostream &operator<<(ostream &os, const Graph &graph)
    {
        if (graph.getAdjacencyMatrix().empty())throw std::invalid_argument("The graph cannot be empty");
        size_t n = graph.getAdjacencyMatrix().size();
        for (size_t i = 0; i < n; ++i)
        {
            os << "[";
            for (size_t j = 0; j < n; ++j)
            {
                os << graph.getA(i, j);
                if (j < n - 1)
                {
                    os << ", ";
                }
            }
            os << "]";
            os << "\n";
            // if (i < n - 1)
            // {
            //     os << ", ";
            // }
        }
        return os;
    }

}