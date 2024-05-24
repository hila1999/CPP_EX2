// hila.shamir99@gmail.com  314906983
#include "Graph.hpp"
#include <iostream>
#include <stdexcept>

using namespace std;

namespace ariel {

    // Loads an adjacency matrix into the graph
    void Graph::loadGraph(const std::vector<std::vector<int>> &matrix) {
        if (matrix.empty()) {
            throw std::invalid_argument("The graph cannot be empty");
        }
        // Check if the matrix is square
        size_t size = matrix.size();
        for (const auto &row : matrix) {
            if (row.size() != size) {
                throw std::invalid_argument("Invalid graph: The graph is not a square matrix.");
            }
        }
        // If the matrix is square, load it into the adjacencyMatrix
        this->adjacencyMatrix = matrix;
    }

    // Prints the graph's adjacency matrix and basic information
    void Graph::printGraph() const {
        unsigned int numVertices = static_cast<unsigned int>(adjacencyMatrix.size());
        int numEdges = 0;
        // Count the number of edges in the graph
        for (size_t i = 0; i < adjacencyMatrix.size(); ++i) {
            for (size_t j = 0; j < adjacencyMatrix.size(); ++j) {
                if (adjacencyMatrix[i][j] != 0) {
                    numEdges += 1;
                }
            }
        }
        // If the graph is undirected, each edge is counted twice
        if (!isDirected()) {
            numEdges /= 2;
        }
        std::cout << "Graph with " << numVertices << " vertices and " << numEdges << " edges." << std::endl;
    }

    // Checks if the graph is directed
    bool Graph::isDirected() const {
        for (size_t i = 0; i < adjacencyMatrix.size(); ++i) {
            for (size_t j = 0; j < adjacencyMatrix[i].size(); ++j) {
                if (adjacencyMatrix[i][j] != adjacencyMatrix[j][i]) {
                    return true;
                }
            }
        }
        return false;
    }

    // Returns the neighbors of a given node
    std::vector<size_t> Graph::getNeighbors(size_t node) const {
        std::vector<size_t> neighbors;
        if (node >= adjacencyMatrix.size()) {
            throw std::out_of_range("Node index out of range");
        }
        for (size_t i = 0; i < adjacencyMatrix[node].size(); ++i) {
            if (adjacencyMatrix[node][i] != 0) {
                neighbors.push_back(i);
            }
        }
        return neighbors;
    }

    // Returns the number of vertices in the graph
    size_t Graph::size() const {
        return adjacencyMatrix.size();
    }

    // Returns the adjacency matrix of the graph
    std::vector<std::vector<int>> Graph::getAdjMatrix() {
        return this->adjacencyMatrix;
    }

    // Checks if two graphs have the same size
    void Graph::check_same_size(const Graph &other) const {
        if (adjacencyMatrix.size() != other.adjacencyMatrix.size()) {
            throw std::invalid_argument("Graphs must be of the same size");
        }
    }

    // Adds two graphs
    Graph Graph::operator+(const Graph &other) const {
        if (other.adjacencyMatrix.empty()) {
            throw std::invalid_argument("The graph cannot be empty");
        }
        check_same_size(other);
        size_t n = adjacencyMatrix.size();
        std::vector<std::vector<int>> result(n, std::vector<int>(n, 0));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                result[i][j] = adjacencyMatrix[i][j] + other.adjacencyMatrix[i][j];
            }
        }
        return Graph(result);
    }

    // Adds another graph to this graph
    Graph &Graph::operator+=(const Graph &other) {
        if (other.adjacencyMatrix.empty()) {
            throw std::invalid_argument("The graph cannot be empty");
        }
        check_same_size(other);
        size_t n = adjacencyMatrix.size();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                adjacencyMatrix[i][j] += other.adjacencyMatrix[i][j];
            }
        }
        return *this;
    }

    // Unary plus operator
    Graph Graph::operator+() const {
        return *this;
    }

    // Subtracts another graph from this graph
    Graph Graph::operator-(const Graph &other) const {
        check_same_size(other);
        size_t n = adjacencyMatrix.size();
        std::vector<std::vector<int>> result(n, std::vector<int>(n, 0));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                result[i][j] = adjacencyMatrix[i][j] - other.adjacencyMatrix[i][j];
            }
        }
        return Graph(result);
    }

    // Subtracts another graph from this graph and assigns the result to this graph
    Graph &Graph::operator-=(const Graph &other) {
        check_same_size(other);
        size_t n = adjacencyMatrix.size();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                adjacencyMatrix[i][j] -= other.adjacencyMatrix[i][j];
            }
        }
        return *this;
    }

    // Unary minus operator
    Graph Graph::operator-() const {
        size_t n = adjacencyMatrix.size();
        std::vector<std::vector<int>> result(n, std::vector<int>(n, 0));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                result[i][j] = -adjacencyMatrix[i][j];
            }
        }
        return Graph(result);
    }

    // Equality operator
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

    // Inequality operator
    bool Graph::operator!=(const Graph &other) const {
        return !(*this == other);
    }

    // Counts the number of edges in the graph
    size_t Graph::countEdges() const {
        size_t count = 0;
        for (const auto &row : this->adjacencyMatrix) {
            for (int val : row) {
                if (val != 0) {
                    count++;
                }
            }
        }
        return count / 2; // Assuming undirected graph, each edge is counted twice
    }

    // Less-than operator
    bool Graph::operator<(const Graph &other) const {
        size_t thisEdges = this->countEdges();
        size_t otherEdges = other.countEdges();
        if (thisEdges != otherEdges) {
            return thisEdges < otherEdges;
        }

        // If number of edges is the same, compare number of vertices
        if (this->size() != other.size()) {
            return this->size() < other.size();
        }

        // If number of edges and vertices are the same, perform element-wise comparison
        for (size_t i = 0; i < this->size(); ++i) {
            for (size_t j = 0; j < this->size(); ++j) {
                if (this->adjacencyMatrix[i][j] != other.adjacencyMatrix[i][j]) {
                    return this->adjacencyMatrix[i][j] < other.adjacencyMatrix[i][j];
                }
            }
        }
        return false;
    }

    // Less-than-or-equal operator
    bool Graph::operator<=(const Graph &other) const {
        return *this < other || *this == other;
    }

    // Greater-than operator
    bool Graph::operator>(const Graph &other) const {
        if (other.adjacencyMatrix.empty()) {
            throw std::invalid_argument("The graph to compare with does not exist.");
        }
        return !(*this <= other);
    }

    // Greater-than-or-equal operator
    bool Graph::operator>=(const Graph &other) const {
        return !(*this < other);
    }

    // Prefix increment operator
    Graph &Graph::operator++() {
        size_t n = adjacencyMatrix.size();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                adjacencyMatrix[i][j]++;
            }
        }
        return *this;
    }

    // Postfix increment operator
    Graph Graph::operator++(int) {
        Graph temp = *this;
        ++(*this);
        return temp;
    }

        // Prefix decrement operator
    Graph &Graph::operator--() {
        size_t n = adjacencyMatrix.size();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                adjacencyMatrix[i][j]--;
            }
        }
        return *this;
    }

    // Postfix decrement operator
    Graph Graph::operator--(int) {
        if (adjacencyMatrix.empty()) throw std::invalid_argument("The graph cannot be empty");
        Graph temp = *this;
        --(*this);
        return temp;
    }

    // Divides all elements of the adjacency matrix by a scalar
    Graph &Graph::operator/=(int scalar) {
        if (scalar != 0) {
            for (size_t i = 0; i < adjacencyMatrix.size(); ++i) {
                for (size_t j = 0; j < adjacencyMatrix[i].size(); ++j) {
                    adjacencyMatrix[i][j] /= scalar;
                }
            }
        }
        return *this;
    }

    // Multiplies all elements of the adjacency matrix by a scalar
    Graph &Graph::operator*=(int scalar) {
        for (size_t i = 0; i < adjacencyMatrix.size(); ++i) {
            for (size_t j = 0; j < adjacencyMatrix[i].size(); ++j) {
                adjacencyMatrix[i][j] *= scalar;
            }
        }
        return *this;
    }

    // Multiplies two graphs
    Graph Graph::operator*(const Graph &other) const {
        check_same_size(other);
        size_t n = adjacencyMatrix.size();
        std::vector<std::vector<int>> result(n, std::vector<int>(n, 0));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    for (size_t k = 0; k < n; ++k) {
                        result[i][j] += adjacencyMatrix[i][k] * other.adjacencyMatrix[k][j];
                    }
                }
            }
        }
        return Graph(result);
    }

    // Stream insertion operator for printing the graph
    ostream &operator<<(ostream &os, const Graph &graph) {
        if (graph.getAdjacencyMatrix().empty()) throw std::invalid_argument("The graph cannot be empty");
        size_t n = graph.getAdjacencyMatrix().size();
        for (size_t i = 0; i < n; ++i) {
            os << "[";
            for (size_t j = 0; j < n; ++j) {
                os << graph.getA(i, j);
                if (j < n - 1) {
                    os << ", ";
                }
            }
            os << "]";
            os << "\n";
        }
        return os;
    }

}

 
