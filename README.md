# Graph Operations Library
hila.shamir99@עmail.com 314906983
## Overview
This project extends the functionality of the Graph class implemented in `Graph.cpp`, allowing representation and manipulation of graphs using adjacency matrices. In this task, arithmetic and comparison operators are added to support operations on graphs.

## Files Included
- `algorithms.cpp`
- `graph.cpp`
- `graph.hpp`
- `README.md`

## Functionality Added

### Arithmetic Operators
1. **Addition (+) and Compound Addition (+=):** Performs element-wise addition of two graphs. Both graphs must be of the same size (n x n).
    ```cpp
    Graph result = graph1 + graph2;
    graph1 += graph2;
    ```
2. **Unary Plus (+):** Returns a copy of the original graph.
    ```cpp
    Graph result = +graph;
    ```

### Comparison Operators
1. **Greater Than (>) and Greater Than or Equal To (>=):** Compares two graphs based on their size and edges.
    ```cpp
    if (graph1 > graph2) { /* do something */ }
    if (graph1 >= graph2) { /* do something */ }
    ```
2. **Less Than (<) and Less Than or Equal To (<=):** Compares two graphs based on their size and edges.
    ```cpp
    if (graph1 < graph2) { /* do something */ }
    if (graph1 <= graph2) { /* do something */ }
    ```
3. **Equal To (==) and Not Equal To (!=):** Compares two graphs based on their size and edges.
    ```cpp
    if (graph1 == graph2) { /* do something */ }
    if (graph1 != graph2) { /* do something */ }
    ```

### Other Operators
1. **Increment (++) and Decrement (--):** Increments or decrements the weights of all edges in the graph by 1.
    ```cpp
    ++graph;
    --graph;
    ```
2. **Scalar Multiplication (*=):** Multiplies the weight of all edges in the graph by a scalar integer.
    ```cpp
    graph *= 2;
    ```
3. **Graph Multiplication (*):** Multiplies two graphs by performing matrix multiplication.
    ```cpp
    Graph result = graph1 * graph2;
    ```

### Output Operator
1. **Output Operator (<<):** Prints a meaningful representation of the graph.
    ```cpp
    cout << graph;
    ```

## Implementation Details
The arithmetic and comparison operators are implemented based on the rules defined for matrices in linear algebra. Graph multiplication is performed using matrix multiplication. Comparison operators consider the size and edges of the graphs.

## How to Use
1. Include `graph.hpp` and `algorithms.cpp` in your project.
2. Create instances of the `Graph` class.
3. Perform desired operations using the implemented operators.
4. Print or use the graphs as required.
## Change in Algorithms.cpp
There is no difference in the algorithms.cpp, because the operators work with the functions I defined.
