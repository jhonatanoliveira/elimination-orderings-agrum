#ifndef HEURISTICS
#define HEURISTICS

#include <agrum/BN/BayesNet.h>
#include <agrum/graphs/undiGraph.h>


template<typename GUM_SCALAR>
int minNeighbours(gum::NodeId& nodeId, gum::UndiGraph& undirectedGraph, gum::BayesNet<GUM_SCALAR>& bn) {
    return undirectedGraph.neighbours(nodeId).size();
}


template<typename GUM_SCALAR>
int minWeight(gum::NodeId& nodeId, gum::UndiGraph& undirectedGraph, gum::BayesNet<GUM_SCALAR>& bn) {
    int product = 1;
    for (auto neighbour : undirectedGraph.neighbours(nodeId)) {
        product *= bn.variable(neighbour).domainSize();
    }
    return product;
}


template<typename GUM_SCALAR>
int minFill(gum::NodeId& nodeId, gum::UndiGraph& undirectedGraph, gum::BayesNet<GUM_SCALAR>& bn) {
    int numEdges = 0;
    gum::NodeSet neighbours = undirectedGraph.neighbours(nodeId);
    for (auto neighbour1 : neighbours) {
        for (auto neighbour2 : neighbours) {
            if (neighbour1 != neighbour2) {
                if (!undirectedGraph.existsEdge(neighbour1, neighbour2) || !undirectedGraph.existsEdge(neighbour2, neighbour1)) {
                    numEdges++;
                }
            }
        }
    }
    // Divide by 2, since if an edge from A to B is needed
    // than an edge from B to A is also needed
    return numEdges/2;
}



template<typename GUM_SCALAR>
int weightedMinFill(gum::NodeId& nodeId, gum::UndiGraph& undirectedGraph, gum::BayesNet<GUM_SCALAR>& bn) {
    int sumWeights = 0;
    gum::NodeSet neighbours = undirectedGraph.neighbours(nodeId);
    for (auto neighbour1 : neighbours) {
        for (auto neighbour2 : neighbours) {
            if (neighbour1 != neighbour2) {
                if (!undirectedGraph.existsEdge(neighbour1, neighbour2) || !undirectedGraph.existsEdge(neighbour2, neighbour1)) {
                    sumWeights += bn.variable(neighbour1).domainSize() * bn.variable(neighbour2).domainSize();
                }
            }
        }
    }
    // Divide by 2, since if an edge from A to B is needed
    // than an edge from B to A is also needed
    return sumWeights/2;
}


#endif
