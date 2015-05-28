#ifndef EOA
#define EOA

#include <iostream>
#include <vector>

#include <agrum/graphs/undiGraph.h>
#include <agrum/BN/BayesNet.h>

void removeNode(const gum::NodeId& nodeId, gum::UndiGraph& undirectedGraph) {
    gum::NodeSet neighbours = undirectedGraph.neighbours(nodeId);
    for (auto neighbour1 : neighbours) {
        for (auto neighbour2 : neighbours) {
            if (neighbour1 != neighbour2) {
                undirectedGraph.addEdge(neighbour1, neighbour2);
            }
        }
    }
    undirectedGraph.eraseNode(nodeId);
}

template <typename GUM_SCALAR>
std::vector<gum::NodeId> findEliminationOrder(std::vector<gum::NodeId> variables,
    gum::BayesNet<GUM_SCALAR>& bn,
    int (*heuristicFunc) (gum::NodeId&, gum::UndiGraph&, gum::BayesNet<GUM_SCALAR>&)) {

    std::vector<gum::NodeId> order;
    gum::UndiGraph undirectedGraph = bn.moralGraph();

    // Eliminate from undirected graph each variable
    int originalSize = variables.size();
    for (int i = 0; i < originalSize; i++) {

        // Score each variable to be eliminated
        int currentMinimumScore;
        bool firstRun = true;
        gum::NodeId chosenVar;
        for (int i = 0; i < variables.size(); i++) {
            int score = heuristicFunc(variables[i], undirectedGraph, bn);
            if (firstRun) {
                currentMinimumScore = score;
                chosenVar = variables[i];
                firstRun = false;
            } else {
                if (score < currentMinimumScore) {
                    currentMinimumScore = score;
                    chosenVar = variables[i];
                }
            }
            // std::cout << ">> Score for var " << bn.variable(variables[i]).name() << " is " << heuristicFunc(variables[i], undirectedGraph, bn) << std::endl;
        }

        // Remove less scored variable and add it to the return list
        // std::cout << ">> Removing " << bn.variable(chosenVar).name() << std::endl;
        // std::cout << ">> UndiGraph before " << variables << std::endl;
        order.push_back(chosenVar);
        removeNode(chosenVar, undirectedGraph);
        variables.erase(std::remove(variables.begin(), variables.end(), chosenVar), variables.end());
        // std::cout << ">> UndiGraph after " << variables << std::endl;
        // std::cout << ">> Nodes quantity " << undirectedGraph.size() << std::endl;
    }

    return order;
}

#endif
