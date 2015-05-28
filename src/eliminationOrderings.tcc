#include <iostream>
#include <vector>

#include <agrum/graphs/undiGraph.h>
#include <agrum/BN/BayesNet.h>

// std::vector<NodeId>
template <typename GUM_SCALAR> void
findEliminationOrder(std::vector<gum::NodeId>& variables,
    gum::BayesNet<GUM_SCALAR> bn,
    int (*heuristicFunc) (gum::NodeId, gum::BayesNet<GUM_SCALAR>)) {

    gum::UndiGraph undirectedGraph = bn.moralGraph();

    // Eliminate from undirected graph each variable
    for (int i = 0; i < variables.size(); i++) {
        // Score each variable to be eliminated
        std::vector<int> scores;
        for (int i = 0; i < variables.size(); i++) {
            // std::cout << heuristicFunc(variables[i], bn) << std::endl;
            scores.push_back( heuristicFunc(variables[i], bn) );
        }
    }

}
