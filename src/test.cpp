#include <iostream>
#include "eliminationOrderings.tcc"
#include "heuristics.tcc"

#include <agrum/BN/BayesNet.h>
#include <agrum/graphs/undiGraph.h>
#include <agrum/BN/io/BIF/BIFReader.h>

int main(void) {

    std::cout << ">>> Loading BN..."  << std::endl;
    gum::BayesNet<double> bn;
    gum::BIFReader<double> bifReader(&bn, "../datatest/asia.bif");
    bifReader.proceed();
    std::cout << ">>> BN loaded with " << bn.size() << " nodes."  << std::endl;

    std::cout << ">>> Finding elimination ordering..."  << std::endl;

    std::vector<gum::NodeId> variables;
    variables.push_back( bn.idFromName("asia") );
    variables.push_back( bn.idFromName("tub") );
    variables.push_back( bn.idFromName("lung") );
    variables.push_back( bn.idFromName("either") );
    variables.push_back( bn.idFromName("dysp") );

    int (*funcPointer) (gum::NodeId&, gum::UndiGraph&, gum::BayesNet<double>&);

    std::cout << ">> Min Neighbours with " << bn.size() << " variables." << std::endl;
    funcPointer = &minNeighbours;
    for (auto varId : findEliminationOrder(variables, bn, funcPointer)) {
        std::cout << bn.variable(varId).name() << " ";
    }
    std::cout << std::endl;

    std::cout << ">> Min Weight with " << bn.size() << " variables." << std::endl;
    funcPointer = &minWeight;
    for (auto varId : findEliminationOrder(variables, bn, funcPointer)) {
        std::cout << bn.variable(varId).name() << " ";
    }
    std::cout << std::endl;

    std::cout << ">> Min Fill with " << bn.size() << " variables." << std::endl;
    funcPointer = &minFill;
    for (auto varId : findEliminationOrder(variables, bn, funcPointer)) {
        std::cout << bn.variable(varId).name() << " ";
    }
    std::cout << std::endl;

    std::cout << ">> WeightedMinFill with " << bn.size() << " variables." << std::endl;
    funcPointer = &weightedMinFill;
    for (auto varId : findEliminationOrder(variables, bn, funcPointer)) {
        std::cout << bn.variable(varId).name() << " ";
    }
    std::cout << std::endl;
}
