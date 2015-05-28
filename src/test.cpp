#include <iostream>
#include "eliminationOrderings.tcc"

#include <agrum/BN/BayesNet.h>
#include <agrum/BN/io/BIF/BIFReader.h>
#include <agrum/graphs/undiGraph.h>

template <typename GUM_SCALAR>
int test(gum::NodeId, gum::BayesNet<GUM_SCALAR>) {
    return 10;
}

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
    variables.push_back( bn.idFromName("dysp") );

    int (*funcPointer) (gum::NodeId, gum::BayesNet<double>);
    funcPointer = &test;

    findEliminationOrder(variables, bn, funcPointer);
}
