#include <iostream>
#include <vector>
#include <set>
#include <string>

void removeNode(gum::NodeId& nodeId, std::vector<std::set<gum::NodeId>>& potentials);
template <typename GUM_SCALAR>
int populationEnergyScore(gum::NodeId nodeId, std::vector<std::set<gum::NodeId>> potentials, gum::BayesNet<GUM_SCALAR>& bn);

template <typename GUM_SCALAR>
std::vector<gum::NodeId> populationEnergyOrder(
    std::vector<gum::NodeId> variables,
    gum::BayesNet<GUM_SCALAR>& bn) {

    std::vector<gum::NodeId> order;

    // Build a vector with sets of potential's header
    std::vector<std::set<gum::NodeId>> potentials;
    for (auto node : bn.nodes()) {
        std::set<gum::NodeId> potentialVars;
        // std::cout << ">>> Vars for potential:" << std::endl;
        for(auto var : bn.cpt(node).variablesSequence()) {
            potentialVars.insert( bn.idFromName(var->name()) );
            // std::cout << var->name() << " ";
        }
        potentials.push_back(potentialVars);
        // std::cout << std::endl;
    }
    int originalSize = variables.size();
    // Score and remove one variable per time
    for (int i = 0; i < originalSize; i++) {
        int currentMinimumScore;
        bool firstRun = true;
        gum::NodeId chosenVar;
        for (int j = 0; j < variables.size(); j++) {
            // Perform Population energy to score the variable.
            // std::cout << ">> Pop energy..." << std::endl;
            int score = populationEnergyScore(variables[j], potentials, bn);
            // std::cout << ">> finished Pop energy" << std::endl;
            // std::cout << ">>> Score for var " << bn.variable(variables[j]).name() << ": " << score << std::endl;
            // Decide if score is minimum
            if (firstRun) {
                currentMinimumScore = score;
                chosenVar = variables[j];
                firstRun = false;
            } else {
                if (score < currentMinimumScore) {
                    currentMinimumScore = score;
                    chosenVar = variables[j];
                }
            }
        }
        // std::cout << ">> i: " << i << std::endl;
        // std::cout << ">> chosenVar: " << chosenVar << std::endl;
        // std::cout << ">> variables: " << variables.size() << std::endl;
        // remove var
        order.push_back(chosenVar);
        // std::cout << ">> Push back" << std::endl;
        removeNode(chosenVar, potentials);
        // std::cout << ">> Remove Node" << std::endl;
        variables.erase(std::remove(variables.begin(), variables.end(), chosenVar), variables.end());
        // std::cout << ">> Erase node. " << variables.size() << std::endl;
    }

    return order;
}


std::set<gum::NodeId> potentialsUnion(std::set<gum::NodeId>& potential1, std::set<gum::NodeId>& potential2) {
    std::set<gum::NodeId> product;
    for (auto var : potential1) {
        product.insert(var);
    }
    for (auto var : potential2) {
        product.insert(var);
    }
    return product;
}


void potentialsCombination(std::set<gum::NodeId>& potential1, std::set<gum::NodeId>& potential2, std::vector<std::set<gum::NodeId>>& potentials) {
    std::set<gum::NodeId> productPotential = potentialsUnion(potential1, potential2);
    potentials.push_back(productPotential);
    potentials.erase(std::remove(potentials.begin(), potentials.end(), potential1), potentials.end());
    potentials.erase(std::remove(potentials.begin(), potentials.end(), potential2), potentials.end());
}


template <typename GUM_SCALAR>
int numOfMult(std::set<gum::NodeId> potential1, std::set<gum::NodeId> potential2, gum::BayesNet<GUM_SCALAR>& bn) {
    std::set<gum::NodeId> product = potentialsUnion(potential1, potential2);
    int numOfMult = 1;
    for (auto var : product) {
        numOfMult *= bn.variable(var).domainSize();
    }
    return numOfMult;
}


std::vector<std::set<gum::NodeId>> findRelevantPotentials(gum::NodeId& nodeId, std::vector<std::set<gum::NodeId>>& potentials) {
    std::vector<std::set<gum::NodeId>> relevantPotentials;
    // Collect relevant potentials (ones involving given node id)
    for (auto potential : potentials) {
        if (potential.find(nodeId) != potential.end()) {
            relevantPotentials.push_back(potential);
            // ### DEBUG
            // for (auto var : potential) {
            //     std::cout << var << " ";
            // }
            // std::cout << std::endl;
            // --- DEBUG
        }
    }
    return relevantPotentials;
}


template <typename GUM_SCALAR>
int populationEnergyScore(gum::NodeId nodeId, std::vector<std::set<gum::NodeId>> potentials, gum::BayesNet<GUM_SCALAR>& bn) {
    std::vector<std::set<gum::NodeId>> relevantPotentials = findRelevantPotentials(nodeId, potentials);
    int sumOfProducts = 0;
    // std::cout << "--> Relevant potentials: " << relevantPotentials.size() << std::endl;
    while (relevantPotentials.size() > 1) {
        std::set<gum::NodeId> firstPotential;
        std::set<gum::NodeId> secondPotential;
        int product;
        int tempProduct;
        bool firstRun = true;
        // Discover which is better order for combination
        for (auto potential1 : relevantPotentials) {
            for (auto potential2 : relevantPotentials) {
                if (potential1 != potential2) {
                    if (firstRun) {
                        firstRun = false;
                        // std::cout << "--> Finished 1st run." << std::endl;
                        firstPotential = potential1;
                        secondPotential = potential2;
                        product = numOfMult(potential1, potential2, bn);
                        // ### DEBUG
                        // std::cout << ">>> VARs 1:" << std::endl;
                        // for (auto var : potential1) {
                        //     std::cout << var << " ";
                        // }
                        // std::cout << ">>> VARs 2:" << std::endl;
                        // for (auto var : potential2) {
                        //     std::cout << var << " ";
                        // }
                        // std::cout << std::endl;
                        // --- DEBUG
                        // std::cout << product << std::endl;
                    } else {
                        // std::cout << "--> 2nd run." << std::endl;
                        tempProduct = numOfMult(potential1, potential2, bn);
                        if (tempProduct < product) {
                            firstPotential = potential1;
                            secondPotential = potential2;
                            product = tempProduct;
                        }
                    }
                }
            }
        }
        // In case of a tie with only 2 potentials
        if (firstPotential.size() == 0 && secondPotential.size() == 0) {
            firstPotential = relevantPotentials[0];
            secondPotential = relevantPotentials[1];
        }
        // std::cout << "--> relevantPotentials: " << relevantPotentials.size() << std::endl;
        // std::cout << "--> tempProduct: " << tempProduct << std::endl;
        // std::cout << "--> product: " << product << std::endl;
        // Do the combination, remove potentials involved from relevant potential list
        // and add the combination to the relevant potential list
        sumOfProducts += product;
        potentialsCombination(firstPotential, secondPotential, relevantPotentials);
    }
    return sumOfProducts;
}


void removeNode(gum::NodeId& nodeId, std::vector<std::set<gum::NodeId>>& potentials) {
    std::vector<std::set<gum::NodeId>> relevantPotentials = findRelevantPotentials(nodeId, potentials);
    while (relevantPotentials.size() > 1) {
        potentialsCombination(relevantPotentials[0], relevantPotentials[1], relevantPotentials);
    }
    potentials.push_back(relevantPotentials[0]);
}
