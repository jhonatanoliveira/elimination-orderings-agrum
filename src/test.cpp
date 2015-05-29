#include <iostream>
#include <fstream>
#include <ctime>
#include <stdlib.h>
#include <string>

#include "eliminationOrderings.tcc"
#include "heuristics.tcc"
#include "populationEnergy.tcc"

#include <agrum/BN/BayesNet.h>
#include <agrum/graphs/undiGraph.h>
#include <agrum/BN/inference/variableElimination.h>
#include <agrum/BN/io/BIF/BIFReader.h>


gum::NodeId chooseVariable(std::vector<gum::NodeId>& allVariables, std::set<gum::NodeId>& chosenVariables);

int main(void) {

    // ######### CONFIGURATIONS
    const std::string databaseNames[] = {"asia", "alarm", "barley", "insurance", "mildew", "hepar2", "diabetes", "pathfinder", "munin"};
    int runQuantities = 1000;

    for (auto databaseName : databaseNames) {

        std::ofstream resultsFile;
        resultsFile.open (databaseName + ".csv");

        for (int runNumber = 0; runNumber < runQuantities; runNumber++) {

            std::cout << ">>> Loading BN " << databaseName << std::endl;
            gum::BayesNet<double> bn;
            gum::BIFReader<double> bifReader(&bn, "../datatest/" + databaseName + ".bif");
            bifReader.proceed();
            std::cout << ">>> BN loaded with " << bn.size() << " nodes."  << std::endl;

            std::cout << ">>> Finding elimination ordering..."  << std::endl;

            // Randomly choose the query
            std::vector<gum::NodeId> allVariables;
            for (auto node : bn.nodes()) {
                allVariables.push_back(node);
            }
            gum::NodeId query;
            std::set<gum::NodeId> chosenVariables;
            std::vector<gum::NodeId> variables;
            bool firstRun = true;
            for (int k = 0; k < bn.size(); k++) {
                if (firstRun) {
                    firstRun = false;
                    query = chooseVariable(allVariables, chosenVariables);
                } else {
                    variables.push_back(chooseVariable(allVariables, chosenVariables));
                }
            }

            std::cout << ">> Chosen QUERY var: " << bn.variable(query).name() << std::endl;
            std::cout << ">> Chosen ELIMIN vars: " << std::endl;
            for (auto var : variables) {
                std::cout << bn.variable(var).name() << " ";
            }
            std::cout << std::endl;

            int (*funcPointer) (gum::NodeId&, gum::UndiGraph&, gum::BayesNet<double>&);
            
            gum::VariableElimination<double> ve(bn);

            // *** MN ***
            std::cout << ">> Min Neighbours with " << bn.size() << " variables." << std::endl;
            funcPointer = &minNeighbours;
            std::vector<gum::NodeId> eliminationOrderMN = findEliminationOrder(variables, bn, funcPointer);
            for (auto varId : eliminationOrderMN) {
                std::cout << bn.variable(varId).name() << " ";
            }
            std::cout << std::endl;
            // do inference
            gum::Potential<double> answerMN;
            ve.setEliminiationOrder(eliminationOrderMN);
            std::clock_t startMN;
            double durationMN;
            startMN = std::clock();
            ve.makeInference();
            answerMN = ve.posterior(query);
            durationMN = ( std::clock() - startMN ) / (double) CLOCKS_PER_SEC;
            std::cout << answerMN << std::endl;
            std::cout << ">>> Time: " << durationMN << std::endl;
            // Record result
            resultsFile << runNumber << "," << databaseName << "," << "MN" << "," << durationMN << "\n";

            // *** MW ***
            std::cout << ">> Min Weight with " << bn.size() << " variables." << std::endl;
            funcPointer = &minWeight;
            std::vector<gum::NodeId> eliminationOrderMW = findEliminationOrder(variables, bn, funcPointer);
            for (auto varId : eliminationOrderMW) {
                std::cout << bn.variable(varId).name() << " ";
            }
            std::cout << std::endl;
            // do inference
            gum::Potential<double> answerMW;
            ve.setEliminiationOrder(eliminationOrderMW);
            std::clock_t startMW;
            double durationMW;
            startMW = std::clock();
            ve.makeInference();
            answerMW = ve.posterior(query);
            durationMW = ( std::clock() - startMW ) / (double) CLOCKS_PER_SEC;
            std::cout << answerMW << std::endl;
            std::cout << ">>> Time: " << durationMW << std::endl;
            // Record result
            resultsFile << runNumber << "," << databaseName << "," << "MW" << "," << durationMW << "\n";

            // *** MF ***
            std::cout << ">> Min Fill with " << bn.size() << " variables." << std::endl;
            funcPointer = &minFill;
            std::vector<gum::NodeId> eliminationOrderMF = findEliminationOrder(variables, bn, funcPointer);
            for (auto varId : eliminationOrderMF) {
                std::cout << bn.variable(varId).name() << " ";
            }
            std::cout << std::endl;
            // do inference
            gum::Potential<double> answerMF;
            ve.setEliminiationOrder(eliminationOrderMF);
            std::clock_t startMF;
            double durationMF;
            startMF = std::clock();
            ve.makeInference();
            answerMF = ve.posterior(query);
            durationMF = ( std::clock() - startMF ) / (double) CLOCKS_PER_SEC;
            std::cout << answerMF << std::endl;
            std::cout << ">>> Time: " << durationMF << std::endl;
            // Record result
            resultsFile << runNumber << "," << databaseName << "," << "MF" << "," << durationMF << "\n";

            // *** WMF ***
            std::cout << ">> WeightedMinFill with " << bn.size() << " variables." << std::endl;
            funcPointer = &weightedMinFill;
            std::vector<gum::NodeId> eliminationOrderWMF = findEliminationOrder(variables, bn, funcPointer);
            for (auto varId : eliminationOrderWMF) {
                std::cout << bn.variable(varId).name() << " ";
            }
            std::cout << std::endl;
            // do inference
            gum::Potential<double> answerWMF;
            ve.setEliminiationOrder(eliminationOrderWMF);
            std::clock_t startWMF;
            double durationWMF;
            startWMF = std::clock();
            ve.makeInference();
            answerWMF = ve.posterior(query);
            durationWMF = ( std::clock() - startWMF ) / (double) CLOCKS_PER_SEC;
            std::cout << answerWMF << std::endl;
            std::cout << ">>> Time: " << durationWMF << std::endl;
            // Record result
            resultsFile << runNumber << "," << databaseName << "," << "WMF" << "," << durationWMF << "\n";

            // *** PE ***
            std::cout << ">> PopulationEnergy with " << bn.size() << " variables." << std::endl;
            std::vector<gum::NodeId> eliminationOrderPE = populationEnergyOrder(variables, bn);
            for (auto varId : eliminationOrderPE) {
                std::cout << bn.variable(varId).name() << " ";
            }
            std::cout << std::endl;
            // do inference
            gum::Potential<double> answerPE;
            ve.setEliminiationOrder(eliminationOrderPE);
            std::clock_t startPE;
            double durationPE;
            startPE = std::clock();
            ve.makeInference();
            answerPE = ve.posterior(query);
            durationPE = ( std::clock() - startPE ) / (double) CLOCKS_PER_SEC;
            std::cout << answerPE << std::endl;
            std::cout << ">>> Time: " << durationPE << std::endl;
            // Record result
            resultsFile << runNumber << "," << databaseName << "," << "PE" << "," << durationPE << "\n";

        }  // finish run
        // Close file
        resultsFile.close();
    } // finish database
}



gum::NodeId chooseVariable(std::vector<gum::NodeId>& allVariables, std::set<gum::NodeId>& chosenVariables) {
    gum::NodeId chosenVar;
    while(true) {
        srand(time(NULL));
        int rnd = rand();
        int randomIdx = rnd % allVariables.size();
        chosenVar = allVariables[randomIdx];
        // std::cout << "------" << std::endl;
        // std::cout << "Random: " << rnd << std::endl;
        // std::cout << "allVariables.size(): " << allVariables.size() << std::endl;
        // std::cout << "Random index: " << randomIdx << std::endl;
        // std::cout << "Chosen var: " << chosenVar << std::endl;
        if (chosenVariables.find(chosenVar) == chosenVariables.end()) {
            // std::cout << "Found!: " << chosenVar << std::endl;
            chosenVariables.insert(chosenVar);
            allVariables.erase(std::remove(allVariables.begin(), allVariables.end(), chosenVar), allVariables.end());
            break;
        }
    }
    return chosenVar;
}
