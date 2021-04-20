#ifndef SENIORDESIGN_GPUDECODER_H
#define SENIORDESIGN_GPUDECODER_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

class gpuDecoder {
private:
    int numSimulations; //exclusive to GPU implementation
    int numberOfReactions; //exclusive to GPU implementation
    std::vector<std::string> listOfSpecies;
    std::vector<int> populationSizes; //indices correspond to the indices in listOfSpecies
    std::vector<std::vector<std::pair<int, int>>> stateChangeVector; //a vector of vectors of pairs
    //each element x of stateChangeArray corresponds to one defined reaction
    //each vector element y of a given element x is a pair of integer elements
    //the first integer is the index of the reactant/product,
    // and the second integer is the net change in the population size of that reactant/product after the reaction occurs
    std::vector<std::vector<std::pair<int, int>>> reactantVector; //a vector of vectors of pairs
    //each element x of reactantArray corresponds to one defined reaction
    //each vector element y of a given element x is a pair of integer elements
    //the first integer is the index of the reactant
    //the second integer is how many copies of that reactant molecule are needed for the reaction to occur
    std::vector<double> kValueVector; //each element is the k value of the corresponding reaction
    std::vector<int> configuration_matrix;
    std::vector<int> state_change_matrix;
    std::vector<int> reactants_table;
    std::vector<double> propensity_matrix;
    std::vector<int> header_indices;
    std::vector<std::string> custom_header;
public:
    void decode(std::string iFile);
    std::string chopOffComments(std::string line);
    void parseReactionSlice(std::string reactionSlice, bool isReversible, bool fencepost, int reactionNumber, bool isReactant);
    void parseReverseReactionSlice(std::string reactionSlice, bool fencepost, int reactionNumber, bool isReactant);
    void updateReactantsVector(int reactionNumber, std::string reactionSlice, bool isReactant);
    void updateReactantsVectorReverse(int reactionNumber, std::string reactionSlice, bool isReactant);
    void updateStateChangeVector(int reactionNumber, std::string reactionSlice, bool isReactant);
    void updateStateChangeVectorReverse(int reactionNumber, std::string reactionSlice, bool isReactant);
    std::vector<std::string> getListOfSpecies();
    std::vector<int> getPopulationSizes();
    std::vector<std::vector<std::pair<int, int>>> getStateChangeVector();
    std::vector<std::vector<std::pair<int, int>>> getReactantVector();
    std::vector<double> getkValueVector();

    //the above methods and fields are identical to the CPU decoder
    //the following methods are used for the GPU implementation

    double calcPropensity(double reactionRate, std::vector<int> moleculeCounts, std::vector<std::pair<int,int>> reactants);
    std::vector<int> getStateChangeMatrix();
    std::vector<int> getConfigurationMatrix();
    std::vector<double> getRRCVector();
    std::vector<int> getReactantsTableVector();
    std::vector<double> getPropensityMatrix();
    std::vector<int> getHeaderIndices();
    std::vector<std::string> getCustomHeader();
    void printVectors();
};

#endif //SENIORDESIGN_GPUDECODER_H
