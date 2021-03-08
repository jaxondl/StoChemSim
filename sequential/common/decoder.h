#ifndef DECODER_H
#define DECODER_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

class decoder {
private:
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
};

#endif //DECODER_H
