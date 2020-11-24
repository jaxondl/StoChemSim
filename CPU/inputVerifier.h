#ifndef INPUTVERIFIER_H
#define INPUTVERIFIER_H

#include <string>
#include <iostream>
#include <fstream>

using namespace std;


class inputVerifier {
public:
    bool verifyFile(string iFile);
    string chopOffComments(string line);
    bool checkReactionSlice(string reactionSlice, int lineNumber, bool errorExists);
    bool containsNonDigitNonDecimal(string reactionSlice);
    bool isReactionRate(string reactionSlice);
    bool isValidReactionRate(string reactionSlice, int lineNumber);
    bool checkReactantsOrProducts(string reactionDefLine, int lineNumber, bool errorExists);
    bool checkSingleMolecule(string molNumAndName, int lineNumber, bool errorExists);
    int checkReactionDefLine(ifstream inputFile, string reactionDefLine);
};


#endif //INPUTVERIFIER_H
