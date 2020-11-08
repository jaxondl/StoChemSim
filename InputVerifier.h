//
// Created by user on 10/30/2020.
//

#ifndef SSA_IMPLEMENTATION_INPUTVERIFIER_H
#define SSA_IMPLEMENTATION_INPUTVERIFIER_H

#include <string>

using namespace std;


class InputVerifier {
public:
    bool verifyFile(string iFile);
    bool checkReactantsOrProducts(string reactionDefLine, int lineNumber, bool errorExists);
    bool checkSingleMolecule(string molNumAndName, int lineNumber, bool errorExists);
    int checkReactionDefLine(ifstream inputFile, string reactionDefLine);
};


#endif //SSA_IMPLEMENTATION_INPUTVERIFIER_H
