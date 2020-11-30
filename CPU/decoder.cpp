#include "decoder.h"

using namespace std;

//the methods in this class assume that the input file is formatted correctly according to the documentation

void decoder::testFunction() {
    cout << "This is a test function for decoder.cpp" << endl;
}

void decoder::decode(string iFile) {
    ifstream inputFile;
    inputFile.open(iFile);
    if (!inputFile) {
        cerr << "Unable to open file";
        exit(1);   // call system to stop
    }

    string inputLine;
    int numReactions;
    string reactionSlice = "";
    int reactionNumber = -1; //this is used for indexing the output vectors/arrays

    if (inputFile.peek() != EOF) {
        //out << "got here" << endl;
        getline(inputFile, inputLine);
        inputLine = chopOffComments(inputLine);

        if (inputLine.at(inputLine.length() - 1) == ' ') { //if there is a space at the end of the line, remove it
            inputLine = inputLine.substr(0, inputLine.length() - 1);
        }

        numReactions = stoi(inputLine);
        //cout << numReactions << endl;
    }

    for (int r = 0; r < numReactions; r++) {
        reactionNumber++;
        bool isReversible = false;
        bool isReactant = true; //set to false one we reach the "->" or "<->"
        bool fencepost = false; //set to true once we reach the fencepost after the while loop

        //cout << "got here" << endl;
        if (inputFile.peek() != EOF) {
            getline(inputFile, inputLine);

            //chop off the comments, if any (denoted by the first instance of #)
            inputLine = chopOffComments(inputLine);

            if (!inputLine.empty()) {
                //if there is a space at the end of the line, remove it
                if (inputLine.at(inputLine.length() - 1) == ' ') {
                    inputLine = inputLine.substr(0, inputLine.length() - 1);
                }

                if (inputLine.find("<->") != std::string::npos) { //if the reaction definition line contains "<->"
                    isReversible = true;
                }

                int j = 0;
                //cout << "inputLine is " << inputLine << endl;
                while (inputLine.find(" ", j) !=
                       std::string::npos) { // while there is still a space in the line on or after index j
                    //cout << "entering while loop on line 53" << endl;
                    reactionSlice = "";
                    //find the next space
                    int spaceIndex = inputLine.find(" ", j);
                    for (int i = j; i < spaceIndex; i++) {
                        reactionSlice += inputLine.at(i); // copy the string up to the next space
                        // (and don't include that space)
                    }
                    //cout << reactionSlice << " is about to be checked" << endl;

                    parseReactionSlice(reactionSlice, isReversible, fencepost, reactionNumber, isReactant);
                    if (isReversible) parseReverseReactionSlice(reactionSlice, fencepost, reactionNumber, isReactant);

                    //increase j until it's at the index after the space
                    j = spaceIndex + 1;
                    //cout << j << endl;
                }
                //now we have reached the fencepost
                fencepost = true;
                //everything after the last space should be a reaction rate (the 2nd if the line contains "<->", or the first if it contains "->")
                reactionSlice = "";
                for (int i = j; i < inputLine.length(); i++) {
                    reactionSlice += inputLine.at(i); // copy the remaining characters
                }
                //cout << "now checking " << reactionSlice << endl;
                parseReactionSlice(reactionSlice, isReversible, fencepost, reactionNumber, isReactant);
                if (isReversible) parseReverseReactionSlice(reactionSlice, fencepost, reactionNumber, isReactant);

                //now all species names and reaction rates from this reaction definition have been stored
                //next, we need to go through the reaction definition again, and store the reactant amounts and state change vectors
                j = 0;
                while (inputLine.find(" ", j) !=
                       std::string::npos) { // while there is still a space in the line on or after index j
                    //cout << "entering while loop on line 105" << endl;
                    reactionSlice = "";
                    //find the next space
                    int spaceIndex = inputLine.find(" ", j);
                    for (int i = j; i < spaceIndex; i++) {
                        reactionSlice += inputLine.at(i); // copy the string up to the next space
                        // (and don't include that space)
                    }
                    //cout << reactionSlice << " is about to be checked" << endl;
                    if (reactionSlice == "->" || reactionSlice == "<->") { //products from here on out for this line
                        isReactant = false;
                    }

                    vector<vector<int>> stateChangesForReactions;
                    this->stateChangeVector.push_back(
                            stateChangesForReactions); //this will be pushed to index reactionNumber
                    vector<int> stateChange{0, 0};
                    this->stateChangeVector[reactionNumber].push_back(
                            stateChange); //push an empty 2-element vector that will be modified by the below method calls

                    vector<vector<int>> reactantsForReactions;
                    this->reactantVector.push_back(
                            reactantsForReactions); //this will be pushed to index reactionNumber
                    vector<int> reactantAmount{0, 0}; //first integer is the index of the reactant species in listOfSpecies; second integer is the number of that reactant needed for the reaction to occur
                    this->reactantVector[reactionNumber].push_back(
                            reactantAmount); //push an empty 2-element vector that will be modified by the below method calls

                    //check this reaction slice; if it's a reactant, update the reactants vector and the state change vector; if it's a product, update the state change vector
                    updateReactantsVector(reactionNumber, reactionSlice, isReactant);
                    if(isReversible) updateReactantsVectorReverse(reactionNumber, reactionSlice, isReactant);

                    updateStateChangeVector(reactionNumber, reactionSlice, isReactant);
                    if(isReversible) updateStateChangeVectorReverse(reactionNumber, reactionSlice, isReactant);



                    //increase j until it's at the index after the space
                    j = spaceIndex + 1;
                    //cout << j << endl;
                }
            }
        }

    } // done checking reaction definitions now

    cout << "List of species:" << endl;
    for (string x : this->listOfSpecies)
        cout << x << " ";
    cout << endl;
    cout << "List of reaction rates for each reaction:" << endl;
    for (float x : this->kValueVector)
        cout << x << " ";
    cout << endl;


}

string decoder::chopOffComments(string line) {
    if (line.find("#") != std::string::npos) {
        string temp = "";
        for (int i = 0; i < line.find("#"); i++) {
            temp += line.at(i);
        }
        line = temp;
    }
    return line;
}

//stores reaction rates and molecule names; does not store state changes or reactant amounts
void decoder::parseReactionSlice(string reactionSlice, bool isReversible, bool fencepost, int reactionNumber, bool isReactant) { //check the forward direction of the reaction
    //cout << "Now checking " << reactionSlice << endl;

    //if the reaction slice is "->" or "<->", ignore it
    if (reactionSlice == "->" || reactionSlice == "<->") {
        return;
    }
    string moleculeName = "";
    bool containsLetter = false;
    for (int i = 0; i < reactionSlice.length(); i++) {
        if (isalpha(reactionSlice.at(i))) {
            containsLetter = true;
        }
    }
    if(!containsLetter && reactionSlice != "0") { //if the reaction slice contains no letter, then it can't be a reactant or product, so it must be a reaction rate
        if(!fencepost) { //if fencepost is false, this must be the forward direction reaction rate of the reaction (and isReversible must be true, otherwise, something is wrong)
            //store the reaction rate
            //cout << std::stof(reactionSlice) << endl;
            this->kValueVector.push_back(std::stof(reactionSlice));
        } else { //if fencepost is true, then this is either the forward reaction rate (if not reversible) or the reverse reaction rate (if reversible)
            if (!isReversible) {
                //cout << std::stof(reactionSlice) << endl;
                this->kValueVector.push_back(std::stof(reactionSlice));
            } else {
                //this is the reverse reaction rate, so let parseReverseReactionSlice handle it
            }
        }
    } else if (reactionSlice == "0") { //there are no reactants or no products for this reaction, so act accordingly

    } else { //the slice contains a letter, so it must be a reactant or a product
        //parse until first non-digit character
        string moleculeCount = "";
        int i = 0;
        while (i < reactionSlice.length()) {
            //cout << i << endl;
            if (isdigit(reactionSlice.at(i))) {
                moleculeCount += reactionSlice.at(i); //copy the digits into a new string
                //cout << "moleculeCount is now " << moleculeCount << endl;
            } else break;
            i++;
        } //moleculeCount now contains the defined number of reactant/product molecules
        //i now points to the first letter; the remaining part of the string is the molecule name
        while (i < reactionSlice.length()) {
            //cout << i << endl;
            moleculeName += reactionSlice.at(i); //copy the name into a new string
            //cout << "moleculeName is now " << moleculeName << endl;
            i++;
        }

        //push the molecule name into listOfSpecies if it is not already there
        if (std::count(this->listOfSpecies.begin(), this->listOfSpecies.end(), moleculeName) <= 0) {
            //if the vector does not already contain this molecule name
            this->listOfSpecies.push_back(moleculeName);
        }
    }

}

//only gets called if the reaction is reversible
void decoder::parseReverseReactionSlice(string reactionSlice, bool fencepost, int reactionNumber, bool isReactant) {
    //forward direction has already been checked, so we won't see any new molecule names here
    //just store the reverse reaction rate

    //if the reaction slice is "->" or "<->", ignore it
    if (reactionSlice == "->" || reactionSlice == "<->") {
        return;
    }
    string moleculeName = "";
    bool containsLetter = false;
    for (int i = 0; i < reactionSlice.length(); i++) {
        if (isalpha(reactionSlice.at(i))) {
            containsLetter = true;
        }
    }

    if(!containsLetter && reactionSlice != "0") { //it's a reaction rate
        if(fencepost) { //if it's not the fencepost, it's not the reverse reaction rate, and vice versa
            //cout << std::stof(reactionSlice) << endl;
            this->kValueVector.push_back(std::stof(reactionSlice)); //store the reverse reaction rate
        }
    }
}

void decoder::updateReactantsVector(int reactionNumber, string reactionSlice, bool isReactant) {

}
void decoder::updateReactantsVectorReverse(int reactionNumber, string reactionSlice, bool isReactant) {

}
void decoder::updateStateChangeVector(int reactionNumber, std::string reactionSlice, bool isReactant) {

}
void decoder::updateStateChangeVectorReverse(int reactionNumber, std::string reactionSlice, bool isReactant) {

}
