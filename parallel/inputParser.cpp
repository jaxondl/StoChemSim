//
// Created by user on 3/5/2021.
//

#include "inputParser.h"

using namespace std;

//the methods in this class assume that the input file is formatted correctly according to the documentation

//modified to return whether or not the number of simulations was specified in the input file
bool gpuDecoderPrototype::decode(string iFile) {
    bool ret = false;
    ifstream inputFile;
    inputFile.open(iFile);
    if (!inputFile) {
        cout << "Unable to open file" << endl;
        exit(1);   // call system to stop
    }

    //cout << "Beginning decoding." << endl;

    string inputLine;
    int numReactions;
    string reactionSlice = "";
    int reactionNumber = -1; //this is used for indexing the output vectors/arrays
    bool firstEntry = true;

    if (inputFile.peek() != EOF) {
        //cout << "got here" << endl;
        getline(inputFile, inputLine);
        inputLine = chopOffComments(inputLine);

        //cout << "got here" << endl;

        if (inputLine.empty()) { //if the line is only comments, keep checking the next line until you get to a non-comment-only line
            while(inputLine.empty()) {
                getline(inputFile,inputLine);
                inputLine = chopOffComments(inputLine);
            }
        }

        int spaceIndex = inputLine.find(" ");
        numReactions = stoi(inputLine.substr(0, spaceIndex));

        // Try to read number of simulations from input file, do not hang if not present (may be in command line)
        this->numSimulations = -1;
        for (int idx = spaceIndex; idx < inputLine.length(); idx++) {
            if (inputLine[idx] != ' ') {
                this->numSimulations = stoi(inputLine.substr(spaceIndex + 1, inputLine.length()));
                ret = true;
                break;
            }
        }

        //numReactions = stoi(inputLine);
        //cout << "Space " << spaceIndex << endl;
        //cout << "numReactions is " << numReactions << endl;
    }

    this->numberOfReactions = 0;
    for (int r = 0; r < numReactions; r++) {
        reactionNumber++;
        this->numberOfReactions++;
        //cout << "working on reaction " << reactionNumber << endl;
        bool isReversible = false;
        bool isReactant = true; //set to false one we reach the "->" or "<->"
        bool fencepost = false; //set to true once we reach the fencepost after the while loop

        bool createdNewVectorForThisReaction = false;

        //cout << "got here" << endl;
        if (inputFile.peek() != EOF) {
            getline(inputFile, inputLine);

            //chop off the comments, if any (denoted by the first instance of #)
            inputLine = chopOffComments(inputLine);

            if (!inputLine.empty()) {
                if (inputLine.find("<->") != std::string::npos) { //if the reaction definition line contains "<->"
                    isReversible = true;
                    this->numberOfReactions++;
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
                        reactionSlice += inputLine.at(i); // copy the string up to the next space (and don't include that space)
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

                    if(!createdNewVectorForThisReaction) {
                        vector <pair<int, int>> stateChangesForReactions;
                        this->stateChangeVector.push_back(stateChangesForReactions); //this will be pushed to index reactionNumber

                        vector <pair<int, int>> reactantsForReactions;
                        this->reactantVector.push_back(reactantsForReactions); //this will be pushed to index reactionNumber

                        if(isReversible) { //do it again for the reverse direction reaction
                            vector <pair<int, int>> stateChangesForReactions;
                            this->stateChangeVector.push_back(stateChangesForReactions);

                            vector <pair<int, int>> reactantsForReactions;
                            this->reactantVector.push_back(reactantsForReactions);
                        }

                        createdNewVectorForThisReaction = true;
                    }

                    //check this reaction slice; if it's a reactant, update the reactants vector and the state change vector; if it's a product, update the state change vector
                    updateReactantsVector(reactionNumber, reactionSlice, isReactant);
                    if(isReversible) updateReactantsVectorReverse(reactionNumber + 1, reactionSlice, isReactant);

                    updateStateChangeVector(reactionNumber, reactionSlice, isReactant);
                    if(isReversible) updateStateChangeVectorReverse(reactionNumber + 1, reactionSlice, isReactant);



                    //increase j until it's at the index after the space
                    j = spaceIndex + 1;
                    //cout << j << endl;
                }
            }
            if (isReversible) {
                reactionNumber++;
            }
        }

    } // done checking reaction definitions now

    //next check the initial population sizes
    for (string x : this->listOfSpecies) {
        this->populationSizes.push_back(0); //these 2 vectors should be the same length,
        //and any species not mentioned in the input file will have a default initial population size of 0
    }

    // keep taking in lines until one is not blank

    while (true) {
        //cout << "Now checking line " << lineNumber << endl;
        if (inputFile.peek() != EOF) {
            getline(inputFile, inputLine); //break if this line is not blank after removing comments
            inputLine = chopOffComments(inputLine);
            if (!inputLine.empty()) break;
        } else { // they didn't give the initial molecule populations
            // (this should never happen because we called inputVerifier->verifyFile before running any of the code in this file)
            cout << "The initial population sizes are not included in the input file, which means inputVerifier.cpp is bugged." << endl;
        }
    }

    // The remaining lines should define molecule populations
    while (inputFile.peek() != EOF) { // check all remaining lines
        //cout << "Now checking line " << lineNumber << endl;
        if (!firstEntry) getline(inputFile, inputLine); //the first time we enter, the line has already been read
        firstEntry = false;
        inputLine = chopOffComments(inputLine);

        //get the molecule name
        string moleculeName = "";
        int j = 0;
        while (inputLine.at(j) != ' ') {
            moleculeName += inputLine.at(j);
            j++;
        }
        j++; //j now points to the start of the number

        //now get the starting population size
        string moleculeCount = "";
        while (j < inputLine.length()) {
            moleculeCount += inputLine.at(j);
            j++;
        }

        //get the index of moleculeName in this->listOfSpecies
        int moleculeIndex = -1;
        bool found = false;
        for (int i = 0; i < this->listOfSpecies.size(); i++) {
            if (this->listOfSpecies[i] == moleculeName) {
                moleculeIndex = i;
                found = true;
                //cout << "moleculeIndex is now " << moleculeIndex << endl;
                break;
            }
        }

        //store moleculeCount in that index of this->populationSizes
        if (found) {
            this->populationSizes[moleculeIndex] = stoi(moleculeCount);
        }
        else { //there exists a molecule in the system that is not a reactant or product for any reaction
            //add the molecule to the list of species, and give it a default population size of 0
            this->listOfSpecies.push_back(moleculeName);
            this->populationSizes.push_back(0);
            //then add its population size
            for (int i = 0; i < this->listOfSpecies.size(); i++) {
                if (this->listOfSpecies[i] == moleculeName) {
                    moleculeIndex = i;
                    //cout << "moleculeIndex is now " << moleculeIndex << endl;
                    break;
                }
            }
            this->populationSizes[moleculeIndex] = stoi(moleculeCount);
        }

    }

    inputFile.close(); //done parsing and creating the needed data structures

//    cout << "List of species:" << endl;
//    for (string x : this->listOfSpecies)
//        cout << x << " ";
//    cout << endl;
//    cout << "List of reaction rates for each reaction:" << endl;
//    for (float x : this->kValueVector)
//        cout << x << " ";
//    cout << endl;
//    cout << "reactantVector:" << endl;
//    cout << "[ ";
//    for (vector<pair<int, int>> x : this->reactantVector) {
//        cout << "[ ";
//        for (pair<int, int> y : x) {
//            cout << "[" << y.first << "," << y.second << "]" << ", ";
//        }
//        cout << "]," << endl;
//    }
//    cout << "]" << endl;
//    cout << "stateChangeVector:" << endl;
//    cout << "[ ";
//    for (vector<pair<int, int>> x : this->stateChangeVector) {
//        cout << "[ ";
//        for (pair<int, int> y : x) {
//            cout << "[" << y.first << "," << y.second << "]" << ", ";
//        }
//        cout << "]," << endl;
//    }
//    cout << "]" << endl;
//    cout << "List of initial population sizes for each species:" << endl;
//    for (int x : this->populationSizes)
//        cout << x << " ";
//    cout << endl;

    return ret;
}

string gpuDecoderPrototype::chopOffComments(string line) {
    if (line.find("#") != std::string::npos) {
        string temp = "";
        for (int i = 0; i < line.find("#"); i++) {
            temp += line.at(i);
        }
        line = temp;
    }

    //if there is a space at the end of the line, remove it
    if (!line.empty()) {
        if (line.at(line.length() - 1) == ' ') {
            line = line.substr(0, line.length() - 1); //remove the very last character
        }
    }

    return line;
}

int gpuDecoderPrototype::getNumSimulations() {
    return this->numSimulations;
}

int gpuDecoderPrototype::getNumReactions() {
    return this->numberOfReactions;
}

int gpuDecoderPrototype::getNumSpecies() {
    return this->listOfSpecies.size();
}

//stores reaction rates and molecule names; does not store state changes or reactant amounts
void gpuDecoderPrototype::parseReactionSlice(string reactionSlice, bool isReversible, bool fencepost, int reactionNumber, bool isReactant) { //check the forward direction of the reaction
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
            this->kValueVector.push_back(std::stod(reactionSlice));
        } else { //if fencepost is true, then this is either the forward reaction rate (if not reversible) or the reverse reaction rate (if reversible)
            if (!isReversible) {
                //cout << std::stof(reactionSlice) << endl;
                this->kValueVector.push_back(std::stod(reactionSlice));
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
void gpuDecoderPrototype::parseReverseReactionSlice(string reactionSlice, bool fencepost, int reactionNumber, bool isReactant) {
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
            this->kValueVector.push_back(std::stod(reactionSlice)); //store the reverse reaction rate
        }
    }
}

//access this->reactantVector[reactionNumber] and store the pair {this->listOfSpecies.indexof(moleculeName), moleculeCount} if isReactant is true
void gpuDecoderPrototype::updateReactantsVector(int reactionNumber, string reactionSlice, bool isReactant) {
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
//    bool containsDigit = false;
//    for (int i = 0; i < reactionSlice.length(); i++) {
//        if (isdigit(reactionSlice.at(i))) {
//            containsDigit = true; //1x is the same as x, so if there is a letter but no digit, pretend there is a 1
//        }
//    }
    if(!containsLetter && reactionSlice != "0") { //if the reaction slice contains no letter, then it can't be a reactant or product, so it must be a reaction rate
        //it's a reaction rate, so we can just ignore it
    } else if (reactionSlice == "0") { //there are no reactants or no products for this reaction, so act accordingly
        //nothing to update
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

        if(moleculeCount.empty()) {
            moleculeCount = "1"; //1x is the same as x, so if there is a letter but no digit, pretend there is a 1
        }

        //i now points to the first letter; the remaining part of the string is the molecule name
        while (i < reactionSlice.length()) {
            //cout << i << endl;
            moleculeName += reactionSlice.at(i); //copy the name into a new string
            //cout << "moleculeName is now " << moleculeName << endl;
            i++;
        }
        //cout << "moleculeName is now " << moleculeName << endl;
        //cout << "moleculeCount is now " << moleculeCount << endl;
        //cout << "stoi(moleculeCount) is " << stoi(moleculeCount) << endl;

        if (isReactant) {
            pair<int, int> pairing;
            pairing.first = -1;

            int moleculeIndex = -1;
            bool found = false;
            for (int i = 0; i < this->listOfSpecies.size(); i++) {
                if (this->listOfSpecies[i] == moleculeName) {
                    moleculeIndex = i;
                    found = true;
                    //cout << "moleculeIndex is now " << moleculeIndex << endl;
                    break;
                }
            }

            if (found) {
                pairing.first = moleculeIndex;
                pairing.second = stoi(moleculeCount);
            }

            if(pairing.first != -1) { //this should always get executed if isReactant is true; if it doesn't, something went wrong
                //cout << "pushing pair [" << pairing.first << "," << pairing.second << "]" << endl;
                this->reactantVector[reactionNumber].push_back(pairing);
            }

        }
    }
}

void gpuDecoderPrototype::updateReactantsVectorReverse(int reactionNumber, string reactionSlice, bool isReactant) {
    //cout << "updating reactantVector for " << reactionSlice << endl;
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
        //it's a reaction rate, so we can just ignore it
    } else if (reactionSlice == "0") { //there are no reactants or no products for this reaction, so act accordingly
        //nothing to update
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

        if(moleculeCount.empty()) {
            moleculeCount = "1"; //1x is the same as x, so if there is a letter but no digit, pretend there is a 1
        }

        //i now points to the first letter; the remaining part of the string is the molecule name
        while (i < reactionSlice.length()) {
            //cout << i << endl;
            moleculeName += reactionSlice.at(i); //copy the name into a new string
            //cout << "moleculeName is now " << moleculeName << endl;
            i++;
        }
        //cout << "moleculeName is now " << moleculeName << endl;
        //cout << "moleculeCount is now " << moleculeCount << endl;

        if (!isReactant) { //if it's on the right side of the "<->" symbol
            pair<int, int> pairing;
            pairing.first = -1;

            int moleculeIndex = -1;
            bool found = false;
            for (int i = 0; i < this->listOfSpecies.size(); i++) {
                if (this->listOfSpecies[i] == moleculeName) {
                    //std::cout << "Element present at index " << i;
                    moleculeIndex = i;
                    found = true;
                    //cout << "moleculeIndex is now " << moleculeIndex << endl;
                    break;
                }
            }

            if (found) {
                pairing.first = moleculeIndex;
                pairing.second = stoi(moleculeCount);
            }

            if(pairing.first != -1) { //this should always get executed if isReactant is true; if it doesn't, something went wrong
                //cout << "pushing pair [" << pairing.first << "," << pairing.second << "]" << endl;
                this->reactantVector[reactionNumber].push_back(pairing);
            }

        }
    }
}

void gpuDecoderPrototype::updateStateChangeVector(int reactionNumber, std::string reactionSlice, bool isReactant) {
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
        //it's a reaction rate, so we can just ignore it
    } else if (reactionSlice == "0") { //there are no reactants or no products for this reaction, so act accordingly
        //nothing to update
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

        if(moleculeCount.empty()) {
            moleculeCount = "1"; //1x is the same as x, so if there is a letter but no digit, pretend there is a 1
        }

        //i now points to the first letter; the remaining part of the string is the molecule name
        while (i < reactionSlice.length()) {
            //cout << i << endl;
            moleculeName += reactionSlice.at(i); //copy the name into a new string
            //cout << "moleculeName is now " << moleculeName << endl;
            i++;
        }
        //cout << "moleculeName is now " << moleculeName << endl;
        //cout << "moleculeCount is now " << moleculeCount << endl;
        //cout << "stoi(moleculeCount) is " << stoi(moleculeCount) << endl;

        //get the index of the molecule
        int moleculeIndex = -1;
        for (int i = 0; i < this->listOfSpecies.size(); i++) {
            if (this->listOfSpecies[i] == moleculeName) {
                moleculeIndex = i;
                //cout << "moleculeIndex is now " << moleculeIndex << endl;
                break;
            }
        }

        bool found = false;
        //check all existing pairs in stateChangeVector[reactionNumber] and set found to true if any of them have moleculeIndex at .first
        for (pair<int, int> y : this->stateChangeVector[reactionNumber]) {
            if (y.first == moleculeIndex) found = true;
        }
        //push a [0,0] pair that will be modified shortly ONLY IF there is not already a pair pushed for that molecule
        if(!found) {
            pair<int, int> pairing;
            pairing.first = moleculeIndex;
            pairing.second = 0; //holds the net change in the species population caused by this reaction
            this->stateChangeVector[reactionNumber].push_back(pairing);
        }

        int pairIndex;
        if (isReactant) { //decrement the corresponding pairing.second by moleculeCount
            for (int i = 0; i < this->stateChangeVector[reactionNumber].size(); i++) {
                if (this->stateChangeVector[reactionNumber][i].first == moleculeIndex) {
                    pairIndex = i;
                    break;
                }
            }
            //cout << this->stateChangeVector[reactionNumber][pairIndex].first << " " << this->stateChangeVector[reactionNumber][pairIndex].second << endl;
            //cout << "decrementing for forward state change vector" << endl;
            this->stateChangeVector[reactionNumber][pairIndex].second -= stoi(moleculeCount);
            //cout << this->stateChangeVector[reactionNumber][pairIndex].first << " " << this->stateChangeVector[reactionNumber][pairIndex].second << endl;
        } else { //it's a product, so increment the corresponding pairing.second by moleculeCount
            for (int i = 0; i < this->stateChangeVector[reactionNumber].size(); i++) {
                if (this->stateChangeVector[reactionNumber][i].first == moleculeIndex) {
                    pairIndex = i;
                    break;
                }
            }
            //cout << this->stateChangeVector[reactionNumber][pairIndex].first << " " << this->stateChangeVector[reactionNumber][pairIndex].second << endl;
            //cout << "incrementing for forward state change vector" << endl;
            this->stateChangeVector[reactionNumber][pairIndex].second += stoi(moleculeCount);
            //cout << this->stateChangeVector[reactionNumber][pairIndex].first << " " << this->stateChangeVector[reactionNumber][pairIndex].second << endl;
        }
    }
}

void gpuDecoderPrototype::updateStateChangeVectorReverse(int reactionNumber, std::string reactionSlice, bool isReactant) {
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
        //it's a reaction rate, so we can just ignore it
    } else if (reactionSlice == "0") { //there are no reactants or no products for this reaction, so act accordingly
        //nothing to update
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

        if(moleculeCount.empty()) {
            moleculeCount = "1"; //1x is the same as x, so if there is a letter but no digit, pretend there is a 1
        }

        //i now points to the first letter; the remaining part of the string is the molecule name
        while (i < reactionSlice.length()) {
            //cout << i << endl;
            moleculeName += reactionSlice.at(i); //copy the name into a new string
            //cout << "moleculeName is now " << moleculeName << endl;
            i++;
        }
        //cout << "moleculeName is now " << moleculeName << endl;
        //cout << "moleculeCount is now " << moleculeCount << endl;
        //cout << "stoi(moleculeCount) is " << stoi(moleculeCount) << endl;

        //get the index of the molecule
        int moleculeIndex = -1;
        for (int i = 0; i < this->listOfSpecies.size(); i++) {
            if (this->listOfSpecies[i] == moleculeName) {
                moleculeIndex = i;
                //cout << "moleculeIndex is now " << moleculeIndex << endl;
                break;
            }
        }

        bool found = false;
        //check all existing pairs in stateChangeVector[reactionNumber] and set found to true if any of them have moleculeIndex at .first
        for (pair<int, int> y : this->stateChangeVector[reactionNumber]) {
            if (y.first == moleculeIndex) found = true;
        }

        //push a [0,0] pair that will be modified shortly ONLY IF there is not already a pair pushed for that molecule
        if(!found) {
            pair<int, int> pairing;
            pairing.first = moleculeIndex;
            pairing.second = 0; //holds the net change in the species population caused by this reaction
            this->stateChangeVector[reactionNumber].push_back(pairing);
        }

        int pairIndex;
        if (!isReactant) { //decrement the corresponding pairing.second by moleculeCount
            for (int i = 0; i < this->stateChangeVector[reactionNumber].size(); i++) {
                if (this->stateChangeVector[reactionNumber][i].first == moleculeIndex) {
                    pairIndex = i;
                    break;
                }
            }
            //cout << this->stateChangeVector[reactionNumber][pairIndex].first << " " << this->stateChangeVector[reactionNumber][pairIndex].second << endl;
            //cout << "decrementing for forward state change vector" << endl;
            this->stateChangeVector[reactionNumber][pairIndex].second -= stoi(moleculeCount);
            //cout << this->stateChangeVector[reactionNumber][pairIndex].first << " " << this->stateChangeVector[reactionNumber][pairIndex].second << endl;
        } else { //it's a product, so increment the corresponding pairing.second by moleculeCount
            for (int i = 0; i < this->stateChangeVector[reactionNumber].size(); i++) {
                if (this->stateChangeVector[reactionNumber][i].first == moleculeIndex) {
                    pairIndex = i;
                    break;
                }
            }
            //cout << this->stateChangeVector[reactionNumber][pairIndex].first << " " << this->stateChangeVector[reactionNumber][pairIndex].second << endl;
            //cout << "incrementing for forward state change vector" << endl;
            this->stateChangeVector[reactionNumber][pairIndex].second += stoi(moleculeCount);
            //cout << this->stateChangeVector[reactionNumber][pairIndex].first << " " << this->stateChangeVector[reactionNumber][pairIndex].second << endl;
        }
    }
}

vector<string> gpuDecoderPrototype::getListOfSpecies() {
    return this->listOfSpecies;
}
vector<int> gpuDecoderPrototype::getPopulationSizes() {
    return this->populationSizes;
}
vector<vector<pair<int, int>>> gpuDecoderPrototype::getStateChangeVector() {
    return this->stateChangeVector;
}
vector<vector<pair<int, int>>> gpuDecoderPrototype::getReactantVector() {
    return this->reactantVector;
}
vector<double> gpuDecoderPrototype::getkValueVector() {
    return this->kValueVector;
}
vector<double> gpuDecoderPrototype::getRRCVector() { //done
    return this->kValueVector;
}
vector<int> gpuDecoderPrototype::getConfigurationMatrix() { //done
    //this is a flattened vector of ints
    //for each row, add the initial population size of each species, in order of index
    //should have numSpecies*numSimulations elements in total
    for (int i = 0; i < this->numSimulations; i++ ) {
        for (int speciesIndex = 0; speciesIndex < listOfSpecies.size(); speciesIndex++) {
            this->configuration_matrix.push_back(populationSizes[speciesIndex]);
        }
    }
    return this->configuration_matrix;
}
vector<int>gpuDecoderPrototype::getConfigurationMatrix(int s) {
    for (int i = 0; i < s; i++) {
        for (int speciesIndex = 0; speciesIndex < listOfSpecies.size(); speciesIndex++) {
            this->configuration_matrix.push_back(populationSizes[speciesIndex]);
        }
    }
    return this->configuration_matrix;
}
vector<int> gpuDecoderPrototype::getReactantsTableVector() { //done
    int maxNumReactantsInAReaction = 0;
    for (int i = 0; i < this->reactantVector.size(); i++) {
        if (this->reactantVector[i].size() > maxNumReactantsInAReaction) maxNumReactantsInAReaction = this->reactantVector[i].size();
    }
    //fill/pad reactants_table with the correct number of zeroes
    for (int i = 0; i < this->numberOfReactions; i++) {
        for (int j = 0; j < maxNumReactantsInAReaction; j++) {
            this->reactants_table.push_back(0); //this slot in the pair represents the index of the species
            this->reactants_table.push_back(0); //this slot in the pair represents the number of molecules of that species needed as reactants
        }
    }

    //go through each reactant/moleculeCount pair and update reactants_table accordingly
    int index;
    int reactionNumber = 0;
    for (vector<pair<int, int>> x : this->reactantVector) {
        index = 2*maxNumReactantsInAReaction*reactionNumber;
        for (pair<int, int> y : x) {
            this->reactants_table[index] = y.first;
            index++;
            this->reactants_table[index] = y.second;
            index++;
        }
        reactionNumber++;
    }
    return this->reactants_table;
}
vector<int> gpuDecoderPrototype::getStateChangeMatrix() { //done
    int numSpecies = listOfSpecies.size();
    //fill/pad state_change_matrix with the correct number of zeroes
    for (int i = 0; i < this->numberOfReactions; i++) {
        for (int j = 0; j < numSpecies; j++) {
            this->state_change_matrix.push_back(0);
        }
    }

    int reactionNumber = 0;
    for (vector<pair<int, int>> x : this->stateChangeVector) {
        for (pair<int, int> y : x) {
            int speciesIndex = y.first;
            int stateChange = y.second;
            this->state_change_matrix.at(speciesIndex+(reactionNumber*numSpecies)) = stateChange;
        }
        reactionNumber++;
    }
    return this->state_change_matrix;
}
vector<double> gpuDecoderPrototype::getPropensityMatrix() { //done
    for (int i = 0; i < this->numSimulations; i++) { //each row is initially identical, so i is not used
        for (int j = 0; j < this->numberOfReactions; j++) {
            double propValue = calcPropensity(this->kValueVector[j], this->populationSizes, this->reactantVector[j]);
            propensity_matrix.push_back(propValue);
        }
    }
    return this->propensity_matrix;
}

vector<double> gpuDecoderPrototype::getPropensityMatrix(int s) {
    return vector<double>(s * this->numberOfReactions, 0);
    //for (int i = 0; i < s; i++) { //each row is initially identical, so i is not used
    //    for (int j = 0; j < this->numberOfReactions; j++) {
    //        double propValue = calcPropensity(this->kValueVector[j], this->populationSizes, this->reactantVector[j]);
    //        propensity_matrix.push_back(propValue);
    //    }
    //}
    //return this->propensity_matrix;
}

double gpuDecoderPrototype::calcPropensity(double reactionRate, vector<int> moleculeAmounts, vector<pair<int, int>> reactants){
    double propensity = reactionRate;
    for (pair<int, int> reactant: reactants) {
        for (int i=0; i<reactant.second; i++) {
            propensity *= (moleculeAmounts[reactant.first] - i);
        }
    }
    return propensity;
}

void gpuDecoderPrototype::printVectors(){
    vector<double> rrc_vector = this->getRRCVector();
    vector<int> state_change_matrix = this->getStateChangeMatrix();
    vector<int> configuration_matrix = this->getConfigurationMatrix();
    vector<double> propensity_matrix = this->getPropensityMatrix();
    vector<int> reactants_table = this->getReactantsTableVector();

    cout << "List of species:" << endl;
    for (string x : this->listOfSpecies)
        cout << x << " ";
    cout << endl;
    cout << "rrc_vector: " << endl;
    cout << "[ ";
    for (float x : rrc_vector)
        cout << x << " ";
    cout << "]" << endl;
    cout << "reactants_table:" << endl;
    cout << "[ ";
    for (int x : reactants_table) {
        cout << x << " ";
    }
    cout << "]" << endl;
    cout << "state_change_matrix:" << endl;
    cout << "[ ";
    for (int x : state_change_matrix) {
        cout << x << " ";
    }
    cout << "]" << endl;
    cout << "configuration_matrix:" << endl;
    cout << "[ ";
    for (int x : configuration_matrix)
        cout << x << " ";
    cout << "]" << endl;
    cout << "propensity_matrix:" << endl;
    cout << "[ ";
    for (double x : propensity_matrix)
        cout << x << " ";
    cout << "]" << endl;
}