#include "inputVerifier.h"
#include <iostream>
#include <fstream>


//returns true if the file has no formatting errors, false otherwise
bool inputVerifier::verifyFile(string iFile) {
    cout << "Beginning input verification." << endl;
    string fullReactionDefLine;
    string reactionDefLine = "";
    string reactionRate = "";
    string reactionSlice = "";
    string moleculeDefLine;
    int lineNumber = 0; //used for warning messages
    int numReactions;

    bool errorExists = false;

    bool firstEntry = true;

    ifstream inputFile;
    inputFile.open(iFile);
    if (!inputFile) {
        cerr << "Unable to open file";
        exit(1);   // call system to stop
    }

    // things to verify:
    // 1) first line contains a single number and nothing else
    if (inputFile.peek() != EOF) {
        //out << "got here" << endl;
        getline(inputFile,
                fullReactionDefLine); //fullReactionDefLine is used here instead of creating a different string that would only be used once
        bool valid = true;
        lineNumber++;
        //cout << fullReactionDefLine << endl;
        fullReactionDefLine = chopOffComments(fullReactionDefLine); //get rid of the comments, if any
        //cout << fullReactionDefLine << endl;
        //if there is a space at the end of the line, remove it
        if (fullReactionDefLine.at(fullReactionDefLine.length() - 1) == ' ') {
            fullReactionDefLine = fullReactionDefLine.substr(0, fullReactionDefLine.length() -
                                                                1); //remove the very last character
        }
        //cout << "Now checking line " << lineNumber << endl;
        //cout << reactionDefLine.length() << endl;
        for (int i = 0; i < fullReactionDefLine.length(); i++) {
            if (!(isdigit(fullReactionDefLine.at(i))) && !errorExists) {
                //first line contains a character that is not a number!
                cout << "Warning: Your first line contains a character that isn't a number (it might be a space)." << endl;
                valid = false;
                errorExists = true;
            }
        }
        //cout << (fullReactionDefLine.at(fullReactionDefLine.length() - 1)) << endl;
        if (valid) {
            //cout << "number of reactions is: " << fullReactionDefLine << endl;
            numReactions = stoi(fullReactionDefLine);
            //cout << numReactions << endl;
        }
    } else {
        errorExists = true;
        cout << "Your file is empty!";
    }

    if (errorExists) {
        inputFile.close();
        //cout << "got here" << endl;
        if(errorExists) {
            cout << "Input file contains errors, so the SSA was not executed." << endl;
        } else {
            cout << "Input file successfully verified." << endl;
        }
        cout << endl;
        return false; //there's an error, so stop checking
    }

    //cout << endl;

    // 2) next lines should start with a reaction definition; number of lines like this should equal the number in the firs tline
    // formatting: should contain either -> or <->
    // molecule names are separated by spaces
    // last thing on the line is the reaction rate, which should be a real number, with up to 1 decimal
    bool notEnoughReactionsDefined = false;
    bool alreadyWarnedAboutReactionCount = false;
    for (int r = 0; r < numReactions; r++) {
        if (!notEnoughReactionsDefined) {
            //cout << "got here" << endl;
            if (inputFile.peek() != EOF) {
                getline(inputFile, fullReactionDefLine);
                lineNumber++;
                //cout << "Now checking line " << lineNumber << endl;

                //chop off the comments, if any (denoted by the first instance of #)
                fullReactionDefLine = chopOffComments(fullReactionDefLine);

                if (!fullReactionDefLine.empty()) {
                    //if there is a space at the end of the line, remove it
                    if (fullReactionDefLine.at(fullReactionDefLine.length() - 1) == ' ') {
                        fullReactionDefLine = fullReactionDefLine.substr(0, fullReactionDefLine.length() -
                                                                            1); //remove the very last character
                    }

                    // number of times "<->" or "->" shows up
                    int numForward = 0;
                    int numReversible = 0;
                    int numReactionRates = 0;

                    int j = 0;
                    //cout << "fullReactionDefLine is " << fullReactionDefLine << endl;
                    while (fullReactionDefLine.find(" ", j) !=
                           std::string::npos) { // while there is still a space in the line on or after index j
                        //cout << "entering while loop on line 188" << endl;
                        reactionSlice = "";
                        //find the next space
                        int spaceIndex = fullReactionDefLine.find(" ", j);
                        for (int i = j; i < spaceIndex; i++) {
                            reactionSlice += fullReactionDefLine.at(
                                    i); // copy the string up to the next space (and don't include that space)
                        }
                        //cout << reactionSlice << " is about to be checked" << endl;
                        if (checkReactionSlice(reactionSlice, lineNumber, errorExists)) { //if there is an error
                            errorExists = true;
                        }

                        if (isReactionRate(reactionSlice)) {
                            if (reactionSlice != "0") numReactionRates++;
                        }

                        if (reactionSlice == "<->") {
                            numReversible++;
                        }
                        if (reactionSlice == "->") {
                            numForward++;
                        }

                        //increase j until it's at the index after the space
                        j = spaceIndex + 1;
                        //cout << j << endl;
                    }
                    //now we have reached the fencepost
                    //everything after the last space should be a reaction rate (the 2nd if the line contains "<->", or the first if it contains "->")
                    reactionSlice = "";
                    for (int i = j; i < fullReactionDefLine.length(); i++) {
                        reactionSlice += fullReactionDefLine.at(i); // copy the remaining characters
                    }
                    //cout << "now checking " << reactionSlice << endl;

                    if (!isValidReactionRate(reactionSlice, lineNumber)) {
                        errorExists = true;
                    } else {
                        if (reactionSlice != "0") numReactionRates++;
                    }
                    if (numReversible + numForward != 1) {
                        cout << "Warning: line " << lineNumber << " doesn't contain exactly one <-> or ->" << endl;
                        cout << numReversible << " " << numForward << endl;
                        errorExists = true;
                    }
                    if (numReversible == 1 && numReactionRates != 2) {
                        cout << "Warning: line " << lineNumber << " contains <-> but contains " << numReactionRates
                             << " valid reaction rate(s)." << endl;
                        errorExists = true;
                    }
                    if (numForward == 1 && numReactionRates != 1) {
                        cout << "Warning: line " << lineNumber << " contains -> but contains " << numReactionRates
                             << " valid reaction rates." << endl;
                        errorExists = true;
                    }

                } else { //line is blank
                    if (!alreadyWarnedAboutReactionCount) {
                        cout << "Warning: Your file has fewer reactions defined than it should, according to line 1."
                             << endl;
                    }
                    errorExists = true;
                    notEnoughReactionsDefined = true;
                    alreadyWarnedAboutReactionCount = true;
                }
            } else {
                if (!alreadyWarnedAboutReactionCount) {
                    cout << "Warning: Your file has fewer reactions defined than it should, according to line 1."
                         << endl;
                }
                errorExists = true;
                notEnoughReactionsDefined = true;
                alreadyWarnedAboutReactionCount = true;
            }
        } else {
            if (!alreadyWarnedAboutReactionCount) {
                cout << "Warning: Your file has fewer reactions defined than it should, according to line 1."
                     << endl;
            }
        }
    } // done checking reaction definitions now

    if (errorExists) {
        inputFile.close();
        //cout << "got here" << endl;
        if(errorExists) {
            cout << "Input file contains errors, so the SSA was not executed." << endl;
        } else {
            cout << "Input file successfully verified." << endl;
        }
        cout << endl;
        return false; //there's an error, so stop checking
    }

    // 3) next lines may or may not be completely blank (excluding comments)
    // keep taking in lines until one is not blank

    while (true) {
        lineNumber++;
        //cout << "Now checking line " << lineNumber << endl;
        if (inputFile.peek() != EOF) {
            getline(inputFile, moleculeDefLine); //break if this line is not blank after removing comments
            moleculeDefLine = chopOffComments(moleculeDefLine);
            if (!moleculeDefLine.empty()) break;
        } else { // they didn't give the initial molecule populations
            cout << "Warning: Your input pile doesn't provide the initial molecule populations." << endl;
            errorExists = true;
        }
    }
    lineNumber--;

    // 4) the remaining lines should define molecule populations
    // formatted as follows: molecule name, then a space, then a number
//    if (inputFile.peek() == EOF) { // they didn't give the initial molecule populations
//        cout << "Warning: Your input pile doesn't provide the initial molecule populations." << endl;
//        errorExists = true;
//    }

//    else {
    while (inputFile.peek() != EOF) { // check all remaining lines
        lineNumber++;
        //cout << "Now checking line " << lineNumber << endl;
        if (!firstEntry) getline(inputFile, moleculeDefLine); //the first time we enter, the line has already been read
        firstEntry = false;
        moleculeDefLine = chopOffComments(moleculeDefLine);

        //if there is a space at the end of the line, remove it
        if (!moleculeDefLine.empty()) {
            if (moleculeDefLine.at(moleculeDefLine.length() - 1) == ' ') {
                moleculeDefLine = moleculeDefLine.substr(0,
                                                         moleculeDefLine.length() - 1); //remove the very last character
            }
        }

        //make sure line contains a space
        if (moleculeDefLine.find(" ") != std::string::npos) {
            //copy the line up until the space
            string moleculeNameOrCount = "";
            int j = 0;
            while (moleculeDefLine.at(j) != ' ') {
                moleculeNameOrCount += moleculeDefLine.at(j);
                j++;
            }
            j++; //j now points to the start of the number
            //moleculeName should start with a letter
            if (!isalpha(moleculeNameOrCount.at(0))) {
                cout << "Warning: Line " << lineNumber
                     << " contains a molecule that starts with a non-alphabet character." << endl;
                errorExists = true;
            }

            //now get the starting population size
            moleculeNameOrCount = "";
            while (j < moleculeDefLine.length()) {
                moleculeNameOrCount += moleculeDefLine.at(j);
                j++;
            }
            //entire number should be digits (obviously)
            for (int i = 0; i < moleculeNameOrCount.length(); i++) {
                if (!isdigit(moleculeNameOrCount.at(i))) {
                    cout << "Warning: Line " << lineNumber
                         << " contains a non-digit character in its population name or size, namely the character '" << moleculeNameOrCount.at(i) << "'." << endl;
                    //cout << "moleculeNameOrCount is '" << moleculeNameOrCount << "'." << endl;
                    errorExists = true;
                }
            }

        } else {
            cout << "Warning: Line " << lineNumber << " appears to be an initial population line but contains no space." << endl;
            errorExists = true;
        }
    }
//}

    // Close the file
    inputFile.close();
    //cout << "got here" << endl;
    if(errorExists) {
        cout << "Input file contains errors, so the SSA was not executed." << endl;
    } else {
        cout << "Input file successfully verified." << endl;
    }
    cout << endl;
    return !errorExists; // returns whether the input file is safe to use
}



//returns true if there is an error, false if there isn't
//only gets called on a slice that starts with a digit and contains a non-digit, non-decimal character
bool inputVerifier::checkSingleMolecule(string molNumAndName, int lineNumber, bool errorExists) {

    //cout << "Now checking " << molNumAndName << " which has length " << molNumAndName.length() << endl;
    bool anyErrors = false;
    string moleculeCount = "";
    int i = 0;
    if (!isdigit(molNumAndName.at(0))) { // if it doesn't start with a number, it must start with an alphabet letter
        if (!isalpha(molNumAndName.at(0))) {
            cout << "Warning: A reactant or product definition on line " << lineNumber << " doesn't start with a number or a letter.";
            anyErrors = true;
            errorExists = true;
        }
    }
//    if (!isalpha(molNumAndName.at(molNumAndName.length() - 1))) {
//        cout << "Warning: A reactant or product on line " << lineNumber << " doesn't end with a letter.";
//        anyErrors = true;
//        errorExists = true;
//    }

    //parse until first non-digit character
    while (i < molNumAndName.length()) {
        //cout << i << endl;
        if (isdigit(molNumAndName.at(i))) {
            moleculeCount += molNumAndName.at(i); //copy the digits into a new string (this string will be returned later)
            //cout << "moleculeCount is now " << moleculeCount << endl;
        }
        else break;
        if (i == molNumAndName.length()) { // if the entire line was digits
            cout << "Warning: Line " << lineNumber << " contains only numbers." << endl;
            anyErrors = true;
            errorExists = true;
            break;
        }
        i++;
    }
    //if (!anyErrors) cout << "number parsed successfully" << endl;
    //i--;
    //cout << "i is " << i << endl;

    string moleculeName = "";
    if (i < molNumAndName.length() && i > 0) { //at least one character remains unchecked so far, and the string started with a digit
        for (int k = i; k < molNumAndName.length(); k++) {
            moleculeName += molNumAndName.at(k); //copy the remaining letters into another string
        }
        //cout << "moleculeName is " << moleculeName << "." << endl;
    }

    if (!isalpha(moleculeName.at(0))) { //molecule name MUST start with an alphabet character
        cout << "Warning: A reactant or product name on line " << lineNumber << " doesn't start with a letter." << endl;
        anyErrors = true;
        errorExists = true;
    }

    //every character in moleculeName should be a letter of the alphabet (this is not true anymore)
//    for (i = 0; i < moleculeName.length(); i++) {
//        if (!isalpha(moleculeName.at(i))) {
//            cout << "Warning: Line " << lineNumber << " has an invalid molecule name." << endl;
//            anyErrors = true;
//            errorExists = true;
//        }
//        else {
//            //cout << "molecule name parsed successfully" << endl;
//        }
//    }

    //if (!anyErrors) cout << molNumAndName << " parsed succesfully" << endl;
    return errorExists;
}

//returns true if there is an error, false if there isn't
bool inputVerifier::checkReactionSlice(string reactionSlice, int lineNumber, bool errorExists) {
    //determine what type of slice this is (reactant/product, ->, <->, or reaction rate)
    //cout << "checking reaction slice " << reactionSlice << " on line " << lineNumber << endl;
    if (reactionSlice.at(0) == 0 && reactionSlice.length() == 1) {
        //if it's just the number 0, this is fine because it's an empty reactant or product
        return false;
    }
    else { //if it's not just the number 0
        if (isdigit(reactionSlice.at(0))) { //if it starts with a digit
            if (containsNonDigitNonDecimal(reactionSlice)) { //it's a (valid?) reactant/product definition
                return checkSingleMolecule(reactionSlice, lineNumber, errorExists);
                //return false;
            } else { //it's a reaction rate, so it should contain up to 1 decimal and the rest digits
                return !isValidReactionRate(reactionSlice, lineNumber); //only gets called if whole string is digits and decimals
            }
        } else { //must either be "<->" or "->"
            if (reactionSlice != "<->" && reactionSlice != "->") {
                cout << "Warning: Line " << lineNumber << " contains something that isn't '<->', '->', or a valid reactant/product, or it contains a reaction rate in the wrong spot." << endl;
                return true; //there is an error
            }
            return false;
        }
    }
}

bool inputVerifier::isReactionRate(string reactionSlice) {
    if (reactionSlice.at(0) == 0 && reactionSlice.length() == 1) {
        return false;
    } else { //if it's not just the number 0
        if (isdigit(reactionSlice.at(0))) { //if it starts with a digit
            if (containsNonDigitNonDecimal(reactionSlice)) { //it's a reactant/product definition
                return false;
            } else { //it's a reaction rate, so it should contain up to 1 decimal and the rest digits
                //cout << reactionSlice << " is a reaction rate" << endl;
                return true;
            }
        }
        else {
            return false;
        }
    }
}

bool inputVerifier::containsNonDigitNonDecimal(string reactionSlice) {
    for (int i = 0; i < reactionSlice.length(); i++) {
        if (!(isdigit(reactionSlice.at(i))) && reactionSlice.at(i) != '.') {
            //cout << reactionSlice.at(i) << " is not a digit or decimal." << endl;
            return true;
        }
    }
    return false;
}


bool inputVerifier::isValidReactionRate(string reactionSlice, int lineNumber){ // should contain up to 1 decimal, and nothing else except digits
    if (reactionSlice.find(".") != std::string::npos) { //if it contains a decimal
        int firstDecimalIndex = reactionSlice.find(".");
        if (reactionSlice.find(".", firstDecimalIndex + 1) != std::string::npos) {
            cout << "Warning: Line " << lineNumber << " contains too many decimals in the reaction rate." << endl;
            return false;
        }
    }
    //now every digit should either be a digit or decimal (I can probably take this out)
    for (int i = 0; i < reactionSlice.length(); i++) {
        if (!(isdigit(reactionSlice.at(i))) && reactionSlice.at(i) != '.') {
            cout << "Warning: Line " << lineNumber << " contains the invalid reaction rate " << reactionSlice << endl;
            return false;
        }
    }
    return true;
}

string inputVerifier::chopOffComments(string line) {
    if (line.find("#") != std::string::npos) {
        string temp = "";
        for (int i = 0; i < line.find("#"); i++) {
            temp += line.at(i);
        }
        line = temp;
    }
    return line;
}

////////////////////////////////////////////////////////////////////////////these methods aren't needed anymore
//int checkReactionDefLine(ifstream inputFile, string reactionDefLine) {
//    return 0;
//}
//
//bool InputVerifier::checkReactantsOrProducts(string reactionDefLine, int lineNumber, bool errorExists) {
//    //fencepost problem
//    string moleculeDecl = "";
//    int j = 0;
//
//    //reactionDefLine must end with a character
//    if (!isalpha(reactionDefLine.at(reactionDefLine.length() - 1))) {
//        cout << "Warning: Line " << lineNumber << " is formatted incorrectly." << endl;
//        errorExists = true;
//    }
//
//    if(!errorExists) {
//        while (reactionDefLine.find(" ", j) != std::string::npos) { // while string contains a space at or after j
//            //cout << "reactionDefLine in checkReactantsOrProducts is " << reactionDefLine << "." << endl;
//
//            while (reactionDefLine.at(j) != ' ') {
//                moleculeDecl += reactionDefLine.at(j); //copy the line up until the first space
//                j++;
//            } //when this while loop gets left, j will point to a space
//            //cout << moleculeDecl << " being checked in while loop" << endl;
//            if(checkSingleMolecule(moleculeDecl, lineNumber, errorExists) == true) {
//                errorExists = true;
//            }
//            moleculeDecl = "";
//            //cout << "got here" << endl;
//            j++; //j now points to the next character after the space, which should be a number
//            if (!isdigit(reactionDefLine.at(j))) {
//                cout << "Warning: The reactants or products in line " << lineNumber << " are formatted incorrectly."
//                     << endl;
//                errorExists = true;
//                break;
//            }
//
//            // while loop will repeat if there is still another space; otherwise, check the remaining characters
//        }
//    }
//    if(!errorExists) {
//        moleculeDecl = "";
//        //cout << j << " " << reactionDefLine.length() << endl;
//        while (j < reactionDefLine.length()) {
//            moleculeDecl += reactionDefLine.at(j);
//            j++;
//        }
//        //cout << "moleculeDecl is " << moleculeDecl << "." << endl;
//
//        if(checkSingleMolecule(moleculeDecl, lineNumber, errorExists) == true) {
//            errorExists = true;
//        } //the final molecule to check
//    }
//    return errorExists;
//}