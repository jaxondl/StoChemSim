//
// Created by user on 10/30/2020.
//

#include "InputVerifier.h"
#include <iostream>
#include <fstream>



bool InputVerifier::verifyFile(string iFile) {
    cout << "Beginning input verification." << endl;
    string fullReactionDefLine;
    string reactionDefLine = "";
    string reactionRate = "";
    string moleculeDefLine;
    int lineNumber = 0; //used for warning messages
    int numReactions;

    bool errorExists = false;

    ifstream inputFile;
    inputFile.open(iFile);
    if (!inputFile) {
        cerr << "Unable to open file";
        exit(1);   // call system to stop
    }

    // things to verify:
    // 1) first line contains a single number and nothing else
    if(inputFile.peek()!=EOF) {
        //out << "got here" << endl;
        getline (inputFile, fullReactionDefLine); //fullReactionDefLine is used here instead of creating a different string that would only be used once
        bool valid = true;
        lineNumber++;
        //cout << "Now checking line " << lineNumber << endl;
        //cout << reactionDefLine.length() << endl;
        for (int i = 0; i < fullReactionDefLine.length(); i++) {
            if (!(isdigit(fullReactionDefLine.at(i)))) {
                //first line contains a character that is not a number!
                cout << "Warning: Your first line contains a character that isn't a number." << endl;
                valid = false;
                errorExists = true;
            }
        }
        //cout << (reactionDefLine.at(reactionDefLine.length() - 1)) << endl;
        if (valid) {
            //cout << "number of reactions is: " << fullReactionDefLine << endl;
            numReactions = stoi(fullReactionDefLine);
            //cout << numReactions << endl;
        }
    }
    else {
        cout << "Your file is empty!";
    }

    //cout << endl;

    // 2) next lines should start with a reaction definition; number of lines like this should equal the number in the firs tline
    // formatting: digits then letters, then a space, then digits, then letters, then a space, etc., then "->"
    // then digits then letters, then a space, then digits, then letters, then a space, etc., then ", k=", then a number
    for (int r = 0; r < numReactions; r++) {
        //cout << "got here" << endl;
        if (inputFile.peek() != EOF) {
            getline(inputFile,
                    fullReactionDefLine); //reactionDefLine is used here instead of creating a different string that would only be used once
            lineNumber++;
            //cout << "Now checking line " << lineNumber << endl;
            if (fullReactionDefLine.find("-") != std::string::npos) {
                //ignore everything everything from "-" onward, and store the rest in reactionDefLine
                int hyphenIndex = fullReactionDefLine.find("-");
                //cout << "fullReactionDefLine is " << fullReactionDefLine << endl;
                for (int i = 0; i < hyphenIndex - 1; i++) { //-1 is so we ignore the space before the hypen
                    reactionDefLine += fullReactionDefLine.at(
                            i); // copy the string up to the space before the hyphen (and don't include that space)
                }
                //cout << "reactionDefLine is " << reactionDefLine << "." << endl;
                if (reactionDefLine.find("-", hyphenIndex + 1) !=
                    std::string::npos) { //if string contains more than one hyphen
                    cout << "Warning: Line " << lineNumber << " contains too many hyphens." << endl;
                    errorExists = true;
                }

                // check the reactants for proper formatting
                errorExists = checkReactantsOrProducts(reactionDefLine, lineNumber, errorExists);

                // now check the products, which start after "-> "
                reactionDefLine = "";
                int commaIndex = fullReactionDefLine.find(",");
                // hyphenIndex must be less than commaIndex
                if (hyphenIndex >= commaIndex) {
                    cout << "Warning: Line " << lineNumber << " has a comma in the wrong spot." << endl;
                } else {
                    for (int i = hyphenIndex + 3; i < commaIndex; i++) {
                        reactionDefLine += fullReactionDefLine.at(i);
                    }
                    //cout << "reactionDefLine is now " << reactionDefLine << "." << endl;
                    errorExists = checkReactantsOrProducts(reactionDefLine, lineNumber, errorExists);

                    if (fullReactionDefLine.find(",", commaIndex + 1) !=
                        std::string::npos) { //if string contains more than one comma
                        cout << "Warning: Line " << lineNumber << " contains too many commas." << endl;
                        errorExists = true;
                    }
                    if (fullReactionDefLine.at(commaIndex + 1) != ' ' ||
                        fullReactionDefLine.at(commaIndex + 2) != 'k' ||
                        fullReactionDefLine.at(commaIndex + 3) != '=') {
                        cout << "Warning: Line " << lineNumber << " is incorrectly formatted after the comma." << endl;
                        errorExists = true;
                    }

                    for (int i = commaIndex + 4; i < fullReactionDefLine.length(); i++) {
                        reactionRate = +fullReactionDefLine.at(i); //copy the remaining characters from the full input line
                    }
                    // make sure every character of reactionRate is a digit
                    for (int i = 0; i < reactionRate.length(); i++) {
                        if (!isdigit(reactionRate.at(i))) {
                            cout << "Warning: Line " << lineNumber << " has an invalid reaction rate listed." << endl;
                            errorExists = true;
                        }
                    }
                }
                //cout << endl;
                // 1 reaction line has now been completely checked
                reactionDefLine = "";
            } else {
                if(!errorExists) {
                    cout << "Warning: Your file has fewer reactions defined than it should, according to line 1."
                         << endl;
                }
                errorExists = true;
            }

        } else {
            if(!errorExists) {
                cout << "Warning: Your file has fewer reactions defined than it should, according to line 1." << endl;
            }
            errorExists = true;
        }
    } // done checking reaction definitions now

    // 3) next line should be completely blank; if it not, they likely define too many reactions
    lineNumber++;
    //cout << "Now checking line " << lineNumber << endl;
    if (inputFile.peek() != EOF) {
        getline (inputFile, reactionDefLine); //this line should be blank
        if (reactionDefLine != "") {
            cout << "Warning: Line " << lineNumber << " should be blank but isn't." << endl;
            errorExists = true;
        }
    }
    else { // they didn't give the initial molecule populations
        cout << "Warning: Your input pile doesn't provide the initial molecule populations." << endl;
        errorExists = true;
    }

    // 4) the remaining lines should define molecule populations
    // formatted as follows: molecule name, then a space, then a number
    if (inputFile.peek() == EOF) { // they didn't give the initial molecule populations
        cout << "Warning: Your input pile doesn't provide the initial molecule populations." << endl;
        errorExists = true;
    }

    else {
        while (inputFile.peek() != EOF) { // check all remaining lines
            lineNumber++;
            //cout << "Now checking line " << lineNumber << endl;
            getline(inputFile, moleculeDefLine);
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
                //entire moleculeName should be letters
                for (int i = 0; i < moleculeNameOrCount.length(); i++) {
                    if (!isalpha(moleculeNameOrCount.at(i))) {
                        cout << "Warning: Line " << lineNumber
                             << " contains a non-alphabet character in its molecule name." << endl;
                        errorExists = true;
                    }
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
                             << " contains a non-digit character in its population name or size." << endl;
                        errorExists = true;
                    }
                }

            } else {
                cout << "Warning: Line " << lineNumber << " contains no space." << endl;
                errorExists = true;
            }
        }
    }

    // Close the file
    inputFile.close();

    if(errorExists) {
        cout << "Input file contains errors, so the SSA was not executed." << endl;
    } else {
        cout << "Input file successfully verified." << endl;
    }
    cout << endl;
    return !errorExists; // returns whether the input file is safe to use
}

int checkReactionDefLine(ifstream inputFile, string reactionDefLine) {
    return 0;
}

bool InputVerifier::checkReactantsOrProducts(string reactionDefLine, int lineNumber, bool errorExists) {
    //fencepost problem
    string moleculeDecl = "";
    int j = 0;

    //reactionDefLine must end with a character
    if (!isalpha(reactionDefLine.at(reactionDefLine.length() - 1))) {
        cout << "Warning: Line " << lineNumber << " is formatted incorrectly." << endl;
        errorExists = true;
    }

    if(!errorExists) {
        while (reactionDefLine.find(" ", j) != std::string::npos) { // while string contains a space at or after j
            //cout << "reactionDefLine in checkReactantsOrProducts is " << reactionDefLine << "." << endl;

            while (reactionDefLine.at(j) != ' ') {
                moleculeDecl += reactionDefLine.at(j); //copy the line up until the first space
                j++;
            } //when this while loop gets left, j will point to a space
            //cout << moleculeDecl << " being checked in while loop" << endl;
            errorExists = checkSingleMolecule(moleculeDecl, lineNumber, errorExists);
            moleculeDecl = "";
            //cout << "got here" << endl;
            j++; //j now points to the next character after the space, which should be a number
            if (!isdigit(reactionDefLine.at(j))) {
                cout << "Warning: The reactants or products in line " << lineNumber << " are formatted incorrectly."
                     << endl;
                errorExists = true;
                break;
            }

            // while loop will repeat if there is still another space; otherwise, check the remaining characters
        }
    }
    if(!errorExists) {
        moleculeDecl = "";
        //cout << j << " " << reactionDefLine.length() << endl;
        while (j < reactionDefLine.length()) {
            moleculeDecl += reactionDefLine.at(j);
            j++;
        }
        //cout << "moleculeDecl is " << moleculeDecl << "." << endl;
        errorExists = checkSingleMolecule(moleculeDecl, lineNumber, errorExists); //the final molecule to check
    }
    return errorExists;
}


bool InputVerifier::checkSingleMolecule(string molNumAndName, int lineNumber, bool errorExists) {
//parse until first non-digit character
    //cout << "Now checking " << molNumAndName << " which has length " << molNumAndName.length() << endl;
    bool anyErrors = false;
    string moleculeCount = "";
    int i = 0;
    if (!isdigit(molNumAndName.at(i))) {
        cout << "Warning: A reactant or product on line " << lineNumber << " doesn't start with a number.";
        anyErrors = true;
        errorExists = true;
    }
    if (!isalpha(molNumAndName.at(molNumAndName.length() - 1))) {
        cout << "Warning: A reactant or product on line " << lineNumber << " doesn't end with a letter.";
        anyErrors = true;
        errorExists = true;
    }
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

    //every character in moleculeName should be a letter of the alphabet
    for (i = 0; i < moleculeName.length(); i++) {
        if (!isalpha(moleculeName.at(i))) {
            cout << "Warning: Line " << lineNumber << " has an invalid molecule name." << endl;
            anyErrors = true;
            errorExists = true;
        }
        else {
            //cout << "molecule name parsed successfully" << endl;
        }
    }

    //if (!anyErrors) cout << molNumAndName << " parsed succesfully" << endl;
    return errorExists;
}