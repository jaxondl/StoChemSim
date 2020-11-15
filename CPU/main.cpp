#include "algorithm.h"
#include "decoder.h"
#include "dependency.h"
#include "tree.h"

#include <iostream>
#include <string>

using namespace std;

int main() {
    std::cout << "Hello, World!" << std::endl;

    algorithm algorithmObject;
    decoder decoderObject;
    dependency dependencyObject;
    tree treeObject;

    algorithmObject.testFunction();
    decoderObject.testFunction();
    dependencyObject.testFunction();
    treeObject.testFunction();

    return 0;
}
