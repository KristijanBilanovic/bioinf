#include <iostream>
#include "spoa/spoa.hpp"           
#include "bioparser/fastq_parser.hpp"
#include "MyClass.hpp"

using namespace std;

int main() {
    MyClass obj("BioProject");
    obj.greet();
    return 0;
}