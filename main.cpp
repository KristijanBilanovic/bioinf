#include <iostream>
#include "spoa/spoa.hpp"           
#include "bioparser/fastq_parser.hpp"
#include "MyClass.hpp"

using namespace std;

struct Sequence{
    public:
        Sequence(
            const char* name, std::uint32_t name_length, 
            const char* data, std::uint32_t data_length,
            const char* quaity, std::uint32_t quality_length
        )
        {
            // TODO
        }
};

int main() {
    MyClass obj("BioProject");
    obj.greet();
    return 0;
}