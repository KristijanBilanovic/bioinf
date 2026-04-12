#include <iostream>
#include <string>
#include "spoa/spoa.hpp"           
#include "bioparser/fastq_parser.hpp"
#include "MyClass.hpp"

using namespace std;

struct Sequence{
    public:
        std::string name;
        std::string data;
        std::string quality;

        Sequence(
            const char* name, std::uint32_t name_length, 
            const char* data, std::uint32_t data_length,
            const char* quaity, std::uint32_t quality_length
        )
        {
            this->name.assign(name, name_length);
            this->data.assign(data, data_length);
            this->quality.assign(quaity, quality_length);
        }
};

int main() {
    MyClass obj("BioProject");
    obj.greet();
    return 0;
}