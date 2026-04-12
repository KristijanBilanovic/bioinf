#include <iostream>
#include <string>
#include <filesystem>
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
    auto parser = bioparser::Parser<Sequence>::Create<bioparser::FastqParser>("../data/fastq/J13_L_CE_IonXpress_033.fastq");

    std::vector<std::unique_ptr<Sequence>> sequences;
    while (true) {
        auto batch = parser->Parse(1ULL << 30);
        if (batch.empty()) {
            break;
        }
        sequences.insert(
            sequences.end(),
            std::make_move_iterator(batch.begin()),
            std::make_move_iterator(batch.end())
        );
    }

    for (const auto& file : std::filesystem::directory_iterator("../data/fastq/")) {
        std::string filename = file.path().filename().string();
        std::cout << filename << std::endl;
    }
}