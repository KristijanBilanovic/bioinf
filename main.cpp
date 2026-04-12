#include <iostream>
#include <string>
#include <filesystem>
#include "spoa/spoa.hpp"           
#include "bioparser/fastq_parser.hpp"

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

std::vector<std::unique_ptr<Sequence>> ParseData(const std::string& path)
{
    std::vector<std::unique_ptr<Sequence>> sequences;

    // Iterate through all files in the directory and parse those starting with 'J'
    for (const auto& file : std::filesystem::directory_iterator(path)) {
        std::string filename = file.path().filename().string();

        if(filename[0] == 'J') {
            // Create a NEW parser for each file
            auto parser = bioparser::Parser<Sequence>::Create<bioparser::FastqParser>(file.path().string());
            
            while (true) {
                // Use smaller buffer: 256MB instead of 1GB
                auto batch = parser->Parse(256ULL << 20);
                if (batch.empty()) {
                    break;
                }
                sequences.insert(
                    sequences.end(),
                    std::make_move_iterator(batch.begin()),
                    std::make_move_iterator(batch.end())
                );
            }
        }
    }
    return sequences;
}

int main() {
    auto sequences = ParseData("../data/fastq/");
    cout << "Parsed " << sequences.size() << " sequences." << endl;

    cout << "First sequence name: " << sequences[0]->name << endl;

    return 0;
}