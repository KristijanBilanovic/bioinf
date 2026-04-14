#include <iostream>
#include <string>
#include <filesystem>
#include "spoa/spoa.hpp"           
#include "bioparser/fastq_parser.hpp"
#include <map>
#include <tuple>

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
        }
    }
    return sequences;
}

// finding most common length
int mode_length(const std::vector<std::unique_ptr<Sequence>>& seqs) {
    std::map<int, int> freq;

    for (const auto& s : seqs) {
        freq[s->data.size()]++;
    }

    int best_length = 0;
    int best_count = 0;

    for (const auto& [length, count] : freq) {
        if (count > best_count) {
            best_count = count;
            best_length = length;
        }
    }
    return best_length;

}
// filtering by most common length 
// keep ones that are within [mode - 5, mode + 5]
std::vector<std::unique_ptr<Sequence>> filter_by_length(std::vector<std::unique_ptr<Sequence>>& seqs) {
    int mode = mode_length(seqs);
    std::cout << "Most common length: " << mode << endl;

    std::vector<std::unique_ptr<Sequence>> filtered;
    for (auto& s : seqs) {
        int len = s->data.size();

        if (len >= mode - 5 && len <= mode + 5) {
            filtered.push_back(std::move(s));
        }

    }
    std::cout << "Filtered down to " << filtered.size() << " sequences." << endl;
    return filtered;

}

int main() {
    auto sequences = ParseData("../data/fastq/");
    cout << "Parsed " << sequences.size() << " sequences." << endl;

    cout << "First sequence name: " << sequences[0]->name << endl;

    return 0;
}