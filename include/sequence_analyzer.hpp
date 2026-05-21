#pragma once
#include <map>
#include <tuple>
#include <cstdint>
#include <vector>
#include <climits>
#include <iostream>
#include <memory>
#include <algorithm>

/*
    Sequence structure to hold the name, data, and quality of a sequence read from a FASTQ file.
    The constructor takes raw character pointers and their lengths to initialize the string members.
    @property name: the identifier of the sequence read
    @property data: the nucleotide sequence
    @property quality: the quality scores associated with the sequence read
*/
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

struct SequenceFA { 
 public:
    std::string name;
    std::string data;
    int name_length;
    int data_length;

    SequenceFA(  
        const char* name, std::uint32_t name_length,
        const char* data, std::uint32_t data_length
    ) 
    {
        this->name.assign(name, name_length);
        this->data.assign(data, data_length);
        this->name_length = name_length;
        this->data_length = data_length;   
    }
};

class SequenceAnalyzer {
public:
    explicit SequenceAnalyzer(const std::vector<std::unique_ptr<Sequence>>& sequences)
        : seqs(sequences) {}

    std::map<int, std::tuple<int, int>> find_nearest_neighbors();

private:
    const std::vector<std::unique_ptr<Sequence>>& seqs;

    int edit_distance(const std::string& a, const std::string& b);
};