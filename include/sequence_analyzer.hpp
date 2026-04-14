#pragma once
#include <map>
#include <tuple>
#include <vector>
#include <climits>
#include <iostream>
#include <memory>

struct Sequence;

class SequenceAnalyzer {
public:
    explicit SequenceAnalyzer(const std::vector<std::unique_ptr<Sequence>>& sequences)
        : seqs(sequences) {}

    std::map<int, std::tuple<int, int>> find_nearest_neighbors();

private:
    const std::vector<std::unique_ptr<Sequence>>& seqs;

    int edit_distance(const std::string& a, const std::string& b);
};