#include <iostream>
#include <string>
#include <filesystem>
#include "spoa/spoa.hpp"           
#include "bioparser/fastq_parser.hpp"
#include <map>
#include <tuple>
#include "sequence_analyzer.hpp"

using namespace std;

// Parses a single FASTQ file and returns a vector of Sequences
std::vector<std::unique_ptr<Sequence>> ParseData(const std::string& filepath) {
    std::vector<std::unique_ptr<Sequence>> sequences;
    auto parser = bioparser::Parser<Sequence>::Create<bioparser::FastqParser>(filepath);
    while (true) {
        auto batch = parser->Parse(1ULL << 30);
        if (batch.empty()) break;
        sequences.insert(
            sequences.end(),
            std::make_move_iterator(batch.begin()),
            std::make_move_iterator(batch.end())
        );
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

    std::vector<std::unique_ptr<Sequence>> filtered;
    for (auto& s : seqs) {
        int len = s->data.size();

        if (len >= mode - 5 && len <= mode + 5) {
            filtered.push_back(std::move(s));
        }

    }
    return filtered;

}
// generates MSA using spoa
std::vector<std::string> generate_msa(
    const std::vector<std::unique_ptr<Sequence>>& seqs)
{
    auto engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW,
        0,   // match
        -1,  // mismatch
        -1   // gap
    );
    spoa::Graph graph{};
    for (const auto& s : seqs) {
        auto alignment = engine->Align(s->data, graph);
        graph.AddAlignment(alignment, s->data);
    }
    return graph.GenerateMultipleSequenceAlignment();
}

//Clusters aligned sequences using greedy approach and Hamming distance
std::vector<std::vector<int>> cluster_sequences(
    const std::vector<std::string>& msa,
    int k = 12) // k is the maximum Hamming distance to cluster together
{
    std::vector<std::vector<int>> clusters;
    std::vector<int> representatives; 

    for (int i = 0; i < (int)msa.size(); i++) {
        bool added = false;

        // compare current sequence with representatives of existing clusters
        for (int c = 0; c < (int)clusters.size(); c++) {
            int rep = representatives[c];

            //hamming distance of current sequence and representative
            int dist = 0;
            for (int p = 0; p < (int)msa[i].size(); p++) {
                if (msa[i][p] != msa[rep][p]) dist++;
            }

            // if close enough to representative, add to cluster and stop looking
            if (dist < k) {
                clusters[c].push_back(i);
                added = true;
                break;
            }
        }

        // if not close to any representative, create new cluster with this sequence as representative
        if (!added) {
            clusters.push_back({i});
            representatives.push_back(i);
        }
    }

    //  sort clusters by size, largest first
    std::sort(clusters.begin(), clusters.end(),
        [](const std::vector<int>& a, const std::vector<int>& b) {
            return a.size() > b.size();
        });

    return clusters;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: ./jelen_analiza <file.fastq>\n";
        return 1;
    }

    auto sequences = ParseData(argv[1]);
    cout << "Parsed: " << sequences.size() << " sequences.\n";

    auto filtered = filter_by_length(sequences);
    cout << "Filtered: " << filtered.size() << " sequences.\n";

    auto msa = generate_msa(filtered);
    cout << "MSA done: " << msa.size() << " sequences.\n";
    cout << "Aligned length: " << msa[0].size() << "\n";

    auto clusters = cluster_sequences(msa);
    cout << "Clusters found: " << clusters.size() << "\n";
    cout << "Largest cluster:  " << clusters[0].size() << " sequences\n";
    cout << "Second largest:   " << clusters[1].size() << " sequences\n";
    // TODO: uncomment later for centroid analysis
    // SequenceAnalyzer analyzer(filtered);
    // auto neighbors = analyzer.find_nearest_neighbors();

    return 0;
}