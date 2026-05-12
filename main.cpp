#include <iostream>
#include <string>
#include <filesystem>
#include "spoa/spoa.hpp"           
#include "bioparser/fastq_parser.hpp"
#include <map>
#include <tuple>
#include "sequence_analyzer.hpp"
#include <vector>
#include <unordered_map>
#include <cmath>
#include <bits/stdc++.h>

/*
    Structure to represent a minimizer (k-mer and its position in the sequence).
*/
struct Minimizer {
    std::string kmer;
    int pos; // position in sequence
};

/*
    Structure to represent a match between two minimizers.
    @property kmer: the k-mer that matches
    @property pos1: position of the k-mer in the first sequence
    @property pos2: position of the k-mer in the second sequence
*/
struct Match {
    std::string kmer;
    int pos1;
    int pos2;
};

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

/*
    Function to generate k-mers from a sequence.
    @param seq: input sequence
    @param k: length of k-mer
    @return: vector of k-mers
*/
vector<string> get_kmers(const string& seq, int k) {
    vector<string> kmers;
    for (int i = 0; i + k <= (int)seq.size(); i++) {
        kmers.push_back(seq.substr(i, k));
    }
    return kmers;
}

/*
    Function to generate minimizers from a sequence.
    @param seq: input sequence
    @param k: length of k-mer
    @param w: window size for minimizer selection
    @return: vector of minimizers (k-mer, pos1, pos2)
*/
std::vector<Minimizer> get_minimizers(const std::string& seq, int k, int w) {
    std::vector<std::string> kmers = get_kmers(seq, k);
    std::vector<Minimizer> mins;

    int n = kmers.size();

    for (int i = 0; i + w <= n; i++) {
        std::string best = kmers[i];
        int best_pos = i;

        // find lexicographically smallest k-mer in the window
        for (int j = 1; j < w; j++) {
            if (kmers[i + j] < best) {
                best = kmers[i + j];
                best_pos = i + j;
            }
        }

        // avoid duplicates - only add if different from last minimizer
        if (mins.empty() || mins.back().kmer != best) {
            mins.push_back({best, best_pos});
        }
    }

    return mins;
}

/*
    Function to find matches between two sets of minimizers.
    @param a: first set of minimizers
    @param b: second set of minimizers
    @return: vector of matches (k-mer, pos in a, pos in b)
*/
vector<Match> get_matches(
    const vector<Minimizer>& a,
    const vector<Minimizer>& b)
{
    unordered_map<string, vector<int>> index;

    // STEP 1: index seq2 minimizers
    for (auto &m : b) {
        index[m.kmer].push_back(m.pos);
    }

    vector<Match> matches;
    matches.reserve(a.size());

    // STEP 2: stream seq1 and lookup
    for (auto &m : a) {
        auto it = index.find(m.kmer);
        if (it == index.end()) continue;

        for (int pos2 : it->second) {
            matches.push_back({m.kmer, m.pos, pos2});
        }
    }

    return matches;
}

/*
    Function to compute the length of the longest increasing subsequence (LIS) of matches.
    This is used to estimate the length of the longest common subsequence (LCS) between two sequences.
    @param seq: vector of positions of matches in one sequence
    @return: length of LIS
*/
int LIS(const std::vector<int>& seq) {
    std::vector<int> dp;

    for (int x : seq) {
        auto it = std::lower_bound(dp.begin(), dp.end(), x);
        if (it == dp.end()) dp.push_back(x);
        else *it = x;
    }

    return dp.size();
}

/*
    Function to compute a distance between two sequences based on their minimizer matches.
    The distance is defined as 1 - (LIS of match positions / number of matches), which estimates how well the sequences align.
    @param s1: first sequence
    @param s2: second sequence
    @param k: length of k-mer for minimizer generation
    @param w: window size for minimizer selection
    @return: distance between 0 and 1, where 0 means identical and 1 means completely different
*/
double minimizer_distance(const string& s1, const string& s2,
                         int k, int w)
{
    auto m1 = get_minimizers(s1, k, w);
    auto m2 = get_minimizers(s2, k, w);

    auto matches = get_matches(m1, m2);

    if (matches.empty()) return 1.0;

    // sort by position in first sequence or second if first is equal
    sort(matches.begin(), matches.end(),
     [](const Match& a, const Match& b) {
         if (a.pos1 == b.pos1)
             return a.pos2 < b.pos2;
         return a.pos1 < b.pos1;
     });

    // extract second positions
    vector<int> seq;
    seq.reserve(matches.size());

    for (auto& m : matches)
        seq.push_back(m.pos2);

    int lis = LIS(seq);

    // normalize → similarity → distance
    double similarity = (double)lis / matches.size();
    return 1.0 - similarity;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: ./jelen_analiza <file.fastq>\n";
        return 1;
    }

    // Step 1: Load ground truth data (for evaluation)
    auto ground_truth_29 = ParseData("../data/fastq/J29_B_CE_IonXpress_005.fastq");
    auto ground_truth_30 = ParseData("../data/fastq/J30_B_CE_IonXpress_006.fastq");

    // Step 2: Parse FASTQ file (sample)
    auto sequences = ParseData(argv[1]);
    cout << "Parsed: " << sequences.size() << " sequences.\n";

    // Step 3: Filter sequences by length (keep those close to mode length)
    auto filtered = filter_by_length(sequences);
    cout << "Filtered: " << filtered.size() << " sequences.\n";

    // Step 4: Cluster sequances before aligning

    // Step 5: Generate MSA using spoa for each cluster and gez its consensus sequence, 
    
    // Step 6: Compare cluster consensus sequences to ground truth using Hamming distance and report results

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