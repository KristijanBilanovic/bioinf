#include <iostream>
#include <string>
#include <filesystem>
#include "spoa/spoa.hpp"           
#include "bioparser/fastq_parser.hpp"
#include "bioparser/fasta_parser.hpp"
#include <map>
#include <tuple>
#include "sequence_analyzer.hpp"
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <climits>
#include <unordered_map>

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


int hamming_distance(const std::string& a,
                     const std::string& b)
{
    if (a.size() != b.size()) {
        return -1;
    }

    int dist = 0;

    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] != b[i]) {
            dist++;
        }
    }

    return dist;
}

std::pair<int, int> best_hamming_match(
    const std::string& consensus,
    const std::string& truth)
{
    if (consensus.size() < truth.size()) {
        return {-1, -1};
    }

    int best_dist = INT_MAX;
    int best_pos = -1;

    for (size_t i = 0;
         i + truth.size() <= consensus.size();
         i++)
    {
        std::string window =
            consensus.substr(i, truth.size());

        int dist = hamming_distance(window, truth);

        if (dist < best_dist) {
            best_dist = dist;
            best_pos = i;
        }
    }

    return {best_dist, best_pos};
}

// Parses a single FASTQ file and returns a vector of Sequences
std::vector<std::unique_ptr<Sequence>> ParseDataFQ(const std::string& filepath) {
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

// Parses a single FASTA file and returns a vector of Sequences
std::vector<std::unique_ptr<SequenceFA>> ParseDataFA(const std::string& filepath) {
    auto parser = bioparser::Parser<SequenceFA>::Create<bioparser::FastaParser>(filepath);
    auto s = parser->Parse(-1);
    return s;
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

// generates SPOA graph
spoa::Graph generate_spoa_graph(
    const std::vector<std::unique_ptr<Sequence>>& seqs)
{
    auto engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kOV,
        3,   // match
        -5,  // mismatch
        -3   // gap
    );
    spoa::Graph graph{};
    for (const auto& s : seqs) {
        auto alignment = engine->Align(s->data, graph);
        graph.AddAlignment(alignment, s->data);
    }
    return graph;
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

    // normalize: how far away are the sequances
    double similarity = (double)lis / min(m1.size(), m2.size());
    return 1.0 - similarity;
}


/*
    Clusters sequences based on their minimizer similarity.
    @param seqs: vector of sequences to cluster
    @return: vector of clusters, where each cluster is a vector of sequence indices
*/
std::vector<std::vector<int>> cluster(
    const std::vector<std::unique_ptr<Sequence>>& seqs)
{
    // parameters for minimizer generation provided in lecture slides
    int k = 11;
    int w = 5;
    double threshold = 0.33;
    
    std::vector<std::vector<int>> clusters;
    
    for (int i = 0; i < (int)seqs.size(); i++) {
        bool assigned = false;
        
        // Try to assign to existing cluster
        for (int c = 0; c < (int)clusters.size(); c++) {
            // Compare to multiple sequences in cluster, not just representative
            int comparisons = 0;
            int similar_count = 0;
            
            for (int idx : clusters[c]) {
                if (comparisons >= 5) break; // Limit comparisons for efficiency
                
                double distance = minimizer_distance(
                    seqs[i]->data, 
                    seqs[idx]->data, 
                    k, w
                );
                
                if (distance <= threshold) {
                    similar_count++;
                }
                comparisons++;
            }
            
            // Assign if similar to majority of compared sequences
            if (similar_count >= (comparisons / 2 + 1)) {
                clusters[c].push_back(i);
                assigned = true;
                break;
            }
        }
        
        if (!assigned) {
            clusters.push_back({i});
        }
    }

    // sort clusters by size (largest first)
    sort(clusters.begin(), clusters.end(),
        [](const vector<int>& a, const vector<int>& b) {
            return a.size() > b.size();
        });
    
    return clusters;
}


/*
    Function to generate consensus sequences for each cluster using the SPOA graph.
    @param seqs: vector of sequences
    @param clusters: vector of clusters, where each cluster is a vector of sequence indices
    @return: vector of consensus sequences, one for each cluster
*/
std::vector<std::string> get_cluster_consensus(
    const std::vector<std::unique_ptr<Sequence>>& seqs,
    const std::vector<std::vector<int>>& clusters)
{
    std::vector<std::string> consensus_sequences;

    for (const auto& cluster_indices : clusters) {
        vector<std::unique_ptr<Sequence>> cluster_seqs;

        // gather sequences for this cluster
        for (int idx : cluster_indices) {
            cluster_seqs.push_back(std::make_unique<Sequence>(
                seqs[idx]->name.c_str(), seqs[idx]->name.size(),
                seqs[idx]->data.c_str(), seqs[idx]->data.size(),
                seqs[idx]->quality.c_str(), seqs[idx]->quality.size()
            ));
        }

        auto graph = generate_spoa_graph(cluster_seqs);
        auto consensus = graph.GenerateConsensus();

        if (!consensus.empty()) {
            consensus_sequences.push_back(consensus); 
        }
    }
    return consensus_sequences;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: ./jelen_analiza <file.fastq>\n";
        return 1;
    }

    // Load ground truth data (for evaluation)
    auto ground_truth_29 = ParseDataFA("../data/J29B_expected.fasta");
    auto ground_truth_30 = ParseDataFA("../data/J30B_expected.fasta");

    // Load FASTQ data from the provided file path
    auto sequences = ParseDataFQ(argv[1]);
    cout << "Parsed: " << sequences.size() << " sequences.\n";

    // Filter sequences by length (keep those close to mode length)
    auto filtered = filter_by_length(sequences);
    cout << "Filtered: " << filtered.size() << " sequences.\n";

    // Cluster sequances before aligning
    auto clusters = cluster(filtered);
    cout << "Clusters: " << clusters.size() << "\n";

    // Generate consensus sequences for each cluster
    auto consensus_sequences = get_cluster_consensus(filtered, clusters); 

    // print consensus sequences
    for (size_t i = 0; i < consensus_sequences.size(); i++) {

        if (clusters[i].size() < 4) continue; // skip small clusters

        cout << "================== consensus_" << i << " (" << clusters[i].size() << " seqs) ================== \n";

        for (size_t j = 0; j < ground_truth_29.size(); j++) {

            auto [dist, pos] =
                best_hamming_match(
                    consensus_sequences[i],
                    ground_truth_29[j]->data
                );

            cout << "vs J29B-" << j + 1
                << " | best Hamming distance = "
                << std::setw(3) << dist
                << " | position = "
                << std::setw(3) <<pos
                << "\n";
        }

        cout << "----------------------------------------------------------\n";

        for (size_t j = 0; j < ground_truth_30.size(); j++) {
            auto [dist, pos] = best_hamming_match(
                consensus_sequences[i], ground_truth_30[j]->data);
            cout << "vs J30B-" << j + 1
                 << " | best Hamming distance = " 
                 << std::setw(3) << dist
                 << " | position = " 
                 << std::setw(3) << pos << "\n";
        }

        cout << "\n";
    }

    return 0;
}