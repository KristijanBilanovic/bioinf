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

// generates SPOA graph
spoa::Graph generate_spoa_graph(
    const std::vector<std::unique_ptr<Sequence>>& seqs)
{
    auto engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW,
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

    std::vector<std::vector<int>> clusters;
    std::vector<int> reps;

    double threshold = 0.4;

    for (int i = 0; i < (int)seqs.size(); i++) {

        bool assigned = false;

        for (int c = 0; c < (int)clusters.size(); c++) {

            double distance = minimizer_distance(
               seqs[i]->data, 
               seqs[reps[c]]->data, 
               k, 
               w
            );

            if (distance <= threshold) {
                clusters[c].push_back(i);
                assigned = true;
                break;
            }
        }

        if (!assigned) {
            clusters.push_back({i});
            reps.push_back(i);
        }
    }

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
    auto ground_truth_29 = ParseData("../data/fastq/J29_B_CE_IonXpress_005.fastq");
    auto ground_truth_30 = ParseData("../data/fastq/J30_B_CE_IonXpress_006.fastq");

    // Filter ground truth sequences by length to get a clean set of reference sequences
    auto gt29_filtered = filter_by_length(ground_truth_29);
    auto gt30_filtered = filter_by_length(ground_truth_30);

    // Cluster ground truth sequences
    auto gt29_clusters = cluster(gt29_filtered);
    auto gt30_clusters = cluster(gt30_filtered);

    // Find index of the largest cluster
    int gt29_largest_cluster_idx = 0;
    int gt30_largest_cluster_idx = 0;
    for (int i = 1; i < (int)gt29_clusters.size(); i++) {
        if (gt29_clusters[i].size() > gt29_clusters[gt29_largest_cluster_idx].size()) {
            gt29_largest_cluster_idx = i;
        }
    }

    for (int i = 1; i < (int)gt30_clusters.size(); i++) {
        if (gt30_clusters[i].size() > gt30_clusters[gt30_largest_cluster_idx].size()) {
            gt30_largest_cluster_idx = i;
        }
    }

    // Generate consensus sequences for ground truth clusters
    auto gt29_cluster_consensus = get_cluster_consensus(gt29_filtered, gt29_clusters);
    auto gt30_cluster_consensus = get_cluster_consensus(gt30_filtered, gt30_clusters);

    // Use the consensus sequence of the largest cluster as the representative ground truth sequence for each file
    string gt29_consensus = gt29_cluster_consensus[gt29_largest_cluster_idx];
    string gt30_consensus = gt30_cluster_consensus[gt30_largest_cluster_idx];

    // Load FASTQ data from the provided file path
    auto sequences = ParseData(argv[1]);
    cout << "Parsed: " << sequences.size() << " sequences.\n";

    // Filter sequences by length (keep those close to mode length)
    auto filtered = filter_by_length(sequences);
    cout << "Filtered: " << filtered.size() << " sequences.\n";

    // Cluster sequances before aligning
    auto clusters = cluster(filtered);
    cout << "Clusters: " << clusters.size() << "\n";

    // Generate consensus sequences for each cluster
    auto consensus_sequences = get_cluster_consensus(filtered, clusters); 
    
    // Compare cluster consensus sequences to ground truth using Hamming distance and report results
    for (size_t i = 0; i < consensus_sequences.size(); i++) {
        const string& consensus = consensus_sequences[i];

        std::vector<std::unique_ptr<Sequence>> concensus_1_and_gt29;
        std::vector<std::unique_ptr<Sequence>> concensus_1_and_gt30;

        concensus_1_and_gt29.push_back(std::make_unique<Sequence>(
            "consensus1", strlen("consensus1"),
            consensus.c_str(), consensus.size(),
            "", 0
        ));

        concensus_1_and_gt29.push_back(std::make_unique<Sequence>(
            "gt29_consensus", strlen("gt29_consensus"),
            gt29_consensus.c_str(), gt29_consensus.size(),
            "", 0
        ));

        concensus_1_and_gt30.push_back(std::make_unique<Sequence>(
            "consensus1", strlen("consensus1"),
            consensus.c_str(), consensus.size(),
            "", 0
        ));

        concensus_1_and_gt30.push_back(std::make_unique<Sequence>(
            "gt30_consensus", strlen("gt30_consensus"),
            gt30_consensus.c_str(), gt30_consensus.size(),
            "", 0
        ));

        
        auto graph_29 = generate_spoa_graph(concensus_1_and_gt29);
        auto graph_30 = generate_spoa_graph(concensus_1_and_gt30);

        // align consensus to ground truth and generate multiple sequence alignment
        auto msa_29 = graph_29.GenerateMultipleSequenceAlignment();
        auto msa_30 = graph_30.GenerateMultipleSequenceAlignment();

        // compute Hamming distance between consensus and ground truth for both files
        int hamming_distance_29 = 0;
        int hamming_distance_30 = 0;

        const std::string& aligned_consensus_29 = msa_29[0];
        const std::string& aligned_gt29 = msa_29[1];

        for (size_t j = 0; j < aligned_consensus_29.size(); j++) {
            if (aligned_consensus_29[j] != aligned_gt29[j]) {
                hamming_distance_29++;
            }
        }

        const std::string& aligned_consensus_30 = msa_30[0];
        const std::string& aligned_gt30 = msa_30[1];

        for (size_t j = 0; j < aligned_consensus_30.size(); j++) {
            if (aligned_consensus_30[j] != aligned_gt30[j]) {
                hamming_distance_30++;
            }
        }

        std::cout << "Cluster " << i << ": Hamming distance to GT29 = " << hamming_distance_29
                  << ", Hamming distance to GT30 = " << hamming_distance_30 << std::endl;
    }


    return 0;
}