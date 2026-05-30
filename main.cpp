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
    Clusters sequences using greedy approach with centroid-based distance.
    Each cluster is represented by its first sequence (centroid).
    New sequence is compared only to the centroid of each existing cluster.
    If distance to centroid is within threshold, sequence is added to that cluster.
    Otherwise, a new cluster is created with this sequence as centroid.
    @param seqs: vector of sequences to cluster
    @return: vector of clusters, where each cluster is a vector of sequence indices
*/
std::vector<std::vector<int>> cluster(
    const std::vector<std::unique_ptr<Sequence>>& seqs)
{
    // parameters from lecture slides
    int k = 11;  // k-mer length for minimizer generation
    int w = 5;   // window size for minimizer selection
    double threshold = 0.33; // maximum distance to join a cluster

    std::vector<std::vector<int>> clusters;
    std::vector<int> centroids; // index of centroid (first sequence) of each cluster

    for (int i = 0; i < (int)seqs.size(); i++) {
        bool assigned = false;

        // compare current sequence to centroid of each existing cluster
        for (int c = 0; c < (int)clusters.size(); c++) {

            // compute minimizer distance between current sequence and centroid
            double distance = minimizer_distance(
                seqs[i]->data,            // current sequence
                seqs[centroids[c]]->data, // centroid of cluster c
                k, w
            );

            // if close enough to centroid, add to this cluster and stop
            if (distance <= threshold) {
                clusters[c].push_back(i);
                assigned = true;
                break;
            }
        }

        // if not close to any centroid, create new cluster
        // current sequence becomes the centroid of the new cluster
        if (!assigned) {
            clusters.push_back({i});
            centroids.push_back(i);
        }
    }

    // sort clusters by size, largest first
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

// generates MSA using spoa
std::vector<std::string> generate_msa(
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
    return graph.GenerateMultipleSequenceAlignment();
}

/*
    Function to save cluster membership information.
    Records which original sequences belong to each cluster.
    @param seqs: vector of sequences
    @param clusters: vector of clusters, where each cluster is a vector of sequence indices
    @param output_filepath: path to output file
*/
void save_cluster_members(
    const std::vector<std::unique_ptr<Sequence>>& seqs,
    const std::vector<std::vector<int>>& clusters,
    const std::string& output_filepath)
{
    std::ofstream outfile(output_filepath);
    
    for (size_t c = 0; c < clusters.size(); c++) {
        const auto& cluster = clusters[c];
        
        // Write cluster header
        outfile << "Cluster_" << c << " (size: " << cluster.size() << "):\n";
        
        // Write each sequence name in the cluster
        for (int idx : cluster) {
            outfile << "  " << seqs[idx]->name << "\n";
        }
        
        outfile << "\n";
    }
    
    outfile.close();
}

/*
    Function to save consensus sequences to a FASTA file.
    @param consensus_sequences: vector of consensus sequences
    @param clusters: vector of clusters (for sizing)
    @param output_filepath: path to output file
*/
void save_consensus_fastq(
    const std::vector<std::string>& consensus_sequences,
    const std::vector<std::vector<int>>& clusters,
    const std::string& output_filepath)
{
    std::ofstream outfile(output_filepath);
    
    for (size_t i = 0; i < consensus_sequences.size(); i++) {
        // Create sequence name based on cluster size
        std::string seq_name = "cluster_size_" + std::to_string(clusters[i].size());
        
        // Write FASTA format: >name, sequence
        outfile << ">" << seq_name << "\n";
        outfile << consensus_sequences[i] << "\n";
    }
    
    outfile.close();
}

/*
    Function to load consensus sequences from an existing FASTA file.
    Returns pairs of (sequence, cluster_size) parsed from the sequence name.
    @param filepath: path to the FASTA file
    @return: vector of pairs (sequence, cluster_size)
*/
std::vector<std::pair<std::string, int>> load_consensus_fastq(const std::string& filepath)
{
    std::vector<std::pair<std::string, int>> consensus_data;
    std::ifstream infile(filepath);
    std::string line;
    
    while (std::getline(infile, line)) {
        if (line[0] == '>') {
            // Parse cluster size from header
            // Format: >cluster_size_<SIZE>
            std::string size_str = line.substr(14); // Skip ">cluster_size_"
            int cluster_size = std::stoi(size_str);
            
            // Read sequence
            std::getline(infile, line);
            std::string sequence = line;
            
            consensus_data.push_back({sequence, cluster_size});
        }
    }
    
    infile.close();
    return consensus_data;
}

/*
    Structure to store consensus sequence metadata for cross-sample analysis.
*/
struct ConsensusMetadata {
    std::string sequence;
    std::string sample_name;
    int cluster_size;
};

/*
    Function to align two sequences using SPOA and compute Hamming distance.
    @param seq1: first sequence
    @param seq2: second sequence
    @return: Hamming distance between aligned sequences
*/
int hamming_distance_aligned(const std::string& seq1, const std::string& seq2)
{
    // If sequences are already the same length, reuse your existing hamming_distance function!
    if (seq1.size() == seq2.size()) {
        return hamming_distance(seq1, seq2);
    }
    
    // Wrap strings into temporary Sequence structures to reuse generate_spoa_graph
    std::vector<std::unique_ptr<Sequence>> pair_seqs;
    pair_seqs.push_back(std::make_unique<Sequence>("s1", 2, seq1.c_str(), seq1.size(), string(seq1.size(), 'A').c_str(), seq1.size()));
    pair_seqs.push_back(std::make_unique<Sequence>("s2", 2, seq2.c_str(), seq2.size(), string(seq2.size(), 'A').c_str(), seq2.size()));
    
    // Reuse your graph generation logic
    auto graph = generate_spoa_graph(pair_seqs);
    auto msa = graph.GenerateMultipleSequenceAlignment();
    
    if (msa.size() < 2) {
        return std::max(seq1.size(), seq2.size());
    }
    
    // Calculate Hamming distance directly on the aligned rows
    int dist = 0;
    for (size_t i = 0; i < msa[0].size(); i++) {
        if (msa[0][i] != msa[1][i]) {
            dist++;
        }
    }
    return dist;
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
    // Load ground truth data (for evaluation)
    auto ground_truth_29 = ParseDataFA("../data/J29B_expected.fasta");
    auto ground_truth_30 = ParseDataFA("../data/J30B_expected.fasta");

    // Ask user if they want to analyze one file or all files
    cout << "\n========== File Selection ==========\n";
    cout << "1) Analyze a single file\n";
    cout << "2) Analyze all files in ../data/fastq/\n";
    cout << "Enter choice (1 or 2): ";
    
    int file_choice;
    cin >> file_choice;
    cin.ignore(); // ignore newline after cin
    
    if (file_choice != 1 && file_choice != 2) {
        cout << "Invalid choice. Using single file analysis.\n";
        file_choice = 1;
    }

    // Ask user which pipeline to use
    cout << "\n========== Pipeline Selection ==========\n";
    cout << "1) Standard (parse -> filter -> msa -> cluster -> consensus)\n";
    cout << "2) With Minimizers (parse -> filter -> cluster -> consensus)\n";
    cout << "Enter choice (1 or 2): ";
    
    int pipeline_choice;
    cin >> pipeline_choice;
    cin.ignore(); // ignore newline after cin
    
    if (pipeline_choice != 1 && pipeline_choice != 2) {
        cout << "Invalid choice. Using standard pipeline.\n";
        pipeline_choice = 1;
    }
    cout << "========================================\n\n";

    // Create clusters output directory if it doesn't exist
    std::filesystem::path clusters_dir("../data/clusters");
    if (!std::filesystem::exists(clusters_dir)) {
        std::filesystem::create_directories(clusters_dir);
        cout << "Created directory: ../data/clusters/\n\n";
    }

    // Create cluster_members output directory if it doesn't exist
    std::filesystem::path cluster_members_dir("../data/cluster_members");
    if (!std::filesystem::exists(cluster_members_dir)) {
        std::filesystem::create_directories(cluster_members_dir);
        cout << "Created directory: ../data/cluster_members/\n\n";
    }

    // Determine which files to process
    std::vector<std::string> files_to_process;
    
    if (file_choice == 1) {
        cout << "Enter FASTQ filename (e.g., J30_B_CE_IonXpress_006.fastq): ";
        std::string filename;
        std::getline(cin, filename);

        if (filename.size() < 6 || filename.substr(filename.size() - 6) != ".fastq") {
            filename += ".fastq";
        }

        files_to_process.push_back("../data/fastq/" + filename);
    } else {
        // Collect all .fastq files from ../data/fastq/ that start with 'J'
        std::filesystem::path fastq_dir("../data/fastq");
        if (std::filesystem::exists(fastq_dir)) {
            for (const auto& entry : std::filesystem::directory_iterator(fastq_dir)) {
                auto fname = entry.path().filename().string();
                if (entry.path().extension() == ".fastq" && !fname.empty() && fname[0] == 'J') {
                    files_to_process.push_back(entry.path().string());
                }
            }
            std::sort(files_to_process.begin(), files_to_process.end());
            cout << "Found " << files_to_process.size() << " FASTQ files to process (starting with 'J').\n\n";
        } else {
            cout << "Error: ../data/fastq/ directory not found.\n";
            return 1;
        }
    }

    // Store consensus metadata for all samples
    std::vector<ConsensusMetadata> all_consensuses;

    // Process each file
    for (const auto& filepath : files_to_process) {
        cout << "\n========== Processing: " << std::filesystem::path(filepath).filename().string() << " ==========\n";

        // Get sequence name and output file path
        std::string seq_name = std::filesystem::path(filepath).stem().string();
        std::string output_file = "../data/clusters/" + seq_name + "_CLUSTERS.fasta";

        // Check if output file already exists
        if (std::filesystem::exists(output_file)) {
            cout << "Output file already exists: " << output_file << "\n";
            cout << "Loading consensus sequences from: " << output_file << "\n";
            
            auto loaded_consensuses = load_consensus_fastq(output_file);
            for (const auto& [seq, size] : loaded_consensuses) {
                all_consensuses.push_back({seq, seq_name, size});
            }
            cout << "Loaded " << loaded_consensuses.size() << " consensus sequences.\n";
            continue;
        }

        // Load FASTQ data
        auto sequences = ParseDataFQ(filepath);
        cout << "Parsed: " << sequences.size() << " sequences.\n";

        // Filter sequences by length
        auto filtered = filter_by_length(sequences);
        cout << "Filtered: " << filtered.size() << " sequences.\n";

        std::vector<std::string> consensus_sequences;
        std::vector<std::vector<int>> clusters;

        if (pipeline_choice == 1) {
            // Standard pipeline
            cout << "Using STANDARD pipeline (MSA + clustering)...\n";
            auto msa = generate_msa(filtered);

            clusters = cluster_sequences(msa);
            cout << "Clusters: " << clusters.size() << "\n";

            consensus_sequences = get_cluster_consensus(filtered, clusters);
        } else {
            // Minimizers pipeline
            cout << "Using MINIMIZERS pipeline...\n";
            
            clusters = cluster(filtered);
            cout << "Clusters: " << clusters.size() << "\n";

            consensus_sequences = get_cluster_consensus(filtered, clusters);
        }

        // Save consensus sequences to FASTQ file
        save_consensus_fastq(consensus_sequences, clusters, output_file);
        cout << "Saved consensus sequences to: " << output_file << "\n";

        // Save cluster membership information
        std::string members_file = "../data/cluster_members/" + seq_name + "_CLUSTER_MEMBERS.txt";
        save_cluster_members(filtered, clusters, members_file);
        cout << "Saved cluster members to: " << members_file << "\n";

        // Store metadata for cross-sample analysis
        for (size_t i = 0; i < consensus_sequences.size(); i++) {
            all_consensuses.push_back({consensus_sequences[i], seq_name, (int)clusters[i].size()});
        }

        cout << "Finished with file: " << seq_name << "\n";
    }

    // Single file analysis: compare to ground truth
    if (file_choice == 1) {
        cout << "\n========== Ground Truth Comparison ==========\n";
        cout << "Which ground truth to compare against?\n";
        cout << "1) J29B\n";
        cout << "2) J30B\n";
        cout << "Enter choice (1 or 2): ";
        
        int truth_choice;
        cin >> truth_choice;
        cin.ignore();
        
        auto& ground_truth = (truth_choice == 1) ? ground_truth_29 : ground_truth_30;
        std::string truth_name = (truth_choice == 1) ? "J29B" : "J30B";
        
        cout << "\n========== Results: " << truth_name << " Comparison ==========\n\n";
        
        for (size_t i = 0; i < all_consensuses.size(); i++) {
            if (all_consensuses[i].cluster_size < 4) continue; // skip small clusters

            cout << "================== " << all_consensuses[i].sample_name 
                 << " - Cluster size: " << all_consensuses[i].cluster_size 
                 << " ================== \n";

            for (size_t j = 0; j < ground_truth.size(); j++) {
                auto [dist, pos] = best_hamming_match(all_consensuses[i].sequence,
                    ground_truth[j]->data);

                cout << "vs " << truth_name << "-" << j + 1
                    << " | best Hamming distance = "
                    << std::setw(3) << dist
                    << " | position = "
                    << std::setw(3) << pos << "\n";
            }

            cout << "\n";
        }
    } 
    // Multi-file analysis: cross-sample validation
    else {
        // Open output file for cross-sample analysis
        std::string analysis_output = "../data/clusters/cross_sample_analysis.txt";
        std::ofstream analysis_file(analysis_output);
        
        cout << "\n========== Cross-Sample Analysis ==========\n";
        analysis_file << "========== Cross-Sample Analysis ==========\n";
        cout << "Comparing small clusters against large clusters for validation...\n\n";
        analysis_file << "Comparing small clusters against large clusters for validation...\n\n";
        
        // Separate consensuses into large and small clusters
        const int LARGE_CLUSTER_THRESHOLD = 100; // clusters >= 100 reads are considered "large"
        const int HAMMING_THRESHOLD = 5; // consensuses within this distance are considered similar
        
        std::vector<size_t> large_cluster_indices;
        std::vector<size_t> small_cluster_indices;
        
        for (size_t i = 0; i < all_consensuses.size(); i++) {
            if (all_consensuses[i].cluster_size >= LARGE_CLUSTER_THRESHOLD) {
                large_cluster_indices.push_back(i);
            } else {
                small_cluster_indices.push_back(i);
            }
        }
        
        cout << "Found " << large_cluster_indices.size() << " large cluster(s) (>= " << LARGE_CLUSTER_THRESHOLD << " reads)\n";
        analysis_file << "Found " << large_cluster_indices.size() << " large cluster(s) (>= " << LARGE_CLUSTER_THRESHOLD << " reads)\n";
        cout << "Found " << small_cluster_indices.size() << " small cluster(s) (< " << LARGE_CLUSTER_THRESHOLD << " reads)\n\n";
        analysis_file << "Found " << small_cluster_indices.size() << " small cluster(s) (< " << LARGE_CLUSTER_THRESHOLD << " reads)\n\n";
        
        // Track validated and unvalidated small clusters
        std::vector<std::pair<size_t, std::vector<size_t>>> validated_small_clusters; // small cluster index -> matching large cluster indices
        std::vector<size_t> unvalidated_small_clusters;
        
        // Compare each small cluster against all large clusters
        for (size_t small_idx : small_cluster_indices) {
            std::vector<size_t> matching_large_clusters;
            const auto& small_cons = all_consensuses[small_idx];
            
            for (size_t large_idx : large_cluster_indices) {
                const auto& large_cons = all_consensuses[large_idx];
                
                // PERFORMANCE FIX 1: Skip if they come from the same FASTQ source
                if (small_cons.sample_name == large_cons.sample_name) {
                    continue; 
                }
                
                // PERFORMANCE FIX 2: Fast heuristic pre-check using Minimizers
                // If their alignment-free k-mer distance is high, they cannot pass a tight Hamming threshold
                double min_dist = minimizer_distance(small_cons.sequence, large_cons.sequence, 11, 5);
                if (min_dist > 0.45) { 
                    continue; // Skip slow graph alignment completely!
                }
                
                // If they pass the pre-check, compute precise alignment distance
                int dist = hamming_distance_aligned(small_cons.sequence, large_cons.sequence);
                if (dist <= HAMMING_THRESHOLD) {
                    matching_large_clusters.push_back(large_idx);
                }
            }
            
            if (!matching_large_clusters.empty()) {
                validated_small_clusters.push_back({small_idx, matching_large_clusters});
            } else {
                unvalidated_small_clusters.push_back(small_idx);
            }
        }
        
        // Report validated small clusters
        cout << "========== VALIDATED Small Clusters (match large clusters) ==========\n";
        analysis_file << "========== VALIDATED Small Clusters (match large clusters) ==========\n";
        cout << "Count: " << validated_small_clusters.size() << "\n\n";
        analysis_file << "Count: " << validated_small_clusters.size() << "\n\n";
        
        for (size_t v = 0; v < validated_small_clusters.size(); v++) {
            auto [small_idx, large_matches] = validated_small_clusters[v];
            const auto& small_cons = all_consensuses[small_idx];
            
            cout << "Small Cluster " << v + 1 << ":\n";
            analysis_file << "Small Cluster " << v + 1 << ":\n";
            cout << "  Sample: " << small_cons.sample_name << "\n";
            analysis_file << "  Sample: " << small_cons.sample_name << "\n";
            cout << "  Size: " << small_cons.cluster_size << " reads\n";
            analysis_file << "  Size: " << small_cons.cluster_size << " reads\n";
            cout << "  Matches with " << large_matches.size() << " large cluster(s):\n";
            analysis_file << "  Matches with " << large_matches.size() << " large cluster(s):\n";
            
            for (size_t large_idx : large_matches) {
                const auto& large_cons = all_consensuses[large_idx];
                int dist = hamming_distance_aligned(small_cons.sequence, large_cons.sequence);
                cout << "    - " << large_cons.sample_name << " (size: " << large_cons.cluster_size 
                     << " reads, distance: " << dist << ")\n";
                analysis_file << "    - " << large_cons.sample_name << " (size: " << large_cons.cluster_size 
                              << " reads, distance: " << dist << ")\n";
            }
            
            cout << "  Sequence: " << small_cons.sequence << "\n";
            analysis_file << "  Sequence: " << small_cons.sequence << "\n";
            cout << "  POSSIBLE VARIANT: Confirmed by appearance in large clusters\n\n\n\n";
            analysis_file << "  POSSIBLE VARIANT: Confirmed by appearance in large clusters\n\n\n\n";
        }
        
        // Report unvalidated small clusters
        cout << "========== UNVALIDATED Small Clusters (no match with large clusters) ==========\n";
        analysis_file << "========== UNVALIDATED Small Clusters (no match with large clusters) ==========\n";
        cout << "Count: " << unvalidated_small_clusters.size() << "\n\n";
        analysis_file << "Count: " << unvalidated_small_clusters.size() << "\n\n";
        
        for (size_t u = 0; u < unvalidated_small_clusters.size(); u++) {
            size_t small_idx = unvalidated_small_clusters[u];
            const auto& small_cons = all_consensuses[small_idx];
            
            cout << "Small Cluster " << u + 1 << ":\n";
            analysis_file << "Small Cluster " << u + 1 << ":\n";
            cout << "  Sample: " << small_cons.sample_name << "\n";
            analysis_file << "  Sample: " << small_cons.sample_name << "\n";
            cout << "  Size: " << small_cons.cluster_size << " reads\n";
            analysis_file << "  Size: " << small_cons.cluster_size << " reads\n";
            cout << "  Sequence: " << small_cons.sequence << "\n";
            analysis_file << "  Sequence: " << small_cons.sequence << "\n";
            cout << "  POSSIBLE NOISE: No match with large clusters\n\n\n\n";
            analysis_file << "  POSSIBLE NOISE: No match with large clusters\n\n\n\n";
        }
        
        // Summary
        cout << "========== Summary ==========\n";
        analysis_file << "========== Summary ==========\n";
        cout << "Total small clusters: " << small_cluster_indices.size() << "\n";
        analysis_file << "Total small clusters: " << small_cluster_indices.size() << "\n";
        cout << "Validated (possible variants): " << validated_small_clusters.size() << "\n";
        analysis_file << "Validated (possible variants): " << validated_small_clusters.size() << "\n";
        cout << "Unvalidated (possible noise): " << unvalidated_small_clusters.size() << "\n";
        analysis_file << "Unvalidated (possible noise): " << unvalidated_small_clusters.size() << "\n";
        
        analysis_file.close();
        cout << "\nCross-sample analysis saved to: " << analysis_output << "\n";
    }

    return 0;
}