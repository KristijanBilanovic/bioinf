#include "sequence_analyzer.hpp"

std::map<int, std::tuple<int,int>> SequenceAnalyzer::find_nearest_neighbors(){
    std::map<int, std::tuple<int,int>> result;
    int n = seqs.size();

    // we find each sequence Si
    for(int i = 0; i < n; i++){
        int min_dist = INT_MAX;
        int best_j = -1;
         // we compare Si to each other sequence Sj, j != i
        for(int j=0; j < n; j++){
            if(i == j) continue;

            int dist = edit_distance(seqs[i]->data, seqs[j]->data);
            if(dist < min_dist){
                min_dist = dist;
                best_j = j;
            }
        }
        // storing result for each sequence i
        result[i] = std::make_tuple(best_j, min_dist);
    }


    return result;
}

int SequenceAnalyzer::edit_distance(const std::string& a, const std::string& b) {
    int n = a.size();
    int m = b.size();
    
    // table [n+1][m+1]
    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1));

    // filling first row and column
    for (int i = 0; i <= n; i++) dp[i][0] = i;
    for (int j = 0; j <= m; j++) dp[0][j] = j;

    // rest of table
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            // if characters match, no cost
            if (a[i - 1] == b[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1]; // no cost
            
            // if they don't match 
            } else {
                // we choose the best among deletion, insertion, substitution
                dp[i][j] = 1 + std::min({dp[i - 1][j],    // deletion
                                         dp[i][j - 1],    // insertion
                                         dp[i - 1][j - 1] // substitution
                                        });
            }
        }
    }
    return dp[n][m];
}
