#include "sequence_analyzer.hpp"
#include "../main.cpp"

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

