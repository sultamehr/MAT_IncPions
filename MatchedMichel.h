#ifndef MatchedMichel_h
#define MatchedMichel_h

#include "Michel.h"
#include "CVUniverse.h"

struct MatchedMichel {
    int m_idx; // Index for matched Michel 
    int m_match_type; // 0 = NULL, 1 = vertex match, 2 = cluster match
        // Write A cut that caluclates the 3D distance 
     std::vector<int> m_best_cluster_idx; // index of the clusters that this michel matches to in each view
};
#endif
