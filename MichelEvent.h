#ifndef MichelEvent_h
#define MichelEvent_h

#include "Michel.h"
#include "CVUniverse.h"
#include "MatchedMichel.h"

struct MichelEvent {
    int m_idx; // Index for Best Michel in nmichels
    double m_bestdist; // in mm 
    std::vector<double> m_best2D; //0: XZ, 1: UZ, 2:VZ   
    double m_best_XZ;
    double m_best_UZ;
    double m_best_VZ;
    std::vector<Michel*> m_nmichels; //nmatched michels
};
#endif
