#ifndef MichelEvent_h
#define MichelEvent_h

#include "Michel.h"
#include "CVUniverse.h"
#include "MatchedMichel.h"

struct MichelEvent {
    int m_idx; // Index for Best Michel in nmichels
    double m_bestdist; // in mm 
   
    std::vector<Michel*> m_nmichels; //nmatched michels
};
#endif
