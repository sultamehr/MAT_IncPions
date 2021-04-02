//PlotUtils includes
#include "PlotUtils/Cut.h"

////Package includes
#include "CVUniverse.h"

template <class CVUniverse, class EVENT>

class Q3RangeReco: public PlotUtils::Cut<CVUniverse, EVENT>

{
  public:
    Q3RangeReco(const double q3min,const double q3max): PlotUtils::Cut<CVUniverse, EVENT>(std::to_string(q3min)+"q3 reco < " + std::to_string(q3max)), 
    fQ3min(q3min),
    fQ3max(q3max)
    {
    }
  private:
    double fQ3min; //minimum q3 allowed in GeV/c;
    double fQ3max; //maximum q3 allowed in GeV/c
    bool checkCut(const CVUniverse& univ, EVENT&) const
    {
      const double q3reco = univ.Getq3();
      return q3reco > fQ3min && q3reco < fQ3max; 
    }

};
