//File: BestMichelDistance.h
//Brief: Require the the distance to vertex for the best Michel electron is below some value.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//PlotUtils includes
#include "PlotUtils/Cut.h"

template <class UNIVERSE, class EVENT>
class BestMichelDistance: public PlotUtils::Cut<UNIVERSE, EVENT>
{
  public:
    BestMichelDistance(const double maxDistance): PlotUtils::Cut<UNIVERSE, EVENT>("Best Michel Distance is < " + std::to_string(maxDistance) + "mm"), m_maxDistance(maxDistance)
    {
    }

  private:
    double m_maxDistance; //Maximum distance from the vertex that the best Michel can have in mm

    bool checkCut(const UNIVERSE& univ, EVENT& evt) const
    {
      return evt.m_bestdist < m_maxDistance;
    }
};
