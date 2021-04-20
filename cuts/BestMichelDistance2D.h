//File: BestMichelDistance2D.h
//Brief: Require the the distance to vertex for the best Michel electron is below some value.
//Author: Andrew Olivier aolivier@ur.rochester.edu, Mehreen Sultana msultana@ur.rochester.edu

//PlotUtils includes
#include "PlotUtils/Cut.h"

template <class UNIVERSE, class EVENT>
class BestMichelDistance2D: public PlotUtils::Cut<UNIVERSE, EVENT>
{
  public:
    BestMichelDistance2D(const double maxDistance): PlotUtils::Cut<UNIVERSE, EVENT>("Best Michel Distance is < " + std::to_string(maxDistance) + "mm"), m_maxDistance(maxDistance)
    {
    }

  private:
    double m_maxDistance; //Maximum distance from the vertex that the best Michel can have in mm

    bool checkCut(const UNIVERSE& univ, EVENT& evt) const
    {
      int best_idx = evt.m_idx;
      double XZ = evt.m_best_XZ;
      double UZ = evt.m_best_UZ;
      double VZ = evt.m_best_VZ;
      bool pass = false;
      if (XZ < m_maxDistance && (UZ < m_maxDistance || VZ < m_maxDistance)) pass = true;
      else if (UZ < m_maxDistance && (XZ < m_maxDistance || VZ < m_maxDistance)) pass = true;
      else if (VZ < m_maxDistance && (XZ < m_maxDistance || UZ < m_maxDistance)) pass = true;
      else pass = false;
      return pass;
    }
};
