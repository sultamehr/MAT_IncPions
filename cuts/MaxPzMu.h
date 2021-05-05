//File: MaxPzMu
//Brief: Don't allow any CC neutrino interactions where the primary lepton's z momentum went above some value.
//       Important when projecting Dan's 2D inclusive analysis efficiency into 1D in pT.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PZMUMAX_H
#define PZMUMAX_H

#include "PlotUtils/Cut.h"

template <class UNIVERSE>
class MaxPzMu: public PlotUtils::SignalConstraint<UNIVERSE>
{
  public:
    MaxPzMu(const double max): PlotUtils::SignalConstraint<UNIVERSE>(std::string("PzMu < ") + std::to_string(max)), fMax(max)
    {
    }

  private:
    bool checkConstraint(const UNIVERSE& univ) const override
    {
      return univ.GetPlepTrue() * cos(univ.GetThetalepTrue()) <= fMax;
    }

    const double fMax;
};

#endif //PZMUMAX_H
