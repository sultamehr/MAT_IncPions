#ifndef EXCLUSIVEVARIABLE_H
#define EXCLUSIVEVARIABLE_H

#include <algorithm>
#include <string>
#include <vector>
#include <functional>

//==============================================================================
//
// EXCLUSIVE VARIABLE (EVENT) BASE CLASS
//
//==============================================================================
#ifndef __CINT__  // CINT doesn't know about std::function
template <class UNIVERSE, class EVENT>
class ExclusiveVariable : public VARIABLE {
 private:
  typedef std::function<double(const UNIVERSE&, const EVENT&)>
      PointerToCVUniverseEventFunction;

 public:
  //============================================================================
  // CTORS
  //============================================================================
  // uniform binning
  ExclusiveVariable(
      const std::string name, const std::string xaxis_label, const int nbins,
      const double xmin, const double xmax,
      PointerToCVUniverseEventFunction reco_func,
      PointerToCVUniverseEventFunction true_func);

  // variable binning
  ExclusiveVariable(
      const std::string name, const std::string xaxis_label,
      const std::vector<double> binning,
      PointerToCVUniverseEventFunction reco_func,
      PointerToCVUniverseEventFunction true_func);

  //Getters
  std::string GetName() const;
  int GetNBins() const;
  std::vector<double> GetBinVec() const;
  void PrintBinning() const;

  //============================================================================
  // GetValue
  //============================================================================
  virtual double GetRecoValue(const UNIVERSE& universe, const EVENT& evt) const
  {
    return m_pointer_to_GetRecoValue(universe, evt);
  }

  //The Truth loop never needs to know about EVENT as of 3/24/2021
  virtual double GetTrueValue(const UNIVERSE& universe) const
  {
    return m_pointer_to_GetTrueValue(universe);
  }

 protected:
  std::string m_xaxis_label;

 private:
  //============================================================================
  // Data members
  //============================================================================
  PointerToCVUniverseEventFunction m_pointer_to_GetRecoValue;
  PointerToCVUniverseEventFunction m_pointer_to_GetTrueValue;

  std::string m_name;
  std::vector<double> m_binning;

  ExclusiveVariable();  // Default is off-limits for the time being

  // Helper functions
  std::vector<double> MakeUniformBinning(const int nbins, const double min,
                                         const double max);
  std::vector<double> GetSortedVector(const std::vector<double>& vin);
};
#endif  // __CINT__

#include "ExclusiveVariable.cxx"

#endif  // EXCLUSIVEVARIABLE_H
