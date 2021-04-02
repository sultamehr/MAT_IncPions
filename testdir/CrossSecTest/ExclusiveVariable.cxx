#ifndef EXCLUSIVEVARIABLE_CXX
#define EXCLUSIVEVARIABLE_CXX

#include "ExclusiveVariable.h" //TODO: This is technically a circular include which is very confusing
#include <iostream>
using namespace PlotUtils;

//==============================================================================
//
// Variable Base
//
//==============================================================================
//==============================================================================
// CTORs
//==============================================================================
// default (off-limits at the moment)
template <class UNIVERSE, class EVENT>
ExclusiveVariable<UNIVERSE, EVENT>::ExclusiveVariable()
    : m_name(),
      m_xaxis_label(),
      m_binning() {}

// variable binning
template <class UNIVERSE, class EVENT>
ExclusiveVariable<UNIVERSE, EVENT>::ExclusiveVariable(const std::string name,
                                     const std::string xaxis_label,
                                     const std::vector<double> binning,
                                     PointerToCVUniverseFunction reco_func,
                                     PointerToCVUniverseFunction true_func)
    : m_name(name),
      m_xaxis_label(xaxis_label),
      m_binning(GetSortedVector(binning)),
      m_pointer_to_GetRecoValue(reco_func),
      m_pointer_to_GetTrueValue(true_func) {}

// uniform binning
template <class UNIVERSE, class EVENT>
ExclusiveVariable<UNIVERSE, EVENT>::ExclusiveVariable(const std::string name,
                                     const std::string xaxis_label,
                                     const int nbins, const double xmin,
                                     const double xmax,
                                     PointerToCVUniverseFunction reco_func,
                                     PointerToCVUniverseFunction true_func)
    : m_name(name),
      m_xaxis_label(xaxis_label),
      m_binning(MakeUniformBinning(nbins, xmin, xmax)),
      m_pointer_to_GetRecoValue(reco_func),
      m_pointer_to_GetTrueValue(true_func) {}

//==============================================================================
// Set/Get
//==============================================================================
template <class UNIVERSE, class EVENT>
std::string ExclusiveVariable<UNIVERSE, EVENT>::GetName() const {
  return m_name;
}

template <class UNIVERSE, class EVENT>
std::string ExclusiveVariable<UNIVERSE, EVENT>::GetAxisLabel() const {
  return m_xaxis_label;
}

template <class UNIVERSE, class EVENT>
int ExclusiveVariable<UNIVERSE, EVENT>::GetNBins() const {
  return m_binning.size() - 1;
}

template <class UNIVERSE, class EVENT>
void ExclusiveVariable<UNIVERSE, EVENT>::PrintBinning() const {
  std::cout << GetName() << " binning: ";
  for (const auto b : m_binning) std::cout << b << " ";
  std::cout << "\n";
}

template <class UNIVERSE, class EVENT>
std::vector<double> ExclusiveVariable<UNIVERSE, EVENT>::GetBinVec() const {
  return m_binning;
}

//==============================================================================
// Get Reco and True Values
//==============================================================================
template <class UNIVERSE, class EVENT>
double ExclusiveVariable<UNIVERSE, EVENT>::GetRecoValue(const UNIVERSE& universe,
                                            const int idx1,
                                            const int idx2) const {
  return m_pointer_to_GetRecoValue(universe);
}

template <class UNIVERSE, class EVENT>
double ExclusiveVariable<UNIVERSE, EVENT>::GetTrueValue(const UNIVERSE& universe,
                                            const int idx1,
                                            const int idx2) const {
  return m_pointer_to_GetTrueValue(universe);
}

//==============================================================================
// Binning Helpers
//==============================================================================
template <class UNIVERSE, class EVENT>
std::vector<double> ExclusiveVariable<UNIVERSE, EVENT>::MakeUniformBinning(
    const int nbins, const double min, const double max) {
  double step_size = (max - min) / nbins;
  double arr[nbins + 1];  // +1 because binning arrays include top edge.
  for (int i = 0; i <= nbins; ++i) arr[i] = min + i * step_size;
  const int size = sizeof(arr) / sizeof(*arr);
  std::vector<double> ret(arr, arr + size);
  return ret;
}

template <class UNIVERSE, class EVENT>
std::vector<double> ExclusiveVariable<UNIVERSE, EVENT>::GetSortedVector(
    const std::vector<double>& vin) {
  std::vector<double> vout = vin;
  std::sort(vout.begin(), vout.end());
  return vout;
}

#endif  // EXCLUSIVEVARIABLE_CXX
