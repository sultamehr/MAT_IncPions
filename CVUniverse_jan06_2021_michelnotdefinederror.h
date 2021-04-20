// =============================================================================
// Base class for an un-systematically shifted (i.e. CV) universe. Implement
// "Get" functions for all the quantities that you need for your analysis.
//
// This class inherits from PU::MinervaUniverse, which in turn inherits from
// PU::BaseUniverse. PU::BU defines the interface with anatuples.
// 
// Within the class, "WeightFunctions" and "MuonFunctions" are included to gain
// access to standardized weight and muon variable getters. See:
// https://cdcvs.fnal.gov/redmine/projects/minerva-sw/wiki/MinervaUniverse_Structure_
// for a full list of standardized functions you can use. In general, if a
// standard version of a function is available, you should be using it.
// =============================================================================
#ifndef CVUNIVERSE_H
#define CVUNIVERSE_H


#include <iostream>

#include "PlotUtils/MinervaUniverse.h"
#include "Cluster.h"
#include "Michel.h"


class CVUniverse : public PlotUtils::MinervaUniverse {
 public:
  #include "PlotUtils/SystCalcs/WeightFunctions.h" // Get*Weight
  #include "PlotUtils/SystCalcs/MuonFunctions.h" // GetMinosEfficiencyWeight

  // ========================================================================
  // Constructor/Destructor
  // ========================================================================
  CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma = 0)
      : PlotUtils::MinervaUniverse(chw, nsigma) {}

  virtual ~CVUniverse() {};

  // ========================================================================
  // Get Weight
  // This provides you with some of MnvGENIEv1 weights.
  // Eventually you'll want to validate which weights YOU want to use.
  // And eventually, there will be a common PlotUtils function that does this
  // for you.
  // ========================================================================
  virtual double GetWeight() const {
    double wgt_flux_and_cv = 1., wgt_genie = 1., wgt_2p2h = 1.;
    double wgt_rpa = 1., wgt_mueff = 1.;

    // flux + cv
    wgt_flux_and_cv = GetFluxAndCVWeight();

    // genie
    wgt_genie = GetGenieWeight();

    // 2p2h
    wgt_2p2h = GetLowRecoil2p2hWeight();

    // rpa
    wgt_rpa = GetRPAWeight();

    // MINOS muon tracking efficiency
    if (!IsTruth() && IsMinosMatchMuon())
      wgt_mueff = GetMinosEfficiencyWeight();

    return wgt_flux_and_cv * wgt_genie * wgt_2p2h * wgt_rpa * wgt_mueff;
  }

  // ========================================================================
  // Write a "Get" function for all quantities access by your analysis.
  // For composite quantities (e.g. Enu) use a calculator function.
  //
  // In order to properly calculate muon variables and systematics use the
  // various functions defined in MinervaUniverse.
  // E.g. GetPmu, GetEmu, etc.
  // ========================================================================

  // Quantities only needed for cuts
  // Although unlikely, in principle these quanties could be shifted by a
  // systematic. And when they are, they'll only be shifted correctly if we
  // write these accessor functions.
  virtual bool IsMinosMatchMuon() const {
    return GetInt("minos_track_match");
  }

  virtual double GetTrueQ2() const {
    return GetDouble("mc_Q2");
  }

  virtual int GetTDead() const {
    return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj");;
  }
    
  virtual int GetNMichels() const{
      return GetInt("FittedMichel_michel_fitPass_sz");
  }
  

  virtual double GetWgenie() const { return GetDouble("mc_w"); }
  virtual int GetInteractionType() const { return GetInt("mc_intType"); }
  virtual int GetTruthNuPDG() const { return GetInt("mc_incoming"); }
  virtual int GetCurrent() const { return GetInt("mc_current"); }


  virtual double GetTpi(const int hadron) const {
    double TLA = GetVecElem("hadron_track_length_area", hadron);
    return 2.3112 * TLA + 37.03;
  }
  


};

  /*
  std::vector<Michel*> CreateMichels(){
  std::vector<Michel*> return_michels;
  unsigned int nmichels = GetNMichels();
  for (unsigned int i = 0; i < nmichels; ++i)
  {
    Michel* current_michel = new Michel(univ, i);
    if (current_michel->is_fitted != 1) continue;
    current_michel->DoesMichelMatchVtx(univ, current_michel);
    current_michel->DoesMichelMatchClus(univ, current_michel);
    return_michels.push_back(current_michel);
  }
  return return_michels;
  }
  */

#endif
