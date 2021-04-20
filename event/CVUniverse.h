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

class CVUniverse : public PlotUtils::MinervaUniverse {
 public:
  #include "PlotUtils/SystCalcs/WeightFunctions.h" // Get*Weight
  #include "PlotUtils/SystCalcs/MuonFunctions.h" // GetMinosEfficiencyWeight
  #include "PlotUtils/SystCalcs/TruthFunctions.h" //Getq3True
  // ========================================================================
  // Constructor/Destructor
  // ========================================================================
  CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma = 0)
      : PlotUtils::MinervaUniverse(chw, nsigma) {}

  virtual ~CVUniverse() {}

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
  
  /*
  std::vector<Michel*> CreateMichels(){
  std::vector<Michel*> return_michels;
  unsigned int nmichels = GetNMichels();
  CVUniverse* univ;
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

  //Muon kinematics
  double GetMuonPT() const //GeV/c
  {
    return GetPmu()/1000. * sin(GetThetamu());
  }

  double GetMuonPTTrue() const //GeV/c
  {
    return GetPlepTrue()/1000. * sin(GetThetalepTrue());
  }

  int GetInteractionType() const {
    return GetInt("mc_intType");
  }
  
  
  virtual bool IsMinosMatchMuon() const {
    return GetInt("has_interaction_vertex") == 1;
  }
  
  ROOT::Math::XYZTVector GetVertex() const
  {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("vtx").data());
    return result;
  }

  ROOT::Math::XYZTVector GetTrueVertex() const
  {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("mc_vtx").data());
    return result;
  }
  
  virtual int GetNMichels() const{
      return GetInt("FittedMichel_michel_fitPass_sz");
  }
  
  virtual double GetMichelEnergy(int i){
      return GetVecElem("FittedMichel_michel_energy", i);
  }
  
  virtual int GetMichelFit(int i){
      return GetVecElem("FittedMichel_michel_fitPass", i);
  }
  
  virtual double GetMichelTime(int i){
      return GetVecElem("FittedMichel_michel_time", i);
  }
  
  virtual double GetMichelX1(int i){
      return GetVecElem("FittedMichel_michel_x1", i);
  }
  
  virtual double GetMichelX2(int i){
      return GetVecElem("FittedMichel_michel_x2", i);
  }
  
  virtual double GetMichelU1(int i){
      return GetVecElem("FittedMichel_michel_u1", i);
  }
  
  virtual double GetMichelU2(int i){
      return GetVecElem("FittedMichel_michel_u2", i);
  }
  
  
  virtual int GetNTruePions() const{
      return GetInt("FittedMichel_all_piontrajectory_trackID_sz");
  }
  
  virtual int GetPionParentID(int i) const {
     return GetVecElem("FittedMichel_all_piontrajectory_ParentID", i);
  }
  
  virtual int GetPionPDG(int i) const{
     return GetVecElem("FittedMichel_all_piontrajectory_pdg", i);
  }
  
  virtual double GetPionE(int i) const{
     return GetVecElem("FittedMichel_all_piontrajectory_energy",i)/pow(10,3);
  }
  
  virtual double GetPionP(int i) const{
    return GetVecElem("FittedMichel_all_piontrajectory_momentum", i)/pow(10,3);
  }
  
  virtual double GetPionMass(int i) const{
    double pionmass = sqrt(pow(GetPionE(i), 2) - pow(GetPionP(i), 2));
    return pionmass;
  }
  
  virtual double GetPionKE(int i) const{
    return GetPionE(i) - GetPionMass(i);
  }

  virtual double GetLowTpi() const {
     double lowesttpi = 9999.;
     int nFSpi = GetNTruePions();
     for (int i = 0; i < nFSpi; i++){
         int pdg = GetPionPDG(i);
         int trackid = GetPionParentID(i);
         if (trackid != 0 ) continue;
         if (pdg != 211) continue;
         double pionKE = GetPionKE(i);
         if (lowesttpi > pionKE) lowesttpi = pionKE;
      }
      //std::cout << "lowest energy pi << " << lowesttpi << std::endl;
      return lowesttpi;
  }
  virtual double GetTrueQ2() const {
    return GetDouble("mc_Q2");
  }

  virtual int GetTDead() const {
    return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj");;
  }
  
  virtual double GetEavail() const {
    return GetDouble("recoilE_SplineCorrected");
  }
  
  virtual double GetQ2Reco() const{
    return GetDouble("qsquared_recoil");
  }
  
  virtual double Getq3() const{
    double eavail = GetEavail()/pow(10,3);
    double q2 = GetQ2Reco() / pow(10,6);
    double q3mec = sqrt(eavail*eavail + q2);
    return q3mec;
  }
   
  virtual int GetCurrent() const { return GetInt("mc_current"); }

  virtual int GetTruthNuPDG() const { return GetInt("mc_incoming"); }

  virtual double GetMuonQP() const {
    return GetDouble("CCQENu_minos_trk_qp");
  }

  //virtual double GetLowTpi() const{

  //}
 
 // virual double GetBestMichelDist(int i) const {
     
 //  }
  
};

#endif
