//==============================================================================
// MINERvA Analysis Toolkit "Minimum Adoption" Event Loop Example
//
// "Minimum adoption" means that it only uses the two essential MAT tools:
// Universes and HistWrappers. For an "full adoption" example that additionally
// makes use of Cuts, MacroUtil, and Variable, refer to the example in
// ../MAT_Tutorial/.
//
// This script follows the canonical event-looping structure:
// Setup (I/O, variables, histograms, systematics)
// Loop events
//   loop universes
//     make cuts
//     loop variables
//       fill histograms
// Plot and Save
//==============================================================================


//PlotUtils includes
//No junk from PlotUtils please!  I already
//know that MnvH1D does horrible horrible things.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

#include "CVUniverse.h"
#include "Systematics.h"

#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#ifndef __CINT__
#include "PlotUtils/MnvPlotter.h"
#include "Variable.h"
#endif

#include "MichelEvent.h"
#include "MatchedMichel.h"

#include "Michel.h"
#include "Cluster.h"
#pragma GCC diagnostic pop

#include "SafeROOTName.cpp" //TODO: This is a very bad idea
#include "Categorized.h"
#include <iostream>
// Histogram binning constants
const int nbins = 30;
const double xmin = 0.;
const double xmax = 20.e3;

//==============================================================================
// Cuts
//==============================================================================
bool PassesCuts(CVUniverse& univ) {
  return //univ.IsMinosMatchMuon() &&
         //univ.GetMuonQP() < 0.0 &&
         //univ.GetTDead() <= 1&&
         univ.GetEmu() > 2e3; // GetEmu is in MeV
}

//==============================================================================
// Plot
//==============================================================================
void PlotErrorSummary(PlotUtils::MnvH1D* h, std::string label) {
  PlotUtils::MnvH1D* hist = (PlotUtils::MnvH1D*)h->Clone("hist");
  PlotUtils::MnvPlotter mnv_plotter( kCompactStyle);
  TCanvas cE("c1", "c1");

  // Group GENIE bands
    mnv_plotter.error_summary_group_map.clear();
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrAbs_N");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrAbs_pi");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrCEx_N");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrCEx_pi");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrElas_N");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrElas_pi");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrInel_N");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrPiProd_N");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back(
        "GENIE_FrPiProd_pi");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_MFP_N");
    mnv_plotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_MFP_pi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_AGKYxF1pi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_AhtBY");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_BhtBY");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_CCQEPauliSupViaKF");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_CV1uBY");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_CV2uBY");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_EtaNCEL");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_MaNCEL");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_MaRES");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_MvRES");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_NormDISCC");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_NormNCRES");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_RDecBR1gamma");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvn1pi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvn2pi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvn3pi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvp1pi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvp2pi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Theta_Delta2Npi");
    mnv_plotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_VecFFCCQEshape");

    mnv_plotter.error_summary_group_map["RPA"].push_back("RPA_HighQ2");
    mnv_plotter.error_summary_group_map["RPA"].push_back("RPA_LowQ2");

  const bool do_fractional_uncertainty = true;
  const bool do_cov_area_norm = false;
  const bool include_stat_error = false;
  const std::string do_fractional_uncertainty_str =
      do_fractional_uncertainty ? std::string("Frac") : std::string("Abs");
  const std::string do_cov_area_norm_str =
      do_cov_area_norm ? std::string("CovAreaNorm") : std::string("");

  mnv_plotter.DrawErrorSummary(hist, "TR", include_stat_error, true,
                              0.0, do_cov_area_norm, "",
                              do_fractional_uncertainty);
  mnv_plotter.AddHistoTitle("Event Selection");
  std::string plotname =
      Form("ErrorSummary_%s_%s_%s", do_fractional_uncertainty_str.c_str(),
           do_cov_area_norm_str.c_str(), label.c_str());
  mnv_plotter.MultiPrint(&cE, plotname, "png");
}

void PlotCVAndError(PlotUtils::MnvH1D* h, std::string label) {
  PlotUtils::MnvH1D* hist = (PlotUtils::MnvH1D*)h->Clone("hist");
  PlotUtils::MnvPlotter mnv_plotter(PlotUtils::kCCNuPionIncStyle);
  TCanvas cE("c1", "c1");
  if (!gPad)
    throw std::runtime_error("Need a TCanvas. Please make one first.");
  PlotUtils::MnvH1D* datahist = new PlotUtils::MnvH1D(
      "adsf", "", nbins, xmin, xmax);
  bool statPlusSys = true;
  int mcScale = 1.;
  bool useHistTitles = false;
  const PlotUtils::MnvH1D* bkgdHist = NULL;
  const PlotUtils::MnvH1D* dataBkgdHist = NULL;
  mnv_plotter.DrawDataMCWithErrorBand(datahist, hist, mcScale, "TL",
                                     useHistTitles, NULL, NULL, false,
                                     statPlusSys);
  mnv_plotter.AddHistoTitle("Event Selection");
  std::string plotname = Form("CV_w_err_%s", label.c_str());
  mnv_plotter.MultiPrint(&cE, plotname, "png");
  delete datahist;
}


//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillEventSelection(
    PlotUtils::ChainWrapper* chain,
    std::map<std::string, std::vector<CVUniverse*> > error_bands,
    util::Categorized<PlotUtils::HistWrapper<CVUniverse>, int> hw_bestpirange) //,PlotUtils::HistWrapper<CVUniverse> hw_emu)
    {
    for (int i=0; i<chain->GetEntries(); ++i) {
    if(i%1000==0) std::cout << (i/1000) << "k " << std::flush;


    int isSignal = 0;
    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : error_bands) {
      std::vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes) {

        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);
        
        // THis is where you would Access/create a Michel
        if (universe->Getq3() > 1.2) continue;
         //std::cout<< "  q3 for the event is " << universe->Getq3() << std::fl;

          int nmichels = universe->GetNMichels();
          
          //std::vector<Michel*> eventmichels = universe->CreateMichels();
         
          //Michel* m = new Michel(*universe, 0);
          
          //This is the function for filling all the Michel info.
          
          double min_dist = 9999.;
          std::vector<Michel*> return_michels;
          for (int i = 0; i < nmichels; ++i)
          {
            Michel* current_michel = new Michel(*universe, i);
            if (current_michel->is_fitted != 1) continue;
            
            //std::cout << "Getting Michel Vtx Match " << std::endl;
            current_michel->DoesMichelMatchVtx(*universe, current_michel);
            //std::cout << "Getting Michel Cluster Match " << std::endl;
            current_michel->DoesMichelMatchClus(*universe, current_michel);
            //std::cout << "Getting Best Michel Match " << std::endl;

            current_michel->GetBestMatch(current_michel);
            //std::cout << "Best Match for this michel is " << current_michel->BestMatch << std::endl;
            
            //std::cout << "INITIALIZED CURRENT MICHEL " << i << std::endl;
            //std::cout << "Filling container of Michel for the event." << std::endl;

            double dist = current_michel->Best3Ddist;
              if (dist <= min_dist) {
              min_dist = dist;
              }
              
            return_michels.push_back(current_michel);
            //delete current_michel;
          }

          if (return_michels.size() == 0) continue; 

          //std::cout << "Best Michel Distance is " << min_dist << std::endl;
          
          /*hw_bestpirange.univHist(universe)->Fill(min_dist,
                                                  universe->GetWeight());*/
          hw_bestpirange[universe->GetInteractionType()].FillUniverse(universe, min_dist, universe->GetWeight());
        //=========================================
        // Cuts in each universe
        //=========================================
        //if(PassesCuts(*universe)) {
        //  hw_emu.univHist(universe)->Fill(universe->GetEmu(), 
                                        //  universe->GetWeight());
          
        //} // End if passed cuts
      } // End band's universe loop
    } // End Band loop
    
    //if (isSignal == 1) std::cout << "THIS EVENT HAS A SIGNAL MICHEL " << std::endl;

  } //End entries loop
}

//==============================================================================
// Main
//==============================================================================
void runEventLoop() {
  TH1::AddDirectory(false);
  // Make a chain of events
  PlotUtils::ChainWrapper* chain = makeChainWrapperPtr("playlist_mc.txt", 
                                                       "CCQENu");

  // You're required to make some decisions
  PlotUtils::MinervaUniverse::SetNuEConstraint(false);
  PlotUtils::MinervaUniverse::SetPlaylist("minervame1A");
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
  PlotUtils::MinervaUniverse::SetNonResPiReweight(false);
  PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false);

  // Make a map of systematic universes
  std::map< std::string, std::vector<CVUniverse*> > error_bands = 
      GetStandardSystematics(chain);

  // HistWrapper
  //PlotUtils::HistWrapper<CVUniverse> hw_emu("hw_emu", "E_{#mu} (MeV)",
                                          //  nbins, xmin, xmax, error_bands);
                                            
  //std::map<int, std::string> signalOrBackground = {{1, "signal"}};
  //For example only.  Don't actually use GENIELabels as your backgrounds!
  //You'd end up with a very model-dependent result, and Luke Pickering
  //would frown on your paper ;)
  std::map<int, std::string> GENIELabels = {{1, "QE"},
                                            {8, "2p2h"},
                                            {2, "RES"},
                                            {3, "DIS"}};
  util::Categorized<PlotUtils::HistWrapper<CVUniverse>, int> hw_bestpirange("hw_bestpirange", "Pion Range (mm)",
                                                                            GENIELabels, 100, 0., 1000., error_bands);

  // Loop entries and fill
  LoopAndFillEventSelection(chain, error_bands, hw_bestpirange);
  //LoopAndFillEventSelection(chain, error_bands,);

  // You must always sync your HistWrappers after filling them
  //hw_emu.SyncCVHistos();
  #ifndef __CINT__ //For lambda function (thing with [](){} is an inline "lambda function" to save me time typing)
  hw_bestpirange.visit([](PlotUtils::HistWrapper<CVUniverse>& wrapper) { wrapper.SyncCVHistos(); });
  #endif //__CINT__

  // Plot stuff
  //hw_emu.hist->GetXaxis()->SetTitle("E_{#mu} (MeV)");
  //PlotCVAndError(hw_emu.hist, "Emu");
  //PlotErrorSummary(hw_emu.hist, "Emu");
  
  //hw_bestpirange.hist->GetXaxis()->SetTitle("Best Michel-Pion Range (mm)");
  #ifndef __CINT__ //For "auto" c++11 feature because Andrew is lazy
  for(const auto& label: GENIELabels)
  {
    PlotCVAndError(hw_bestpirange[label.first].hist, label.second);
    PlotErrorSummary(hw_bestpirange[label.first].hist, label.second);
  }
  #endif //__CINT__

  std::cout << "Success" << std::endl;
}
/*
void CreateMichels(CVUnivere &univ,std::vector<Michel*> &return_michels){
  unsigned int nmichels = univ->GetNMichels();
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
