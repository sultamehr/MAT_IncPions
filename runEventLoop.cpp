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

#include "event/CVUniverse.h"
#include "systematics/Systematics.h"

#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/MacroUtil.h"
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/cuts/CCInclusiveCuts.h"
#include "PlotUtils/cuts/CCInclusiveSignal.h"
#include "util/Categorized.h"
#include "PlotUtils/Cutter.h"
#include "util/Variable.h"
#include "cuts/BestMichelDistance.h"
#include "cuts/BestMichelDistance2D.h"
#include "cuts/SignalDefinition.h"
#include "cuts/q3RecoCut.h"
#include "cuts/hasMichel.h"
#include "studies/Study.h"
#include "studies/PerMichelVarByGENIELabel.h"
#include "studies/PerMichelEventVarByGENIELabel.h" 
#include "studies/CreateDataHistPerMichelEvent.h"
#include "studies/CreateDataHistPerMichel.h"
#include "event/Michel.h"
#include "event/Cluster.h"
#include "event/MichelEvent.h"
#include "event/MatchedMichel.h"
//#include "Binning.h" //TODO: Fix me
#pragma GCC diagnostic pop
#include <iostream>

class Variable;

// Histogram binning constants
const int nbins = 30;
const double xmin = 0.;
const double xmax = 20.e3;

//==============================================================================
// Plot
//==============================================================================
void PlotErrorSummary(PlotUtils::MnvH1D* h, std::string label)
{
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

void PlotCVAndError(PlotUtils::MnvH1D* h, std::string label)
{
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
    std::vector<Variable*> vars,
    std::vector<Study*> studies,
    PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts)
{
  for (int i=0; i<chain->GetEntries(); ++i)
  {
    if(i%1000==0) std::cout << (i/1000) << "k " << std::flush;

    int isSignal = 0;
    
    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : error_bands)
    {
      std::vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes)
      {
         MichelEvent myevent; // make sure your event is inside the error band loop. 
    
        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);
        
        // THis is where you would Access/create a Michel

        //std::cout << "Best Michel Distance is " << min_dist << std::endl;
        
        if (!michelcuts.isMCSelected(*universe, myevent, 1).all()) continue; //all is another function that will later help me with sidebands

        //For each Variable
        //if(best_michel != -1) //Equivalent to if (return_michels.size() == 0) continue;
      
        for(auto& var: vars)
        {
          for(auto& study: studies) study->SelectedSignal(*universe, myevent, 1); //TODO: Last argument is weight
          (*var->m_bestPionByGENIELabel)[universe->GetInteractionType()].FillUniverse(universe, var->GetRecoValue(*universe, myevent.m_idx), universe->GetWeight());

          //Fill other per-Variable histograms here
          
        }
      } // End band's universe loop
    } // End Band loop
  } //End entries loop
}


void LoopAndFillData( PlotUtils::ChainWrapper* data,
			        std::vector<CVUniverse*> data_band,
				std::vector<Variable*> vars,
                                std::vector<Study*> studies,
				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts)

{
  for (auto universe : data_band) {
    for (int i=0; i<data->GetEntries(); ++i) {
      universe->SetEntry(i);
      if(i%1000==0) std::cout << (i/1000) << "k " << std::flush;
      MichelEvent myevent; 
      if (!michelcuts.isDataSelected(*universe, myevent).all()) continue;
      for(auto& study: studies) study->Selected(*universe, myevent, 1); 
      for(auto& var: vars)
      {
        (var->dataHist)->FillUniverse(universe, var->GetRecoValue(*universe, myevent.m_idx), 1);
      }
    }
  }
}



void LoopAndFillEffDenom( PlotUtils::ChainWrapper* truth,
    				std::map<std::string, std::vector<CVUniverse*> > truth_bands,
    				std::vector<Variable*> vars,
    				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts)
{
  for (int i=0; i<truth->GetEntries(); ++i)
  {
    if(i%1000==0) std::cout << (i/1000) << "k " << std::flush;

    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : truth_bands)
    {
      std::vector<CVUniverse*> truth_band_universes = band.second;
      for (auto universe : truth_band_universes)
      {
        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);

        if (!michelcuts.isEfficiencyDenom(*universe, 1)) continue; 
        //Fill efficiency denominator now: 
        //TODO ADD PLOTS
      }
    }
  }
}

//==============================================================================
// Main
//==============================================================================
int main(const int /*argc*/, const char** /*argv*/)
{
  TH1::AddDirectory(false);
  // Make a chain of events
  PlotUtils::ChainWrapper* chain = makeChainWrapperPtr("playlist_mc.txt", 
                                                       "CCQENu");
  PlotUtils::ChainWrapper* truth = makeChainWrapperPtr("playlist_mc.txt",
						       "Truth");
  PlotUtils::ChainWrapper* data = makeChainWrapperPtr("playlist_data.txt",
						      "CCQENu");
  const std::string mc_file_list("playlist_mc.txt");
  const std::string data_file_list("playlist_data.txt");
  const std::string plist_string("minervame1a");
  const std::string reco_tree_name("CCQENu");
  const bool do_truth = false;
  const bool is_grid = false;
  // You're required to make some decisions
  PlotUtils::MinervaUniverse::SetNuEConstraint(false);
  PlotUtils::MinervaUniverse::SetPlaylist("minervame1A");
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
  PlotUtils::MinervaUniverse::SetNonResPiReweight(false);
  PlotUtils::MinervaUniverse::SetDeuteriumGeniePiTune(false);

  // Make a map of systematic universes
  std::map< std::string, std::vector<CVUniverse*> > error_bands = 
      GetStandardSystematics(chain);

  std::map< std::string, std::vector<CVUniverse*> > truth_bands =
      GetStandardSystematics(truth);

  std::vector<Variable*> vars = {
  new Variable("tpi", "T#pi [MeV]", 100, 0., .5, &CVUniverse::GetLowTpi, &CVUniverse::GetLowTpi),
  new Variable("q3", "q3 (GeV)", 100, 0.0, 1.3, &CVUniverse::Getq3, &CVUniverse::Getq3)
  };
  std::vector<Study*> studies;

  std::function<double(const CVUniverse&, const MichelEvent&, const int)> oneVar = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                 {
                                   return evt.m_nmichels[whichMichel]->Best3Ddist;
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> michelE = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
                                 {
                                   return evt.m_nmichels[whichMichel]->energy;
                                 };
  std::function<double(const CVUniverse&, const MichelEvent&, const int)> delta_t = [](const CVUniverse& univ, const MichelEvent& evt, const int whichMichel)
				 {
                                   double micheltime = evt.m_nmichels[whichMichel]->time;
                                   double vtxtime = univ.GetVertex().t();
                                   double deltat = (micheltime - vtxtime/1000.); //hopefully this is in microseconds (mus)
				   return deltat;
				 };

  std::function<double(const CVUniverse&)> low_tpi = [](const CVUniverse& univ)
				{
				   double lowesttpi = 9999.;
                                   int nFSpi = univ.GetNTruePions();
                                   for (int i = 0; i < nFSpi; i++){
         			       int pdg = univ.GetPionPDG(i);
				       int trackid = univ.GetPionParentID(i);
				       if (trackid != 0 ) continue;
  				       if (pdg != 211) continue;
 				       double pionKE = univ.GetPionKE(i);
                                       if (lowesttpi > pionKE) lowesttpi = pionKE;
				   }
                                   return lowesttpi;
				};
  std::function<double(const CVUniverse&, const MichelEvent&)> bestdistvar = [](const CVUniverse& univ, const MichelEvent& evt)
                                 {
                                   return evt.m_bestdist;
                                 }; 
  studies.push_back(new PerMichelVarByGENIELabel(oneVar, "Best Distance", "mm", 100, 0., 1000., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(michelE, "Michel energy", "MeV", 20, 0., 100., error_bands));
  studies.push_back(new PerMichelVarByGENIELabel(delta_t, "Michel Time Diff", "#mus", 30, 0.0, 9.0, error_bands));
  studies.push_back(new PerMichelEventVarByGENIELabel(bestdistvar, "Best Distance", "mm", 100, 0., 1000., error_bands));

  //Creating the single Data universe 
  PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list,
                    plist_string, do_truth, is_grid);
  CVUniverse* data_universe = new CVUniverse(util.m_data);
  std::vector<CVUniverse*> data_band = {data_universe};
  std::map<std::string, std::vector<CVUniverse*> > data_error_bands;
  data_error_bands["cv_data"] = data_band;
  
  std::vector<Study*> data_studies;
  // data_studies.push_back(new CreateDataHistPerMichel(oneVar, "Best Distance", "mm", 100, 0, 1000., data_band));
  // data_studies.push_back(new CreateDataHistPerMichel(michelE, "Michel energy", "MeV", 20, 0., 100., data_band));
  // data_studies.push_back(new CreateDataHistPerMichel(delta_t, "Michel Time Diff", "#mus", 30, 0.0, 9.0, data_band));
  // data_studies.push_back(new CreateDataHistPerMichelEvent(bestdistvar, "Best Distance", "mm", 100, 0, 1000., data_band));
   
  data_studies.push_back(new PerMichelVarByGENIELabel(oneVar, "Best Distance", "mm", 100, 0., 1000., data_error_bands));
  data_studies.push_back(new PerMichelVarByGENIELabel(michelE, "Michel energy", "MeV", 20, 0., 100., data_error_bands));
  data_studies.push_back(new PerMichelVarByGENIELabel(delta_t, "Michel Time Diff", "#mus", 30, 0.0, 9.0, data_error_bands));
  data_studies.push_back(new PerMichelEventVarByGENIELabel(bestdistvar, "Best Distance", "mm", 100, 0., 1000., data_error_bands));

  for(auto& var: vars) var->InitializeMCHists(error_bands);
  for(auto& var: vars) var->InitializeDATAHists(data_band);
  PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t sidebands;
  auto precuts = reco::GetCCInclusive2DCuts<CVUniverse, MichelEvent>();
  precuts.emplace_back(new Q3RangeReco<CVUniverse, MichelEvent>(0.0, 1.20));
  precuts.emplace_back(new hasMichel<CVUniverse, MichelEvent>());
  precuts.emplace_back(new BestMichelDistance2D<CVUniverse, MichelEvent>(102.));
  auto signalDefinition = truth::GetCCInclusive2DSignal<CVUniverse>();
  signalDefinition.emplace_back(new Q3Limit<CVUniverse>(1.2));

  PlotUtils::Cutter<CVUniverse, MichelEvent> mycuts(std::move(precuts), std::move(sidebands) , std::move(signalDefinition),std::move(truth::GetCCInclusive2DPhaseSpace<CVUniverse>()));
  // Loop entries and fill
  LoopAndFillEventSelection(chain, error_bands, vars, studies, mycuts);
  LoopAndFillEffDenom(truth, truth_bands, vars, mycuts);
  TFile* outDir = TFile::Open("OutputMichelHists.root", "RECREATE");
  LoopAndFillData(data, data_band,vars, data_studies, mycuts);

  for(auto& var: vars)
  {
    // You must always sync your HistWrappers before plotting them
    var->SyncCVHistos();

    //Categorized makes sure GetTitle() is the same as the labels you were looping over before
    var->m_bestPionByGENIELabel->visit([](PlotUtils::HistWrapper<CVUniverse>& categ)
                                       {
                                         PlotCVAndError(categ.hist, categ.hist->GetTitle());
                                         PlotErrorSummary(categ.hist, categ.hist->GetTitle());
                                       });
   //var->Write(*outDir); 
  }

  //TFile* outDir = TFile::Open("PerMichelHists.root", "RECREATE");
  for(auto& study: studies) study->SaveOrDraw(*outDir);
  for(auto& var: vars) var->Write(*outDir);
  outDir->Write();
  TFile* dataDir = TFile::Open("DataStudyHists.root", "RECREATE");
  for(auto& study: data_studies) study->SaveOrDraw(*dataDir);  
  std::cout << mycuts << std::endl; 
  std::cout << "Success" << std::endl;
}
