#ifndef Systematics_h
#define Systematics_h

//==============================================================================
// Get Several standard MINERvA systematics
//==============================================================================

#include "event/CVUniverse.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"

typedef std::map<std::string, std::vector<CVUniverse*>> UniverseMap;

UniverseMap GetStandardSystematics(PlotUtils::ChainWrapper* chain,
                                   bool is_truth = false) {
  // return map
  UniverseMap error_bands;

  // CV
  error_bands[std::string("cv")].push_back(new CVUniverse(chain));

  //========================================================================
  // FLUX
  //========================================================================
  UniverseMap bands_flux =
      PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain, CVUniverse::GetNFluxUniverses());
  error_bands.insert(bands_flux.begin(), bands_flux.end());

  //========================================================================
  // GENIE
  //========================================================================
  // Standard
  UniverseMap bands_genie =
      PlotUtils::GetGenieSystematicsMap<CVUniverse>(chain); //PlotUtils::GetStandardGenieSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_genie.begin(), bands_genie.end());

  //========================================================================
  // MnvTunes
  //========================================================================
  // RPA
  UniverseMap bands_rpa = PlotUtils::GetRPASystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_rpa.begin(), bands_rpa.end());

  // 2P2H
  UniverseMap bands_2p2h = PlotUtils::Get2p2hSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_2p2h.begin(), bands_2p2h.end());

    /*
  //========================================================================
  // Muons
  //========================================================================
  // Muon reco in MINERvA -- Catchall systematic for pmu reco in minerva.
  // Lateral-only. Shifts pmu.
  UniverseMap bands_muon_minerva =
      PlotUtils::GetMinervaMuonSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_muon_minerva.begin(), bands_muon_minerva.end());

  // Muons in MINOS -- Catchall systematic for wiggle solution -- correlates
  // flux universes and minos muon momentum reco.
  // Lateral AND Vertical systematic. Shifts Pmu and GetFluxAndCVUniverse.
  //
  // Expect a non-zero systematic even when no pmu involved.
  UniverseMap bands_muon_minos =
     PlotUtils::GetMinosMuonSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_muon_minos.begin(), bands_muon_minos.end());

  // Vertical only
  UniverseMap bands_minoseff =
      PlotUtils::GetMinosEfficiencySystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_minoseff.begin(), bands_minoseff.end());

  for (auto band : error_bands) {
    std::vector<CVUniverse*> universes = band.second;
    for (auto universe : universes) universe->SetTruth(is_truth);
  }
     */
  return error_bands;
}

#endif  // Systematics_h
