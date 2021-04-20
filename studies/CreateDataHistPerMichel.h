#include "studies/Study.h"
#include "event/MichelEvent.h"
#include "util/Categorized.h"
#include "event/CVUniverse.h"

#include <functional> //for std::function

class CreateDataHistPerMichel: public Study

{
  public:
    using reco_t = std::function<double(const CVUniverse&, const MichelEvent&, const int)>;
    CreateDataHistPerMichel(reco_t reco, const std::string& varName, const std::string& varUnits, const int nBins, const double minBin, const double maxBin, std::vector<CVUniverse*>& univs): Study(), fReco(reco)
    {
       std::string name = varName+"_data_";
       std::string title = varName+ " [" + varUnits + "]";
       dataHist = new HIST(varName.c_str(),varName.c_str(),nBins,minBin, maxBin, univs);
    } 



    void SaveOrDraw(TDirectory& outDir)
    {
       outDir.cd();
       dataHist->SyncCVHistos();
       dataHist->hist->Write();
       //TODO: You could do plotting here
    }

  private:
    using HIST = PlotUtils::HistWrapper<CVUniverse>;
    reco_t fReco;
    HIST* dataHist;  
    void fillSelected(const CVUniverse& univ, const MichelEvent& evt, const double weight) {
        for(size_t whichMichel = 0; whichMichel < evt.m_nmichels.size(); ++whichMichel)
      {
         (*dataHist).FillUniverse(&univ, fReco(univ, evt, whichMichel), weight);
      }
    }

};
