#ifndef VARIABLE_H
#define VARIABLE_H

#include "util/SafeROOTName.h"
#include "PlotUtils/VariableBase.h"
#include "util/Categorized.h"

class Variable: public PlotUtils::VariableBase<CVUniverse>
{
  private:
    typedef PlotUtils::HistWrapper<CVUniverse> Hist;
  public:
    template <class ...ARGS>
    Variable(ARGS... args): PlotUtils::VariableBase<CVUniverse>(args...)
    {
    }

    //TODO: It's really silly to have to make 2 sets of error bands just because they point to different trees.
    //      I'd rather the physics of the error bands remain the same and just change which tree they point to.
    void InitializeMCHists(std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
                           std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands)
    {
      //For example only.  Don't actually use GENIELabels as your backgrounds!
      //You'd end up with a very model-dependent result, and Luke Pickering
      //would frown on your paper ;)
      std::map<int, std::string> GENIELabels = {{1, "QE"},
                                                {8, "2p2h"},
                                                {2, "RES"},
                                                {3, "DIS"}};
 
      m_bestPionByGENIELabel = new util::Categorized<Hist, int>((GetName() + "_by_GENIE_Label").c_str(),
                                                                GetName(), GENIELabels,
                                                                GetNBins(), GetBinVec(), mc_error_bands);

      efficiencyNumerator = new Hist((GetName() + "_efficiency_numerator").c_str(), GetName().c_str(), GetNBins(), GetBinVec(), mc_error_bands);
      efficiencyDenominator = new Hist((GetName() + "_efficiency_denominator").c_str(), GetName().c_str(), GetNBins(), GetBinVec(), truth_error_bands);
    }

    //Histograms to be filled
    util::Categorized<Hist, int>* m_bestPionByGENIELabel;
    Hist* dataHist;  
    Hist* efficiencyNumerator;
    Hist* efficiencyDenominator;

    void InitializeDATAHists(std::vector<CVUniverse*>& data_error_bands)
    {
	std::vector<double> bins = GetBinVec();
        const char* name = GetName().c_str();
  	dataHist = new Hist(Form("_data_%s", name), name, GetNBins(), bins, data_error_bands);
 
    }

    void Write(TFile& file)
    {
      SyncCVHistos();

      file.cd();

      m_bestPionByGENIELabel->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?
                                    });
     
      if (dataHist->hist) {
		dataHist->hist->SetDirectory(&file);
		dataHist->hist->Write();
      }

      if(efficiencyNumerator)
      {
        efficiencyNumerator->hist->SetDirectory(&file); //TODO: Can I get around having to call SetDirectory() this many times somehow?
        efficiencyNumerator->hist->Write();
      }

      if(efficiencyDenominator)
      {
        efficiencyDenominator->hist->SetDirectory(&file);
        efficiencyDenominator->hist->Write();
      }
    }

    //Only call this manually if you Draw(), Add(), or Divide() plots in this
    //program.
    //Makes sure that all error bands know about the CV.  In the Old Systematics
    //Framework, this was implicitly done by the event loop.
    void SyncCVHistos()
    {
      m_bestPionByGENIELabel->visit([](Hist& categ) { categ.SyncCVHistos(); });
      if(dataHist) dataHist->SyncCVHistos();
      if(efficiencyNumerator) efficiencyNumerator->SyncCVHistos();
      if(efficiencyDenominator) efficiencyDenominator->SyncCVHistos();
    }
};

#endif //VARIABLE_H
