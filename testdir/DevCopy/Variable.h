#ifndef VARIABLE_H
#define VARIABLE_H

#ifndef __CINT__ //PlotUtils/VariableBase.h uses std::function which is from c++11
#include "SafeROOTName.cpp" //TODO: This is a very bad idea
#endif //__CINT__
#include "PlotUtils/VariableBase.h"
#include "Categorized.h"
namespace dont
{
  //For example only.  Don't actually use GENIELabels as your backgrounds!
  //You'd end up with a very model-dependent result, and Luke Pickering
  //would frown on your paper ;)
  std::map<int, std::string> GENIELabels = {{1, "QE"},
                                            {8, "2p2h"},
                                            {2, "RES"},
                                            {3, "DIS"}};
}

class Variable: public PlotUtils::VariableBase<CVUniverse>
{
  private:
    typedef PlotUtils::HistWrapper<CVUniverse> Hist;
  public:
    #ifndef __CINT__ //For variadic template parameter ARGS
    template <class ...ARGS>
    Variable(ARGS... args): PlotUtils::VariableBase<CVUniverse>(args...)
    {
    }

    #endif //__CINT__
    //Very dissappointed in Ben for not allowing me to pass universes to the constructor.
    //This is not a solution I like :(
    void InitializeMCHists(const std::map<std::string, std::vector<CVUniverse*>>& error_bands)
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
                                                                GetNBins(), GetBinVec(), error_bands);
    }
    #ifndef __CINT__
    util::Categorized<Hist, int>* m_bestPionByGENIELabel;
    #endif //__CINT__
    Hist* dataHist;  
    void InitializeDATAHists(const std::map<std::string, std::vector<CVUniverse*>>& data_error_bands)
    {
	std::vector<double> bins = GetBinVec();
        const char* name = GetName().c_str();
        //dataHist = new Hist(Form("_data_%s", name), name, GetNBins(), bins, universe);
  	dataHist = new Hist(Form("_data_%s", name), name, GetNBins(), bins);
 
    }

    void Write(TFile& file)
    {
      SyncCVHistos();

      file.cd();

      #ifndef __CINT__ //For labmda function [](auto& categ){}
      m_bestPionByGENIELabel->visit([&file](Hist& categ)
                                    {
                                      categ.hist->SetDirectory(&file);
                                      categ.hist->Write(); //TODO: Or let the TFile destructor do this the "normal" way?
                                    });
     
      if (dataHist->hist) dataHist->hist->Write();
      #endif //__CINT__
    }

    //Only call this manually if you Draw(), Add(), or Divide() plots in this
    //program.
    //Makes sure that all error bands know about the CV.  In the Old Systematics
    //Framework, this was implicitly done by the event loop.
    void SyncCVHistos()
    {
      #ifndef __CINT__ //For labmda function [](auto& categ){}
      m_bestPionByGENIELabel->visit([](Hist& categ) { categ.SyncCVHistos(); });
      //dataHist->SyncCVHistos();
      #endif //__CINT__
    }

//Write a destructor
//  ~Variable()
//  {
//    delete dataHist;
//    dataHist = NULL;
//   }
};

#endif //VARIABLE_H
