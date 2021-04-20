#include "TString.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "Cintex/Cintex.h"

#include <iostream>

void load()
{
#ifdef __CINT__
  cout << "Run load.C+" << endl;
  exit(1);
#else
  gSystem->SetAclicMode(TSystem::kDebug);

  // MnvH1D hides approximately everything, so just turn off the pages
  // of compiler warnings. It would have been easier to do this by
  // using SetFlagsDebug(), but those flags get put before the default
  // settings in the compile line, and so the default settings win
  TString makeSharedLib(gSystem->GetMakeSharedLib());
  makeSharedLib.ReplaceAll("-Woverloaded-virtual", "-Wno-overloaded-virtual");
  gSystem->SetMakeSharedLib(makeSharedLib);

#endif
}
