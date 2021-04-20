# CC Inclusive Cross Section Example for MINERvA 101 2021
An example of a charged current (CC) inclusive cross section analysis using the MINERvA Analysis Toolkit (MAT).  Starts from the official MINERvA data preservation "anaTuple" .root files, produced by MasterAnaDev (MAD), and produces all histograms needed to extract a cross section.  The basic tutorial does not attempt to constrain backgrounds with sidebands.  Based on Ben Messerly's examples of using the MINERvA Analysis Toolkit and work by Mehreen Sultana on her thesis analysis.

## Installation Instructions
1. Install [dependencies](#Dependencies)
2. Make a working directory.  Mine is called `MINERvA101_2021/`
3. Download the complete source from [MinervaExpt](https://github.com/MinervaExpt/MAT_IncPions): `git clone https://github.com/MinervaExpt/MAT_IncPions.git src #Makes a src subdirectory for you`
4. Make a build directory: `mkdir opt && cd opt && mkdir build && cd build #opt for optimized build as opposed to debug build`
5. Run cmake to generate a build system: ``cmake ../../src -DCMAKE_INSTALL_PREFIX=`pwd`/.. -DCMAKE_BUILD_TYPE=Release``
6. Compile and install: `make install`

For debugging, produce a `debug` build by creating a second build directory and replacing step 5 with: ``cmake ../../src -DCMAKE_INSTALL_PREFIX=`pwd`/.. -DCMAKE_BUILD_TYPE=Debug``.  This will turn on `assert()`s and include debug symbols that `gdb` and `valgrind` use to tell you line numbers of memory errors.  You probably want to rebuild each dependency in debug mode too to maximize helpfulness of the debugger.

## Usage
TODO

## Dependencies
0. [git](https://git-scm.com/downloads): version control system.  **You probably already have this**
1. [CMake 2.8](https://cmake.org/install/): build system generator for lots of operating systems.  **You probably already have this**
2. [ROOT](https://root.cern/install/): object-oriented toolkit for high energy physics analysis.  Make sure to enable at least xrootd, kerberos, and Minuit for this tutorial.  **Already installed on MINERvA GPVMs**
3. [PlotUtils](https://cdcvs.fnal.gov/redmine/projects/minerva-sw/repository/show/AnalysisFramework/Ana/PlotUtils): the MINERvA Analysis Toolkit and associated systematics tools (TODO: via github)
4. [UnfoldUtils](https://cdcvs.fnal.gov/redmine/projects/minerva-sw/repository/show/AnalysisFramework/Ana/UnfoldUtils): MINERvA's fork of [RooUnfold](https://gitlab.cern.ch/RooUnfold/RooUnfold) with a compatibility layer for embedding systematics in histograms and MINERvA-specific systematics tweaks (TODO: via github)
5. [MParamFiles](https://cdcvs.fnal.gov/redmine/projects/minerva-sw/repository/show/AnalysisFramework/MParamFiles): additional reweighting and calibration parameters for MINERvA systematics

TODO: Write an `ExternalProject_add()`-based installer
