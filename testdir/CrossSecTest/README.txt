MINERvA Analysis Toolkit "Minimum Adoption" Example

"Minimum adoption" means that it only uses the two essential MAT tools:
Universes and HistWrappers. For an "full adoption" example that additionally
makes use of Cuts, MacroUtil, and Variable, refer to the example in
../MAT_Tutorial/.

The runEventLoop.C is the main script here. It loops through a single MC CCNuPi
tuple, makes a few very basic CC inclusive cuts, and fills a neutrino energy
histogram + several standard systematic error bands. At the end, it plots the
enu distribution with error bands along with a separate fractional error
summary.

CVUniverse.h defines the interface with the anatuple and the event-by-event
weight. Most of its core branch accessors are defined in its grandparent class
PlotUtils/BaseUniverse.*.

The load.C file just hides some benign PlotUtils compiler warnings.

Run it like this:
root -l -b -q load.C+ runEventLoop.C+

