# ClusterECorr
ROOT-based GUI, [1,2], to extract corrections to cluster energies
It allows to open a root file and read a TH2 histogram in order to step through x-slices and fit the resulting y-projections
with a gaussian, a landau or a simplified Crystal Ball function. The mean/MPV (correction) is stored in a TH1 histogram.
When saving the extracted corrections, the TH1 distribution is smoothed and transferred to a second TH1 histogram with
binwidths spanning 1 MeV (expecting the x-axis of the original TH2 to be given in MeV). The empty bins are filled by 
interpolation.

## Usage
1. In the main function, add the path and name of the root file as well as the name of the TH2 histogram
   (Optionally: modify the global parameter "slicewidth")
2. Load the macro and run the main function in a root session
### GUI options
Skip(n) - Go to next slice without storing current mean/MPV
Save(m) - Go to next slice and store current mean/MPV
FitLandau(f) - Fit current y-projection with a Landau
FitGaus(g) - Fit current y-projection with a Gaus
FitCB(h) - Fit current y-projection with a Crystal Ball function
Number entry - This number, multiplied with the standard deviation of the current y-projection, determines the fit range
Save(s) - Creates the smoothed and interpolated final correction distributions and saves it into a root file
Exit(e) - Exists the program

## References
[1] https://root.cern
[2] https://root.cern.ch/root/htmldoc/guides/users-guide/WritingGUI.html
