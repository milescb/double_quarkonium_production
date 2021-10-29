# Uncorrelated Production of Double Quarkonium

The code presented here was used in the analysis of double quarkonium yields as part of an ongoing research project at the University of California, Davis as a part of the CMS collaboration. This research would not be possible without the UCD REU program and was funded through the NSF. 

## Paper and Poster

Both a paper (`quarkonium_production_paper.pdf`) containing a full summary of the research and a poster (`quarkonium_poster.pdf`) summarizing the research findings are included in the repo. The poster was presented as part of the CEU program at the 2021 DNP conference. 

## Summary of code

All code is written in `C++` and uses Cern's `ROOT` library. 

1. `Glauber_model` contains code for a standard Monte-Carlo Glauber model. Parameters were used to match data at center-of-mass energy 5.02 TeV. 
2. `Pythia_accept_eff` contains code for a `PYTHIA` simulation and analysis. The code is broken into two folders. `pythia_model` which contains the the pythia code and excecutables as well as appropriate `Makefiles`. This code is written to interface with `ROOT` and thus appropriate compliation flags must be used when compiling (see `Makefile`). The `analysis` folder contains analysis of the `PYTHIA` results done in `ROOT` to find values of acceptance and efficiency of the CMS detector at 5.02 TeV, integrate luminosity of 1.7/nb. 
3. Finally, `Quarkonium_production` contains the files used to create final estimates of uncorrelated yields of the Upsilon and J/psi mesons.

