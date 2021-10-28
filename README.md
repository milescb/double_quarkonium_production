# Uncorrelated Production of Double Quarkonium

The code presented here was used in the analysis of double quarkonium yields as part of an ongoing research project at the University of California, Davis as a part of the CMS collaboration. This research would not be possible without the UCD REU program and was funded through the NSF. 

## Summary of code

All code is written in `C++` and uses Cern's `ROOT` library. 

1. `Glauber_model` contains code for a standard Monte-Carlo Glauber model. Parameters were used to match data at <img src="http://www.sciweavers.org/tex2img.php?eq=%5Csqrt%7Bs%7D%20%3D%205.02%20%5C%2C%20%5Ctextrm%7BTeV%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="\sqrt{s} = 5.02 \, \textrm{TeV}" width="118" height="19" />. 
2. `Pythia_accept_eff` contains code for a `PYTHIA` simulation and analysis. The code is broken into two folders. `pythia_model` which contains the the pythia code and excecutables as well as appropriate `Makefiles`. This code is written to interface with `ROOT` and thus appropriate compliation flags must be used when compiling (see `Makefile`). The `analysis` folder contains analysis of the `PYTHIA` results done in `ROOT` to find values of acceptance and efficiency of the CMS detector at <img src="https://render.githubusercontent.com/render/math?math=\sqrt{s} = 5.02 TeV">, <img src="https://render.githubusercontent.com/render/math?math=\mathcal{L}_{int} = 1.7 nb^{-1}"> . 
3. Finally, `Quarkonium_production` contains the files used to create final estimates of uncorrelated yields of the <img src="https://render.githubusercontent.com/render/math?math=\Upsilon"> and <img src="https://render.githubusercontent.com/render/math?math=J/\psi"> mesons.
