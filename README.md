# RCF_Analysis
This repo holds the Matlab code used to convert RCF dose data into Protons/MeV.
The code is organized as follows:
### RCFanalysis.m
This is the main file that generates the response function of a designated RCF stack, reads the mean dose on the RCF stack, and converts that dose (Gy/sr) into Protons/MeV/sr. This code must be vastly editted to meet your needs, this simply reflects the proper way to use the code to generate what you want. The RCF stack configurtion must be changed to match your own individual needs. The dose data importing must also be matched to your own individual needs. The plotting must be changed to what you want to plot.

### packmaker.m
This code is a function that will create an output struct file to be used as the Result in the RCFanalysis.m code. This generates the response function for the input pack over the energy range that you set

### filmpack.m
This is found in packmaker.m
It is used to determine how far a proton of each energy step specified in packmaker.m will travel and in what layer it will deposit most of its energy according to the Bragg peak. 

### ExtractSRIM.m
This is found in filmpack.m
It is used to extract the dE/dx information from SRIM calculations for each layer of material. This must be done prior to running this code and the directory to your SRIM data must be set.

### deltaE.m
This code is a function that determines the range of proton energies that are deposited into each active layer of your film pack. This is used as &Delta;E in Protons/MeV/sr.
