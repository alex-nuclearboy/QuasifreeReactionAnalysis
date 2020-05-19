# Quasi-free reaction analysis

In order to determine the luminosity dependence on the excess energy Q, we used the quasi-elastic proton–proton scattering process.
In the proton–deuteron collisions, the protons from the beam scatter on the protons in the deuteron target and the neutrons from the deuteron play a role of spectators.

All files are distributed under the terms of the GNU General Public Licence Version 3.

## Selection criteria for quasi-free proton-proton reaction studies

Each of simulated event was weighted by the differential cross section for the quasi-free reaction which is a function of the scattering angle and the effective proton beam momentum. 
We computed it based on the SAID partial-wave program (http://gwdac.phys.gwu.edu/).

In order to select quasi-free proton-proton events the condition of exactly one charged particle in the Forward Detector (FD) and one particle in the Central Detector (CD) was applied. 
Then, the dominating background processes were subtracted. 
The events corresponding to the charged pions registered in the CD were subtracted based on energy deposited in the electromagnetic calorimeter (SEC) and thin plastic scintillator (PSB), while the background coming from elastic proton-deuteron scattering was subtracted applying the cut in polar angle.

## Software

### ROOT 
ROOT is an object-oriented framework for large scale data analyses based on C++.
It was developed for experiments at CERN, in which an impressive amount of data has to be processed. It was originally designed for particle physics data analysis.
ROOT has an extensive library with modules in physics, mathematics, statistics, histograms and graphics. It contains statistics tools, which can be used for data analysis.

### RootSorter
RootSorter is based on the ROOT data analysis framework. The framework is organized in different parts to handle different tasks like, decoding of data, calibration of the individual detector, track and $
The experimental data which come directly from the electronics are stored in a Hit Bank Raw. The simulated data are stored in Hit Bank MC. The experimental data goes through the calibration whereas the s$
Then the track reconstruction, energy reconstruction and particle identification are done taking calibrated data or the filtered data. The reconstructed kinematic information are written to a ROOT format$

## Needed environment variables:

path where ROOT is installed

    ROOTSYS

path where RootSorter is instaled

    ROOTSORTERSYS

path where preselected data from experiment are located

    PRESEL

path where the MC simulation results are stored

    WMC_DATA

path to store analysis results for raw data

    OUTPUT_DATA

path to store analysis results for MC data

    OUTPUT_MC

## Running analysis

    ./runAnalysis-mc.sh <reaction>

for analysis of simulation results or

    ./runAnalysis-data.sh

for the analysis of raw data.
