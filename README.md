## Compton FastFit
This repository contains the code used in the paper "Computational method for the optimization of quasimonoenergetic Laser Compton X-ray sources for imaging applications" as well as code to reproduce all figures in the paper.

## About this code
This repository includes code that takes csv-formatted data of the scattering angles and energies of laser-Compton scattered (LCS) X-rays and then produces an interpolation map based on polynomial fitting to rapidly produce new spectra at different endpoint energies without any loss of local energy distribution information at each scattering angle.

All code is written using MATLAB R2021a and is intended to be used in the MATLAB GUI.

## How to use
This code is by no means optimized for user-friendliness and alterations, so please keep that in mind.

Properly running this code requires use of the "Run Section" option in MATLAB as well as many (too many) workspace variables.

We suggest starting with "getfigures.m," which tells you exactly which sections (SECS.) of the workhorse script, comptonfastfit.m, need to be run to produce each figure in the paper.

Sections 1-4 of comptonfastfit.m are ultimately the heart of Compton FastFit.
  Section 1: Load LCS X-ray spectral data
  Section 2: Convert the data into polynomial matrices
  Section 3: Create an interpolation map
  Section 4: Produce an interpolated LCS X-ray spectrum at a new energy
