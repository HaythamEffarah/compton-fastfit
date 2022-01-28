## Compton FastFit
This repository contains the code used in the paper "Computational method for the optimization of quasimonoenergetic Laser Compton X-ray sources for imaging applications" as well as code to reproduce all figures in the paper.

## About this code
This repository includes code that takes csv-formatted data of the scattering angles and energies of laser-Compton scattered (LCS) X-rays and then produces an interpolation map based on polynomial fitting to rapidly produce new spectra at different endpoint energies without any loss of local energy distribution information at each scattering angle.

All code is written using MATLAB R2021a and is intended to be used in the MATLAB GUI.

## How to use
Properly running this code requires use of the "Run Section" option in MATLAB as well as workspace variables.

We suggest starting with "getfigures.m," which tells you exactly which sections (SECS.) of the workhorse script, comptonfastfit.m, need to be run to produce each figure in the paper.

Sections 1-4 of comptonfastfit.m are ultimately the heart of Compton FastFit <br>

Section 1: Load LCS X-ray spectral data <br>
Section 2: Convert the data into polynomial matrices <br>
Section 3: Create an interpolation map <br>
Section 4: Produce an interpolated LCS X-ray spectrum at a new energy <br>

## Data

This code is intended to be run with data that was generated using laser-Compton scattering simulation software and can be found here with instructions on which directories to store the data: https://doi.org/10.6084/m9.figshare.19083344.v1
