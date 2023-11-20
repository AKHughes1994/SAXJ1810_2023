#!/bin/bash

################################################################################################
# For this to work you need both HEASOFT and the Swift-XRT Python API installed and functional #
################################################################################################

# Fit each spectra individually
python3 fit_spectra_individual.py --emin 0.5 --emax 10.0 
python3 fit_spectra_individual.py --emin 1.0 --emax 10.0 

# Fit the joint spectra
python3 fit_spectra_full.py --emin 0.5 --emax 10.0 
python3 fit_spectra_full.py --emin 1.0 --emax 10.0 

# Fit the low-count spectra cash statistic epochs
python3 fit_spectra_lowcount.py --emin 0.5 --emax 10.0 
python3 fit_spectra_lowcount.py --emin 1.0 --emax 10.0 

# Fit the spectra with the 3-component fits
python3 fit_spectra_3comp.py
