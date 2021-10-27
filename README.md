# Glitches
Programs to calculate glitches in single crystal x-ray optics


process_spectrum_norm.cpp - processing text logs from ESRF

nk_fit_rev_conf_par_v1.cpp - finding initial orientation of the optical element

nk_v1.py - simple calulation of glitches spectra

nk_analit_refine_and_error.py - refining the orientation of the optical element and calculating resulting error

conf_mul_scan_chi0_ang - example configuration file for nk_fit_rev_conf_par_v1

nk_fit_rev_conf_par_v2.cpp - new version of the orientation determination (indexing), much faster (for JSR'21 paper)

nk_fit_rev2.py - the same but in Python (much slower, just for educational purpose)

nk_calc_error_rev2.py - just calculate error

nk_analit_refine_and_error2.py - new refinement (for JSR'21 paper)

Glitches, plot at ω=0° and φ=0°

![alt text](plot_om0_phi0.gif?raw=true "") 

Glitches, reciprocal space at ω=0° and φ=0°

![alt text](reci_om0_phi0.gif?raw=true " ")

Glitches, plot at ω=4° and φ=4°

![alt text](plot_om4_phi4.gif?raw=true " ") 

Glitches, reciprocal space at ω=4° and φ=4°

![alt text](reci_om4_phi4.gif?raw=true " ")

