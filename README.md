# Full MonteCarlo method icecube muon flux simulation

The file random_generaotr.py generates one muon with random position, direction, and energy according to the atmospheric muon flux. 

F_B_test.py accepts this random muon and propagates it till hit or energy below threshold. 

WMM.COF and geomag.py provides the simulated geomagnetic field model from https://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml

HTCondor.py submits this job to UW-Madison HTCondor cluster.
