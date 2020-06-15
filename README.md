# Correction
Correction of the droop and undershoot behavior on pulses.
This scripts run on Python3.

#reconstruction.py file
Aims to correct pulses using real data from a High Voltage Board
 connected to a function generator (5us, 50 mV pulse)

#reconstruction_square_pulse.py file
Creates an ideal pulse that gets deformed by an arbitrary time constant tau.
Avaliable pulses are: 
-square, 
-gaussian, 
-2 gaussians, 
-poisson and 
-noisy gaussian.

The pulse amplitude and time cosntant tau, can be choosen (default 1.5).

The script gives back png and svg files.
The files correspond to:
-The original pulse shape
-The contribution of both, the current droop behavior and the cummulative undershoot from previous bins.
-The original pulse and the correction
-Error of the correction
