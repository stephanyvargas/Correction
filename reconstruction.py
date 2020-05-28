from matplotlib import pyplot as plt
from random import randint
import numpy as np
import glob
import sys


'''This code uses the single tau reconstruction method.
First is to load the constant table database. The important
constants for this algorithm are p0, p1, p2. The other constants
are also loaded for verification afterwards'''

droop_under_file = '/home/stephy/ICECUBE/undershoot/20191203/constant_table/DroopUnder.txt'
droopunder = np.loadtxt(fname=droop_under_file)
board = droopunder[:,0]
batch = droopunder[:,6]
alpha = droopunder[:,1]
beta = droopunder[:,3]

droop_temp_file = '/home/stephy/ICECUBE/undershoot/20191203/constant_table/TauDroop_Temp.txt'
drooptemp =np.loadtxt(fname=droop_temp_file)
board_dt = drooptemp[:,0]
batch_dt = drooptemp[:,8]
p0 = drooptemp[:,1]
p1 = drooptemp[:,3]
p2 = drooptemp[:,5]

under_temp_file = '/home/stephy/ICECUBE/undershoot/20191203/constant_table/TauUnder_Temp.txt'
undertemp = np.loadtxt(fname=under_temp_file)
board_ut = undertemp[:,0]
batch_ut = undertemp[:,8]
p3 = undertemp[:,1]
p4 = undertemp[:,3]
p5 = undertemp[:,5]


'''Load a random Waveform for a given temperature
for a given HVB, each measurement has 100 waveforms
This reconstruction is done on an individual measurement.
The correction begins at 90% of the pulse height.'''
Batch = int(sys.argv[1])
Directory_board = int(sys.argv[2])
Directory_temp = sys.argv[3]
temp = float(sys.argv[4])
random_wf = randint(0,99)
waveform = '/home/stephy/ICECUBE/database/B%s_HVB_%s/%s/*.CSV' %(Batch, Directory_board, Directory_temp)
filename = sorted(glob.glob(waveform))
random_waveform = filename[random_wf]
data = np.loadtxt(fname=random_waveform, delimiter=',', skiprows=25)
time = data[:,0]
volt = data[:,1]
volts = volt - np.mean(volt[0:100])
dt = np.diff(time)[0]
mask = volts >= 0.9*max(volts) #begin the reconstruction at 90% of the maximum ->trigger
droop_trigger = time[mask]
time_trigger = droop_trigger[0]
#print(time_trigger, dt)


'''Recover the droop time constant from the temperature
and the database values. The time and the tau are in
seconds.'''
mask1 = batch_dt == Batch
mask2 = board_dt == Directory_board
mask = np.logical_and(mask1, mask2)
p_0 = p0[mask]
p_1 = p1[mask]
p_2 = p2[mask]
tau = (p_0 + p_1/(1+np.exp(-temp/p_2)))*1E-6
A = tau * (1-np.exp(-dt/tau))
print(A)
A = A[0]
print(tau, A)
Mask = time >= time_trigger
Y = volts[Mask]
T = time[Mask]
X = [1/A * Y[0]]
Sj = [0]

for j in range(1, len(Y)):
    S = X[j-1] + np.exp(-dt/tau) * Sj[j-1]
    Xj = 1/A *Y[j] + (A*A*dt)/tau * (X[j-1] + np.exp(-dt/tau)*Sj[j-1])
    X.append(Xj)
    Sj.append(S)

X = np.asarray(X)
print(Sj[0], Sj[1], X[0])
plt.plot(time, volts)
#plt.plot(T,Sj)
plt.show()

#fj = -(1/tau)*A^2*np.exp(-(j-1)/tau)
