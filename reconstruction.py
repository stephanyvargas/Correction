from matplotlib import pyplot as plt
from random import randint
import numpy as np
import glob
import sys
plt.rcParams.update({'font.size': 13})

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

def random_wfcorrection(Batch, Directory_board, Directory_temp):
    #random_wf = randint(0,99)
    waveform = '/home/stephy/ICECUBE/database/B%s_HVB_%s/%s/*.CSV' %(Batch, Directory_board, Directory_temp)
    filename = sorted(glob.glob(waveform))
    random_wf = randint(0,len(filename))
    random_waveform = filename[random_wf]
    data = np.loadtxt(fname=random_waveform, delimiter=',', skiprows=25)
    time = data[:,0]
    volt = data[:,1]
    volts = volt - np.mean(volt[0:100])
    dt = np.diff(time)[0]
    mask = volts >= 0.9*max(volts) #begin the reconstruction at 90% of the maximum ->trigger
    droop_trigger = time[mask]
    time_trigger = droop_trigger[0]
    return time, volts, dt, time_trigger

def mean_wfcorrection(Batch, Directory_board, Directory_temp):
    waveform = '/home/stephy/ICECUBE/database/B%s_HVB_%s/%s/*.CSV' %(Batch, Directory_board, Directory_temp)
    filename = sorted(glob.glob(waveform))
    NumF = len(filename)
    Times = []
    Voltages = []
    for i in range(0, NumF):
        data = np.loadtxt(fname=filename[i], delimiter=',', skiprows=25)
        time = data[:,0]
        volt = data[:,1]
        volts = volt - np.mean(volt[0:100])
        Times.append(time)
        Voltages.append(volts)

    Times = np.asarray(Times)
    Voltages = np.asarray(Voltages)
    TIME = sum(Times)/len(Times)
    VOLT = sum(Voltages)/len(Voltages)
    TIME = np.asarray(TIME)
    VOLT = np.asarray(VOLT)

    dt = np.diff(TIME)[0]
    mask = VOLT >= 0.9*max(VOLT) #begin the reconstruction at 90% of the maximum ->trigger
    droop_trigger = time[mask]
    time_trigger = droop_trigger[0]
    return TIME, VOLT, dt, time_trigger

#time, volts, dt, time_trigger = random_wfcorrection(Batch, Directory_board, Directory_temp)
time, volts, dt, time_trigger = mean_wfcorrection(Batch, Directory_board, Directory_temp)

'''Recover the droop time constant from the temperature
and the database values. The time and the tau are in
seconds.'''
mask1 = batch_dt == Batch
mask2 = board_dt == Directory_board
mask = np.logical_and(mask1, mask2)
p_0 = p0[mask]
p_1 = p1[mask]
p_2 = p2[mask]
p_3 = p3[mask]
p_4 = p4[mask]
p_5 = p5[mask]


tau_2 = (p_0 + p_1/(1+np.exp(-temp/p_2)))*1E-6
tau = (p_3 + p_4/(1+np.exp(-temp/p_5)))*1E-6
mask_trigger = time >= time_trigger
Y = volts[mask_trigger]
T = time[mask_trigger]


def chris_correction(tau, dt, Y, T):
    A = (tau/dt) * (1-np.exp(-dt/tau))
    S = 0
    X0 = (1/A * Y[0])
    X = [X0]
    for j in range(1, len(T)):
        sj = X[j-1] + S*np.exp(-dt/tau)
        xj = (1/A) * Y[j] + (A*dt)/tau * sj
        S = sj
        X.append(xj)
    return X

def stephy_correction(tau, dt, Y, T):
    A = np.max(Y)
    B = np.min(Y)
    mask = Y == np.min(Y)
    width = T[mask]

    X = []
    time = T[0]
    j = 0
    #Droop Correction
    while time < width[0]:
        xj = Y[j] + A * (1 - np.exp(-(T[j]-T[0])/tau))
        X.append(xj)
        time = T[j]
        j += 1
    #Undershoot correction
    for i in range(j, len(T)):
        xi = Y[i] - B * (np.exp(-(T[i]-width[0])/tau))
        X.append(xi)
        time = T[i]
    return X

Xc = chris_correction(tau, dt, Y, T)
Xs = stephy_correction(tau, dt, Y, T)

TAU = tau[0]*1E6
Xc = np.asarray(Xc)
Xs = np.asarray(Xs)
#plt.title(f'Correction on HV Board Number {Directory_board}, Batch {Batch}, Temperature {temp}')
plt.plot(time*1E6, volts*1000, label=fr'Drooped pulse with $\tau=$ {TAU:.2f} $\mu$s')
plt.plot(T*1E6, Xc*1000, label='Correction (Chris)')
#plt.plot(T*1E6, Xs*1000, label='Correction (Amplitude dependent)')
#plt.plot(T*1E6, (Xc/Xs))#, label='Ratio (Chris/Steph)')
#plt.xlim(T[0]*1E6,10)
#plt.ylim(0.9,1.1)
plt.xlabel(r'Time [$\mu$s]')
plt.ylabel('Voltage [mV]')
#plt.ylabel('Ratio (Chris/Steph)')
plt.legend(loc='best')
plt.grid(linestyle='dotted')
plt.show()

#fj = -(1/tau)*A^2*np.exp(-(j-1)/tau)
