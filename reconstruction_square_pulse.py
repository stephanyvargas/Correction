from matplotlib import pyplot as plt
import numpy as np
plt.rcParams.update({'font.size': 13})


tau = 1.5
width = 1


def calculate_integral(dy, uy):
    I_dy = np.sum(dy*(1/101))
    I_uy = np.sum(uy*(9/900))
    print('Integral for deformed pulses:\n',
          'Droop integral:', I_dy, '\n', 'Undershoot integral:', I_uy)


#artificial square pulse algorithm
up_y = np.ones(101)
low_y = np.zeros(900)
square_pulse = np.concatenate((up_y, low_y))


#upper part of the pulse for the droop, lower for the undershoot
up_x = np.linspace(0,1.,101)
low_x = np.linspace(1.01,10.,900)


#deformed the pulse by a tau constant
dy = np.exp(-up_x/tau)
uy = (np.exp(-width/tau)-1)*np.exp(-(low_x-width)/tau)
x = np.linspace(0,10,1001)
y = np.concatenate((dy, uy))


#inital check for both integrals (droop and undershoot)
'''#between 0 to 10 the integrals are: (0 to 1) 0.787103603082587, (1 to 10) -0.7793227436228078
#between 0 to 100 the integrals are: (0 to 1) 0.787103603082587, (1 to 100) -0.7888281781410758'''
#calculate_integral(dy, uy)


#Algorithm to reconstruct by slices the pulse using single tau approximation
#width of time slice
dt = np.diff(x)[0]
print(dt)
A = (tau/dt) * (1-np.exp(-dt/tau))
#print(A) Amplitude at this point A= 0.996670739874859
S = 0
X0 = (1/A)*y[0]
X = [X0]

for i in range(1, len(x)):
    sj = X[i-1] + S*np.exp(-dt/tau)
    xj = (1/A)*y[i] + (dt*A/tau)*sj
    S = sj
    X.append(xj)


plt.plot(x, square_pulse, label='Square pulse')
#plt.plot(x,y, label=r'drooped pulse with $\tau=$ {}'.format(tau))
#plt.plot(x, X, label='Droop correction')
#plt.axhline(1/A)
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.legend(loc='best')
plt.grid(linestyle='dotted')
plt.show()


plt.plot(x, square_pulse-X)
plt.title('Difference: Original square pulse - Droop corrected pulse')
plt.xlabel('Time')
plt.ylabel('Difference (Original - Corrected pulse)')
plt.grid(linestyle='dotted')
#plt.show()
