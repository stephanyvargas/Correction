from matplotlib import pyplot as plt
import numpy as np

tau = 1.5
delta = 0.01
width = 1


def calculate_integral(dy, uy):
    I_dy = np.sum(dy*(1/100))
    I_uy = np.sum(uy*(9/900))
    print('Integral for deformed pulses:\n',
          'Droop integral:', I_dy, '\n', 'Undershoot integral:', I_uy)


#artificial square pulse algorithm
bl_y = [0]
up_y = np.ones(100)
low_y = np.zeros(900)
square_pulse = np.concatenate((bl_y, up_y, low_y))

#upper part of the pulse for the droop, lower for the undershoot
bl_x = [-0.01]
up_x = np.linspace(0,1,100)
low_x = np.linspace(1,10,900)

#deformed the pulse by a tau constant
dy = np.exp(-up_x/tau)
uy = (np.exp(-width/tau)-1)*np.exp(-(low_x-1)/tau)
x = np.linspace(-0.01,10,1001)
y = np.concatenate((bl_y, dy, uy))


#inital check for both integrals (droop and undershoot)
'''#between 0 to 10 the integrals are: (0 to 1) 0.787103603082587, (1 to 10) -0.7793227436228078
#between 0 to 100 the integrals are: (0 to 1) 0.787103603082587, (1 to 100) -0.7888281781410758'''
#calculate_integral(dy, uy)

S1 = y[0]*np.exp(-delta/tau)

#Algorithm to reconstruct by slices the pulse using single tau approximation
A = tau * (1-np.exp(-width/tau))
f0 = A
Sj = [0]
X0 = (1/A)*bl_y[0]
X = [X0]


for i in range(1, len(x)):
    #t0 =

    sj = X[i-1] + np.exp(-delta/tau)*Sj[i-1]
    Sj.append(sj)
    xj = (1/A)*x[i] + (A/tau)*sj
    X.append(xj*0.01)

#print('value of sj', Sj)
#print('value of xj', X, '\n')

plt.plot(x, square_pulse, label='original square pulse')
plt.plot(x,y, label=r'drooped pulse with $\tau=$ {}'.format(tau))
plt.plot(x, X, label='Droop correction')
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.legend(loc='best')
plt.grid(linestyle='dotted')
plt.show()
