from matplotlib import pyplot as plt
import numpy as np

tau = 1
delta = 0.01
#upper part of the pulse for the droop, lower for the undershoot
up_x = np.linspace(0,1,100)
low_x = np.linspace(1,10,900)


#deformed pulse
dy = np.exp(-up_x/tau)
uy = (np.exp(-delta/tau)-1)*np.exp(-(low_x-1)/tau)
x = np.linspace(0,10,1000)
y = np.concatenate((dy,uy))


#integral of both parts
#between 0 to 10 the integrals are: 0.787103603082587 -0.7793227436228078
#between 0 to 100 the integrals are: 0.787103603082587 -0.7888281781410758
I_dy = np.sum(dy*(1/100))
I_uy = np.sum(uy*(9/900))


#Algorithm to reconstruct by slices the pulse
A = tau * (1-np.exp(-delta/tau))
f0 = A
Sj = [0]
X0 = (1/A)*uy[0]
X = [X0]

for i in range(1, len(low_x)):
    sj = X[i-1] + np.exp(-delta/tau)*Sj[i-1]
    Sj.append(sj)
    xj = (1/A)*low_x[0] + (A/tau)*sj
    X.append(xj*0.01)

#print('value of sj', Sj)
#print('value of xj', X, '\n')

plt.plot(low_x, X)
plt.plot(x,y)
plt.grid(linestyle='dotted')
plt.show()
