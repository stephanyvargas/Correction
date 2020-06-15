from matplotlib import pyplot as plt
from scipy.special import factorial
import numpy as np
import random
plt.rcParams.update({'font.size': 13})


def calculate_integral(dy, uy):
    I_dy = np.sum(dy*(1/101))
    I_uy = np.sum(uy*(9/900))
    print('Integral for deformed pulses:\n',
          'Droop integral:', I_dy, '\n', 'Undershoot integral:', I_uy)

def plot_pulse(x,y,name):
    plt.plot(x, y, label='Original pulse')
    plt.xlabel('Time')
    plt.ylabel('Voltage')
    plt.title(f'Ideal {name} pulse')
    plt.legend(loc='best')
    plt.grid(linestyle='dotted')
    plt.savefig(f'/home/stephy/ICECUBE/undershoot/20200519/ideal_pulse_correction/png/pulse_plot_{name}.png')
    plt.savefig(f'/home/stephy/ICECUBE/undershoot/20200519/ideal_pulse_correction/svg/pulse_plot_{name}.svg')
    #plt.show()
    plt.clf()
    plt.cla()
    plt.close()

def ideal_square_pulse():#artificial square pulse algorithm
    x = np.linspace(0,10,1001)
    up_y = Amp * np.ones(101)
    low_y = np.zeros(900)
    y = np.concatenate((up_y, low_y))
    plot_pulse(x,y,'Square')
    return y

def ideal_gaussian_pulse():#artificial gaussian pulse algorithm
    x = np.linspace(0,10,1001)
    y = Amp * np.exp(-(x-2)**2/(0.2))
    plot_pulse(x,y,'Gaussian')
    return y


def ideal_2gaussian_pulse():#artificial gaussian pulse algorithm
    x = np.linspace(0,10,1001)
    y = Amp*0.05*np.exp(-(x-1.5)**2/(0.2)) + Amp * np.exp(-(x-3.5)**2/(0.2))
    plot_pulse(x,y,'Gaussian')
    return y

def ideal_poisson_pulse():#artificial gaussian pulse algorithm
    x = np.linspace(0,10,1001)
    y = Amp*np.exp(-5)*np.power(5, x)/factorial(x)
    plot_pulse(x,y,'Poisson')
    return y

def noisy_gaussian_pulse():#artificial gaussian pulse algorithm
    x = np.linspace(0,10,1001)
    y = []
    for i in range(0, len(x)):
        noise = random.randrange(-10, 10, 1)
        yi = Amp * (np.exp(-(x[i]-2)**2/(0.2))+ 0.007*noise)
        y.append(yi)
    y = np.asarray(y)
    plot_pulse(x,y,'Noisy Gaussian')
    return y

def droop_pulse(int_check):#manually deforme a square pulse by a tau constant
    up_x = np.linspace(0,1.,101)#upper part of the pulse for the droop, lower for the undershoot
    low_x = np.linspace(1.01,10.,900)
    dy = Amp * np.exp(-up_x/tau)
    uy = Amp * (np.exp(-width/tau)-1)*np.exp(-(low_x-width)/tau)
    x = np.linspace(0,10,1001)
    y = np.concatenate((dy, uy))
    return x, y

def droop(int_check, ideal_pulse):
    x = np.linspace(0,10,1001)
    y = ideal_pulse
    dt = np.diff(x)[0] #width of time slice
    A = (tau/dt) * (1-np.exp(-dt/tau))
    S = 0
    X = [A*y[0]]
    for i in range(1, len(x)):
        sj = y[i-1] + S*np.exp(-dt/tau)
        xj = A*y[i] - (dt*A*A/tau)*sj
        S = sj
        X.append(xj)
    if int_check == 'True':
        #inital check for both integrals (droop and undershoot)
        '''#between 0 to 10 the integrals are: (0 to 1) 0.787103603082587, (1 to 10) -0.7793227436228078
        #between 0 to 100 the integrals are: (0 to 1) 0.787103603082587, (1 to 100) -0.7888281781410758'''
        calculate_integral(dy, uy)
    return x, X

def correct(int_check, ideal_pulse, name):
    #Algorithm to reconstruct by slices the pulse using single tau approximation
    x, y = droop(int_check, ideal_pulse)
    dt = np.diff(x)[0] #width of time slice
    A = (tau/dt) * (1-np.exp(-dt/tau))
    S = 0
    X0 = ideal_pulse[0]
    X = [X0]
    #get the contribution of each term X[i-1] and S*np.exp(-dt/tau)
    T1 = [(1/A)*y[0]]
    T2 = [0]
    for i in range(1, len(x)):
        sj = X[i-1] + S*np.exp(-dt/tau)
        xj = (1/A)*y[i] + (dt*A/(tau))*sj
        t1 = (1/A)*y[i]
        t2 = (dt*A/tau)*sj
        S = sj
        X.append(xj)
        T1.append(t1)
        T2.append(t2)

    #Plot of the drooped pulse and the correction
    plt.plot(x[1:],y[1:], label=r'drooped pulse with $\tau=$ {}'.format(tau))
    plt.plot(x[1:], X[1:], label='Droop correction')
    plt.xlabel('Time')
    plt.ylabel('Voltage')
    plt.title('Ideal pulse correction')
    plt.legend(loc='best')
    plt.grid(linestyle='dotted')
    plt.savefig(f'/home/stephy/ICECUBE/undershoot/20200519/ideal_pulse_correction/png/pulse_correction_{name}.png')
    plt.savefig(f'/home/stephy/ICECUBE/undershoot/20200519/ideal_pulse_correction/svg/pulse_correction_{name}.svg')
    #plt.show()
    plt.clf()
    plt.cla()
    plt.close()

    #Plot of the contributions of each term for the correction
    plt.plot(x,T1, label=r'Term $\frac{1}{A}\,V_j$')
    plt.plot(x, T2, label=r'Term $\frac{A\,dt}{\tau}\,S_j$')
    plt.title('Contribution of each term to the correction')
    plt.xlabel('Time')
    plt.ylabel('Voltage')
    plt.legend(loc='best')
    plt.grid(linestyle='dotted')
    plt.savefig(f'/home/stephy/ICECUBE/undershoot/20200519/ideal_pulse_correction/png/pulse_each_term_{name}.png')
    plt.savefig(f'/home/stephy/ICECUBE/undershoot/20200519/ideal_pulse_correction/svg/pulse_each_term_{name}.svg')
    #plt.show()
    plt.clf()
    plt.cla()
    plt.close()

    #plot of the difference between the original pulse and the corrected
    plt.plot(x[1:], ideal_pulse[1:]-X[1:])
    plt.title('Original-Correction')
    plt.xlabel('Time')
    plt.ylabel('Difference (Original - Corrected pulse)')
    plt.grid(linestyle='dotted')
    plt.savefig(f'/home/stephy/ICECUBE/undershoot/20200519/ideal_pulse_correction/png/diff_plot_{name}.png')
    plt.savefig(f'/home/stephy/ICECUBE/undershoot/20200519/ideal_pulse_correction/svg/diff_plot_{name}.svg')
    #plt.show()
    plt.clf()
    plt.cla()
    plt.close()


def main():
    check_integral = 'False'
    type_pulse_sq = ideal_square_pulse()
    correct(check_integral, type_pulse_sq, 'Square')
    type_pulse_gaus = ideal_gaussian_pulse()
    correct(check_integral, type_pulse_gaus, 'Gaussian')
    type_pulse_2gaus = ideal_2gaussian_pulse()
    correct(check_integral, type_pulse_2gaus, '2Gaussian')
    #type_pulse = ideal_poisson_pulse()
    #correct(check_integral, type_pulse, 'Poisson')
    type_pulse_noise = noisy_gaussian_pulse()
    correct(check_integral, type_pulse_noise, 'Noisy_Gaussian')



if __name__ == "__main__":
    global tau
    global width
    global Amp
    tau = 1.5
    width = 1 #only for drooping manually an square pulse
    Amp = 1.5
    main()
