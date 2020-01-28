# import all the necessary functions
import numpy as np
import random
import scipy.stats as sps
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import astropy.constants as const
import scipy.constants as const2

#make a binary class
class Binary:

    #define underlying distribution of semi-major axis(log-flat) and eccentricity(uniform) here
    def __init__(self):
        self.a = self.log_a()
        self.e = np.random.uniform(0.0,1.0)
        self.pi = const2.pi

    #function creating log-flat semi-major axis distribution
    def log_a(self):
        rand = np.random.uniform(0.0,1.0)
        power_value = 1+(3*rand)
        a = 10**(power_value)
        return a

    #randomly generate a number between 0 and 2*pi, with equal possibility
    def random_generator(self):
        guess = np.random.uniform(0.0,2*(self.pi))
        return guess

    
    #Newton-Raphson method
    def guess_E(self):
        M = self.random_generator()
        e = self.e
        Eguess = M
        old_guess = M
        while True:  
            old_guess = Eguess
            Eguess = (Eguess) - ((Eguess - (e*(np.sin(Eguess))) - M)/(1 - (e*(np.cos(Eguess)))))  
        
            if np.absolute(old_guess - Eguess) < (1.48*(10**-8)):
                E = Eguess
                break
    
            else:
                pass
        
        return E

    def true_anomaly(self):
        e = self.e
        E = self.guess_E()
        theta = 2*(np.arctan((np.sqrt((1+e)/(1-e)))*(np.tan(E/2))))
        return theta

    def calc_r(self):
        a = self.a
        e = self.e
        theta = self.true_anomaly()
        r = ((a*(1-(e**2)))/(1+(e*(np.cos(theta)))))
        return theta,r


    def pol_to_cart(self):
        theta,r = self.calc_r()
        x = r*(np.cos(theta))
        y = r*(np.sin(theta))
        z = 0.0    
        return x,y,z

#another binary class with log-normal semi-major axis distribution (all other parameters inherited from Binary class)
class AnotherBinary(Binary):
    def __init__(self):
        self.a = self.log_a()
        self.e = np.random.uniform(0.0,1.0)
        self.pi = const2.pi
    
    def log_a(self):
        rand =  np.random.normal(2.5,1.0)
        a = 10**(rand)
        return a

#change x and y paramters when changing observer angle
def change_observer_angle(x,y):
    phi = binary.random_generator()
    x2 = (x*(np.cos(phi))) - (y*(np.sin(phi)))
    y2 = (x*(np.sin(phi))) + (y*(np.cos(phi)))
    return x2,y2,phi
        
#further change x and z when changing the inclination 
def change_inclination(x2,z):
    sin_i = np.random.uniform(0.0,1.0)
    i = np.arcsin(sin_i)
    x3 = (x2*(np.cos(i))) - (z*(np.sin(i)))
    z2 = (x2*(np.sin(i))) + (z*(np.cos(i)))
    return x3,z2,i

#calculate the separation after putting in a, e, and all other random parameters
def separation(x3,z2):
    s = np.sqrt((x3**2)+(z2**2))
    return s

#run Monte Carlo simulations on both binary of different semi-major axis distribution and collect the statistical results
def collect_rand_obs(iters):
    count = 0
    sep = []
    sep2 = []
    a = []
    a2 = []

    while count < iters:
        binary = Binary()
        binary2 = AnotherBinary()
        aa = binary.a
        ab = binary2.a
        a.append(aa)
        a2.append(ab)
        x,y,z = binary.pol_to_cart()
        xx,yy,zz = binary2.pol_to_cart()
        x_primed,y_primed,phi = change_observer_angle(x,y)
        xx_primed,yy_primed,phi2 = change_observer_angle(xx,yy)
        x_double_primed,z_primed,i = change_inclination(x_primed,z)
        xx_double_primed,zz_primed,ii = change_inclination(xx_primed,zz)
        s = separation(x_double_primed,z_primed)
        s2 = separation(xx_double_primed,zz_primed)
        sep.append(s)
        sep2.append(s2)
        count += 1

    return sep,sep2,a,a2

#plot the statistical result obtained by Monte-Carlo simulation
def plot_CDF(parameter1,parameter2,iters):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    X1 = np.sort(parameter1)
    F1 = (np.array(range(iters))+1)/float(iters)
    X2 = np.sort(parameter2)
    F2 = (np.array(range(iters))+1)/float(iters)
    ax.plot(X1,F1,label = 'log flat distribution',linestyle = '-.')
    ax.plot(X2,F2,label = 'log normal distribution',linestyle = '-')
    #ax.scatter(X1, F1,label = 'log flat distribution', marker = '.',s = 1, color = 'red', alpha = 1)
    #ax.scatter(X2, F2,label = 'log normal distribution', marker = '.',s = 1, color = 'royalblue', alpha = 1)
    ax.set_xscale('log')
    plt.legend()
    plt.ylabel('Cumulative distribution')
    return X1, X2

#call the necessary functions defined above (output)
binary = Binary()
binary2 = AnotherBinary()
s_flat_distribution,s_normal_distribution,a,a2 = collect_rand_obs(100)
X1, X2 = plot_CDF(s_flat_distribution,s_normal_distribution,100)
#D, p = sps.ks_2samp(X1, X2)
#print(D)
#print(p)
#plt.axvline((10**2.5), linewidth = 1, color = 'black', linestyle = '-')
#plt.axvline((10**1.5), linewidth = 1, color = 'grey', linestyle = '--')
#plt.axvline((10**3.5), linewidth = 1, color = 'grey', linestyle = '--')
#plt.axvline((10**(0.5)), linewidth = 1, color = 'silver', linestyle = ':')
#plt.axvline((10**4.5), linewidth = 1, color = 'silver', linestyle = ':')
plt.xlabel('Separation(AU)')
plt.show()

