from scipy.integrate import *
from pandas import *
from numpy import *
from pylab import *
from scipy import integrate
import math
from scipy.interpolate import InterpolatedUnivariateSpline

#parameter values
phi1 = 0.001 #interaction strength for p1
phi2 = 0.001 #interaction strength for p2
eps1 = 0.1 #transfer efficiency for p1
eps2 = 0.1 #transfer efficiency for p2
delta_c = 0.1 #consumer mortality
k_w = 0.02 #light attenuation by just water
k_p = 0.001 #attenuation due to producers
delta1 = 0.12 #death rate of producer 1
delta2 = 0.12 #death rate of producer 2
mu1 = 1 #resource affinity parameter for p1
mu2 = 1 #resource affinity parameter for p2
alpha_n = 1 #saturation of nutrients
alpha_i = 1 #saturation of light
S_R = 0.4 #supply rate

# time change (discrete time)
T = 36500.0 #days
delt = 1.0/240.0
t = np.linspace(0, T, int(T/delt)) #time array
#print(t)

# array  to store  the  solution if you were using Euler's method
N = np.zeros(len(t))
P1 = np.zeros(len(t))
P2 = np.zeros(len(t))
C = np.zeros(len(t))
N[0] = 500
P1[0] = 100
P2[0] = 100
C[0] = 100

#create an array with your starting values for odeint
y0 = np.array([N[0], P1[0],P2[0],C[0]])
print(y0)

zetas = np.linspace(0,100,50) #depth array
print(len(zetas))

#initial irradiance equation (seasonality)
I_naughts = (np.sin(((t)/365*2*math.pi-math.pi/2))+1)/2
plt.plot(t,I_naughts)

#only plot 1 year
plt.xlim([0,365])
plt.show()

# creating irradiance array
Is = [] #empty list
for I_naught in I_naughts:
    # accounting for initial values of phytoplankton
    I = I_naught * np.exp((-k_w*zetas)-(k_p*(P1[0]+P2[0])))
    Is.append(I)
light_att = np.array(Is) # make list into an array
lights = light_att.T #transpose

print(zetas.shape)
print(lights.shape)


