import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import math
import sys

#parameter values
phi1 = 0.2 #interaction strength for p1
phi2 = 0.1 #interaction strength for p2
eps1 = 0.4 #transfer efficiency for p1
eps2 = 0.3 #transfer efficiency for p2
delta_c = 0.1 #consumer mortality
k_w = 0.02 #light attenuation by just water
k_p = 0.001 #attenuation due to producers
delta1 = 0.1 #death rate of producer 1
delta2 = 0.1 #death rate of producer 2
mu1 = 3.0 #resource affinity parameter for p1
mu2 = 2.0 #resource affinity parameter for p2
alpha_n = 0.01 #saturation of nutrients
alpha_i = 0.01 #saturation of light
S_R = 0.4 #supply rate

# time change (discrete time)
T = 365.0 #days
dt = 0.04
#time array
t = np.linspace(0, T, int(T/dt)) 

#depth array
zetas = np.linspace(0,100,50) 

# array  to store  the  solution if you were using Euler's method
N = np.zeros((len(t),len(zetas)))
P1 = np.zeros((len(t),len(zetas)))
P2 = np.zeros((len(t),len(zetas)))
C = np.zeros((len(t),len(zetas)))
N[0,:] = 0.4
P1[0,:] = 0.5
P2[0,:] = 0.5
C[0,:] = 0.5
#print(N.shape)

#light function
def irradiance_function(zetas, P1,P2, k, l, dep_size):
    I = I_naughts[k] * np.exp(-(k_w*zetas[l]+np.sum(k_p*(P1[:l]+P2[:l])*dep_size)))
    return I

#light array
l = 0
k = 25
dep_size = 2 #depth (meters) in steps
I_naughts = (np.sin(((t)/365*2*math.pi-math.pi/2))+1)/2
light_list = []
for I_naught in I_naughts:
    loop = irradiance_function(zetas, P1, P2, k, l, dep_size)
    light_list.append(loop)
light_array = np.array(light_list)
print('light array shape from I function', light_array.shape)

# creating irradiance array for figure
Is = [] #empty list
for I_naught in I_naughts:
    # accounting for initial values of phytoplankton
    I = I_naught * np.exp((-k_w*zetas)-(k_p*(P1[0]+P2[0])))
    Is.append(I)
light_att = np.array(Is) # make list into an array
lights = light_att.T #transpose
print('irradiance array from contour figure', lights.shape)
 
#building 2 dimensional arrays for nutrients, producer, and consumer groups
for i in  range(1, len(t)):
    for j in range(0,len(zetas)):
        I = I_naughts[i] * np.exp(-(k_w*zetas[j]+np.sum(k_p*(P1[0,:j]+P2[0,:j])*j)))
        mu1_growth = np.minimum((N[i-1,j]/(N[i-1,j]+(mu1/alpha_n))), (I/(I+(mu1/alpha_i))))
        mu2_growth = np.minimum((N[i-1,j]/(N[i-1,j]+(mu2/alpha_n))), (I/(I+(mu2/alpha_i))))
        dNdt = S_R - (mu1*mu1_growth*P1[i-1,j]) - (mu2*mu2_growth*P2[i-1,j])
        dP1dt = (mu1*mu1_growth*P1[i-1,j]) - (phi1*P1[i-1,j]*C[i-1,j]) - (delta1*P1[i-1,j])
        dP2dt = (mu2*mu2_growth*P2[i-1,j]) - (phi2*P2[i-1,j]*C[i-1,j]) - (delta2*P2[i-1,j])
        dCdt = (eps1*phi1*P1[i-1,j]*C[i-1,j]) + (eps2*phi2*P2[i-1,j]*C[i-1,j]) - (delta_c*C[i-1,j])
        N[i] = N[i-1,j-1] + dNdt * dt
        P1[i] = P1[i-1,j-1] + dP1dt * dt
        P2[i] = P2[i-1,j-1] + dP2dt * dt
        C[i] = C[i-1,j-1] + dCdt * dt

# print shape of arrays for nutrients, producer, and consumer groups
print('N shape = ', N.shape)
print('P1 shape = ', P1.shape)
print('P2 shape = ', P2.shape)
print('C shape = ', C.shape)
            
fig, axs = plt.subplots(1,3,figsize=(16,6))

#contour plot of seasonal irradiance
cb = axs[0].contourf(t,zetas,lights,cmap = 'hot')
axs[0].invert_yaxis()
axs[0].set_xlabel('Time (Days)', color = 'k')
axs[0].set_ylabel('Depth (Meters)', color = 'k')
axs[0].set_title('Irradiance')
divider = make_axes_locatable(axs[0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cb, cax=cax, orientation='vertical')

#plotting nutrients
axs[1].plot(t,N, color = 'red')
axs[1].set_xlabel('Time (days)', color = 'k')
axs[1].set_ylabel('Concentration (mmol C/ m$^{3}$)', color = 'k')
axs[1].set_title('Nutrients')
axs[1].legend(['Nutrients'], loc = 'best')

#plotting biomass
pro1 = axs[2].plot(t,P1, '--', color = 'orange')
pro2 = axs[2].plot(t,P2, 'blue')
con = axs[2].plot(t, C, color = 'grey', linestyle = 'dotted')
axs[2].set_xlabel('Time (days)', color = 'k')
axs[2].set_ylabel('Biomass (mmol C/ m$^{3}$ day)', color = 'k')
axs[2].set_title('Biomass')
# axs[2].legend([pro1,pro2,con], ['Producer 1','Producer 2', 'Consumer'])
custom_lines = [Line2D([0], [0], color='orange', linestyle = 'dashed', lw=1),
                Line2D([0], [0], color='blue',   lw=1),
                Line2D([0], [0], color='grey', linestyle = 'dotted',  lw=1)]
axs[2].legend(custom_lines, ['Producer 1', 'Producer 2', 'Consumer'])

plt.tight_layout()
plt.show()
