import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import math
import sys

# time change (discrete time)
T = 5475.0 #days
dt = 1
#time array
t = np.linspace(0, T, int(T/dt))

#first 5 years = 0-1825 days
#second 5 years = 1826 - 3650 days
#third 5 years = 3651 - 5475 days

supply_rate = np.empty(3650)
supply_rate.fill(0.04)

year11 = np.empty(365)
year11.fill(0.08)
supply_rate = np.append(supply_rate,year11)

year12 = np.empty(365)
year12.fill(0.11)
supply_rate = np.append(supply_rate,year12)

year13 = np.empty(365)
year13.fill(0.14)
supply_rate = np.append(supply_rate,year13)

year14 = np.empty(365)
year14.fill(0.17)
supply_rate = np.append(supply_rate,year14)

year15 = np.empty(365)
year15.fill(0.20)
supply_rate = np.append(supply_rate,year15)

k_w = 0.02 #light attenuation by just water
k_p = 0.001 #attenuation due to producers

#depth array
deltaz = 0.5
depth = 100
zetas = np.linspace(0,depth,int(depth/deltaz)) 

# array  to store  the  solution if you were using Euler's method
N = np.zeros((len(t),len(zetas)))
P1 = np.zeros((len(t),len(zetas)))
P2 = np.zeros((len(t),len(zetas)))
C = np.zeros((len(t),len(zetas)))
light = np.zeros((len(t),len(zetas)))
N[0,:] = 0.4
P1[0,:] = 0.5
P2[0,:] = 0.5
C[0,:] = 0.5

#initial light irradiance (constant)
I_naughts = np.empty(len(t))

#inital light irradiance (seasonal)
I_max = 1000
I_min = 300

i = 0
for a in t:
    if a<=1825.0:
        I_naughts[i] = I_min
    else:
        I_naughts[i] = I_min + (I_max-I_min)/2.*(np.sin(((a)/365*2*np.pi-np.pi/2))+1)
    i = i + 1
        
#parameter values
phi1 = 0.1 #interaction strength for p1
phi2 = 0.1 #interaction strength for p2
eps1 = 0.3 #transfer efficiency for p1
eps2 = 0.3 #transfer efficiency for p2
delta_c = 0.001 #consumer mortality
k_w = 0.02 #light attenuation by just water
k_p = 0.001 #attenuation due to producers
delta1 = 0.001 #death rate of producer 1
delta2 = 0.001 #death rate of producer 2
mu1 = 1.0 #resource affinity parameter for p1
mu2 = 0.8 #resource affinity parameter for p2
alpha_n1 = 0.14 #saturation of nutrients for producer 1
alpha_i1 = 0.017 #saturation of light for producer 1
alpha_n2 = 0.14 #saturation of nutrients for producer 2
alpha_i2 = 0.014 #saturation of light for producer 2
 
#building 2 dimensional arrays for nutrients, producer, and consumer groups
total_time_arr1 = []
total_zeta_arr1 = []
total_color_arr1 = []
total_time_arr2 = []
total_zeta_arr2 = []
total_color_arr2 = []
for i in range(1, len(t)):
    #print(t[i])
    for j in range(0,len(zetas)):
        I = I_naughts[i] * np.exp(-(k_w*zetas[j]*(deltaz)+np.sum(k_p*(P1[0,:j]+P2[0,:j])*(deltaz))))

        if (N[i-1,j]/(N[i-1,j]+(mu1/alpha_n1)) > I/(I+(mu1/alpha_i1))):
            total_time_arr1.append(i)
            total_zeta_arr1.append(zetas[j])
            total_color_arr1.append('tab:red')
        else:
            total_time_arr1.append(i)
            total_zeta_arr1.append(zetas[j])
            total_color_arr1.append('tab:blue')
            
        if (N[i-1,j]/(N[i-1,j]+(mu2/alpha_n2)) > I/(I+(mu2/alpha_i2))):
            total_time_arr2.append(i)
            total_zeta_arr2.append(zetas[j])
            total_color_arr2.append('tab:red')
        else:
            total_time_arr2.append(i)
            total_zeta_arr2.append(zetas[j])
            total_color_arr2.append('tab:blue')

        mu1_growth = np.minimum((N[i-1,j]/(N[i-1,j]+(mu1/alpha_n1))), (I/(I+(mu1/alpha_i1))))
        mu2_growth = np.minimum((N[i-1,j]/(N[i-1,j]+(mu2/alpha_n2))), (I/(I+(mu2/alpha_i2))))
        dNdt = supply_rate[i] - (mu1*mu1_growth*P1[i-1,j]) - (mu2*mu2_growth*P2[i-1,j]) + (delta1*P1[i-1,j]) +(delta2*P2[i-1,j])
        dP1dt = (mu1*mu1_growth*P1[i-1,j]) - (phi1*P1[i-1,j]*C[i-1,j]) - (delta1*P1[i-1,j])
        dP2dt = (mu2*mu2_growth*P2[i-1,j]) - (phi2*P2[i-1,j]*C[i-1,j]) - (delta2*P2[i-1,j])
        dCdt = (eps1*phi1*P1[i-1,j]*C[i-1,j]) + (eps2*phi2*P2[i-1,j]*C[i-1,j]) - (delta_c*C[i-1,j])
        N[i,j] = N[i-1,j] + dNdt * dt
        P1[i,j] = P1[i-1,j] + dP1dt * dt
        P2[i,j] = P2[i-1,j] + dP2dt * dt
        C[i,j] = C[i-1,j] + dCdt * dt
        light[i,j] = I

# print shape of arrays for nutrients, producer, and consumer groups
print('N shape = ', N.shape)
print('P1 shape = ', P1.shape)
print('P2 shape = ', P2.shape)
print('C shape = ', C.shape)

N = N.T
P1 = P1.T
P2 = P2.T
C = C.T

fig, axs = plt.subplots(3,4,figsize=(30,12))

axs[0,0].plot(t,np.log(N[0, :]), color = 'red')
axs[0,0].set_title('Nutrients at surface')
axs[0,0].set_xlabel('time (days)')
axs[0,0].set_ylabel('Nutrient Concentration')

axs[1,0].plot(t,np.log(N[50, :]), color = 'red')
axs[1,0].set_title('Nutrients at half depth')
axs[1,0].set_xlabel('time (days)')
axs[1,0].set_ylabel('Nutrient Concentration')

axs[2,0].plot(t,np.log(N[100, :]), color = 'red')
axs[2,0].set_title('Nutrients at full depth')
axs[2,0].set_xlabel('time (days)')
axs[2,0].set_ylabel('Nutrient Concentration')

axs[0,1].plot(t,np.log(P1[0, :]), color = 'orange', linestyle = 'dashed')
axs[0,1].set_title('Producer 1 at surface')
axs[0,1].set_xlabel('time (days)')
axs[0,1].set_ylabel('Producer Biomass')

axs[1,1].plot(t,np.log(P1[50, :]), color = 'orange', linestyle = 'dashed')
axs[1,1].set_title('Producer 1 at half depth')
axs[1,1].set_xlabel('time (days)')
axs[1,1].set_ylabel('Producer Biomass')

axs[2,1].plot(t,np.log(P1[100, :]), color = 'orange', linestyle = 'dashed')
axs[2,1].set_title('Producer 1 at full depth')
axs[2,1].set_xlabel('time (days)')
axs[2,1].set_ylabel('Producer Biomass')

axs[0,2].plot(t,np.log(P2[0, :]), color = 'blue')
axs[0,2].set_title('Producer 2 at surface')
axs[0,2].set_xlabel('time (days)')
axs[0,2].set_ylabel('Producer Biomass')

axs[1,2].plot(t,np.log(P2[50, :]), color = 'blue')
axs[1,2].set_title('Producer 2 at half depth')
axs[1,2].set_xlabel('time (days)')
axs[1,2].set_ylabel('Producer Biomass')

axs[2,2].plot(t,np.log(P2[100, :]), color = 'blue')
axs[2,2].set_title('Producer 2 at full depth')
axs[2,2].set_xlabel('time (days)')
axs[2,2].set_ylabel('Producer Biomass')

axs[0,3].plot(t,np.log(C[0, :]), color = 'grey', linestyle = 'dotted')
axs[0,3].set_title('Consumer at surface')
axs[0,3].set_xlabel('time (days)')
axs[0,3].set_ylabel('Consumer Biomass')

axs[1,3].plot(t,np.log(C[50, :]), color = 'grey', linestyle = 'dotted')
axs[1,3].set_title('Consumer at half depth')
axs[1,3].set_xlabel('time (days)')
axs[1,3].set_ylabel('Consumer Biomass')

axs[2,3].plot(t,np.log(C[100, :]), color = 'grey', linestyle = 'dotted')
axs[2,3].set_title('Consumer at full depth')
axs[2,3].set_xlabel('time (days)')
axs[2,3].set_ylabel('Consumer Biomass')


plt.show()





