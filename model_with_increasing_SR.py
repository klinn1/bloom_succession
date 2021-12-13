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
phi1 = 0.2 #interaction strength for p1
phi2 = 0.2 #interaction strength for p2
eps1 = 0.4 #transfer efficiency for p1
eps2 = 0.4 #transfer efficiency for p2
delta_c = 0.001 #consumer mortality
k_w = 0.02 #light attenuation by just water
k_p = 0.001 #attenuation due to producers
delta1 = 0.001 #death rate of producer 1
delta2 = 0.001 #death rate of producer 2
mu1 = 1.0 #resource affinity parameter for p1
mu2 = 1.0 #resource affinity parameter for p2
alpha_n1 = 0.17 #saturation of nutrients for producer 1
alpha_i1 = 0.017 #saturation of light for producer 1
alpha_n2 = 0.17 #saturation of nutrients for producer 2
alpha_i2 = 0.017 #saturation of light for producer 2
 
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
            
fig, axs = plt.subplots(2,4,figsize=(30,12))

#plotting supply rate over time
axs[0,0].plot(t,supply_rate, color = 'black', linestyle = 'dotted')
axs[0,0].set_title('Supply Rate')
axs[0,0].set_xlabel('time (days)')
axs[0,0].set_ylabel('$S_R$')

#contour plot of seasonal irradiance
cb = axs[0,1].pcolor(light.T,cmap = 'hot')
axs[0,1].invert_yaxis()
axs[0,1].set_xlabel('Time (Days)', color = 'k')
axs[0,1].set_ylabel('Depth (Meters)', color = 'k')
axs[0,1].set_title('Irradiance')
divider = make_axes_locatable(axs[0,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cb, cax=cax, orientation='vertical')

#showing min as either nutrients or light for producer 1
im = axs[0,2].scatter(total_time_arr1, total_zeta_arr1, c=total_color_arr1)
axs[0,2].set_ylim(axs[0,2].get_ylim()[::-1])
axs[0,2].set_ylabel('Depth (meters)')
axs[0,2].set_xlabel('Time (days)')
axs[0,2].set_title('Limiting Resource for Producer 1')

#showing min as either nutrients or light for producer 2
im = axs[0,3].scatter(total_time_arr2, total_zeta_arr2, c=total_color_arr2)
axs[0,3].set_ylim(axs[0,3].get_ylim()[::-1])
axs[0,3].set_ylabel('Depth (meters)')
axs[0,3].set_xlabel('Time (meters)')
axs[0,3].set_title('Limiting Resource for Producer 2')

#plotting nutrients
axs[1,0].pcolor(np.log(N.T),cmap = 'hot')
axs[1,0].set_xlabel('Time (days)', color = 'k')
axs[1,0].set_ylabel('Concentration (mmol C/ m$^{3}$)', color = 'k')
axs[1,0].set_title('Nutrient Concentration')
divider = make_axes_locatable(axs[1,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cb, cax=cax, orientation='vertical')

#plotting biomass
pro1 = axs[1,1].pcolor(np.log(P1.T),cmap = 'hot')
axs[1,1].invert_yaxis()

axs[1,1].set_xlabel('Time (days)', color = 'k')
axs[1,1].set_ylabel('Biomass (mmol C/ m$^{3}$ day)', color = 'k')
axs[1,1].set_title('Biomass of Producer 1')
divider = make_axes_locatable(axs[1,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cb, cax=cax, orientation='vertical')

pro2 = axs[1,2].pcolor(np.log(P2.T), cmap = 'hot')
axs[1,2].invert_yaxis()
axs[1,2].set_xlabel('Time (days)', color = 'k')
axs[1,2].set_ylabel('Biomass (mmol C/ m$^{3}$ day)', color = 'k')
axs[1,2].set_title('Biomass of Producer 2')
divider = make_axes_locatable(axs[1,2])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cb, cax=cax, orientation='vertical')

con = axs[1,3].pcolor(np.log(C.T),cmap = 'hot')
axs[1,3].invert_yaxis()
axs[1,3].set_xlabel('Time (days)', color = 'k')
axs[1,3].set_ylabel('Biomass (mmol C/ m$^{3}$ day)', color = 'k')
axs[1,3].set_title('Biomass of Consumer')
divider = make_axes_locatable(axs[1,3])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cb, cax=cax, orientation='vertical')

#supply rate versus producers (succesional pattern)
#axs[4].plot(supply_rate,P1.T[0], '--', color = 'orange')
# #axs[4].plot(supply_rate,P2.T[0], color = 'blue')
# #axs[4].set_title('Producer Biomass as $S_R$ Increases')
# #axs[4].set_ylabel('Biomass (mmol C/ m$^{3}$ day)', color = 'k')
# #axs[4].set_xlabel('$S_R$')


plt.tight_layout()
plt.show()




