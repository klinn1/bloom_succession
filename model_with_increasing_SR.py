import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import math
import sys

# time change (discrete time)
T = 3650.0 #days
dt = 1
#time array
t = np.linspace(0, T, int(T/dt))

#supply rate array for 10 years, increasing stepwise
supply_rate = np.empty(365)
supply_rate.fill(0.04)

year2 = np.empty(365)
year2.fill(0.11)
supply_rate = np.append(supply_rate,year2)

year3 = np.empty(365)
year3.fill(0.18)
supply_rate = np.append(supply_rate,year3)

year4 = np.empty(365)
year4.fill(0.26)
supply_rate = np.append(supply_rate,year4)

year5 = np.empty(365)
year5.fill(0.33)
supply_rate = np.append(supply_rate,year5)

year6 = np.empty(365)
year6.fill(0.40)
supply_rate = np.append(supply_rate,year6)

year7 = np.empty(365)
year7.fill(0.48)
supply_rate = np.append(supply_rate,year7)

year8 = np.empty(365)
year8.fill(0.55)
supply_rate = np.append(supply_rate,year8)

year9 = np.empty(365)
year9.fill(0.62)
supply_rate = np.append(supply_rate,year9)

year10 = np.empty(365)
year10.fill(0.70)
supply_rate = np.append(supply_rate,year10)

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

#depth array
deltaz = 1
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
#print(N.shape)

#inital light irradiance (seasonal)
#I_max = 1000
#I_naughts = (np.sin(((t)/365*2*math.pi-math.pi/2))+1)*I_max/2

#initial light irradiance (constant)
I_naughts = np.empty(int(T/dt))
I_naughts.fill(2)

 
#building 2 dimensional arrays for nutrients, producer, and consumer groups
for i in  range(1, len(t)):
    #print(t[i])
    for j in range(0,len(zetas)):
        I = I_naughts[i] * np.exp(-(k_w*zetas[j]*(deltaz)+np.sum(k_p*(P1[0,:j]+P2[0,:j])*(deltaz))))
        mu1_growth = np.minimum((N[i-1,j]/(N[i-1,j]+(mu1/alpha_n))), (I/(I+(mu1/alpha_i))))
        mu2_growth = np.minimum((N[i-1,j]/(N[i-1,j]+(mu2/alpha_n))), (I/(I+(mu2/alpha_i))))
        dNdt = supply_rate[i] - (mu1*mu1_growth*P1[i-1,j]) - (mu2*mu2_growth*P2[i-1,j])
        dP1dt = (mu1*mu1_growth*P1[i-1,j]) - (phi1*P1[i-1,j]*C[i-1,j]) - (delta1*P1[i-1,j])
        dP2dt = (mu2*mu2_growth*P2[i-1,j]) - (phi2*P2[i-1,j]*C[i-1,j]) - (delta2*P2[i-1,j])
        dCdt = (eps1*phi1*P1[i-1,j]*C[i-1,j]) + (eps2*phi2*P2[i-1,j]*C[i-1,j]) - (delta_c*C[i-1,j])
        N[i,j] = N[i-1,j] + dNdt * dt
        P1[i,j] = P1[i-1,j] + dP1dt * dt
        P2[i,j] = P2[i-1,j] + dP2dt * dt
        C[i,j] = C[i-1,j] + dCdt * dt
        light[i,j] = I

print('light shape = ',light.shape)
# print shape of arrays for nutrients, producer, and consumer groups
print('N shape = ', N.shape)
print('P1 shape = ', P1.shape)
print('P2 shape = ', P2.shape)
print('C shape = ', C.shape)
            
fig, axs = plt.subplots(2,3,figsize=(30,12))

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

#plotting nutrients
#axs[2].plot(t,N.T[0], color = 'red')
axs[0,2].pcolor(N,cmap = 'hot')
axs[0,2].set_xlabel('Time (days)', color = 'k')
axs[0,2].set_ylabel('Concentration (mmol C/ m$^{3}$)', color = 'k')
axs[0,2].set_title('Nutrient Concentration')

#plotting biomass
pro1 = axs[1,0].pcolor(P1,cmap = 'hot')
axs[1,0].invert_yaxis()
axs[1,0].set_xlabel('Time (Days)', color = 'k')
axs[1,0].set_ylabel('Depth (Meters)', color = 'k')
axs[1,0].set_title('Irradiance')
axs[1,0].set_xlabel('Time (days)', color = 'k')
axs[1,0].set_ylabel('Biomass (mmol C/ m$^{3}$ day)', color = 'k')
axs[1,0].set_title('Biomass of Producer 1')
divider = make_axes_locatable(axs[1,0])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cb, cax=cax, orientation='vertical')

pro2 = axs[1,1].pcolor(P2, cmap = 'hot')
axs[1,1].invert_yaxis()
axs[1,1].set_xlabel('Time (Days)', color = 'k')
axs[1,1].set_ylabel('Depth (Meters)', color = 'k')
axs[1,1].set_title('Irradiance')
axs[1,1].set_xlabel('Time (days)', color = 'k')
axs[1,1].set_ylabel('Biomass (mmol C/ m$^{3}$ day)', color = 'k')
axs[1,1].set_title('Biomass of Producer 2')
divider = make_axes_locatable(axs[1,1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cb, cax=cax, orientation='vertical')

con = axs[1,2].pcolor(C,cmap = 'hot')
axs[1,2].invert_yaxis()
axs[1,2].set_xlabel('Time (Days)', color = 'k')
axs[1,2].set_ylabel('Depth (Meters)', color = 'k')
axs[1,2].set_title('Irradiance')
axs[1,2].set_xlabel('Time (days)', color = 'k')
axs[1,2].set_ylabel('Biomass (mmol C/ m$^{3}$ day)', color = 'k')
axs[1,2].set_title('Biomass of Consumer')
divider = make_axes_locatable(axs[1,2])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cb, cax=cax, orientation='vertical')

#supply rate versus producers (succesional pattern)
#axs[4].plot(supply_rate,P1.T[0], '--', color = 'orange')
#axs[4].plot(supply_rate,P2.T[0], color = 'blue')
#axs[4].set_title('Producer Biomass as $S_R$ Increases')
#axs[4].set_ylabel('Biomass (mmol C/ m$^{3}$ day)', color = 'k')
#axs[4].set_xlabel('$S_R$')

#add_axis ... supply bounding box for each one

plt.tight_layout()
plt.show()




