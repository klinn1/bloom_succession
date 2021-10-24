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
year11.fill(0.11)
supply_rate = np.append(supply_rate,year11)

year12 = np.empty(365)
year12.fill(0.18)
supply_rate = np.append(supply_rate,year12)

year13 = np.empty(365)
year13.fill(0.26)
supply_rate = np.append(supply_rate,year13)

year14 = np.empty(365)
year14.fill(0.33)
supply_rate = np.append(supply_rate,year14)

year15 = np.empty(365)
year15.fill(0.40)
supply_rate = np.append(supply_rate,year15)

k_w = 0.02 #light attenuation by just water
k_p = 0.001 #attenuation due to producers

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

#initial light irradiance (constant)
I_naughts = np.empty(len(t))
I_naughts.fill(2)

#inital light irradiance (seasonal)
I_max = 1000
I_min = 500
# I_naughts2 = (np.sin(((t)/365*2*math.pi-math.pi/2))+1)*I_max/2

#combine I_naught arrays 
# I_naughts = np.append(I_naughts,I_naughts2)

i = 0
for a in t:
    if a<=1825.0:
        I_naughts[i] = I_min
    else:
        I_naughts[i] = I_min + (I_max-I_min)/2.*(np.sin(((a)/365*2*np.pi-np.pi/2))+1)
    i = i + 1
        

#building 2 dimensional arrays for nutrients, producer, and consumer groups
for i in  range(1, len(t)):
    for j in range(0,len(zetas)):
        I = I_naughts[i] * np.exp(-(k_w*zetas[j]*(deltaz)+np.sum(k_p*(P1[0,:j]+P2[0,:j])*(deltaz))))
        light[i,j] = I

print('light shape = ',light.shape)


fig, axs = plt.subplots(1,2,figsize=(30,6))

#plotting supply rate over time
axs[0].plot(t,supply_rate, color = 'black', linestyle = 'dotted')
axs[0].set_title('Supply Rate')
axs[0].set_xlabel('time (days)')
axs[0].set_ylabel('$S_R$')

#contour plot of seasonal irradiance
cb = axs[1].pcolor(light.T,cmap = 'hot')
axs[1].invert_yaxis()
axs[1].set_xlabel('Time (Days)', color = 'k')
axs[1].set_ylabel('Depth (Meters)', color = 'k')
axs[1].set_title('Irradiance')
divider = make_axes_locatable(axs[1])
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(cb, cax=cax, orientation='vertical')

plt.show()
