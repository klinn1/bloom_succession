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

#plotting supply rate over time
plt.plot(t,supply_rate, color = 'black', linestyle = 'dotted')
plt.title('Supply Rate')
plt.xlabel('time (days)')
plt.ylabel('$S_R$')
plt.show()

