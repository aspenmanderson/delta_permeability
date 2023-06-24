# Aspen Anderson
# Delft3D - Create Sediment Distribution  
# Last Modified: 11/9/2018


from numpy import *
import matplotlib.pyplot as plt
import math


##############################################################################
# User Input: Set Parameters
##############################################################################
dmm50 = 0.5 #mm
sigma_mm = 0.10 #mm


#phi50 = 1 #mediam grain size
#sigma_phi = 2 #standard deviation
phi50 = -1*math.log(dmm50,2.0)
sigma_phi = -1*math.log(sigma_mm,2.0)

low_phi = -1
high_phi = 9

low = low_phi-1
high = high_phi+1

# Calculate discharge of sediment
qw = 1890 #total river discharge [m^3/s]
ps = 2650 #sediment density [kg/m^3]
sed_con = 0.5 #concentration [kg/m^3]  
    #this value is the same as 0.0377 m^3/s sediment for 1000 m^3/s water

qs = sed_con*qw/ps
print("sed discharge = ",qs) #total sediment discharge [m^3/s]
    #qs = 0.0377 is what Caldwell and Edmounds used (total SUSPENDED sediment discharge [m^3/s])
    # mine at 0.34 m^3/s is about an order or magnitude higher


#############################################################################
# Calculate
##############################################################################

# Create normal distribution 
#particles = random.normal(dmm50, sigma_mm, 1000) 
particles = random.normal(phi50, sigma_phi, 1000) 

# Calculate in mm
#D = particles
D = 2**(-1*particles)
D50 = 2**(-1*phi50)
print("D50 = ", D50)

k = 0
for i in arange(0,len(D)):
    if D[i] < 0.064:
        k += 1
    #print(k)
        
if k > 0:
    percent_cohesion = float(k)/float(len(D))*100
else: 
    percent_cohesion = 0

# Create histogram 
    # bins from -1phi 10 10pih, where -1phi is bin 0
edges = arange(low,high+2,2)
#centers =
#frequency = digitize(particles, edges , right=False)
hist, bin_edges = histogram(particles, bins=edges, density=1)
    #phi values are the center of each bin

# Calculate concentration of each sediment
concentration = multiply(sed_con,hist*diff(bin_edges)) #[kg/m^3]


print(concentration)
print("sum concentration = ", sum(concentration)) #[kg/m^3]
print("% cohesion = ", percent_cohesion)

# Phi 5 std 2 concentrations: [0.00050302 0.00513078 0.02464789 0.03822938 0.02545272 0.00603622]. Percent cohesion: 70.3

# Phi 5 std 1 concentrations: [0     0.0002 0.0149 0.0692 0.0156 0.0001].Percent cohesion: 85.3

# Phi 1 std 2 concentrations: [0.02448759 0.04066882 0.02664509 0.00733549 0.000863   0     ].Percent cohesion: 7.7


##############################################################################
# Plot
##############################################################################
f, (ax1) = plt.subplots(1, 1, sharex='col', figsize=(8,6))
#first plot
#ax1.hist(frequency,bin_edges,color='darkgrey',edgecolor='black',linewidth=1.2)
#ax1.hist(particles,bin_edges,density=True,color='darkgrey',edgecolor='black',linewidth=1.2)
#ax1.plot(arange(-1,11,0.01),1/(sigma_phi*sqrt(2*pi))*exp(-(arange(-1,11,0.01)-phi50)**2/(2*sigma_phi**2)),linewidth=2, color='r')

#second plot
#ax1.plot(arange(low,high,0.01),1/(sigma_phi*sqrt(2*pi))*exp(-(arange(low,high,0.01)-1)**2/(2*sigma_phi**2)),linewidth=2,color='black',linestyle='--',label='$\phi_{50}=1$')
ax1.plot(arange(low,high,0.01),1/(sigma_phi*sqrt(2*pi))*exp(-(arange(low,high,0.01)-3)**2/(2*sigma_phi**2)),linewidth=2,color='black',label='$\phi_{50}=3$')
#ax1.plot(arange(low,high,0.01),1/(sigma_phi*sqrt(2*pi))*exp(-(arange(low,high,0.01)-5)**2/(2*sigma_phi**2)),linewidth=2,color='dimgrey',linestyle='--',label='$\phi_{50}=5$')
ax1.hist(particles,bin_edges,density=True,color='darkgrey',edgecolor='black',linewidth=1.2)

ax1.set_ylabel('percent',size=16) #Mean Flow Velocity 
ax1.tick_params(axis='y', labelsize=14)
ax1.set_xlabel(r'$\phi_{50} $ [mm]',size=16) #Natural Log of Depth
ax1.tick_params(axis='x', labelsize=14)
ax1.set_xlim([-2,10])
#ax1.legend(bbox_to_anchor=(1, -0.15),prop={'size': 14},ncol=4)
#ax1.set_ylim([0,34.9])


"""