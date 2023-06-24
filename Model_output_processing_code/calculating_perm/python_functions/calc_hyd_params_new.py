# Aspen Anderson
# Calculate permeability function
# Last Modified: 11/27/2019

import numpy as np
from scipy.interpolate import interp1d

import matplotlib.pylab as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

# colormaps
virids_cmap = cm.get_cmap('viridis', 256)
blue_cmap = cm.get_cmap('Blues', 256)
gnbu_cmap = cm.get_cmap('GnBu', 256)
ylgnbu_cmap = cm.get_cmap('YlGnBu', 256)
terrain_cmap = cm.get_cmap('terrain_r', 256) #_r reversed
ocean_cmap = cm.get_cmap('ocean', 256)
gist_cmap = cm.get_cmap('gist_earth', 256)
#my_cmap = cmap(np.arange(cmap.N)) # set alpha
#my_cmap[:,-1] = np.concatenate([np.linspace(0.75, 0.4, cmap.N/2), np.linspace(0.4, 0.75, cmap.N/2)]) # create new colormap
#my_cmap = ListedColormap(my_cmap) 

###############################################################################
# get data
###############################################################################

import sys
sys.path.append('C:\Users\labuser10\Desktop') 
from functions import import_delft3d as impor

class data:
    def __init__(self, model):
        self.model = impor.model.get_sed_params()

namcon, sednum, msed_cup2, _, _, _, lyrfrac, _ = model.get_sed_params() 

###############################################################################
# calculate permeability
###############################################################################











def perm(msed,low_phi,high_phi,sed_num,M,N,step,layer):
        
    # calculate total mass in each cell
    for i in np.arange(0,M):
        for j in np.arange(0,N):
            msed_rev = msed[step-1,:,layer-1,i,j].data.astype('Float64')
    
    msed_rev = msed_rev[::-1]
    
    totmass = np.zeros((M,N))
    for i in np.arange(0,M):
        for j in np.arange(0,N):
            totmass[i,j] = np.sum(msed[step-1,:-1,layer-1,i,j].data.astype('Float64')) #[kg]
            
            
    # plot mass percent fractions
    '''
    fig, axes = plt.subplots(2,3, figsize=(18,10))
    im = axes[0,0].imshow(percent_mass[:,:,0].T, cmap = my_cmap, vmin=0, vmax=1, origin='lower')
    axes[0,0].set_title("phi50 = 7")
    fig.colorbar(im, ax=axes[0, 0], orientation='vertical',label='%')
    
    im =axes[0,1].imshow(percent_mass[:,:,1].T, cmap = my_cmap, vmin=0, vmax=1, origin='lower')
    axes[0,1].set_title("phi50 = 5")
    fig.colorbar(im, ax=axes[0, 1], orientation='vertical',label='%')
    
    im = axes[0,2].imshow(percent_mass[:,:,2].T, cmap = my_cmap, vmin=0, vmax=1, origin='lower')
    axes[0,2].set_title("phi50 = 3")
    fig.colorbar(im, ax=axes[0, 2], orientation='vertical',label='%')
    
    im = axes[1,0].imshow(percent_mass[:,:,3].T, cmap = my_cmap, vmin=0, vmax=1, origin='lower')
    axes[1,0].set_title("phi50 = 1")
    fig.colorbar(im, ax=axes[1, 0], orientation='vertical',label='%')
    
    im = axes[1,1].imshow(percent_mass[:,:,4].T, cmap = my_cmap, vmin=0, vmax=1, origin='lower')
    axes[1,1].set_title("phi50 = -1")
    fig.colorbar(im, ax=axes[1, 1], orientation='vertical',label='%')
    
    im = axes[1,2].imshow(percent_mass[:,:,5].T, cmap = my_cmap, origin='lower')
    axes[1,2].set_title("total")
    fig.colorbar(im, ax=axes[1, 2], orientation='vertical',label='mass (kg/m^2)')
    '''
            
    # calculate percent retained (also just the percent of mass for each sediment type)
    percent_mass = np.zeros((M,N,sed_num))
    for k in np.arange(0,sed_num):
        for i in np.arange(0,M):
            for j in np.arange(0,N):
                percent_mass[i,j,k] = np.divide(msed[step-1,sed_num-k-1,layer-1,i,j],totmass[i,j]) #[%] sed_num-k is to reverse (so it is small to large sediment)
    
    # calculate cumulative percent retained
    cumpercent_retained = np.zeros((M,N,sed_num))
    for i in np.arange(0,M):
        for j in np.arange(0,N):
            cumpercent_retained[i,j,:] = np.cumsum(percent_mass[i,j,:]) #[%]
                
    # calculate percent finer than 
    percent_finer = np.zeros((M,N,sed_num))
    for k in np.arange(0,sed_num):
        for i in np.arange(0,M):
            for j in np.arange(0,N):
                percent_finer[i,j,k] = np.subtract(1,cumpercent_retained[i,j,k]) #[%]
    
    
    # interpolate between data points
    low_dmm = 2**(-1*low_phi) # increasing phi scale is decreasing dmm
    high_dmm = 2**(-1*high_phi)
    
    sed_size = [2**(-1*9), 2**(-1*8), 2**(-1*7), 2**(-1*6), 2**(-1*5), 2**(-1*4), 2**(-1*3), 2**(-1*2), 2**(-1*1), 2**(-1*0), 2**(-1*-1)] # size of sediment input into delft
    xnew = np.linspace(low_dmm, high_dmm, num = 100, endpoint=True) # sediment size d [mm] - more refined for interpolation
    
    x50 = np.zeros((M,N))
    for i in np.arange(0,M):
        for j in np.arange(0,N):
            f2 = interp1d(sed_size, percent_finer[i,j,:], kind='linear')  #'linear' 'cubic'
            ind = min(range(len(f2(xnew))), key=lambda l:abs(f2(xnew)[l]-0.50)) #get index of value closest to 0.50
            x50[i,j] = xnew[ind] # find x [d50 in mm] at y ~ 0.50
    
    
    # plot interpolation
    '''
    xind = 150
    yind = 91
    f2 = interp1d(sed_size, percent_finer[xind,yind,:], kind='linear') # flip of plotting percent finer
    ind =  min(range(len(f2(xnew))), key=lambda l:abs(f2(xnew)[l]-0.50)) 
    x50_plt = xnew[ind] 
    
    fig, ax = plt.subplots(1,1, figsize=(5,5))
    ax.plot(sed_size, percent_finer[xind,yind,:], 'o')
    ax.plot(xnew, f2(xnew), '--',color = 'forestgreen')
    ax.plot(x50_plt, 0.50, 'x', markersize=12, color='red')
    ax.invert_xaxis()
    ax.legend(['data', 'interpolated', 'd50'], loc='best', fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_ylabel('percent finer than', fontsize=14)
    ax.set_xlabel('grain size (mm)', fontsize=14)
    
    fig, ax = plt.subplots(1,1, figsize=(5,5))
    ax.semilogx(sed_size, percent_finer[xind,yind,:], 'o')
    ax.semilogx(xnew, f2(xnew), '--',color = 'forestgreen')
    ax.semilogx(x50_plt, 0.50, 'x', markersize=12, color='red')
    ax.invert_xaxis()
    ax.legend(['data', 'interpolated', 'd50'], loc='best', fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_ylabel('percent finer than', fontsize=14)
    ax.set_xlabel('grain size (mm)', fontsize=14)
    '''
    
    # calculate premeability
        #c and j values from Shepard for beach material
            #c = 587.896 #12000 [gal_per_day/ft^2] /7.4805 * 0.36648 -> 587.896 [darcy]
    c = 5.8087E-6 #12000 [gal_per_day/ft^2] /7.4805 * 3.6210E-9 -> 5.8087E-6 [cm^2] or 40.746 for m
    j = 1.65 
        #x50_m = x50/100 #[mm] -> [m]
    perm = np.multiply(c,np.power(x50,j)) # [m^2] create conductivity matrix for a single timestep 
        #saywer got values -3, -5, and -7 m/s 
        
    # plot permeability
    fig, ax1 = plt.subplots(1,1, figsize=(16,10))
    im = ax1.imshow(perm[0:M,0:N].T, cmap = terrain_cmap, origin='lower')
    cbar = fig.colorbar(im, ax=ax1, orientation='vertical')
    cbar.set_label('permeability (darcy)',  size = 12)
    ax1.tick_params(axis = 'both', which = 'major', labelsize = 12)

    return (perm)


###############################################################################
# 
###############################################################################

