# Aspen Anderson
# Calculate permeability function
# Last Modified: 11/27/2019

import numpy as np

def calc_perm(msed,low_phi,high_phi,sed_num,M,N,step,layer):
        
        
    # calculate total mass in each cell
    for i in np.arange(0,M):
        for j in np.arange(0,N):
            msed_rev = msed[step-1,:,layer-1,i,j].data.astype('Float64')
    
    msed_rev = msed_rev[::-1]
    
    totmass = np.zeros((M,N))
    for i in np.arange(0,M):
        for j in np.arange(0,N):
            totmass[i,j] = np.sum(msed[step-1,:-1,layer-1,i,j].data.astype('Float64')) #[kg]
            
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
    num_sed_type = sed_num
    
    sed_size = [2**(-1*9), 2**(-1*8), 2**(-1*7), 2**(-1*6), 2**(-1*5), 2**(-1*4), 2**(-1*3), 2**(-1*2), 2**(-1*1), 2**(-1*0), 2**(-1*-1)] # size of sediment input into delft
     #sed_size_rev = sed_size[::-1]
     #x = np.linspace(low_dmm, high_dmm, num = num_sed_type, endpoint=True) # sediment size d [mm] - same size as cumfrac
    xnew = np.linspace(low_dmm, high_dmm, num = 100, endpoint=True) # sediment size d [mm] - more refined for interpolation
    
    x50 = np.zeros((M,N))
    for i in np.arange(0,M):
        for j in np.arange(0,N):
            f2 = interp1d(sed_size, percent_finer[i,j,:], kind='linear')  #'linear' 'cubic'
            ind = min(range(len(f2(xnew))), key=lambda l:abs(f2(xnew)[l]-0.50)) #get index of value closest to 0.50
            x50[i,j] = xnew[ind] # find x [d50 in mm] at y ~ 0.50
    
    
    # calculate premeability
        # c and j values from Shepard for beach material
        #c = 587.896 #12000 [gal_per_day/ft^2] /7.4805 * 0.36648 -> 587.896 [darcy]
    c = 5.8087E-6 #12000 [gal_per_day/ft^2] /7.4805 * 3.6210E-9 -> 5.8087E-6 [cm^2] or 40.746 for m
    j = 1.65 
        #x50_m = x50/100 #[mm] -> [m]
    perm = np.multiply(c,np.power(x50,j)) # [m^2] create conductivity matrix for a single timestep 
        #saywer got values -3, -5, and -7 m/s 
    

    return (perm)
