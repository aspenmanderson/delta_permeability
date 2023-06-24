# Aspen Anderson
# Calculate permeability function - compatible with matlab
# Last Modified: 4/5/2021

'''
To call from matlab:
   1) change working directory to folder where function resides (addpath is not enough)
   2) set python environment (code below)
   2) load variables from Delft3D output (.dat file)
   3) import function from file and run (code below)
    
 % set python environment 
pe = pyenv;
if pe.Status == "NotLoaded"
    %[~,exepath] = system("where python");
    pypath = 'C:\Users\labuser10\Anaconda2\pythonw.exe'; %if where doesn't work, hard code path
        % found path using > import sys > sys.executable in Spyder
    pe = pyenv('Version',pypath);
end  

% import function from file and run
clear mod
mod = py.importlib.import_module('calc_hyd_params_for_matlab')
py.reload(mod) 
perm = py.calc_hyd_params_for_matlab.calc_perm(msed,9,-1,11,M,N,step,layer)
'''

###############################################################################
# test function is working in matlab
###############################################################################
 
def test(input1):
    print("Test successful!") 
    
    output1 = input1 + input1
    import numpy as np
    
    return output1

###############################################################################
# calculate permeability
###############################################################################

def calc_perm(perm):
        
    #import numpy as np
    #from scipy.interpolate import interp1d
    M = 302
    N = 227
    
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
    
    
    '''
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
    '''  
    # calculate premeability
        #c and j values from Shepard for beach material
            #c = 587.896 #12000 [gal_per_day/ft^2] /7.4805 * 0.36648 -> 587.896 [darcy]
    c = 5.8087E-6 #12000 [gal_per_day/ft^2] /7.4805 * 3.6210E-9 -> 5.8087E-6 [cm^2] or 40.746 for m
    j = 1.65 
        #x50_m = x50/100 #[mm] -> [m]
    perm = np.multiply(c,np.power(x50,j)) # [m^2] create conductivity matrix for a single timestep 
        #saywer got values -3, -5, and -7 m/s 
    ''' 
    
    # plot permeability
    fig, ax1 = plt.subplots(1,1, figsize=(16,10))
    im = ax1.imshow(perm[0:M,0:N].T, cmap = terrain_cmap, origin='lower')
    cbar = fig.colorbar(im, ax=ax1, orientation='vertical')
    cbar.set_label('permeability (darcy)',  size = 12)
    ax1.tick_params(axis = 'both', which = 'major', labelsize = 12)

    #return (perm)


###############################################################################
# 
###############################################################################
