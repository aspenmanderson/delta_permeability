# Aspen Anderson
# Import Delft3D data
# Last Modified: 11/27/2019

import numpy as np
import netCDF4
from netCDF4 import Dataset


###############################################################################
# import data from .nc file
###############################################################################
class nc:
    def __init__(self, file_in_path):
        self.data = Dataset(file_in_path,"r")
        
    def get_grids_params(self): 
        xz = self.data.variables['XZ'][:]
        #yz = self.data.variables['YZ'][:]
        M = len(xz[:,0])
        N = len(xz[0,:])
        cellsize = 25 #[m]
        return M, N, cellsize
   
    def get_time_params(self): 
        times = self.data.variables['time']
        dates = netCDF4.num2date(times[:], times.units)
        return times, dates  
    
    def get_sed_params(self):
        namcon = self.data.variables['NAMCON'][:] #name of constituent quantity (LSTSCI, strlen20)
        msed = self.data.variables['MSED'][:] #mass of sediment in layer [kg/m2] (time, LSEDTOT, nlyr, M, N)
        r1 = self.data.variables['R1'][:] #concentrations per layer in zeta point (time, LSTSCI, KMAXOUT_RESTR, M, N)
        rho = self.data.variables['RHO'][:] #density per layer (time, KMAXOUT_RESTR, M, N)
        rhocon = self.data.variables['RHOCONST'][:] #density per layer (time, KMAXOUT_RESTR, M, N)
        #rsedeq = self.data.variables['RSEDEQ'][:] #user specified concentration (int)
        lyrfrac = self.data.variables['LYRFRAC'][:] #volume fraction of sediment in layer (time, LSEDTOT, nlyr, M, N)
        sednum = len(lyrfrac[0,:,0,0,0].data.astype('int'))
        dp_bedlyr = self.data.variables['DP_BEDLYR'][:] #vertical position of sediment layer interface (time, nlyrp1, M, N)
        #epspor = self.data.variables['EPSPOR'][:] #porosity coefficent [ratio] (time, nlayer, M, N)
        return namcon, sednum, msed, r1, rho, rhocon, lyrfrac, dp_bedlyr
        
    def get_mor_params(self):
        mftavg = self.data.variables['MFTAVG'][:] #morphological time (days since simulation)
        moravg = self.data.variables['MORAVG'][:] #morphological time (days since simulation)
        return mftavg, moravg
        
    def get_flow_params(self):
        kcs = self.data.variables['KCS'][:] #non-active/active water level points (time, KMAXOUT_RESTR, M, N)
        s1 = self.data.variables['S1'][:] #water-level in zeta point (time, M, N)
        return kcs, s1
        
        
    def get_domain_params(self):
        dps = self.data.variables['DPS'][:] # initial bottom depth at zeta points (time, layer, M, N)
        gsqs = self.data.variables['GSQS'][:] #horizontal area of computational cell (M, N)
        return dps, gsqs




