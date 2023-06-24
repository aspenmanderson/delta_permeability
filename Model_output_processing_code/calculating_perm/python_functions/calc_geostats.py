# Aspen Anderson
# Calculate geostats for 2D variables
# Last Modified: 11/27/2019

import numpy as np
import matplotlib.pylab as plt

import gstools as gs
from gstools import Exponential, Spherical, Circular
from scipy.stats import skew as skew

###############################################################################
# functions
###############################################################################
class field:
    def __init__(self, data):
        M = len(data[:,0])
        N = len(data[0,:])
        self.field = data[5:M-5,25:N-5].T #clips field to get rid of shoreline and boundaries
        self.field_logT = np.log10(data[5:M-5,25:N-5].T)

def stat(self):
    # plot distribution of data
    a = np.reshape(self.field,self.field.size)
    fig, ax = plt.subplots(1,1, figsize=(10,7))
    ax.hist(a,bins=100)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_ylabel('frequency', fontsize=14)
    ax.set_xlabel('permeability (Darcy)', fontsize=14)
    
    # find range
    #upper = np.mean(field)+3*np.std(field)
    #lower = np.mean(field)-3*np.std(field)
    
    # calculate general stats
    return np.mean(a), np.std(a), np.min(a), np.max(a), 
    np.median(a), np.var(a), skew(a)
            
def create_semivario(self, cellsize):
    
    # caluclate gamma
    vario_x = gs.variogram.vario_estimate_structured(field, direction='x') 
    vario_y = gs.variogram.vario_estimate_structured(field, direction='y')
    
    # plot semivarience
    fig, ax1 = plt.subplots(1,1, figsize=(11,7))
    ax1.semilogy(np.linspace(0,len(vario_x)/2,num = len(vario_x)/2)*cellsize,vario_x[0:len(vario_x)/2],'.-',label='x direction', color='steelblue', markersize=14)
    ax1.semilogy(np.linspace(0,len(vario_y)/2,num = len(vario_y)/2)*cellsize,vario_y[0:len(vario_y)/2],'.-',label='y direction', color='seagreen', markersize=14)
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.tick_params(axis = 'both', which = 'major', labelsize = 30)
    ax1.legend(loc=4, prop={'size': 24})
    
    #change border color
    for ax, color in zip([ax1], ['dimgrey']):
        for ticks in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
            ticks.set_color(color)
        for pos in ['top', 'bottom', 'right', 'left']:
            ax.spines[pos].set_edgecolor(color)

    
    