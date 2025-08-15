### This file contains important functions for plotting large HDF5 files from COMPAS

#imports!

# let's import things
import h5py as h5 
import pandas as pd
from astropy.table import Table
from astropy import units as u
from astropy import constants as const
import numpy as np
import matplotlib.pyplot as plt
plt.rc('text.latex', preamble=r'\usepackage{textgreek}')
plt.rc('font', family='serif')
import sys
import os
from scipy import stats
import seaborn as sns
import matplotlib as mpl

# Add the subdir to sys.path
sys.path.append('/home/jovyan/home/research_work/useful_py_scripts/')

# Now you can import the module
import useful_fncs 

# import for axes labels 
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif"
})




def triangle_plot_fnc(pathToH5, title):

    """
    Plotting the triangle plot
    pathTOH5 = path to the HDF5 file
    """

# Set the appropriate path to the data file + read in the data

    pathToH5 = pathToH5
    Data  = h5.File(pathToH5, "r")

# To make the triangle plot we need the stellar types, masses, rates, and DCO mask

    DCOs = Data['BSE_Double_Compact_Objects'] # gathering the DCO group

    DCO_mask = Data['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['DCOmask'][()]

    stellar_types_1 = DCOs['Stellar_Type(1)'][()][DCO_mask]
    stellar_types_2 = DCOs['Stellar_Type(2)'][()][DCO_mask]

    mass1 = DCOs['Mass(1)'][()][DCO_mask]
    mass2 = DCOs['Mass(2)'][()][DCO_mask]

    rates_z0_DCO = Data['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['merger_rate_z0'][()]

# Let's add this data to a dataframe so we can mask and manipulate the data more efficently
    data = {
    "Stellar_Type(1)": stellar_types_1,
    "Stellar_Type(2)": stellar_types_2,
    "Mass(1)": mass1,
    "Mass(2)": mass2
    }

    DCOs_masked = pd.DataFrame(data)

# let's analyze how many systems of interest there are in our output

# NSNS
    NSNS_systems_bool = np.logical_and(stellar_types_1==14, stellar_types_2==14)
    print("There are {} NSNS systems from the DCO mask." .format(sum(NSNS_systems_bool)))

# WDWD
    WDWD_bool = np.logical_and(np.isin(stellar_types_1,[10,11,12]),np.isin(stellar_types_2,[10,11,12]))
    print("There are {} WDWD systems from the DCO masks." .format(sum(WDWD_bool)))

# COWD
    HeWD_bool,COWD_bool,ONeWD_bool,HeCOWD_bool,HeONeWD_bool,COHeWD_bool,COONeWD_bool,ONeHeWD_bool,ONeCOWD_bool = useful_fncs.WD_BINARY_BOOLS(DCOs_masked)
    carbon_oxygen_bool = np.logical_or(ONeCOWD_bool,np.logical_or(COONeWD_bool,np.logical_or(COHeWD_bool,np.logical_or(COWD_bool,HeCOWD_bool))))
    print("There are {} COWD systems from the DCO masks." .format(sum(carbon_oxygen_bool)))

# Making sure we are only masking the COWD 
    Mass1 = np.array(mass1[carbon_oxygen_bool]) 
    Mass2 = np.array(mass2[carbon_oxygen_bool])

# Now, let's make sure M1 represents the more massive object and M2 the less massive
    M1 = np.maximum(Mass1, Mass2)
    M2 = np.minimum(Mass1, Mass2)

# Let's plot!

# change figuresize
    fig, ax = plt.subplots(figsize = (6,6))

# let's define a few things first
    WDWD_merger_rate_Z0 = rates_z0_DCO[carbon_oxygen_bool]
    vmin = 10**-3
    vmax = 10**5

# for 2D histogram or 2D plots, we need to use something that will allow us to use a color bar so something like pcolormesh or plt.contour or hexbin
    hb = plt.hexbin(M1,M2,C=WDWD_merger_rate_Z0, gridsize=(15,15), reduce_C_function = np.sum, 
                    cmap=sns.color_palette("Spectral_r",as_cmap=True),norm='log',vmin=vmin,vmax=vmax) # C is how much each point is weighted
                # use symlog when you also want to cover negative values    
# right now we are not dividing by the bin size, so when we chage the bins - it changes the shape of our dist 
    zvalue_array = hb.get_array() # the merger rates of the histogram
    # print(min(zvalue_array),max(zvalue_array)) # helps us detemine what vim and vmax should be and what the bin size should be 

# colorbar
    cb = plt.colorbar()
    cb.set_label(label="Event Rate at $\mathrm{z =0}$, $\mathrm{dNdGpc^{-3}dyr^{-1}}$", fontsize = 15)

    max_mass_lim = 1.4
    plt.ylim(0.15,max_mass_lim)
    plt.xlim(0.15,max_mass_lim)

# let's add the mass restrictions for each case of binary WDs as prompted by Shen 2025

    xlim = max(M1)
    ylim = max(M1)

    linecolors = 'white'
    linewidths = 4

# Setting our line boundaries
    plt.axline((0,0), (max(mass1[carbon_oxygen_bool]),max(mass1[carbon_oxygen_bool])), color=linecolors, ls='--', lw=linewidths, transform=plt.gca().transAxes)

#Helium WD cutoff
    plt.vlines(x=[0.5], ymin=0, ymax=0.5, colors=linecolors, ls='--', lw=linewidths) # vertical line
    plt.plot([0.5,max_mass_lim],[0.5,0.5],color=linecolors,lw=linewidths, ls='--') # horizontal line

#Carbon oxygen WD cutoff
    plt.vlines(x=[1.1], ymin=0, ymax=1.1, colors=linecolors, ls='--', lw=linewidths) # vertical line
    plt.plot([1.1,max_mass_lim],[1.1,1.1],color=linecolors,lw=linewidths, ls='--') # horizontal line

# purple region - 2 star SN Ia
    plt.plot([0.8,1],[0.8,0.5],color='purple') # bottom boundary
    plt.plot([1.0,1.0],[1.0,0.5],color='purple') # side boundary
    plt.plot([0.8,1.0],[0.8,1.0],color='purple',ls='--') # top boundary 

# # chandrasekar mass line
    plt.plot((1.4,0.7),(0,0.7),color='black', lw=3, ls='--', alpha = 0.4)


# red region - hypervelocity WDs
    plt.plot([0.8,1],[0.8,0.5],color='red',ls='--') # overlapping boundary
    plt.plot([0.75,1.0],[0.75,0.375],color='red') # bottom boundary
    plt.plot([1.0,1.0],[0.5,0.375],color='red') # side boundary
    plt.plot([0.75,0.8],[0.75,0.8],color='red',ls='--') # top boundary
    plt.plot([1.1,1.1],[0.7,1.1],color='orange',ls='--') # left side boundary 

# orange region - 2003fg HVS
    plt.plot([1.0,1.1],[0.8,0.7],color='orange') # botton boundary
    plt.plot([1.0,1.0],[1.0,0.8],color='orange',ls='--') # left side overlapping boundary
    plt.plot([1.0,1.1],[1.0,1.1],color='orange',ls='--') # top boundary


    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel("$M_{1}$[$M_{\odot}$]",fontsize=20)
    plt.ylabel("$M_{2}$[$M_{\odot}$]",fontsize=20)
    plt.title(title)

## save figure:
    # plt.savefig("./figures/triangle_plots/triangle_CEalpha1.png",bbox_inches='tight',pad_inches=0.1)

# Finally, let's close the HDF5 file
    Data.close()

