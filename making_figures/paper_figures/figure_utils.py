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
import utils_from_others

# import for axes labels 
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif"
})



def systems_of_interest_counter(pathToH5):
    """
    Printing the number of systems that we choose to mask
    pathToH5 = path to the HDF5 file
    """
# Set the appropriate path to the data file + read in the data

    Data  = h5.File(pathToH5, "r")

# To count our systems of interest we need the stellar types, masses, rates, and DCO mask

    DCOs = Data['BSE_Double_Compact_Objects'] # gathering the DCO group

    DCO_mask = Data['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['DCOmask'][()]

    stellar_types_1_all = DCOs['Stellar_Type(1)'][()]
    stellar_types_1 = stellar_types_1_all[DCO_mask]

    stellar_types_2_all = DCOs['Stellar_Type(2)'][()]
    stellar_types_2 = stellar_types_2_all[DCO_mask]

    mass1_all = DCOs['Mass(1)'][()]
    mass1 = mass1_all[DCO_mask]

    mass2_all = DCOs['Mass(2)'][()]
    mass2 = mass2_all[DCO_mask]

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

# Let's close the HDF5 File
    Data.close()




def triangle_plot_fnc(pathToH5, title, plot_output, filename):

    """
    Plotting the triangle plot
    pathTOH5 = path to the HDF5 file
    title = the name of the plot that corresponds to the file name
    """

# Set the appropriate path to the data file + read in the data

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

# let's make masks for each type of system we care to analyze 

# NSNS
    NSNS_systems_bool = np.logical_and(stellar_types_1==14, stellar_types_2==14)
# WDWD
    WDWD_bool = np.logical_and(np.isin(stellar_types_1,[10,11,12]),np.isin(stellar_types_2,[10,11,12]))
# COWD
    HeWD_bool,COWD_bool,ONeWD_bool,HeCOWD_bool,HeONeWD_bool,COHeWD_bool,COONeWD_bool,ONeHeWD_bool,ONeCOWD_bool = useful_fncs.WD_BINARY_BOOLS(DCOs_masked)
    carbon_oxygen_bool = np.logical_or(ONeCOWD_bool,np.logical_or(COONeWD_bool,np.logical_or(COHeWD_bool,np.logical_or(COWD_bool,HeCOWD_bool))))

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
    plt.title(title, pad=20)

# save figure:
    plt.savefig(plot_output + filename +".png",bbox_inches='tight',pad_inches=0.1)

# Finally, let's close the HDF5 file
    Data.close()





def redshift_rates_plotter(pathToH5, title, plot_output, filename):

    """
    Plotting the redhshift vs. rates plot for NSNS systems and COWD + WD systems
    pathTOH5 = path to the HDF5 file
    title = the name of the plot that corresponds to the file name
    """

# Set the appropriate path to the data file + read in the data

    Data  = h5.File(pathToH5, "r")


# let's gather the data we need for the redshift rates plots

# we want to use information in the double compact object group
    DCOs = Data['BSE_Double_Compact_Objects']
# gathering the double compact objects that we have computed rates for
    DCO_mask = Data['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['DCOmask'][()]

# first we want to investigate how many of each type of DCO there is
    stellar_types_1 = DCOs['Stellar_Type(1)'][()][DCO_mask]
    stellar_types_2 = DCOs['Stellar_Type(2)'][()][DCO_mask]

# we need the masses
    mass1 = DCOs['Mass(1)'][()][DCO_mask]
    mass2 = DCOs['Mass(2)'][()][DCO_mask]

# we also need the rates and redshift data
    rates_DCO = Data['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['merger_rate'][()]
    redshifts = Data['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['redshifts'][()]


# let's add everything to a dataframe so we can analyze things easier

    data = {
    "Stellar_Type(1)": stellar_types_1,
    "Stellar_Type(2)": stellar_types_2,
    "Mass(1)": mass1,
    "Mass(2)": mass2
    }

    DCOs_masked = pd.DataFrame(data)


# let's make bools for each type of system we care about 

# NSNS
    NSNS_systems_bool = np.logical_and(stellar_types_1==14, stellar_types_2==14)

# COWD + WD
    HeWD_bool,COWD_bool,ONeWD_bool,HeCOWD_bool,HeONeWD_bool,COHeWD_bool,COONeWD_bool,ONeHeWD_bool,ONeCOWD_bool = useful_fncs.WD_BINARY_BOOLS(DCOs_masked)
    carbon_oxygen_bool = np.logical_or(ONeCOWD_bool,np.logical_or(COONeWD_bool,np.logical_or(COHeWD_bool,np.logical_or(COWD_bool,HeCOWD_bool))))

# WDWD
    WDWD_bool = np.logical_and(np.isin(stellar_types_1,[10,11,12]),np.isin(stellar_types_2,[10,11,12]))


# let's now calculate the rates of these systems of interest
    cowd_rate = np.sum(rates_DCO[carbon_oxygen_bool], axis=0)
    NSNS_rate = np.sum(rates_DCO[NSNS_systems_bool], axis=0)


# extracting the redshifts and rates from Briel et al
# units in the appendix should be in h^-3 y^-1 Gpc^-3 so we must convert below to get yr^-1 Gpc^-3
    h_little = 0.6766

    redshifts_briel = [
        0, 0.01, 0.03, (0.025+0.050)/2, 0.073, (0.05+0.15)/2, (0.075+0.125)/2, 0.11, 0.11, 0.13, 
        0.15, (0.125+0.175)/2, 0.16, (0.175+0.225)/2, 0.2, 0.25, (0.15+0.35)/2, (0.225+0.275)/2, 
        0.26, 0.3, (0.275+0.325)/2, 0.35, 0.35, 0.42, 0.44, 0.45, 0.45, (0.35+0.55)/2, 0.46, 0.47, 
        0.47, 0.55, 0.55, 0.55, 0.62, 0.65, (0.55+0.75)/2, 0.65, 0.74, 0.75, 0.75, 0.75, 0.8, 0.83, 0.85, 
        0.85, 0.94, 0.95, 0.95, 1.05, 1.1, 1.14, 1.21, 1.23, 1.25, 1.59, 1.61, 1.69, 1.75, 2.25
    ]

    rates_briel = [
        0.77, 0.82, 0.82, 0.81, 0.71, 1.60, 0.76, 1.08, 0.72, 0.58, 0.93, 0.90, 0.41, 1.01, 0.58,
        1.05, 1.14, 1.06, 0.82, 0.99, 1.27, 0.99, 1.05, 1.34, 0.76, 0.90, 1.05, 1.52, 1.40, 1.22, 
        2.33, 0.93, 1.40, 1.52, 3.76, 1.40, 2.01, 1.43, 2.30, 1.49, 1.98, 1.69, 2.45, 3.79, 2.27, 
        1.66, 1.31, 2.22, 2.24, 2.30, 2.16, 2.06, 3.85, 2.45, 1.87, 1.31, 1.22, 2.97, 2.10, 1.43
    ]

# converting the rates to the correct units
    rates_briel = np.array(rates_briel)
    converted_rates_briel = (rates_briel*(10**5))*(h_little**3)

## uncertainties
    lower_limits = [
        -0.10, -0.26, -0.32, -0.24, -0.08, -0.85, -0.13, -0.29, -0.20, -0.18, -0.67, -0.10, -0.26, -0.09, 
        -0.23, -0.76, -0.35, -0.08, -0.20, -0.44, -0.10, -0.55, -0.17, -0.93, -0.39, -0.44, -0.17, -0.38, 
        -0.50, -0.17, -0.79, -0.41, -0.17, -0.26, -1.66, -0.15, -0.52, -0.50, -1.20, -0.55, -0.61, -0.17, 
        -0.54, -0.79, -0.64, -0.15, -0.55, -0.73, -0.23, -0.82, -0.35, -0.53, -0.85, -0.82, -0.64, -0.64, 
        -0.67, -1.08, -0.87, -1.11
    ]

    lower_limits = np.array(lower_limits)
    converted_lower_limits = (lower_limits*(10**5)*(h_little**3))

    upper_limits = [
        0.10, 0.26, 0.32, 0.33, 0.08, 1.46, 0.15, 0.29, 0.08, 0.20, 0.67, 0.11, 0.26, 0.09, 0.23,
        1.75, 0.38, 0.09, 0.20, 0.47, 0.11, 0.55, 0.17, 1.22, 0.67, 0.44, 0.17, 0.32, 0.50, 0.17, 
        1.08, 0.41, 0.17, 0.29, 2.57, 0.15, 0.55, 0.50, 0.96, 0.79, 0.61, 0.17, 0.67, 0.96, 0.64, 
        0.15, 0.64, 0.73, 0.23, 0.82, 0.35, 0.70, 1.05, 0.73, 0.90, 0.99, 1.14, 1.57, 1.31, 2.77
    ]

    upper_limits = np.array(upper_limits)
    converted_upper_limits = (upper_limits*(10**5)*(h_little**3))

# multiplied the lower errors by -1 so make them positive to avoid the plt.errorbar error 
    y_error = [-1*(converted_lower_limits), converted_upper_limits]


#Let's now actually plot!
    plt.figure(figsize=(9,6))
    plt.plot(redshifts[()],cowd_rate,linewidth=2,linestyle='--',color='mediumblue',label=r'$\mathrm{COWD + WD}$') # all COWD

# NSNS Rate
    plt.plot(redshifts[()],NSNS_rate,linewidth=2,color='grey',alpha=0.7,label='NSNS')

## LVK BNS rate
    plt.fill_between([0.1,0.3], 
                    10,
                    1700, 
                    alpha=0.15, 
                    color="grey")#,label=r'LVK BNS Rate $\mathrm{z=0.2}$')



## seeing if this plot matches Max Briel's paper
    plt.errorbar(redshifts_briel,converted_rates_briel,yerr=y_error, fmt='o', color = 'lightsteelblue', alpha=0.7)#,label='Briel et al. 2020')


## axis
    plt.xlim(0,8)
    plt.ylim(10**0,5*10**5)
    plt.yscale('log')
    plt.ylabel(r"Event Rate, $\mathrm{dNdGpc^{-3}dyr^{-1}}$",fontsize=20)
    plt.xlabel(r"Redshift",fontsize=25)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.title(title, pad=20)
    plt.legend()

# save figure:
    plt.savefig(plot_output + filename +".png",bbox_inches='tight',pad_inches=0.1)


# Closing the HDF5 File
    Data.close()





def metallicity_plotter(pathToH5, title, plot_output, filename):

    """
    Plotting the distribution of metallicites in log10 space for NSNS systems and COWD + WD systems
    pathTOH5 = path to the HDF5 file
    title = the name of the plot that corresponds to the file name
    """

# Read in the data
    Data  = h5.File(pathToH5, "r")


# we need the metallicities, stellar types, masses, mixture weights, merges hubble time
    DCOs = Data['BSE_Double_Compact_Objects']
# gathering the double compact objects that we have computed rates for
    DCO_mask = Data['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['DCOmask'][()]

# first we want to investigate how many of each type of DCO there is
    stellar_types_1 = DCOs['Stellar_Type(1)'][()][DCO_mask]
    stellar_types_2 = DCOs['Stellar_Type(2)'][()][DCO_mask]

# masses
    mass1 = DCOs['Mass(1)'][()][DCO_mask]
    mass2 = DCOs['Mass(2)'][()][DCO_mask]

# mixture weights
    mixture_weights = DCOs['mixture_weight'][()][DCO_mask]

# merges hubble time
    merges_compas = DCOs['Merges_Hubble_Time'][()][DCO_mask]

# metallicities
    metallicities = DCOs['Metallicity@ZAMS(1)'][()][DCO_mask]


# Let's add everything to a dataframe so we can analyze things easier

    data = {
    "Stellar_Type(1)": stellar_types_1,
    "Stellar_Type(2)": stellar_types_2,
    "Mass(1)": mass1,
    "Mass(2)": mass2,
    "mixture_weights": mixture_weights,
    "Merges_Hubble_time": merges_compas,
    "Metallicity@ZAMS(1)": metallicities
    }

    DCOs_masked = pd.DataFrame(data)


# let's make evenly spaced metallicity bins in log 
    metallicities = np.array(DCOs_masked['Metallicity@ZAMS(1)'])
    metallicities_log = np.log10(metallicities)
    bins_Z = np.linspace(-4, np.log10(0.03), 20) # making sure we are considereing the bounds of COMPAS


# Defining what constants to use for the star forming mass per binary
    m1min = min(Data['BSE_System_Parameters']['Mass@ZAMS(1)'][()])
    m1max = max(Data['BSE_System_Parameters']['Mass@ZAMS(1)'][()])
    m2min = Data['Run_Details']["minimum-secondary-mass"][0]

# calling the function for star forming mass per binary 
    analytical_star_forming_mass_per_binary = utils_from_others.analytical_star_forming_mass_per_binary_using_kroupa_imf(m1min, m1max, m2min)



# let's use the "Merges_Hubble_Time" flag to flag when binaries merge and produce graviational waves
    merges_comaps_bool = DCOs_masked['Merges_Hubble_time']==1

# let's create masks for our systems of intersts so we can plot the metallicity dist for each
# COWD+WD systems
    HeWD_bool,COWD_bool,ONeWD_bool,HeCOWD_bool,HeONeWD_bool,COHeWD_bool,COONeWD_bool,ONeHeWD_bool,ONeCOWD_bool = useful_fncs.WD_BINARY_BOOLS(DCOs_masked)
    carbon_oxygen_bool = np.logical_or(ONeCOWD_bool,np.logical_or(COONeWD_bool,np.logical_or(COHeWD_bool,np.logical_or(COWD_bool,HeCOWD_bool))))
    COWD_merged = carbon_oxygen_bool*merges_comaps_bool
# NSNS systems
    NSNS_bool = np.logical_and(DCOs_masked["Stellar_Type(1)"]==14, DCOs_masked['Stellar_Type(2)']==14)
    NSNS_merged = NSNS_bool*merges_comaps_bool


#Plotting!
# let's plot our normalized counts vs metallicities!
    fig, ax = plt.subplots(2, figsize=(10, 7))

    counts_bins_mergers_cowd, bin_edges_mergers_cowd = np.histogram(np.log10(DCOs_masked['Metallicity@ZAMS(1)'][COWD_merged]), weights=DCOs_masked['mixture_weights'][COWD_merged], bins=bins_Z)
    # left edges + right edges/2 = center of each bin
    center_bins_mergers_cowd = (bin_edges_mergers_cowd[:-1] + bin_edges_mergers_cowd[1:])/2
    # we need to consider the bin width so that if we change the normalization constant, the shape stays the same

    bin_width_mergers_cowd = np.diff(bin_edges_mergers_cowd)

    # we now need to see how we can make the counts we get mroe realistic and rep the universe not just what we simulated
    numer_of_binaries_simulated_per_bin = 1e5 #this works because we have systems with metallicites that are uniform in log 
    total_SFM_per_bin = analytical_star_forming_mass_per_binary*numer_of_binaries_simulated_per_bin # total SFM that COMPAS simulation represents in each bin - makes hist more realistic 

    ax[0].step(center_bins_mergers_cowd, (counts_bins_mergers_cowd/total_SFM_per_bin)/bin_width_mergers_cowd, where='mid', label='COWD+WD', color='mediumblue')


    # ax[0].set_xlabel(r"$\mathrm{Log_{10}}(\mathrm{Z)}$")
    ax[0].set_ylabel(r"Merged Systems Per Solar Mass Formed")
    ax[0].legend()
    ax[0].set_title(title, pad=20)



# let's do this again but for NSNS systems
    counts_bins_mergers_NSNS, bin_edges_mergers_NSNS = np.histogram(np.log10(DCOs_masked['Metallicity@ZAMS(1)'][NSNS_merged]), weights=DCOs_masked['mixture_weights'][NSNS_merged], bins=bins_Z)
    # left edges + right edges/2 = center of each bin
    center_bins_mergers_NSNS = (bin_edges_mergers_NSNS[:-1] + bin_edges_mergers_NSNS[1:])/2
    # we need to consider the bin width so that if we change the normalization constant, the shape stays the same

    bin_width_mergers_NSNS = np.diff(bin_edges_mergers_NSNS)

    # we now need to see how we can make the counts we get mroe realistic and rep the universe not just what we simulated
    numer_of_binaries_simulated_per_bin = 1e5 #this works because we have systems with metallicites that are uniform in log 
    total_SFM_per_bin = analytical_star_forming_mass_per_binary*numer_of_binaries_simulated_per_bin # total SFM that COMPAS simulation represents in each bin - makes hist more realistic 

    ax[1].step(center_bins_mergers_NSNS, (counts_bins_mergers_NSNS/total_SFM_per_bin)/bin_width_mergers_NSNS, where='mid', label="NSNS", color='grey')
    ax[1].set_xlabel(r"$\mathrm{Log_{10}}(\mathrm{Z)}$")
    ax[1].set_ylabel(r"Merged Systems Per Solar Mass Formed")
    ax[1].legend()


# save figure:
    plt.savefig(plot_output + filename +".png",bbox_inches='tight',pad_inches=0.1)

# Closing the HDF5 File
    Data.close()




def time_dist_plotter(pathToH5, title, plot_output, filename):

    """
    Plotting the time distirbution plot for NSNS systems and COWD + WD systems
    pathTOH5 = path to the HDF5 file
    title = the name of the plot that corresponds to the file name
    """

# Read in the data
    Data  = h5.File(pathToH5, "r")


# we want to use information in the double compact object group
    DCOs = Data['BSE_Double_Compact_Objects']
    # gathering the double compact objects that we have computed rates for
    DCO_mask = Data['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['DCOmask'][()]

    # first we want to investigate how many of each type of DCO there is
    stellar_types_1 = DCOs['Stellar_Type(1)'][()][DCO_mask]
    stellar_types_2 = DCOs['Stellar_Type(2)'][()][DCO_mask]


    # we need the masses, mixture weight, and rate info
    mass1 = DCOs['Mass(1)'][()][DCO_mask]
    mass2 = DCOs['Mass(2)'][()][DCO_mask]
    mixture_weights = DCOs['mixture_weight'][()][DCO_mask]

    # times
    lifetimes = DCOs['Time'][()][DCO_mask]
    col_times = DCOs['Coalescence_Time'][()][DCO_mask]
    delay_times = DCOs['Time'][()][DCO_mask] + DCOs['Coalescence_Time'][()][DCO_mask]


# let's add everything to a dataframe so we can analyze things easier
    data = {
    "Stellar_Type(1)": stellar_types_1,
    "Stellar_Type(2)": stellar_types_2,
    "Mass(1)": mass1,
    "Mass(2)": mass2,
    "mixture_weight": mixture_weights,
    "Time": lifetimes,
    "Coalescence_Time": col_times,
    "Delay_Time": delay_times
    }

    DCOs_masked = pd.DataFrame(data)


# let's gather our masks for our data so we can plot our systems separately
    HeWD_bool,COWD_bool,ONeWD_bool,HeCOWD_bool,HeONeWD_bool,COHeWD_bool,COONeWD_bool,ONeHeWD_bool,ONeCOWD_bool = useful_fncs.WD_BINARY_BOOLS(DCOs_masked)
    carbon_oxygen_bool = np.logical_or(ONeCOWD_bool,np.logical_or(COONeWD_bool,np.logical_or(COHeWD_bool,np.logical_or(COWD_bool,HeCOWD_bool))))

    NSNS_bool = np.logical_and(DCOs_masked["Stellar_Type(1)"]==14, DCOs_masked['Stellar_Type(2)']==14)



# Let's plot the time distributions!
#Let's first start with the COWD+WD
    ## hubble time
    age_universe = 13.7e9

    fig, ax = plt.subplots(2,figsize=(12,12))

    ## all systems w/ COWD + WD life time
    time_life_log_cowd_wd = np.log10((DCOs_masked['Time'][carbon_oxygen_bool]*(1e6)))
    hist_cowd_wd_life, bin_edges_cowd_wd_life = np.histogram(time_life_log_cowd_wd, weights=DCOs_masked['mixture_weight'][carbon_oxygen_bool],bins=np.linspace(7.5,np.log10(age_universe),50))
    center_bins_cowd_wd_life = (bin_edges_cowd_wd_life[:-1] + bin_edges_cowd_wd_life[1:])/2
    bin_width_cowd_wd_life = np.diff(bin_edges_cowd_wd_life)

    ax[0].plot(center_bins_cowd_wd_life,(hist_cowd_wd_life/bin_width_cowd_wd_life)*1e-3,color='lightcoral',lw=2, label=r'$\mathrm{t_{life}}$')

    ## all systems w/ COWD + WD coalescence time
    time_col_log_cowd_wd = np.log10(((DCOs_masked['Coalescence_Time'][carbon_oxygen_bool]*(1e6))))
    hist_cowd_wd_col, bin_edges_cowd_wd_col = np.histogram(time_col_log_cowd_wd, weights=DCOs_masked['mixture_weight'][carbon_oxygen_bool],bins=np.linspace(7.5,np.log10(age_universe),50))
    center_bins_cowd_wd_col = (bin_edges_cowd_wd_col[:-1] + bin_edges_cowd_wd_col[1:])/2
    bin_width_cowd_wd_col = np.diff(bin_edges_cowd_wd_col)

    ax[0].plot(center_bins_cowd_wd_col,(hist_cowd_wd_col/bin_width_cowd_wd_col)*1e-3,color='lightblue',lw=2, label=r'$\mathrm{t_{inspiral}}$')


    ## all systems w/ COWD + WD delay time
    time_delay_log_cowd_wd = np.log10((DCOs_masked['Delay_Time'][carbon_oxygen_bool]*(1e6)))
    hist_cowd_wd, bin_edges_cowd_wd = np.histogram(time_delay_log_cowd_wd, weights=DCOs_masked['mixture_weight'][carbon_oxygen_bool],bins=np.linspace(7.5,np.log10(age_universe),50))
    center_bins_cowd_wd = (bin_edges_cowd_wd[:-1] + bin_edges_cowd_wd[1:])/2
    bin_width_cowd_wd = np.diff(bin_edges_cowd_wd)

    ax[0].plot(center_bins_cowd_wd,(hist_cowd_wd/bin_width_cowd_wd)*1e-3,color='mediumpurple',lw=2, label=r'$\mathrm{t_{delay}}$')




    ax[0].text(0.03, 1.03, r'$\times 10^3$', fontsize = 15,  transform = ax[0].transAxes)

    # ax[0].set_xlabel(r"$\log_{10}$(t/$\mathrm{Yr}$)",fontsize=15)
    # plt.yscale('log')
    # plt.xlim(0,15e9)
    # plt.ylim(1e-1,1e3)
    ax[0].set_ylabel(r"$\mathrm{dN}$/$\mathrm{d}\log_{10}(t$)",fontsize=15) 
    ax[0].tick_params(axis='both', which='major', labelsize=15)

    ax[0].axvline(np.log10(age_universe), color='r', linestyle='--', linewidth=2)#,label='Hubble Time')
    ax[0].legend(fontsize=15)
    ax[0].set_title(title+" COWD", fontsize = 20, pad=20)



    # Let's do this again for NSNS systems

    ## all systems w/ NSNS life time
    time_life_log_NSNS = np.log10((DCOs_masked['Time'][NSNS_bool]*(1e6)))
    hist_cowd_NSNS_life, bin_edges_cowd_NSNS_life = np.histogram(time_life_log_NSNS, weights=DCOs_masked['mixture_weight'][NSNS_bool],bins=np.linspace(7.5,np.log10(age_universe),50))
    center_bins_cowd_NSNS_life = (bin_edges_cowd_NSNS_life[:-1] + bin_edges_cowd_NSNS_life[1:])/2
    bin_width_cowd_NSNS_life = np.diff(bin_edges_cowd_NSNS_life)

    ax[1].plot(center_bins_cowd_NSNS_life,(hist_cowd_NSNS_life/bin_width_cowd_NSNS_life)*1e-3,color='lightcoral',lw=2, label=r'$\mathrm{t_{life}}$')

    ## all systems w/ NSNS coalescence time
    time_col_log_NSNS = np.log10(((DCOs_masked['Coalescence_Time'][NSNS_bool]*(1e6))))
    hist_cowd_NSNS_col, bin_edges_cowd_NSNS_col = np.histogram(time_col_log_NSNS, weights=DCOs_masked['mixture_weight'][NSNS_bool],bins=np.linspace(7.5,np.log10(age_universe),50))
    center_bins_cowd_NSNS_col = (bin_edges_cowd_NSNS_col[:-1] + bin_edges_cowd_NSNS_col[1:])/2
    bin_width_cowd_NSNS_col = np.diff(bin_edges_cowd_NSNS_col)

    ax[1].plot(center_bins_cowd_NSNS_col,(hist_cowd_NSNS_col/bin_width_cowd_NSNS_col)*1e-3,color='lightblue',lw=2, label=r'$\mathrm{t_{inspiral}}$')


    ## all systems w/ NSNS delay time
    time_delay_log_NSNS_wd = np.log10((DCOs_masked['Delay_Time'][NSNS_bool]*(1e6)))
    hist_NSNS, bin_edges_NSNS = np.histogram(time_delay_log_NSNS_wd, weights=DCOs_masked['mixture_weight'][NSNS_bool],bins=np.linspace(7.5,np.log10(age_universe),50))
    center_bins_NSNS = (bin_edges_NSNS[:-1] + bin_edges_NSNS[1:])/2
    bin_width_NSNS = np.diff(bin_edges_NSNS)

    ax[1].plot(center_bins_NSNS,(hist_NSNS/bin_width_NSNS)*1e-3,color='mediumpurple',lw=2, label=r'$\mathrm{t_{delay}}$')

    ax[1].text(0.03, 1.03, r'$\times 10^3$', fontsize = 15,  transform = ax[1].transAxes)

    ax[1].set_xlabel(r"$\log_{10}$(t/$\mathrm{Yr}$)",fontsize=15)
    # plt.yscale('log')
    # plt.xlim(0,15e9)
    # plt.ylim(1e-1,1e3)
    ax[1].set_ylabel(r"$\mathrm{dN}$/$\mathrm{d}\log_{10}(t$)",fontsize=15) 
    ax[1].tick_params(axis='both', which='major', labelsize=15)

    ax[1].axvline(np.log10(age_universe), color='r', linestyle='--', linewidth=2)#,label='Hubble Time')
    ax[1].legend(fontsize=15)
    ax[1].set_title(title+" NSNS", fontsize = 20, pad=20)

# save figure:
    plt.savefig(plot_output + filename +".png",bbox_inches='tight',pad_inches=0.1)

# Close the HDF5 File
    Data.close()