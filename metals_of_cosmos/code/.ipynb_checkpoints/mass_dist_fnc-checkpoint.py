### Py file to gather the mass distirbution information given an HDF5 file

import h5py as h5 
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

import useful_fncs
import utils_from_others
import figure_utils


def mass_dis_info(pathtoh5_NSNS, pathtoh5_WDWD):
    ## NSNS optimized run

    # let's first look at the NSNS_output
    pathToH5_NSNS = pathtoh5_NSNS

    Data_NSNS  = h5.File(pathToH5_NSNS, "r")

    DCOs_NSNS = Data_NSNS['BSE_Double_Compact_Objects'] # getting the DCO objects

    # gathering the double compact objects that we have computed rates for
    DCO_mask_NSNS = Data_NSNS['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['DCOmask'][()]

    # gathering information to mask the data even more
    # merges in a Hubble Time
    # times (these should be in Myr)
    lifetimes_all = DCOs_NSNS['Time'][()]
    lifetimes_DCO = lifetimes_all[DCO_mask_NSNS]

    col_times_all = DCOs_NSNS['Coalescence_Time'][()]
    col_times_DCO = col_times_all[DCO_mask_NSNS]

    delay_times_DCO = lifetimes_DCO + col_times_DCO
    condition_mergers = delay_times_DCO < 13800 # Myr


    # gathering just the DCO objects that merge within a Hubble Time
    stellar_types_all_1 = DCOs_NSNS['Stellar_Type(1)'][()]
    stellar_types_1_DCO = stellar_types_all_1[DCO_mask_NSNS]
    stellar_types_1_merged = stellar_types_1_DCO[condition_mergers]

    stellar_types_all_2 = DCOs_NSNS['Stellar_Type(2)'][()]
    stellar_types_2_DCO = stellar_types_all_2[DCO_mask_NSNS]
    stellar_types_2_merged = stellar_types_2_DCO[condition_mergers]

    # bool for just the NSNS systems
    NSNS_systems_bool = np.logical_and(stellar_types_1_merged==13, stellar_types_2_merged==13)

    # gathering the masses
    mass_1_all_NSNS = DCOs_NSNS['Mass(1)'][()]
    mass_1_DCO_NSNS = mass_1_all_NSNS[DCO_mask_NSNS]
    mass_1_merged_NSNS = mass_1_DCO_NSNS[condition_mergers]

    mass_2_all_NSNS = DCOs_NSNS['Mass(2)'][()]
    mass_2_DCO_NSNS = mass_2_all_NSNS[DCO_mask_NSNS]
    mass_2_merged_NSNS = mass_2_DCO_NSNS[condition_mergers]

    # we are going to conditions that M1>M2 (not considering mass ratio reversal cases)
    M1_NSNS = np.maximum(mass_1_merged_NSNS, mass_2_merged_NSNS)
    M2_NSNS = np.minimum(mass_1_merged_NSNS, mass_2_merged_NSNS)

    M1 = M1_NSNS[NSNS_systems_bool]
    M2 = M2_NSNS[NSNS_systems_bool]

    # gathering the local rates data
    rates_z0_DCO_NSNS = Data_NSNS['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['merger_rate_z0'][()]
    rates_z0_merged_hubble = rates_z0_DCO_NSNS[condition_mergers]
    rates_z0_merged_NSNS = rates_z0_merged_hubble[NSNS_systems_bool]



    # Let's do the same for the WDWD systems

    ## WDWD optimized run
    pathToH5_WDWD = pathtoh5_WDWD

    Data_WDWD  = h5.File(pathToH5_WDWD, "r")

    DCOs_WDWD = Data_WDWD['BSE_Double_Compact_Objects'] # getting the DCO objects

    # gathering the double compact objects that we have computed rates for
    DCO_mask_WDWD = Data_WDWD['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['DCOmask'][()]

    DATA_SPS_WDWD = Data_WDWD['BSE_System_Parameters']

    # gathering information to mask the data even more
    # merges in a Hubble Time
    lifetimes_all_WDWD = DCOs_WDWD['Time'][()]
    lifetimes_DCO_WDWD = lifetimes_all_WDWD[DCO_mask_WDWD]

    col_times_all_WDWD = DCOs_WDWD['Coalescence_Time'][()]
    col_times_DCO_WDWD = col_times_all_WDWD[DCO_mask_WDWD]

    delay_times_DCO_WDWD = lifetimes_DCO_WDWD + col_times_DCO_WDWD
    condition_mergers_WDopt = delay_times_DCO_WDWD < 13800

    # gathering the rates data
    rates_DCO_WDopt = Data_WDWD['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['merger_rate'][()]
    rates_DCO_masked_WDopt = rates_DCO_WDopt[condition_mergers_WDopt]

    redshifts_WDWD = Data_WDWD['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['redshifts'][()]

    stellar_types_all_1_WDopt = DCOs_WDWD['Stellar_Type(1)'][()]
    stellar_types_1_DCO_WDopt = stellar_types_all_1_WDopt[DCO_mask_WDWD]
    stellar_types_1_merged_WDopt = stellar_types_1_DCO_WDopt[condition_mergers_WDopt]

    stellar_types_all_2 = DCOs_WDWD['Stellar_Type(2)'][()]
    stellar_types_2_DCO = stellar_types_all_2[DCO_mask_WDWD]
    stellar_types_2_merged_WDopt = stellar_types_2_DCO[condition_mergers_WDopt]

    # gathering the masses
    mass_1_all = DCOs_WDWD['Mass(1)'][()]
    mass_1_DCO = mass_1_all[DCO_mask_WDWD]
    mass_1_merged = mass_1_DCO[condition_mergers_WDopt]

    mass_2_all = DCOs_WDWD['Mass(2)'][()]
    mass_2_DCO = mass_2_all[DCO_mask_WDWD]
    mass_2_merged = mass_2_DCO[condition_mergers_WDopt]

    # we are going to conditions that M1>M2 (not considering mass ratio reversal cases)
    M1_DWD = np.maximum(mass_1_merged, mass_2_merged)
    M2_DWD = np.minimum(mass_1_merged, mass_2_merged)

    # let's find the bools for each of our progenitor systems

    # WDWD bool with at least one COWD + WD
    HeWD_bool,COWD_bool,ONeWD_bool,HeCOWD_bool,HeONeWD_bool,COHeWD_bool,COONeWD_bool,ONeHeWD_bool,ONeCOWD_bool = useful_fncs.WD_BINARY_BOOLS(stellar_types_1_merged_WDopt, stellar_types_2_merged_WDopt)
    carbon_oxygen_bool_WDWD_merged_WDopt = np.logical_or(ONeCOWD_bool,np.logical_or(COONeWD_bool,np.logical_or(COHeWD_bool,np.logical_or(COWD_bool,HeCOWD_bool))))

    # gather the local rates data
    rates_z0_DCO = Data_WDWD['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['merger_rate_z0'][()]
    rates_z0_merged = rates_z0_DCO[condition_mergers_WDopt]
    rates_z0_merged_COWD = rates_z0_merged[carbon_oxygen_bool_WDWD_merged_WDopt]

    local_rates = [rates_z0_merged_NSNS, rates_z0_merged_COWD]

    M1_COWD = M1_DWD[carbon_oxygen_bool_WDWD_merged_WDopt]
    M2_COWD = M2_DWD[carbon_oxygen_bool_WDWD_merged_WDopt]

    masses = [M1, M2, M1_COWD, M2_COWD]

    Data_NSNS.close()
    Data_WDWD.close()
    
    return(local_rates, masses)




