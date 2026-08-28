### Py file to gather the redshift rates information given an HDF5 file

import h5py as h5 
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

import useful_fncs
import utils_from_others
import figure_utils


def redshift_rates_info(pathtoh5_NSNS, pathtoh5_WDWD):
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
    condition_mergers = delay_times_DCO < 14100 # Myr

    # gathering the rates data
    rates_DCO = Data_NSNS['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['merger_rate'][()]
    rates_DCO_masked = rates_DCO[condition_mergers]

    redshifts_NSNS = Data_NSNS['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['redshifts'][()]

    # gathering just the DCO objects that merge within a Hubble Time
    stellar_types_all_1 = DCOs_NSNS['Stellar_Type(1)'][()]
    stellar_types_1_DCO = stellar_types_all_1[DCO_mask_NSNS]
    stellar_types_1_merged = stellar_types_1_DCO[condition_mergers]

    stellar_types_all_2 = DCOs_NSNS['Stellar_Type(2)'][()]
    stellar_types_2_DCO = stellar_types_all_2[DCO_mask_NSNS]
    stellar_types_2_merged = stellar_types_2_DCO[condition_mergers]

    # bool for just the NSNS systems
    NSNS_systems_bool = np.logical_and(stellar_types_1_merged==13, stellar_types_2_merged==13)

    # gathering the mixture weight info
    mixture_weights_all = DCOs_NSNS['mixture_weight'][()]
    mixtrue_weights_DCO = mixture_weights_all[DCO_mask_NSNS]
    mixture_weights_merged = mixtrue_weights_DCO[condition_mergers]
    mixture_weights_merged_NSNS = mixture_weights_merged[NSNS_systems_bool]



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
    condition_mergers_WDopt = delay_times_DCO_WDWD < 14100

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
    M1 = np.maximum(mass_1_merged, mass_2_merged)
    M2 = np.minimum(mass_1_merged, mass_2_merged)

    # gathering the mixture weight info
    mixture_weights_all_WDopt = DCOs_WDWD['mixture_weight'][()]
    mixtrue_weights_dco_WDopt = mixture_weights_all_WDopt[DCO_mask_WDWD]
    mixtrue_weights_merged_WDopt = mixtrue_weights_dco_WDopt[condition_mergers_WDopt]

    # let's find the bools for each of our progenitor systems

    # WDWD bool with at least one COWD + WD
    HeWD_bool,COWD_bool,ONeWD_bool,HeCOWD_bool,HeONeWD_bool,COHeWD_bool,COONeWD_bool,ONeHeWD_bool,ONeCOWD_bool = useful_fncs.WD_BINARY_BOOLS(stellar_types_1_merged_WDopt, stellar_types_2_merged_WDopt)
    carbon_oxygen_bool_WDWD_merged_WDopt = np.logical_or(ONeCOWD_bool,np.logical_or(COONeWD_bool,np.logical_or(COHeWD_bool,np.logical_or(COWD_bool,HeCOWD_bool))))

    # COWD + COWD/HeWD with Mchan > 1.4
    tot_mass_cond = mass_1_merged + mass_2_merged > 1.4
    super_chan_bool = carbon_oxygen_bool_WDWD_merged_WDopt*tot_mass_cond

    # D6 HVS
    SN_Ia_HVS,two_star_SNIA,Champagne_Supernova = useful_fncs.check_if_SNIA(mass_1_merged, mass_2_merged)

    # Let's compute the merger rates

    NSNS_rate = np.sum(rates_DCO_masked[NSNS_systems_bool], axis=0)
    cowd_rate = np.sum(rates_DCO_masked_WDopt[carbon_oxygen_bool_WDWD_merged_WDopt], axis=0)
    super_chan_rate = np.sum(rates_DCO_masked_WDopt[super_chan_bool], axis=0)
    HVS_rate = np.sum(rates_DCO_masked_WDopt[SN_Ia_HVS], axis=0)

    rates = np.array([NSNS_rate, cowd_rate, super_chan_rate, HVS_rate])

    # Let's boostrap these rates
    NSNS_rate_2D = rates_DCO_masked[NSNS_systems_bool] # NSNS optimized
    percentiles = useful_fncs.bootstrapping_intervals(NSNS_rate_2D, 50, redshifts_NSNS)

    # WDWD_rate_2D = rates_DCO_masked_WDopt[carbon_oxygen_bool_WDWD_merged_WDopt] # WDWD optimized
    # percentiles_WDWD = useful_fncs.bootstrapping_intervals(WDWD_rate_2D, 50, redshifts_WDWD)

    boostrap_percentiles = np.array([percentiles])#, percentiles_WDWD])

    # Let's compute the ratio of each WDWD subpopulation over redshift
    tot_subpop = cowd_rate + super_chan_rate + HVS_rate
    # the ratios
    cowd_rate_ratio = cowd_rate/tot_subpop
    super_chan_rate_ratio = super_chan_rate/tot_subpop
    HVS_rate_ratio = HVS_rate/tot_subpop
    ratios = np.array([cowd_rate_ratio, super_chan_rate_ratio, HVS_rate_ratio])
    
    return(rates, boostrap_percentiles, ratios)




