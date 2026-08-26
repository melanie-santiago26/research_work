"""This py file will serve as a place where all of the functions I have created
for this project will be stored"""

# general imports that are needed for these functions and future ones
import h5py as h5  
import pandas as pd
from astropy.table import Table
from astropy import units as u
from astropy import constants as const
import numpy as np
import matplotlib.pyplot as plt
import utils_from_others

### time coalescence function based on Ilya Mandel "An Accurate Analytical Fit to the Gravitational-wave Inspiral Duration for Eccentric Binaries" (2021)
# Eq 5

# let's now make the definition for the coalescence time (no variations for the very small or very larger eccentiricies)
# this function requires for the untis of each parameter to be inputted

# let's now make the definition for the coalescence time (no variations for the very small or very larger eccentiricies)

def tgw(a,e,Mmoremass,Mlessmass,Data,key,parameter):

    """
    Calcualte the coalescence time (inspiral time) in Myrs
    a = semi major axis (expected in AU or Rsun)
    e = eccentricity
    Mmoremass = mass of the more massive compact object (expected in solar masses!)
    Mlessmass = mass of the less massive compact object (expected in solar masses!)
    """

    SYS = Data[key]
    sep_unit = SYS[parameter].attrs['units']
    if sep_unit == b'Rsol':
        a = (a * u.Rsun).to(u.m)

    elif sep_unit == b'AU':
        a = (a * u.AU).to(u.m)      

    Mmoremass = (Mmoremass * u.Msun).to(u.kg)
    Mlessmass = (Mlessmass * u.Msun).to(u.kg)

    tc = ((((5*((a)**4)*(const.c**5))/(256*(const.G**3)*(Mmoremass)*(Mlessmass)*((Mmoremass)+(Mlessmass))))*(1+(0.27*e**10)+(0.33*e**20)+(0.2*e**1000))*(1-(e**2))**(7/2)))*((3.171e-8)*(u.yr/u.s))*((1e-6)*(u.Myr/u.yr))

    return tc.values


def separations(e,Mmoremass,Mlessmass,t_life):

    """
    This function will give the maximum separation needed for a bianry to merge within a hubble time given the massesand lifetime of the system
    e = eccentricity
    t_hubble = age of the universe (expected in Myr)
    Mmoremass = mass of the more massive compact object (expected in solar masses!)
    Mlessmass = mass of the less massive compact object (expected in solar masses!)
    """

    age_universe = (13.7e9*u.yr).to(u.s)

    if t_life > age_universe.value:
        return print("The age of your binary surpasses that of the age of the universe (ypur binary is still forming).")

    else:

        Mmoremass = (Mmoremass * u.Msun).to(u.kg)
        Mlessmass = (Mlessmass * u.Msun).to(u.kg)

        a_min_num = (1/(5*(const.c**5)))*((age_universe-((t_life*u.Myr).to(u.s)))*(256*(const.G**3)*Mmoremass*Mlessmass*(Mmoremass+Mlessmass)))
        a_min_den = ((1+(0.27*e**10)+(0.33*e**20)+(0.2*e**1000))*(1-(e**2))**(7/2))**(1/4)
        a_min_final = ((a_min_num/a_min_den)**(1/4)).to(u.Rsun)


        return (a_min_final).value


### This functions serves as a tool to select for WD+WD systems
def WD_BINARY_BOOLS(stellar_type_1, stellar_type_2):

# let's first look at if there are only helium white dwarf WD binaries
    HeWD_bool = np.logical_and(stellar_type_1==10,stellar_type_2==10)
# then carbon oxygen WD
    COWD_bool = np.logical_and(stellar_type_1==11,stellar_type_2==11)
# then oxgen neon WD
    ONeWD_bool = np.logical_and(stellar_type_1==12,stellar_type_2==12)

# let's look at the combination of WD binaries

# Helium WD combos
    HeCOWD_bool = np.logical_and(stellar_type_1==10,stellar_type_2==11)
    HeONeWD_bool = np.logical_and(stellar_type_1==10,stellar_type_2==12)

# Carbon Oxygen WD combos
    COHeWD_bool = np.logical_and(stellar_type_1==11,stellar_type_2==10)
    COONeWD_bool = np.logical_and(stellar_type_1==11,stellar_type_2==12)

# Oxygen Neon WD combos
    ONeHeWD_bool = np.logical_and(stellar_type_1==12,stellar_type_2==10)
    ONeCOWD_bool = np.logical_and(stellar_type_1==12,stellar_type_2==11)

    # let's return all of these bools
    return(HeWD_bool,COWD_bool,ONeWD_bool,HeCOWD_bool,HeONeWD_bool,COHeWD_bool,COONeWD_bool,ONeHeWD_bool,ONeCOWD_bool)



def WDWD_bools(dataframe,stellar_type1 = 'Stellar_Type(1)',stellar_type2 = 'Stellar_Type(2)'):
    """
    This function is used to get all of the WD+WD bianries from one dataframe rather than the separate bools
    """
    BWD_BOOL = np.logical_or(np.logical_and(dataframe[stellar_type1]==12,dataframe[stellar_type2]==11),np.logical_or(np.logical_and(dataframe[stellar_type1]==12,dataframe[stellar_type2]==10),np.logical_or(np.logical_and(dataframe[stellar_type1]==11,dataframe[stellar_type2]==12),np.logical_or(np.logical_and(dataframe[stellar_type1]==11,dataframe[stellar_type2]==10),np.logical_or(np.logical_and(dataframe[stellar_type1]==10,dataframe[stellar_type2]==12),np.logical_or(np.logical_and(dataframe[stellar_type1]==10,dataframe[stellar_type2]==11),np.logical_or(np.logical_and(dataframe[stellar_type1]==10,dataframe[stellar_type2]==10),np.logical_or(np.logical_and(dataframe[stellar_type1]==11,dataframe[stellar_type2]==11),np.logical_and(dataframe[stellar_type1]==12,dataframe[stellar_type2]==12)))))))))
    return BWD_BOOL



def WDWD_bools_from_array(stellar_type1,stellar_type2):
    """
    This function is used to get all of the WD+WD bianries from one stellar types array
    """
    BWD_BOOL = np.logical_or(np.logical_and(stellar_type1==12,stellar_type2==11),np.logical_or(np.logical_and(stellar_type1==12,stellar_type2==10),np.logical_or(np.logical_and(stellar_type1==11,stellar_type2==12),np.logical_or(np.logical_and(stellar_type1==11,stellar_type2==10),np.logical_or(np.logical_and(stellar_type1==10,stellar_type2==12),np.logical_or(np.logical_and(stellar_type1==10,stellar_type2==11),np.logical_or(np.logical_and(stellar_type1==10,stellar_type2==10),np.logical_or(np.logical_and(stellar_type1==11,stellar_type2==11),np.logical_and(stellar_type1==12,stellar_type2==12)))))))))
    return BWD_BOOL



def COWD_bool(dataframe,stellar_type1 = 'Stellar_Type(1)',stellar_type2 = 'Stellar_Type(2)'):
    """
    This function is used to make a selection of binary systems with at least a carbon oxygen WD
    """
    cowd_bool = np.logical_or(dataframe[stellar_type1]==11,dataframe[stellar_type2]==11)
    return cowd_bool


def line(x,slope,b):
    """
    This functions is to define the lines used as boundaries in check_if_SNIA.
    """
    y = slope*x +b
    return(np.array(y)) 


# let's test the function for just the red region
def check_if_SNIA(mass1,mass2):

    """
    This functions takes systems and assigns flags to what type of SN Ia it is if any.
    The boundaries for the mass cuts and SN Ia categorizations come from Shen 2025: https://arxiv.org/abs/2502.04451. 
    """
    # let's select the masses from the compas output
    M_more_massive = np.maximum(mass1,mass2)
    M_less_massive = np.minimum(mass1,mass2)

    # # empty flag arrays that we will add to our dataset later
    # SN_Ia_HVS = np.empty_like(M_more_massive)
    # two_star_SNIA = np.empty_like(M_more_massive)
    # Champagne_Supernova = np.empty_like(M_more_massive)

    # let's now make regimes based on Ken Shen 2025
    # red region, we define the border cases to be read from left to right so a region does not take systems that are on the left border but it does take its right border
    # Mass 1 condition
    red_more_massive_bool = np.logical_and(M_more_massive<1.0,
                                           M_more_massive>=M_less_massive)
    # Mass 2 condition
    red_less_massive_bool = np.logical_and(M_less_massive>=line(M_more_massive,-1.5,1.875),
                                       M_less_massive<line(M_more_massive,-1.5,2.0))

    # let's now mask which masses fall within the red region
    SN_Ia_HVS = red_more_massive_bool*red_less_massive_bool


    # purple region
    # Mass 1 condition
    purple_more_massive_bool = np.logical_and(M_more_massive<1.0,
                                           M_more_massive>=M_less_massive)
    # Mass 2 condition
    purple_less_massive_bool = M_less_massive>=line(M_more_massive,-1.5,2.0)
    
    # let's now mask which masses fall within the purple region
    two_star_SNIA = purple_more_massive_bool*purple_less_massive_bool


    # orange region
    # Mass 1 condition
    orange_more_massive_bool = np.logical_and(M_more_massive<1.1,
                                           M_more_massive>=1.0)
    # Mass 2 condition
    orange_less_massive_bool = np.logical_and(M_less_massive>line(M_more_massive,-1.0,1.8),M_less_massive<=M_more_massive)
    
    # let's now mask which masses fall within the purple region
    Champagne_Supernova = orange_more_massive_bool*orange_less_massive_bool

    return(SN_Ia_HVS,two_star_SNIA,Champagne_Supernova)



# making this a function

def bootstrapping_intervals(rate_2D, boostraps_num, redshifts):

    """
    This function creates the condifence intervals on the rates through bootstrapping inspired by Lieke van Son's code
    https://github.com/LiekeVanSon/LowMBH_and_StableChannel/blob/master/Code/Fig1_MassDistributions.ipynb
    rate_2D = the 2 dimensional array of the the rates that has already been masked
    boostraps_num = number of times you would like to bootstrap
    redshifts = you array of redhshift bins
    """
    
    indices = np.arange(len(rate_2D)) # indicies is same size of what we want to bootstrap (aka number of systems)
    rates_DCO_boots = np.zeros((boostraps_num, len(redshifts))) # we want to start with an empty 2D array 

    for b in utils_from_others.progressbar(range(len(rates_DCO_boots)), "Bootstrapping :"): # looping through each redshift bin

        boot_index = np.random.choice(indices, size=len(indices), replace=True)
        boots_rate = rate_2D[boot_index] # taking a random index corresponding to a particular system and getting what the rate is at each redshift
        boots_instance = np.sum(boots_rate, axis=0) # now actually computing the rate at each redshift 

        rates_DCO_boots[b] = boots_instance # adding our rate at each redshift to the previous empty array of zeros 

    percentiles = np.percentile(rates_DCO_boots, [5., 50., 95.], axis=0)

    return(percentiles)



def merger_rate_z0_result_WDWD(pathToH5):

    """
    This function gathers the local rate for each subpopulation in the WDWD optimized output of COMPAS
    pathToH5 = file path to the COMPAS output hdf5 file 
    """
    
    # we first need to gather the information for the file of interest
    # reading in the HDF5 file

    Data  = h5.File(pathToH5, "r")

    # To make the variations plot we need the stellar types, masses, and rate at z=0

    DCOs = Data['BSE_Double_Compact_Objects'] # gathering the DCO group
    DCO_mask = Data['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['DCOmask'][()]


    # merges in a Hubble Time

    lifetimes_all = DCOs['Time'][()]
    lifetimes_DCO = lifetimes_all[DCO_mask]

    col_times_all = DCOs['Coalescence_Time'][()]
    col_times_DCO = col_times_all[DCO_mask]

    delay_times_DCO = lifetimes_DCO + col_times_DCO
    condition_mergers = delay_times_DCO < 14100 # Myr

    
    # HDF5 files are most efficent if you apply the mask after reading in the key of interest

    stellar_types_1_all = DCOs['Stellar_Type(1)'][()]
    stellar_types_1 = stellar_types_1_all[DCO_mask]
    stellar_type_1_merged = stellar_types_1[condition_mergers]

    stellar_types_2_all = DCOs['Stellar_Type(2)'][()]
    stellar_types_2 = stellar_types_2_all[DCO_mask]
    stellar_type_2_merged = stellar_types_2[condition_mergers]
    
    mass1_all = DCOs['Mass(1)'][()]
    mass1 = mass1_all[DCO_mask]
    mass1_merged = mass1[condition_mergers]

    mass2_all = DCOs['Mass(2)'][()]
    mass2 = mass2_all[DCO_mask]
    mass2_merged = mass2[condition_mergers]
    
    rates_z0_DCO = Data['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['merger_rate_z0'][()]
    rates_z0_DCO_merged = rates_z0_DCO[condition_mergers]

    Data.close()

    # let's make sure that at least one of these white dwarfs are COWD
    HeWD_bool,COWD_bool,ONeWD_bool,HeCOWD_bool,HeONeWD_bool,COHeWD_bool,COONeWD_bool,ONeHeWD_bool,ONeCOWD_bool = WD_BINARY_BOOLS(stellar_type_1_merged, stellar_type_2_merged)
    carbon_oxygen_bool = np.logical_or(ONeCOWD_bool,np.logical_or(COONeWD_bool,np.logical_or(COHeWD_bool,np.logical_or(COWD_bool,HeCOWD_bool))))


    # let's now sort our data into the mass regimes we care about
    # let's add the flags for specific calssifications of SN Ia

    SN_Ia_HVS,two_star_SNIA,Champagne_Supernova = check_if_SNIA(mass1_merged[carbon_oxygen_bool],mass2_merged[carbon_oxygen_bool])

    # let's now gather the merger rate at redshift zero
    WDWD_merger_rate_Z0 = rates_z0_DCO_merged[carbon_oxygen_bool]
    cowd_rate = np.sum(WDWD_merger_rate_Z0)

    # # let's now get the values of these merger rates for all of the systems that fall within this specific regime
    two_star_SNIA_z0_rate = np.sum(WDWD_merger_rate_Z0[two_star_SNIA==True])


    #let's get the merger rates for Mtot>mchan
    mtot_chan_more = mass1_merged + mass2_merged > 1.4
    mtot_chan_bool = carbon_oxygen_bool*mtot_chan_more
    mchan_rate = np.sum(rates_z0_DCO_merged[mtot_chan_bool])


    return([cowd_rate,mchan_rate,two_star_SNIA_z0_rate])


def merger_rate_z0_result_NSNS(pathToH5):

    """
    This function gathers the local rate for each subpopulation in the NSNS optimized output of COMPAS
    pathToH5 = file path to the COMPAS output hdf5 file 
    """
    
    # we first need to gather the information for the file of interest
    # reading in the HDF5 file

    Data  = h5.File(pathToH5, "r")

    # To make the variations plot we need the stellar types, masses, and rate at z=0

    DCOs = Data['BSE_Double_Compact_Objects'] # gathering the DCO group
    DCO_mask = Data['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['DCOmask'][()]

    # merges in a Hubble Time

    lifetimes_all = DCOs['Time'][()]
    lifetimes_DCO = lifetimes_all[DCO_mask]

    col_times_all = DCOs['Coalescence_Time'][()]
    col_times_DCO = col_times_all[DCO_mask]

    delay_times_DCO = lifetimes_DCO + col_times_DCO
    condition_mergers = delay_times_DCO < 14100 # Myr

    
    # HDF5 files are most efficent if you apply the mask after reading in the key of interest

    stellar_types_1_all = DCOs['Stellar_Type(1)'][()]
    stellar_types_1 = stellar_types_1_all[DCO_mask]
    stellar_type_1_merged = stellar_types_1[condition_mergers]

    stellar_types_2_all = DCOs['Stellar_Type(2)'][()]
    stellar_types_2 = stellar_types_2_all[DCO_mask]
    stellar_type_2_merged = stellar_types_2[condition_mergers]

    rates_z0_DCO = Data['Rates_mu00.025_muz-0.049_alpha-1.79_sigma01.129_sigmaz0.048']['merger_rate_z0'][()]
    rates_z0_DCO_merged = rates_z0_DCO[condition_mergers]

    Data.close()

    # the NSNS merger rate
    # let's gather the NSNS merger information
    NSNS_bool = np.logical_and(stellar_type_1_merged==13,stellar_type_2_merged==13)
    nsns_rate = np.sum(rates_z0_DCO_merged[NSNS_bool],axis=0)

    return(nsns_rate)

