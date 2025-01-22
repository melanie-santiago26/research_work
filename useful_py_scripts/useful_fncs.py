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

### time coalescence function based on Ilya Mandel "An Accurate Analytical Fit to the Gravitational-wave Inspiral Duration for Eccentric Binaries" (2021)
# Eq 5

# let's now make the definition for the coalescence time (no variations for the very small or very larger eccentiricies)
# this function requires for the untis of each parameter to be inputted

# let's now make the definition for the coalescence time (no variations for the very small or very larger eccentiricies)

def tgw(a,e,Mmoremass,Mlessmass,Data,key,parameter):

    """
    Calcualte the coalescence time (inspiral time)
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

"""
This function will give the maximum separation needed for a bianry to merge within a hubble time given the massesand lifetime of the system
"""
def separations(e,Mmoremass,Mlessmass,t_life):

    """
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
def WD_BINARY_BOOLS(dataframe):

# let's first look at if there are only helium white dwarf WD binaries
    HeWD_bool = np.logical_and(dataframe['Stellar_Type(1)']==10,dataframe['Stellar_Type(2)']==10)
# then carbon oxygen WD
    COWD_bool = np.logical_and(dataframe['Stellar_Type(1)']==11,dataframe['Stellar_Type(2)']==11)
# then oxgen neon WD
    ONeWD_bool = np.logical_and(dataframe['Stellar_Type(1)']==12,dataframe['Stellar_Type(2)']==12)

# let's look at the combination of WD binaries

# Helium WD combos
    HeCOWD_bool = np.logical_and(dataframe['Stellar_Type(1)']==10,dataframe['Stellar_Type(2)']==11)
    HeONeWD_bool = np.logical_and(dataframe['Stellar_Type(1)']==10,dataframe['Stellar_Type(2)']==12)

# Carbon Oxygen WD combos
    COHeWD_bool = np.logical_and(dataframe['Stellar_Type(1)']==11,dataframe['Stellar_Type(2)']==10)
    COONeWD_bool = np.logical_and(dataframe['Stellar_Type(1)']==11,dataframe['Stellar_Type(2)']==12)

# Oxygen Neon WD combos
    ONeHeWD_bool = np.logical_and(dataframe['Stellar_Type(1)']==12,dataframe['Stellar_Type(2)']==10)
    ONeCOWD_bool = np.logical_and(dataframe['Stellar_Type(1)']==12,dataframe['Stellar_Type(2)']==11)

    # let's return all of these bools
    return(HeWD_bool,COWD_bool,ONeWD_bool,HeCOWD_bool,HeONeWD_bool,COHeWD_bool,COONeWD_bool,ONeHeWD_bool,ONeCOWD_bool)


"""
This function is used to get all of the WD+WD bianries from one dataframe rather than the separate bools
"""
def WDWD_bools(dataframe,stellar_type1 = 'Stellar_Type(1)',stellar_type2 = 'Stellar_Type(2)'):

    BWD_BOOL = np.logical_or(np.logical_and(dataframe[stellar_type1]==12,dataframe[stellar_type2]==11),np.logical_or(np.logical_and(dataframe[stellar_type1]==12,dataframe[stellar_type2]==10),np.logical_or(np.logical_and(dataframe[stellar_type1]==11,dataframe[stellar_type2]==12),np.logical_or(np.logical_and(dataframe[stellar_type1]==11,dataframe[stellar_type2]==10),np.logical_or(np.logical_and(dataframe[stellar_type1]==10,dataframe[stellar_type2]==12),np.logical_or(np.logical_and(dataframe[stellar_type1]==10,dataframe[stellar_type2]==11),np.logical_or(np.logical_and(dataframe[stellar_type1]==10,dataframe[stellar_type2]==10),np.logical_or(np.logical_and(dataframe[stellar_type1]==11,dataframe[stellar_type2]==11),np.logical_and(dataframe[stellar_type1]==12,dataframe[stellar_type2]==12)))))))))
    return BWD_BOOL


"""
This function is used to make a selection of binary systems with at least a carbon oxygen WD
"""
def COWD_bool(dataframe,stellar_type1 = 'Stellar_Type(1)',stellar_type2 = 'Stellar_Type(2)'):
    cowd_bool = np.logical_or(dataframe[stellar_type1]==11,dataframe[stellar_type2]==11)
