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




def triangle_plot_fnc(pathToH5):

    """
    Plotting the triangle plot
    pathTOH5 = path to the HDF5 file
    """

# Set the appropriate path to the data file + read in the data

    pathToH5 = pathToH5
    Data  = h5.File(pathToH5, "r")

