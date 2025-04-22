#!/usr/bin/env python3

"""

    Filename :   header.py
    Notes :      Header file for importing packages/modules and declaring global variables required for working with PASSAGE data.
    Author :    Ayan
    Created: 11-06-24
    Last modified: 11-06-24

"""

from __future__ import print_function

import numpy as np
import argparse
import os
import copy
import importlib
import glob
import sys
import re
import json
import shutil
import drizzlepac
import subprocess
import itertools
import h5py
import io
import colorsys

from datetime import datetime, timedelta
from collections import defaultdict
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from importlib import reload
from PIL import Image
from venn import venn, pseudovenn, generate_petal_labels, draw_venn, generate_colors
from reproject import reproject_interp, reproject_adaptive
from scipy.stats import mode as sci_mode
from scipy.interpolate import NearestNDInterpolator
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit

import regions
from regions import Regions, PixCoord
from regions import RectangleSkyRegion, RectanglePixelRegion

import vorbin
from vorbin.voronoi_2d_binning import voronoi_2d_binning

import requests
from bs4 import BeautifulSoup
from urllib.parse import quote as urlencode
from urllib.request import urlretrieve
from mastquery import utils as mastutils
from mastquery import query
from astroquery.mast import Observations
import plotly.express as px
import gzip

from astropy.table import Table, join
from astropy import units as u
from astropy.coordinates import SkyCoord, angular_separation
from astropy.nddata.utils import Cutout2D
from astropy import wcs as pywcs
from astropy.io import fits
from astropy.cosmology import Planck18
from astropy.cosmology import FlatLambdaCDM
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel


# grizli stuff
import grizli
from grizli import utils, multifit, fitting, prep, model, jwst_utils, grismconf
from grizli.pipeline import auto_script
from grizli.multifit import MultiBeam
from grizli.aws import visit_processor
from grizli.pipeline.auto_script import get_yml_parameters
from grizli.utils import GTable

import datashader as dsh
from datashader.mpl_ext import dsshow
from datashader import transfer_functions as dstf

from uncertainties import unumpy as unp
from uncertainties import ufloat

from PIL import Image
import fitz

import warnings
warnings.filterwarnings("ignore")

import pprint
pp = pprint.PrettyPrinter(indent=4)

from os import PathLike
from numpy.typing import ArrayLike
from functools import partial
import bagpipes

#from pywwt.qt import WWTQtClient

from NebulaBayes import NB_Model
import NebulaBayes

import matplotlib
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.image as mpimg
from matplotlib import colors as mplcolors
from matplotlib import colormaps as mplcolormaps
from matplotlib import cm as mpl_cm
from matplotlib.backends.backend_pdf import PdfPages
import mplcyberpunk
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker

from matplotlib import pyplot as plt
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.top'] = True

HOME = Path.home()

rest_wave_dict = {'Lya': 121.6, 'OII': 372.7, 'NeIII-3867':386.9, 'Hd': 434.0, 'OIII-4363': 436.3, 'Hb': 486.1, 'OIII': 500.7,
                  'Ha+NII': 655.5, 'Ha': 656.2, 'SII': 671.7, 'ArIII-7138':713.6,
                  'SIII': 953.0, 'PaD': 1004.6, 'PaG': 1093.5, 'PaB': 1281.4,
                  'PaA': 1874.5}  # approximate wavelengths in nm

# Below: NIRISS filter data taken from Table 1 in https://jwst-docs.stsci.edu/jwst-near-infrared-imager-and-slitless-spectrograph/niriss-instrumentation/niriss-filters#gsc.tab=0
filter_waverange_dict = {'F090W': [0.796, 1.005], 'F115W': [1.013, 1.283], 'F150W': [1.330, 1.671],
                         'F200W': [1.751, 2.226], 'F277W': [2.413, 3.143], 'F356W': [3.140, 4.068],
                         'F444W': [3.880, 5.023], 'F140M': [1.331, 1.480], 'F158M': [1.488, 1.688],
                         'F380M': [3.726, 3.931], 'F430M': [4.182, 4.395], 'F480M': [4.668, 4.971]}  # microns

ids_to_re_extract_dict = {'Par051': [141, 159, 179, 207, 220, 239, 285, 304, 315, 333, 337, 383, 393, 406, 416, 420, 461, 520, 539, 549, 576,
               585, 586, \
               601, 607, 620, 672, 673, 700, 709, 733, 735, 744, 751, 761, 785, 806, 807, 820, 851, 858, 891, 899, 934,
               944, \
               1007, 1025, 1030, 1044, 1104, 1111, 1121, 1153, 1224, 1265, 1282, 1283, 1290, 1293, 1294, 1350, 1389,
               1395, 1396, \
               1411, 1413, 1440, 1441, 1456, 1474, 1479, 1485, 1509, 1510, 1555, 1556, 1559, 1593, 1594, 1641, 1647, 1650,
               1676, 1687, \
               1708, 1709, 1721, 1740, 1748, 1772, 1777, 1796, 1810, 1824, 1825, 1843, 1844, 1857, 1864, 1901, 1902,
               1923, 1946, 1965, \
               2007, 2029, 2045, 2072, 2092, 2093, 2103, 2109, 2137, 2138, 2169, 2170, 2182, 2188, 2190, 2244, 2260,
               2262, 2275, \
               2314, 2320, 2326, 2327, 2364, 2372, 2373, 2374, 2395, 2413, 2420, 2438, 2539, 2478, 2480, 2481, 2504,
               2515, 2516, 2584, \
               2655, 2772,
               2777]} # IDs determined by visual inspection, to re-extract the spectra using wider (0.1 < z < 8) z constraints, and perhaps better choices for other parameters too

fields_with_2PA = ['Par027', 'Par028', 'Par006', 'Par051', 'Par009']

passage_fields_in_cosmos = [3, 5, 6, 17, 20, 23, 24, 25, 26, 28, 29, 46, 47, 48, 49, 51, 52, 53]

fields_uploaded_by_nuo = [17, 23, 24, 25, 26, 49, 53]

# Set up necessary variables for cosmological calculations.
cosmo = FlatLambdaCDM(H0=69.5, Om0=0.285, Ob0=0.0461)
