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
import regions
import shutil
import drizzlepac

from datetime import datetime, timedelta
from collections import defaultdict
from regions import Regions
from pathlib import Path
import pandas as pd
from tqdm import tqdm
from importlib import reload
from PIL import Image

import vorbin
from vorbin.voronoi_2d_binning import voronoi_2d_binning

import requests
from urllib.parse import quote as urlencode
from mastquery import utils as mastutils
from mastquery import query

import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import wcs as pywcs
from astropy.io import fits
from astropy.cosmology import WMAP9 as cosmo

# grizli stuff
import grizli
from grizli import utils, multifit, fitting, prep, model, jwst_utils, grismconf
from grizli.pipeline import auto_script
from grizli.multifit import MultiBeam
from grizli.aws import visit_processor
from grizli.pipeline.auto_script import get_yml_parameters
from grizli.utils import GTable


import warnings
warnings.filterwarnings("ignore")

import pprint
pp = pprint.PrettyPrinter(indent=4)

#from pywwt.qt import WWTQtClient

HOME = Path.home()
