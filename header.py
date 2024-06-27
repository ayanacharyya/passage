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
#import requests
import regions

from urllib.parse import quote as urlencode
from datetime import datetime, timedelta
from collections import defaultdict
from regions import Regions
from pathlib import Path
import pandas as pd

import matplotlib
from matplotlib import pyplot as plt

from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.io import fits

# grizli stuff
import grizli
import grizli.utils
from grizli.pipeline import auto_script
from grizli import utils, fitting, multifit, prep, model
from grizli.multifit import MultiBeam

import warnings
warnings.filterwarnings("ignore")

import pprint
pp = pprint.PrettyPrinter(indent=4)

#from pywwt.qt import WWTQtClient

HOME = Path.home()
