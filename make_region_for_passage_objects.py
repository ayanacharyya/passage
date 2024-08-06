'''
    Filename: make_region_for_passage_objects.py
    Notes: Make region files for all objects in a given PASSAGE field
    This code was provided by Kalina, who received it from Mason
    Author : Mason, Kalina
    Updated by : Ayan
    Updated on: 01-08-24
    Example: run make_region_for_passage_objects.py

'''

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

#Direct image
hdul = fits.open('Par28_speccat.fits')
cat=hdul[1].data
f = open("Par28v2direct.reg",'a')
for i in range(len(cat['ra'])):
    f.write("WCS;circle("+str(cat['ra'][i])+','+str(cat['dec'][i])+',0.5") # color=green text={'+str(cat['id'][i])+' z='+str(round(cat['redshift'][i],3))+'} font="times 10 bold italic" textangle=30\n')
f.close()

#This and subsequent are for the first order beams. Offsets taken from the config files
f = open("Par28F115r_grism.reg",'a')
w = WCS('Par28_228_f115w-gr150r_drz_sci.fits')
for i in range(len(cat['ra'])):
    ra, dec = (cat['ra'][i], cat['dec'][i])
    x, y = w.all_world2pix(ra, dec, 1)
    f.write("box("+str(x-0.6548681566074263)+','+str(y-33.73739138173772)+',10.0,93.54,0.0) # color=red text={'+str(cat['id'][i])+'} font="times 10 bold italic" textangle=30\n')
f.close()


f = open("Par28F115c_grism.reg",'a')
w = WCS('Par28_228_f115w-gr150c_drz_sci.fits')
for i in range(len(cat['ra'])):
    ra, dec = (cat['ra'][i], cat['dec'][i])
    x, y = w.all_world2pix(ra, dec, 1)
    f.write("box("+str(x-31.91107156101387)+','+str(y-1.3922939626209256)+',97.28751330105166,10.0,0.0) # color=red text={'+str(cat['id'][i])+'} font="times 10 bold italic" textangle=30\n')
f.close()


f = open("Par28F150r_grism.reg",'a')
w = WCS('Par28_228_f150w-gr150r_drz_sci.fits')
for i in range(len(cat['ra'])):
    ra, dec = (cat['ra'][i], cat['dec'][i])
    x, y = w.all_world2pix(ra, dec, 1)
    f.write("box("+str(x-0.6548681566074263)+','+str(y-106.79254657227568)+',10.0,93.54,0.0) # color=red text={'+str(cat['id'][i])+'} font="times 10 bold italic" textangle=30\n')
f.close()


f = open("Par28F150c_grism.reg",'a')
w = WCS('Par28_228_f150w-gr150c_drz_sci.fits')
for i in range(len(cat['ra'])):
    ra, dec = (cat['ra'][i], cat['dec'][i])
    x, y = w.all_world2pix(ra, dec, 1)
    f.write("box("+str(x-96.44444)+','+str(y-0.6548681566074263)+',93.54,10.0,0.0) # color=red text={'+str(cat['id'][i])+'} font="times 10 bold italic" textangle=30\n')
f.close()


f = open("Par28F200r_grism.reg",'a')
w = WCS('Par28_228_f200w-gr150r_drz_sci.fits')
for i in range(len(cat['ra'])):
    ra, dec = (cat['ra'][i], cat['dec'][i])
    x, y = w.all_world2pix(ra, dec, 1)
    f.write("box("+str(x-0.6548681566074263)+','+str(y-204.8370874255101)+',10.0,131.78,0.0) # color=red text={'+str(cat['id'][i])+'} font="times 10 bold italic" textangle=30\n')
f.close()


f = open("Par28F200c_grism.reg",'a')
w = WCS('Par28_228_f200w-gr150c_drz_sci.fits')
for i in range(len(cat['ra'])):
    ra, dec = (cat['ra'][i], cat['dec'][i])
    x, y = w.all_world2pix(ra, dec, 1)
    f.write("box("+str(x-200.9228)+','+str(y-0.6548681566074263)+',127.806,10.0,0.0) # color=red text={'+str(cat['id'][i])+'} font="times 10 bold italic" textangle=30\n')
f.close()