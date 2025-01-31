'''
    Filename: narrow_band.py
    Notes: Computes narrow band image from a given 3D IFU data cube (for Ramona)
    Author : RA, AA
    Created: 30-01-25
    Example: run narrow_band.py
'''
from datetime import datetime, timedelta
import numpy as np
from astropy.io import fits
from astropy import wcs as pywcs
from matplotlib import pyplot as plt
from pathlib import Path
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def plot_imhsow(map, ax, label='', cmap='Grays', percentiles=[1, 99], hide_cbar=False):
    '''
    Plots given 2D map on a given axis handle with cmap stretched within given percentiles
    Returns axis handle
    '''

    img = ax.imshow(map, cmap=cmap, vmin=np.percentile(map[np.isfinite(map)], percentiles[0]), vmax=np.percentile(map[np.isfinite(map)], percentiles[1]), origin='lower')
    if not hide_cbar:
        cb = plt.colorbar(img)
        cb.set_label(label)

    return ax
# --------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    #plt.close('all')
    # Variables to change
    mxdf_cube_filename = Path.home() / 'Downloads/DATACUBE_MXDF.fits'
    wavelength_start = 5446.9 # loscat['fitparam_mu1'].min()  # in Angstroms (example)
    wavelength_end = 5465.1 # loscat['fitparam_mu1'].max()  # in Angstroms (example)
    do_weight = True

    # Read in datacube
    mxdf_cube = fits.open(mxdf_cube_filename)
    mxdf_cube_data = mxdf_cube[1].data
    weights_cube = 1. / mxdf_cube_data
    header = mxdf_cube[1].header
    wcs = pywcs.WCS(mxdf_cube[1].header, naxis=2)

    # Determine wavelength array and pixel indices for slicing
    wave_arr = np.array([header['CD3_3'] * (item - header['CRPIX3']) + header['CRVAL3'] for item in range(np.shape(mxdf_cube_data)[0])]) # Wavelength array corresponding to the datacube's [0] axis, in Angstrom
    px_start = np.where(wave_arr >= wavelength_start)[0][0]
    px_end = np.where(wave_arr <= wavelength_end)[0][-1]
    npix_cont = px_end - px_start # no. of pixels on left and right side of the window, for continuum

    # Full integrated data
    full_integrated_map = np.nansum(mxdf_cube_data, axis=0)

    # Slice the cube to extract the narrowband data
    narrowband_map = np.average(mxdf_cube_data[px_start:px_end, :, :], axis=0, weights=weights_cube[px_start:px_end, :, :] if do_weight else None) # Sum the data along the spectral axis, because sum = sum of all flux within that wavelength window
    cont_left_map = np.nanmedian(mxdf_cube_data[px_start - npix_cont:px_start, :, :], axis=0)
    cont_right_map = np.nanmedian( mxdf_cube_data[px_end:px_end + npix_cont, :, :], axis=0)
    cont_mean_map = np.nanmean([cont_left_map, cont_right_map], axis=0) # Mean along the wavelength axis
    narrowband_map_contsub = narrowband_map - cont_mean_map * (px_end - px_start) # Multiplying the mean continuum by delta_pix gives the total continuum within the wavelength window

    # Make the plots
    print(f'Collapsed 3D cube of shape {np.shape(mxdf_cube_data)} from pixels {px_start} ({wavelength_start} A) to {px_end} ({wavelength_end} A) into a shape of {np.shape(narrowband_map)} and subtracted continuum from a {npix_cont}-pixels wide window on either side')
    fig = plt.figure(figsize=(8,6))
    plt.subplot(projection=wcs)

    ax = plt.gca()
    #ax = plot_imhsow(np.log10(full_integrated_map), ax, label='log(full integrated image)', cmap='Greys_r', percentiles=[1, 99])
    ax = plot_imhsow(narrowband_map_contsub, ax, label='narrow band image', cmap='Greys_r', percentiles=[1, 99], hide_cbar=True)

    df = pd.DataFrame({'ra':[53.1625000], 'dec': [-27.7922222], 'vel':[78]})
    scat = ax.scatter(df['ra'], df['dec'], c=df['vel'], vmin=-100, vmax=100, cmap='RdBu', transform=ax.get_transform('world'))
    cb = plt.colorbar(scat)
    cb.set_label('Velocity (km/s)')

    # Save and display figure
    figname = mxdf_cube_filename.parent / Path(mxdf_cube_filename.stem + '_nb_image.png')
    fig.savefig(figname)
    print(f'Saved figure as {figname}')
    plt.show(block=False)
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
