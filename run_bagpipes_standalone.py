'''
    Filename: run_bagpipes_standalone.py
    Notes: Computes stellar masses for a given list of galaxies that have fluxes and redshift, and run SED fitting with BAGPIPES
    Author : Ayan
    Created: 20-05-25
    Example: run run_bagpipes_standalone.py
'''
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from pathlib import Path
import os
import bagpipes
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from numpy.typing import ArrayLike
from functools import partial
from uncertainties import unumpy as unp
from uncertainties import ufloat
from scipy.stats import mode as sci_mode

start_time = datetime.now()
cosmo = FlatLambdaCDM(H0=70., Om0=0.3)

# -------------------------------------------------------------------------------------------------------
def load_photom_bagpipes(str_id, phot_cat, id_colname = 'bin_id', zeropoint = 28.9, cat_hdu_index = 0, extra_frac_err = 0.1):
    '''
    Code written by P J Watson
    Load photometry from a catalogue to bagpipes-formatted data.

    The output fluxes and uncertainties are scaled to microJanskys.

    Parameters
    ----------
    str_id : str
        The ID of the object in the photometric catalogue to fit.
    phot_cat : os.PathLike
        The location of the photometric catalogue.
    id_colname : str, optional
        The name of the column containing ``str_id``, by default
        ``"bin_id"``.
    zeropoint : float, optional
        The AB magnitude zeropoint, by default ``28.9``.
    cat_hdu_index : int | str, optional
        The index or name of the HDU containing the photometric catalogue,
        by default ``0``.
    extra_frac_err : float, optional
        An additional fractional error to be added to the photometric
        uncertainties. By default ``extra_frac_err=0.1``, i.e. 10% of the
        measured flux will be added in quadrature to the estimated
        uncertainty.

    Returns
    -------
    ArrayLike
        An Nx2 array containing the fluxes and their associated
        uncertainties in all photometric bands.
    '''

    if not isinstance(phot_cat, Table):
        try: phot_cat = Table.read(phot_cat, hdu=cat_hdu_index)
        except TypeError: phot_cat = Table.read(phot_cat)
    phot_cat[id_colname] = phot_cat[id_colname].astype(str)

    row_idx = (phot_cat[id_colname] == str(str_id)).nonzero()[0][0]
    fluxes = []
    errs = []

    for c in phot_cat.colnames:
        if '_sci' in c.lower(): fluxes.append(phot_cat[c][row_idx])
        elif '_var' in c.lower(): errs.append(np.sqrt(phot_cat[c][row_idx]))
        elif '_err' in c.lower(): errs.append(phot_cat[c][row_idx])

    #if zeropoint == 28.9: flux_scale = 1e-2
    #else: flux_scale = 10 ** ((8.9 - zeropoint) / 2.5 + 6)
    flux_scale = 1

    flux = np.asarray(fluxes) * flux_scale
    flux_err = np.asarray(errs) * flux_scale
    flux_err = np.sqrt(flux_err**2 + (0.1 * flux) ** 2)
    flux = flux.copy()
    flux_err = flux_err.copy()
    bad_values = (~np.isfinite(flux) | (flux <= 0) | ~np.isfinite(flux_err) | (flux_err <= 0))
    flux[bad_values] = 0.0
    flux_err[bad_values] = 1e30

    return np.c_[flux, flux_err]

# -------------------------------------------------------------------------------------------------------
def generate_fit_params(obj_z, z_range = 0.01, num_age_bins = 5, min_age_bin = 30):
    '''
    Code written by P J Watson
    Generate a default set of fit parameters for ``bagpipes``.

    Parameters
    ----------
    obj_z : float | ArrayLike
        The redshift of the object to fit. If ``ArrayLike``, this
        indicates the maximum range of redshifts to fit to.
    z_range : float, optional
        The maximum redshift range to search over, by default 0.01. To fit
        to a single redshift, pass a single value for ``obj_z``, and set
        ``z_range=0.0``. If ``obj_z`` is ``ArrayLike``, this parameter is
        ignored.
    num_age_bins : int, optional
        The number of age bins to fit, each of which will have a constant
        star formation rate following Leja+19. By default ``5`` bins are
        generated.
    min_age_bin : float, optional
        The minimum age of any bin in Myr, by default 30.

    Returns
    -------
    dict
        A dictionary containing the necessary fit parameters for
        ``bagpipes``.
    '''
    fit_params = {}
    if (z_range == 0.0) or (type(obj_z) is ArrayLike): fit_params['redshift'] = obj_z
    else: fit_params['redshift'] = (obj_z - z_range / 2, obj_z + z_range / 2)

    age_at_z = cosmo.age(np.nanmax(fit_params['redshift'])).value

    age_bins = np.geomspace(min_age_bin, age_at_z * 1e3, num=num_age_bins)
    age_bins = np.insert(age_bins, 0, 0.0)

    continuity = {
        'massformed': (6.0, 11.0),
        'metallicity': (0.0, 3.0),
        'metallicity_prior_mu': 1.0,
        'metallicity_prior_sigma': 0.5,
        'bin_edges': age_bins.tolist(),
    }

    for i in range(1, len(continuity['bin_edges']) - 1):
        continuity['dsfr' + str(i)] = (-10.0, 10.0)
        continuity['dsfr' + str(i) + '_prior'] = 'student_t'
        continuity['dsfr' + str(i) + '_prior_scale'] = (
            0.5  # Defaults to 0.3 (Leja19), we aim for a broader sample
        )
        continuity['dsfr' + str(i) + '_prior_df'] = (
            2  # Defaults to this value as in Leja19, but can be set
        )

    fit_params['continuity'] = continuity

    fit_params['dust'] = {
        'type': 'Cardelli',
        'Av': (0.0, 2.0),
        'eta': 2.0,
    }
    fit_params['nebular'] = {'logU': (-3.5, -2.0)}
    fit_params['t_bc'] = 0.02

    return fit_params

# --------------------------------------------------------------------------------------------------------------------
def run_bagpipes(photcat_filename_sed, filter_dir, output_dir, run, idcol='PASSAGE_ID', start_id=0, ncpus=1):
    '''
    Code written by A Acharyya

    Run bagpipes on each object in given phot catalog photcat_filename, save stellar mass in the catalog and write out the resulting catalog
    Returns resulting catalog, containing stellar mass
    '''
    # -------reading and polishing phot catalog------------------
    df = pd.read_csv(photcat_filename_sed)
    filter_list = [str(filter_dir) + '/' + item[:item.lower().find('_sci')] + '.txt' for item in df.columns if '_sci' in item]
    print(f'Resultant photcat has {len(filter_list)} filters')

    # ---------loading bagpipes-related stuff------------------------
    os.chdir(output_dir)
    load_fn = partial(load_photom_bagpipes, phot_cat=photcat_filename_sed, id_colname=idcol, zeropoint=28.9)

    # --------create columns to store stellar masses---------------
    new_columns_dict = {'log_mass_bgp':('stellar_mass', False), 'z_bgp': ('redshift', False), 'log_sfr_bgp':('sfr', True)} # dict of new_label:(quantity_in_bagpipe, needs_to_be_log)
    new_columns = np.hstack([[item, item + '_u'] for item in list(new_columns_dict.keys())])
    for thiscol in new_columns: df[thiscol] = np.zeros(len(df))

    # ---------Loop over the objects-------------
    n_objects = len(df)
    df = df.iloc[start_id:]
        
    for index, obj in df.iterrows():
        print(f'\nLooping over object {index + 1} of {n_objects}..')
        fit_params = generate_fit_params(obj_z=obj['redshift'], z_range=0.01, num_age_bins=5, min_age_bin=30) # Generate the fit parameters

        galaxy = bagpipes.galaxy(ID=int(obj[idcol]), load_data=load_fn, filt_list=filter_list, spectrum_exists=False) # Load the data for this object
        fit = bagpipes.fit(galaxy=galaxy, fit_instructions=fit_params, run=run) # Fit this galaxy
        fit.fit(verbose=True, sampler='nautilus', pool=ncpus)

        # ---------Make some plots---------
        fig, ax = fit.plot_spectrum_posterior(save=True, show=True)
        fig = fit.plot_sfh_posterior(save=True, show=True)
        fig = fit.plot_corner(save=True, show=True)

        # --------Save the stellar masses---------------
        for thisquant in list(new_columns_dict.keys()):
            peak_value = sci_mode(fit.posterior.samples[new_columns_dict[thisquant][0]]).mode
            low_value, up_value = np.percentile(fit.posterior.samples[new_columns_dict[thisquant][0]], (16, 84))
            err_value = (up_value - low_value) / 2
            quant = ufloat(peak_value, err_value)
            if new_columns_dict[thisquant][1]: quant = unp.log10(quant)
            df.loc[index, thisquant] = unp.nominal_values(quant)
            df.loc[index, thisquant + '_u'] = unp.std_devs(quant)

    # ------writing modified df with stellar masses etc-------------------
    df_filename_sed = Path(str(photcat_filename_sed).replace('.csv', f'_withSED_{run}.csv'))
    df.to_csv(df_filename_sed, index=None)
    print(f'Added SED results to df and saved in {df_filename_sed}.')

    return df

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    
    # --------global variables to change----------
    idcol = 'ID_cweb' # column name of ID in the photometry catalog
    run = 'for_cy5' # string label with which the runs would be saved
    ncpus = 1 # BAGPIPES can parallelise, it will use ncpus number of processors

    # ------the transmission curves should be in this directory------
    filter_directory = Path('/Users/acharyya/Work/astro/passage/passage_data/v0.5') / 'COSMOS' / 'transmission_curves'
    # -------this is supposed to have following columns: photometry fluxes (as FILTERNAME_sci), corresponding errors (FILTERNAME_err), redshift (redshift) and object id (idcol)------
    photcat_filename_sed = Path('/Users/acharyya/Work/astro/passage/passage_data/v0.5/COSMOS') / 'zCOSMOS_WFC3+ACS_Web_2020.csv'
    # -------this is where the outputs would be stored-------------
    output_dir = Path('/Users/acharyya/Work/astro/passage/passage_output/v0.5')

    # ----------reading the photometric catalog for SED fitting-----------------
    df = pd.read_csv(photcat_filename_sed)

    # --------running SED fitting---------------------
    dfm = run_bagpipes(photcat_filename_sed, filter_directory, output_dir, run, idcol=idcol, ncpus=ncpus)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

