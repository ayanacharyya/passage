import numpy as np
from pathlib import Path
import pandas as pd
from functools import partial
import bagpipes
import sys
from astropy.table import Table
from os import PathLike
from numpy.typing import ArrayLike
from compute_stellar_masses import load_photom_bagpipes, generate_fit_params
from astropy.cosmology import FlatLambdaCDM, Planck18

cosmo = FlatLambdaCDM(H0=69.5, Om0=0.285, Ob0=0.0461)

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

    row_idx = (phot_cat[id_colname] == str_id).nonzero()[0][0]
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
    flux_err = np.sqrt(flux_err**2 + (extra_frac_err * flux) ** 2)
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
if __name__ == "__main__":
    filter_dir = Path('/Users/acharyya/Work/astro/passage/passage_data/v0.5/COSMOS/transmission_curves')
    photcat_filename_sed = Path('/Users/acharyya/Work/astro/passage/passage_output/v0.5/catalogs/Par028_v0.5_venn_OII,NeIII-3867,Hb,OIII,SNR>2.0,mass_passage_cosmos_fluxes_wCWebb_for_bagpipe_for_paper_only_st.csv')

    df_sed = pd.read_csv(photcat_filename_sed)
    filter_list = [str(filter_dir) + '/' + item[:item.lower().find('_sci')] + '.txt' for item in df_sed.columns if '_sci' in item]
    obj = df_sed[df_sed['objid'] == 1303].iloc[0]

    load_fn = partial(load_photom_bagpipes, phot_cat=photcat_filename_sed, id_colname='objid', zeropoint=28.9)
    fit_params = generate_fit_params(obj_z=obj['redshift'], z_range=0.0, num_age_bins=5, min_age_bin=30)
    
    galaxy = bagpipes.galaxy(ID=obj['objid'], load_data=load_fn, filt_list=filter_list, spectrum_exists=False) # Load the data for this object

    fit = bagpipes.fit(galaxy=galaxy, fit_instructions=fit_params, run='test3') # Fit this galaxy
    fit.fit(verbose=True, sampler='nautilus', n_live=400, pool=1, n_eff=10000)
    plot_restframe = True

    # ---------converting to restframe-----------------
    if plot_restframe:
        fit.posterior.get_advanced_quantities()
        redshift = np.median(fit.posterior.samples['redshift'])
        fit.galaxy.photometry[:, 0] = fit.galaxy.photometry[:, 0] / (1 + redshift)
        if fit.galaxy.spectrum_exists: fit.galaxy.spectrum[:, 0] = fit.galaxy.spectrum[:, 0] / (1 + redshift)
        fit.galaxy.filter_set.eff_wavs = fit.galaxy.filter_set.eff_wavs / (1 + redshift)
        fit.posterior.model_galaxy.wavelengths = fit.posterior.model_galaxy.wavelengths / (1 + redshift)
    
    # ---------Make some plots---------
    fig, ax = fit.plot_spectrum_posterior(save=True, show=True, log_x=True, xlim=[2.8, 4.5], ylim=[0, 4])
    fig, ax = fit.plot_spectrum_posterior(save=True, show=True, log_x=False, xlim=[500, 30000], ylim=[0, 4])
    fig = fit.plot_sfh_posterior(save=True, show=True, xlim=None, ylim=[0, 10])
    fig = fit.plot_corner(save=True, show=True)
