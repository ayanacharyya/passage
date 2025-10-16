'''
    Filename: get_bagpipes_spectra_for_etc.py
    Notes: Extracts and saves the best-fit Bagpipes spectra model for a given object in a way suitable for JWST ETC
    Author : Ayan
    Created: 08-10-25
    Example: run get_bagpipes_spectra_for_etc.py --id 549577
'''

from header import *
from util import *
from run_bagpipes_standalone import generate_fit_params, load_photom_bagpipes

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
def calc_total_snr(filename=None, dirname=None, z=1, wave_cen=6562.819, outdir=None, id=''):
    '''
    Calculates the integrated SNR from files downloaded from JWST ETC
    '''
    if filename is None: filename = Path(dirname) / 'lineplot' / 'lineplot_sn.fits'
    else: filename = Path(filename)

    df_f = Table.read(filename).to_pandas()
    df_n = Table.read(str(filename).replace('flux', 'noise')).to_pandas()

    df_spec = pd.merge(df_f, df_n, on='WAVELENGTH')
    df_spec.rename(columns={'WAVELENGTH':'obswave'}, inplace=True)
    df_spec['restwave'] = df_spec['obswave'] * 1e4 / (1 + z) # Angstrom

    # -----plotting SNR---------------
    fig, ax = plt.subplots(1)
    ax.plot(df_spec['restwave'], df_spec['extracted_flux'], c='b')
    ax.plot(df_spec['restwave'], df_spec['extracted_noise'], c='gray')
    ax.set_xlabel('restframe wavelength (A)')
    ax.set_ylabel('Flux (el/sec)')
    ax.set_title(f'redshift={z:.1f}')
   
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticklabels([f'{item * 1e-4 * (1 + z):.1f}' for item in list(ax2.get_xticks())])
    ax2.set_xlabel('observed wavelength (micron)')

    # ------calculating SNR------------
    vel_fwhm = 500 # km/s
    lambda_fwhm = vel_fwhm * wave_cen / 3e5 # A 
    left_wave = wave_cen - lambda_fwhm
    right_wave = wave_cen + lambda_fwhm

    ax.fill_betweenx([0, ax.get_ylim()[1]], left_wave, right_wave, color='green', alpha=0.2)
    df_chunk = df_spec[df_spec['restwave'].between(left_wave, right_wave)]
    total_snr = np.sum(df_chunk['extracted_flux']) / np.sqrt(np.sum(df_chunk['extracted_noise'] **2))
    print(f'\nTotal SNR within the a {2*lambda_fwhm:.1f} A window around {wave_cen} A line is {total_snr:.1f}')
    ax.text(0.95, 0.95, f'Total SNR across {2*lambda_fwhm:.1f} A = {total_snr:.1f}', ha='right', va='top', transform=ax.transAxes)

    # --------saving bagpipes spectra plot----------
    if outdir is None: outdir = filename.parent
    figname = outdir / Path(f'{id}_{filename.stem}.png')
    fig.savefig(figname)
    print(f'Saved SNR plot at {figname} and .png')
    plt.show(block=False)

    return df_spec

# -------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # --------global variables to change----------
    idcol = 'ID_c2025' # column name of ID in the photometry catalog
    #idcol = 'id' # column name of ID in the photometry catalog
    run = 'for_cy5_c2025' # string label with which the runs would be saved
    args.id = args.id[0]
    #wave_cen = 6560 # corresponding to halpha
    #wave_cen = 4850 # corresponding to hbeta
    wave_cen = 6566 # corresponding to nii

    # ------directory paths------
    photcat_filename_sed = args.output_dir / 'catalogs/zCOSMOS_2025_cy5_targets_for_bagpipe.csv' # choose any *_for_bagpipe.csv file which lists the desired object
    #photcat_filename_sed = args.output_dir / 'catalogs/vandels_uds_aftercuts_single_191447_for_bagpipe.csv' # choose any *_for_bagpipe.csv file which lists the desired object
    filter_dir = args.input_dir / 'COSMOS/transmission_curves' # the transmission curves should be in this directory
    etc_dir = HOME / 'Documents/writings/telescope_proposals/JWSTCy5/ETC' # this is where the ETC related outputs would be stored

    if args.input_etc_dir is None:
        etc_result_path = Path('/Users/acharyya/Downloads')
        etc_result_dict = {549577: 'c2', 432033: 'c4', 270854: 'c5', 191447: 'c6'}
        possible_folders = glob.glob(str(etc_result_path) + f'/wb264545*{etc_result_dict[args.id]}*')
        if len(possible_folders) > 0: args.input_etc_dir = max(possible_folders, key=os.path.getctime)
        else: args.input_etc_dir = etc_result_path / 'dummy'
        print(f'Reading {args.input_etc_dir}..')
    etc_result_filename = Path(args.input_etc_dir) / 'lineplot/lineplot_extracted_flux.fits'

    # ----------reading the photometric catalog for SED fitting-----------------
    df = pd.read_csv(photcat_filename_sed)
    obj = df.iloc[df[df[idcol] == args.id].index.values[0]]

    if not etc_result_filename.exists() or args.clobber:
        # ---------loading bagpipes-related stuff------------------------
        filter_list = [str(filter_dir) + '/' + item[:item.lower().find('_sci')] + '.txt' for item in df.columns if '_sci' in item]
        load_fn = partial(load_photom_bagpipes, phot_cat=photcat_filename_sed, id_colname=idcol, zeropoint=28.9)

        # --------running SED fitting---------------------
        os.chdir(args.output_dir)
        fit_params = generate_fit_params(obj_z=obj['redshift'], z_range=0.01, num_age_bins=5, min_age_bin=30) # Generate the fit parameters
        galaxy = bagpipes.galaxy(ID=int(obj[idcol]), load_data=load_fn, filt_list=filter_list, spectrum_exists=False) # Load the data for this object
        fit = bagpipes.fit(galaxy=galaxy, fit_instructions=fit_params, run=run) # Fit this galaxy
        fit.fit(verbose=True, sampler='nautilus', pool=1)

        # ----------extract bagpipes spectra----------------------------
        fit.posterior.get_advanced_quantities()
        spec = np.percentile(fit.posterior.samples['spectrum_full'], 50, axis=0).T.astype(float) # spec in ergs/s/cm^2/A
        restwave = fit.posterior.model_galaxy.wavelengths # wave in Angstrom
        df_etc = pd.DataFrame({'restwave':restwave, 'flux':spec})
        os.chdir(args.code_dir)
    
        # ----------tailoring to JWST ETC----------------------------
        approx_ha_flux = df_etc[df_etc['restwave'].between(6500, 6600)]['flux'].sum()
        print(f'Approx Ha flux= {approx_ha_flux:.2e} ergs/s/cm2/A')
        df_etc['obswave'] = df_etc['restwave'] * (1 + obj['redshift'])
        df_etc = df_etc[df_etc['obswave'].between(0.9e4, 5e4)] # including 5 micron, so that F44W bandwidth is included, as normalisation is performed on F444W magnitude from NIRCam
        df_etc['flux_mjy'] = df_etc['flux'] * 3.33564095E+04 * (df_etc['restwave'])**2 #mJy
        df_etc['restwave_mu'] = df_etc['restwave'] / 1e4 # micron
        df_etc = df_etc[['restwave_mu', 'flux_mjy']]
        
        # --------saving bagpipes spectra----------
        filename = etc_dir / f'C2025_{obj[idcol]:.0f}_bagpipes_median_spectra_foretc_highres.txt'
        df_etc.to_csv(filename, sep='\t', header=None, index=None)

        # ----------plot bagpipes spectra----------------------------
        fig, ax = plt.subplots(1)
        ax.plot(df_etc['restwave_mu'] * (1 + obj["redshift"]), df_etc['flux_mjy'])
        ax.set_xlabel('observed wavelength (micron)')
        ax.set_ylabel('flux (mJy)')
        ax.set_title(f'redshift={obj["redshift"]:.2f}')
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticklabels([f'{item * 1e4 / (1 + obj["redshift"]):.0f}' for item in list(ax2.get_xticks())])
        ax2.set_xlabel('restframe wavelength (A)')

        # --------saving bagpipes spectra plot----------
        filename = filename.parent / Path(filename.stem + '.png')
        fig.savefig(filename)
        print(f'Saved bagipes spectra at {filename} and .txt')
    else:
        df_spec = calc_total_snr(filename=etc_result_filename, z=obj['redshift'], wave_cen=wave_cen, outdir=etc_dir, id=args.id)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

