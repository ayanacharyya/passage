'''
    Filename: ifu_analysis.py
    Notes: Perform emission line fitting, make emission line maps, and derived properties' maps from a given IFU datacube
    Author : Ayan
    Created: 29-10-25
    Example: run ifu_analysis.py --input_file /Users/acharyya/Work/astro/jwst_uvit_cospar/project/science_ready_data/5001-GS/f170lp_g235h-f17.fits
'''

from header import *
from util import *
from make_diagnostic_maps import compute_EB_V, compute_SFR, compute_Z_C19, take_safe_log_ratio, take_safe_log_sum, get_linelist, plot_linelist, plot_2D_map, plot_MAPPINGS_lines, get_dereddened_flux

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def save_fig(fig, figname, args, silent=False):
    '''
    Saves a given figure handle as a given output filename
    '''

    if args.fortalk:
        #mplcyberpunk.add_glow_effects()
        #try: mplcyberpunk.make_lines_glow()
        #except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass
        plot_output_dir = args.root_dir / 'zgrad_paper_plots_fortalk'
    else:
        plot_output_dir = args.root_dir / 'zgrad_paper_plots'

    if 'pjw' in args.drv: plot_output_dir = Path(str(plot_output_dir) + '_pjw')

    plot_output_dir.mkdir(exist_ok=True, parents=True)
    figname = plot_output_dir / figname
    fig.savefig(figname, transparent=args.fortalk)
    if not silent: print(f'Saved figure as {figname}')
    plt.show(block=False)

    return

# --------------------------------------------------------------------------------------------------------------------
def continuum_subtract(df_spec_input, wavelengths_to_fit, contmask_lambda_window=200, npix=10):
    '''
    Function for fitting (with a spline) and subtracting the continuum, from an input dataframe of spectra
    Returns the dataframe with continuum subtracted spectra, and an additional column 'cont' carrying the continuum
    '''
    df_spec = df_spec_input.copy()
    # -------mask around lines of interest-------
    for wavelength in wavelengths_to_fit:
        left_window = wavelength - contmask_lambda_window
        right_window = wavelength + contmask_lambda_window
        df_spec = df_spec[(df_spec['restwave'] < left_window) | (df_spec['restwave'] > right_window)]

    # ------spline fit the continuum--------
    f = interp1d(df_spec['restwave'].values, df_spec['flux'].values, fill_value='extrapolate')
    df_spec_input['cont'] = f(df_spec_input['restwave'].values)
    df_spec_input['flux'] = df_spec_input['flux'] / df_spec_input['cont']
    df_spec_input['flux_u'] = df_spec_input['flux_u'] / df_spec_input['cont']

    return df_spec_input

# --------------------------------------------------------------------------------------------------------------------
def linefit_function(x, wavelengths_to_fit, *popt):
    '''
    Function for line fitting: a continuum + a set of emission lines (input as df_linelist in rest)
    Returns the flux array
    '''
    c = 3e5 # km/s
    model = popt[0]
    for index, wave in enumerate(wavelengths_to_fit):
        model += popt[3 + index] * np.exp(-((x - wave * (1. + popt[1]/c)) ** 2) / (2 * (wave * popt[2] / c / 2.35482) ** 2))
    
    return model

# --------------------------------------------------------------------------------------------------------------------
def linefit_1dspec(df_spec, df_linelist, contmask_lambda_window=200, xpix=0, ypix=0, args=None, silent=True):
    '''
    Fits a set of emission lines (input as df_linelist, in restframe) simultaenously to a given spectrum (input as df_spec)
    Returns a dataframe containing line fluxes, systemtic velocities and velocity dispersions, and corresponding uncertainties
    '''
    c = 3e5 # km/s
    # ------continuum subtraction-------------
    df_spec = df_spec[(df_spec['flux'] > 0) & (df_spec['flux_u'] > 0) & (df_spec['flux_u'] != 999) & (~np.isnan(df_spec['flux']))]
    if len(df_spec) > 0:
        wavelengths_to_fit = df_linelist['restwave'].values
        df_spec = continuum_subtract(df_spec, wavelengths_to_fit, contmask_lambda_window=contmask_lambda_window)

        # ------line fitting-------------
        p_init = np.hstack([1, 200, 500, [10] * len(wavelengths_to_fit)]) # [continuum, vel, vdisp, [flux] * nlines]
        p_low = np.hstack([0, 0, 111, [0] * len(wavelengths_to_fit)]) # [continuum, vel, vdisp, [flux] * nlines]
        p_high = np.hstack([1.3, np.inf, np.inf, [np.inf] * len(wavelengths_to_fit)]) # [continuum, vel, vdisp, [flux] * nlines]
        try:
            popt, pcov = curve_fit(lambda x, *popt: linefit_function(x, wavelengths_to_fit, *popt), df_spec['restwave'].values, df_spec['flux'].values, p0=p_init, sigma=df_spec['flux_u'].values, absolute_sigma=True, bounds=(p_low, p_high))
            bad_fit = False
        except:
            print(f'Line fit did not work for this pixel; hence returning NaNs for this pixel.')
            bad_fit = True
        
        # -----------plotting the linefit-----------------
        if args is not None:
            norm_factor = 1
            df_spec['cont'] = 1
            if args.fortalk: col_arr = ['salmon', 'peachpuff', 'lawngreen', 'cyan']
            else: col_arr = ['k', 'gray', 'crimson', 'cornflowerblue']
            fig, ax = plt.subplots(figsize=(14, 6))
            fig.subplots_adjust(left=0.07, right=0.98, top=0.9, bottom=0.13)
            
            # ----plotting the observed spectra in restframe------
            ax.step(df_spec['restwave'], (df_spec['flux'] * df_spec['cont']) / norm_factor, lw=2, c=col_arr[0], label='Data')
            # ax.step(df_spec['restwave'], df_spec['cont'] / norm_factor, lw=0.5, c=col_arr[1], label='Modeled cont.')
            ax.fill_between(df_spec['restwave'], (df_spec['flux'] - df_spec['flux_u']/2) * df_spec['cont'] / norm_factor, (df_spec['flux'] + df_spec['flux_u']/2) * df_spec['cont'] / norm_factor, lw=0, color=col_arr[0], alpha=0.5, step='pre')#, drawstyle='steps')
            
            # -----plotting the fitted spectra---------
            if not bad_fit:
                yarr = linefit_function(df_spec['restwave'], wavelengths_to_fit, *popt)
                ax.plot(df_spec['restwave'], (yarr * df_spec['cont']) / norm_factor, lw=1, c=col_arr[2], label='Modeled spectra')
        
            # ---------annotating plots------
            ax.set_xlabel(r'Rest-frame wavelength ($\AA$)', fontsize=args.fontsize)
            #ax.set_ylabel(r'f$_{\lambda}$ ' + '(%.0e ' % norm_factor + r'ergs/s/cm$^2$/A)', fontsize=args.fontsize)
            ax.set_ylabel(r'f$_{\lambda}$ ' + '(continuum normalised)', fontsize=args.fontsize)
            if args.flam_max is None: args.flam_max = 10
            ax.set_ylim(0, args.flam_max) # flam_max should be in units of 1e-19 ergs/s/cm^2/A
            if args.xmin is None: args.xmin = np.min(df_spec['restwave'])
            if args.xmax is None: args.xmax = np.max(df_spec['restwave'])
            ax.set_xlim(args.xmin, args.xmax)
            ax.legend(fontsize=args.fontsize, loc='upper left')
            ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

            # -----plotting velocity windows for masking lines----------
            if not args.fortalk:
                for wavelength in wavelengths_to_fit:
                    left_window = wavelength - contmask_lambda_window
                    right_window = wavelength + contmask_lambda_window
                    df_spec = df_spec[(df_spec['restwave'] < left_window) | (df_spec['restwave'] > right_window)]
                    ax.fill_betweenx([ax.get_ylim()[0], ax.get_ylim()[1]], left_window, right_window, color=col_arr[1], lw=0, alpha=0.3)
                    #ax.axvline(wavelength, color=col_arr[3], lw=1, ls='dashed', alpha=0.8)
                
            # ---observed wavelength axis-------
            ax2 = ax.twiny()
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticks(ax.get_xticks())
            ax2.set_xticklabels(['%.2F' % (item * (1 + args.z) / 1e4) for item in ax2.get_xticks()], fontsize=args.fontsize)
            ax2.set_xlabel(r'Observed wavelength ($\mu$)', fontsize=args.fontsize)

            # ---vertical lines for emission line wavelengths------
            if args.plot_mappings: ax = plot_MAPPINGS_lines(ax)
            else: ax = plot_linelist(ax, fontsize=args.fontsize)

            bad_fit_text = '_failed_fit' if bad_fit else ''
            out_dir = args.output_fig_dir / 'linefits'
            out_dir.mkdir(parents=True, exist_ok=True)
            figname = out_dir / f'linefit_pixel_{xpix:.0f}_{ypix:.0f}{bad_fit_text}.png'
            save_fig(fig, figname, args, silent=True)
            plt.close()
        
        # -----------collate and convert all fitted parameters--------
        if not bad_fit:
            flux = [np.sqrt(2 * np.pi) * popt[3 + index] * (wave * popt[2] / c / 2.35482) for index, wave in enumerate(wavelengths_to_fit)]
            flux_u = [np.sqrt(2 * np.pi) * np.sqrt(pcov[3 + index][3 + index]) * (wave * np.sqrt(pcov[2][2]) / c / 2.35482) for index, wave in enumerate(wavelengths_to_fit)]
            vel, vel_u = popt[1], np.sqrt(pcov[1][1])
            vdisp, vdisp_u = popt[1], np.sqrt(pcov[2][2])
            df_fit = pd.DataFrame({'LineID': df_linelist['LineID'].values, 'flux': flux, 'flux_u': flux_u, 'vel': vel, 'vel_u': vel_u, 'vdisp': vdisp, 'vdisp_u': vdisp_u})
        else:
            df_fit = pd.DataFrame({'LineID': df_linelist['LineID'].values, 'flux': [np.nan] * len(df_linelist), 'flux_u': [np.nan] * len(df_linelist), 'vel': [np.nan] * len(df_linelist), 'vel_u': [np.nan] * len(df_linelist), 'vdisp': [np.nan] * len(df_linelist), 'vdisp_u': [np.nan] * len(df_linelist)})  
    else:
        df_fit = pd.DataFrame({'LineID': df_linelist['LineID'].values, 'flux': [np.nan] * len(df_linelist), 'flux_u': [np.nan] * len(df_linelist), 'vel': [np.nan] * len(df_linelist), 'vel_u': [np.nan] * len(df_linelist), 'vdisp': [np.nan] * len(df_linelist), 'vdisp_u': [np.nan] * len(df_linelist)})
        print(f'This pixel had {len(df_spec)} spectral elements to fit lines to; hence returning NaNs for this pixel.')
    return df_fit

# --------------------------------------------------------------------------------------------------------------------
def linefit_cube(data_cube, err_cube, rest_wave_arr, df_linelist, contmask_lambda_window=200, args=None):
    '''
    Fits a set of emission lines (input as df_linelist, in restframe) simultaenously to a given IFU datacube (input as data_cube)
    Returns three data cubes (x_coord x y_coord x n_lines) containing line fluxes, systemtic velocities and velocity dispersions
    '''
    nlines = len(df_linelist)
    print(f'\nFitting {nlines} lines to the full datacube of shape {np.shape(data_cube)}..')
    nwave, nrows, ncols = np.shape(data_cube)
    flux_cube, vel_cube, vdisp_cube = np.tile(unp.uarray(0, 0), (nlines, nrows, ncols)), np.tile(unp.uarray(0, 0), (nlines, nrows, ncols)), np.tile(unp.uarray(0, 0), (nlines, nrows, ncols))

    val_cube = data_cube.data
    err_cube = err_cube.data

    # ------looping over spaxels------
    for i in range(nrows):#[26:27]:
        for j in range(ncols):#[22:23]:
            if args is not None and not args.silent: print(f'Fitting spaxel ({i}, {j}) which is {i * ncols + j + 1} out of {nrows * ncols}..')
            df_spec = pd.DataFrame({'restwave': rest_wave_arr, 'flux': val_cube[:, i, j], 'flux_u': err_cube[:, i, j]})
            df_fit = linefit_1dspec(df_spec, df_linelist, contmask_lambda_window=contmask_lambda_window, xpix=i, ypix=j, args=args, silent=False)

            flux_cube[:, i, j] = unp.uarray(df_fit['flux'].values, df_fit['flux_u'].values)
            vel_cube[:, i, j] = unp.uarray(df_fit['vel'].values, df_fit['vel_u'].values)
            vdisp_cube[:, i, j] = unp.uarray(df_fit['vdisp'].values, df_fit['vdisp_u'].values)

    return flux_cube, vel_cube, vdisp_cube

# --------------------------------------------------------------------------------------------------------------------
def save_fitted_cube(flux_cube, vel_cube, vdisp_cube, df_linelist, outfilename, args):
    '''
    Saves the input three data cubes (corresponding to fitted line properties) as a fits cube
    '''
    hdr1, hdr2 = fits.Header(), fits.Header()
    hdr1['filename'] = Path(args.input_file).stem
    hdr1['redshift'] = args.z
    hdr1['NLINES'] = len(df_linelist)
    hdr1['LINES'] = ','.join(df_linelist['LineID'].values)
    hdr2['FLUX_UNIT'] = 'ergs/s/cm^2'
    hdr2['VEL_UNIT'] = 'km/s'
    primary_hdu = fits.PrimaryHDU(header=hdr1)

    hdul = [primary_hdu]
    for index, row in df_linelist.iterrows():
        line_label = row['LineID']
        this_flux_hdu = fits.ImageHDU(data=unp.nominal_values(flux_cube[index, :, :]), name=f'FLUX_{line_label}', header=hdr2)
        this_flux_err_hdu = fits.ImageHDU(data=unp.std_devs(flux_cube[index, :, :]), name=f'FLUXERR_{line_label}', header=hdr2)
        this_vel_hdu = fits.ImageHDU(data=unp.nominal_values(vel_cube[index, :, :]), name=f'VEL_{line_label}', header=hdr2)
        this_vel_err_hdu = fits.ImageHDU(data=unp.std_devs(vel_cube[index, :, :]), name=f'VELERR_{line_label}', header=hdr2)
        this_vdisp_hdu = fits.ImageHDU(data=unp.nominal_values(vdisp_cube[index, :, :]), name=f'VDISP_{line_label}', header=hdr2)
        this_vdisp_err_hdu = fits.ImageHDU(data=unp.std_devs(vdisp_cube[index, :, :]), name=f'VDISPERR_{line_label}', header=hdr2)

        hdul = np.hstack([hdul, [this_flux_hdu, this_flux_err_hdu, this_vel_hdu, this_vel_err_hdu, this_vdisp_hdu, this_vdisp_err_hdu]])

    hdul = fits.HDUList(list(hdul))
    hdul.writeto(outfilename, overwrite=True)

    print(f'Saved metallicity maps in {outfilename}')

# --------------------------------------------------------------------------------------------------------------------
def get_emission_line_map(line_label, full_hdu, EB_V=0, ext='FLUX', dered=True, silent=True):
    '''
    Extract the emission line flux or vel or vdisp (depending on the provided ext) of a given line from an input FITS HDUlist
    Returns 2D unumpy.uarray with flux map and corresponding uncertainty map 
    '''
    val_map = full_hdu[f'{ext}_{line_label}'].data
    err_map = full_hdu[f'{ext}ERR_{line_label}'].data
    flux_map = unp.uarray(val_map, err_map)

    # -----------getting the dereddened flux value-----------------
    if dered:
        line_wavelength = df_linelist[df_linelist['LineID']==line_label]['restwave'].values[0]
        line_map_quant = get_dereddened_flux(flux_map, line_wavelength, EB_V)
        val_map = unp.nominal_values(line_map_quant)
        err_map = unp.std_devs(line_map_quant)
        flux_map = unp.uarray(val_map, err_map)

    return flux_map

# --------------------------------------------------------------------------------------------------------------------
def get_emission_line_maps(full_hdu, line_labels, EB_V=0, dered=True, silent=True):
    '''
    Computes and returns the spatially resolved as well as intregrated metallicity based on a given Curti+2019 calibration, from a given HDU
    '''
    line_map_arr = []
    for line in line_labels:
        line_map = get_emission_line_map(line, full_hdu, EB_V=EB_V, dered=dered, silent=silent)
        line_map_arr.append(line_map)

    return line_map_arr

# --------------------------------------------------------------------------------------------------------------------
def get_EB_V(full_hdu):
    '''
    Computes and returns the spatially resolved as well as integrated dust extinction map from a given HDU
    Based on Eqn 4 of Dominguez+2013 (https://iopscience.iop.org/article/10.1088/0004-637X/763/2/145/pdf)
    '''

    ha_map = get_emission_line_map('H6562', full_hdu, dered=False) # do not need to deredden the lines when we are fetching the flux in order to compute reddening
    hb_map = get_emission_line_map('HBeta', full_hdu, dered=False)

    EB_V_map = compute_EB_V(ha_map, hb_map)

    return EB_V_map

# --------------------------------------------------------------------------------------------------------------------
def get_SFR(full_hdu, distance, EB_V=0.):
    '''
    Computes and returns the spatially resolved as well as intregrated SFR from a given HDU
    '''
    ha_map = get_emission_line_map('H6562', full_hdu, EB_V=EB_V)

    SFR_map = compute_SFR(ha_map, distance)

    return SFR_map

# --------------------------------------------------------------------------------------------------------------------
def get_Z_C19(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated metallicity based on a given Curti+2019 and Cataldi+2025 calibrations, from a given HDU
    '''
    # ------getting appropriate emission lines and calibration coefficients--------------
    if args.Zdiag == 'O3O2':
        line_map_arr = get_emission_line_maps(full_hdu, ['OIII5007', 'OII3727'], EB_V=args.EB_V, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        ratio_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1])
        if args.use_C25: coeff = [0.01124, 0.03072, 0.1251, -0.01470] # from Table 4 of Cataldi+25
        else: coeff = [-0.691, -2.944, -1.308]  # c0-2 parameters from Table 2 of Curti+19 3rd row (O3O2)

    elif args.Zdiag == 'R3':
        line_map_arr = get_emission_line_maps(full_hdu, ['OIII5007', 'HBeta'], EB_V=args.EB_V, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        ratio_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1])
        if args.use_C25: coeff = [-3.587, -4.475, 1.345, -0.08951] # from Table 4 of Cataldi+25
        else: coeff = [-0.277, -3.549, -3.593, -0.981]  # c0-3 parameters from Table 2 of Curti+19 2nd row (R3)

    elif args.Zdiag == 'R23':
        line_map_arr = get_emission_line_maps(full_hdu, ['OIII5007', 'OII3727', 'HBeta'], EB_V=args.EB_V, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        R1_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[2], skip_log=True)
        R2_map = take_safe_log_ratio(line_map_arr[1], line_map_arr[2], skip_log=True)
        ratio_map = take_safe_log_sum(R1_map, R2_map)
        if args.use_C25: coeff = [-2.555, -3.192, 0.9630, -0.06351] # from Table 4 of Cataldi+25
        else: coeff = [0.527, -1.569, -1.652, -0.421]  # c0-3 parameters from Table 2 of Curti+19 4th row (R23)

    elif args.Zdiag == 'O3S2':
        line_map_arr = get_emission_line_maps(full_hdu, ['OIII5007', 'HBeta', 'SII6730', 'H6562'], EB_V=args.EB_V, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        R1_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1], skip_log=True)
        R2_map = take_safe_log_ratio(line_map_arr[2], line_map_arr[3], skip_log=True)
        ratio_map = take_safe_log_ratio(R1_map, R2_map)
        coeff = [0.191, -4.292, -2.538, 0.053, 0.332]  # c0-4 parameters from Table 2 of Curti+19 last row (O3S2)

    elif args.Zdiag == 'S2':
        line_map_arr = get_emission_line_maps(full_hdu, ['SII6730', 'H6562'], EB_V=args.EB_V, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        ratio_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1])
        coeff = [-0.442, -0.360, -6.271, -8.339, -3.559]  # c0-3 parameters from Table 2 of Curti+19 3rd-to-last row (S2)

    elif args.Zdiag == 'RS32':
        line_map_arr = get_emission_line_maps(full_hdu, ['OIII5007', 'HBeta', 'SII6730', 'H6562'], EB_V=args.EB_V, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        R1_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1], skip_log=True)
        R2_map = take_safe_log_ratio(line_map_arr[2], line_map_arr[3], skip_log=True)
        ratio_map = take_safe_log_sum(R1_map, R2_map)
        coeff = [-0.054, -2.546, -1.970, 0.082, 0.222]  # c0-3 parameters from Table 2 of Curti+19 2nd-to-last row (RS32)

    elif args.Zdiag == 'R2':
        line_map_arr = get_emission_line_maps(full_hdu, ['OII3727', 'HBeta'], EB_V=args.EB_V, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        ratio_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1])
        if args.use_C25: coeff = [-6.481, 0.8163] # from Table 4 of Cataldi+25
        else: coeff = [0.435, -1.362, -5.655, -4.851, -0.478, 0.736]  # c0-3 parameters from Table 2 of Curti+19 1st row (R2)

    else:
        print(f'Could not apply any of the metallicity diagnostics, so returning NaN metallicities')
        return None, ufloat(np.nan, np.nan), ufloat(np.nan, np.nan)

    coeff = coeff[::-1] # because in Curti+2019 the coefficients are listed in the reverse order compared to what np.poly1d prefers

    # ----setting up Z diag debugging plots---------
    if args.debug_Zdiag:
        fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
        ax[0].set_xlabel('Possible values of log(O/H)+12', fontsize=args.fontsize)
        if args.Zbranch == 'high': ax[1].set_xlabel('log(O/H)+12 = max(solution)', fontsize=args.fontsize)
        elif args.Zbranch == 'low': ax[1].set_xlabel('log(O/H)+12 = min(solution)', fontsize=args.fontsize)
        ax[0].set_ylabel(f'Observed log {args.Zdiag}', fontsize=args.fontsize)
        allowed_Z_limit = [7.0, 8.4] if args.use_C25 else [7.6, 8.9] # Z limits within which each calibration is valid
        Z_limits = [6.5, 9.5]
        ratio_limits = [-0., 1.5]

        metallicity_offset = 0 if args.use_C25 else 8.69
        xarr_valid = np.linspace(allowed_Z_limit[0], allowed_Z_limit[1], 100)
        ax[0].plot(xarr_valid, np.poly1d(coeff)(xarr_valid - metallicity_offset), lw=3, c='k', ls='solid', label='Valid C25 calibration' if args.use_C25 else 'Valid C20 calibration')
        ax[1].plot(xarr_valid, np.poly1d(coeff)(xarr_valid - metallicity_offset), lw=3, c='k', ls='solid')

        xarr_full = np.linspace(Z_limits[0], Z_limits[1], 100)
        ax[0].plot(xarr_full, np.poly1d(coeff)(xarr_full - metallicity_offset), lw=0.5, c='k', ls='dotted', label='Extrapolated C25 calibration' if args.use_C25 else 'Extrapolated C20 calibration')
        ax[1].plot(xarr_full, np.poly1d(coeff)(xarr_full - metallicity_offset), lw=0.5, c='k', ls='dotted')

        # for ax1 in ax:
        #     ax1.fill_betweenx([-5, 5], Z_limits[0], allowed_Z_limit[0], color='grey', lw=0, alpha=0.2, label='Calibration invalid')
        #     ax1.fill_betweenx([-5, 5], allowed_Z_limit[1], Z_limits[1], color='grey', lw=0, alpha=0.2)

        ax[0].set_ylim(ratio_limits[0], ratio_limits[1])
        ax[0].set_xlim(Z_limits[0], Z_limits[1])
        ax[1].set_xlim(Z_limits[0], Z_limits[1])
    else:
        ax = None

    # -------estimating the metallicities---------------
    logOH_map = compute_Z_C19(ratio_map, coeff, ax=ax, branch=args.Zbranch) if not args.only_integrated else None

    # -------saving the debugging plots---------------
    if ax is not None:
        Zbranch_text = '' if args.Zdiag in ['NB', 'P25', 'Te'] else f'-{args.Zbranch}'
        figname = args.fig_dir / f'{args.field}_{args.id:05d}_metallicity_debug_Zdiag_{args.Zdiag}{Zbranch_text}.png'
        fig.savefig(figname, transparent=args.fortalk, dpi=200)
        print(f'\nSaved figure at {figname}')

    return logOH_map

# --------------------------------------------------------------------------------------------------------------------
def background_subtract(filename, filename_bkg_sub, args, bkg_region=[[1,5], [1, 5]], smooth=True):
    '''
    Extract the median spectrum within the given region, subtract that from the spectra of each pixel of the given filename
    Save the resulting 3D data cube as filename_bkg_sub
    Return the 3D data cube
    '''
    print(f'\nPerforming background subtraction within pixel region of {bkg_region}..')
    # ------reading in original input IFU data cube---------------
    primary_header, flux_cube, err_cube, flux_header, err_header, wave_arr = read_datacube(filename, args, return_components=True)

    # ------getting the bkg spectra---------
    bkg_cube = flux_cube[:, bkg_region[0][0] : bkg_region[0][1], bkg_region[1][0] : bkg_region[1][1]]
    bkg_flux = np.nanmedian(np.nanmedian(unp.nominal_values(bkg_cube.data), axis=2), axis=1)

    # ----------saving the bkg spectra--------
    bkg_filename = filename.parent / 'fits' / 'bg.csv'
    bkg_df = pd.DataFrame({'obswave': wave_arr, 'bg_flux': bkg_flux})
    bkg_df.to_csv(bkg_filename)
    print(f'Saved background spectra as {bkg_filename}')
 
    # --------performing bkg subtraction---------
    nwave, nrow, ncol = np.shape(flux_cube)
    bkg_cube = np.tile(bkg_flux, nrow * ncol).reshape(nwave, nrow, ncol)
    
    flux_cube_bkg_sub = flux_cube - bkg_cube
    err_cube_bkg_sub = err_cube - bkg_cube

    # ------performing kernel smoothing-------------
    if smooth:
        print(f'\nPerforming gaussian smoothing with {args.kernel_size} pixel kernel..')
        smoothing_kernel = Box2DKernel(args.kernel_size, mode=args.kernel_mode)
        for index in range(nwave):
            print(f'Smoothing wavelength slice {index + 1} of {nwave}...')
            flux_cube_bkg_sub[index, :, :] = convolve(flux_cube_bkg_sub[index, :, :], smoothing_kernel)
            err_cube_bkg_sub[index, :, :] = convolve(err_cube_bkg_sub[index, :, :], smoothing_kernel)

    # ------writing the bkg subtracted fits file--------
    primary_hdu = fits.PrimaryHDU(header=primary_header)
    flux_header['BACKGROUND_SUBTRACTED'] = 'True'
    err_header['BACKGROUND_SUBTRACTED'] = 'True'

    if smooth:
        flux_header['SMOOTHED'] = f'Smoothed by {args.kernel_size} pixel kernel'
        err_header['SMOOTHED'] = f'Smoothed by {args.kernel_size} pixel kernel'
       
    flux_hdu = fits.ImageHDU(data=flux_cube_bkg_sub, name='SCI', header=flux_header)
    err_hdu = fits.ImageHDU(data=err_cube_bkg_sub, name='ERR', header=err_header)

    hdul = fits.HDUList([primary_hdu, flux_hdu, err_hdu])
    hdul.writeto(filename_bkg_sub, overwrite=True)

    print(f'Saved metallicity maps in {filename_bkg_sub}')

# --------------------------------------------------------------------------------------------------------------------
def read_datacube(filename, args, return_components=False):
    '''
    Read the 3D IFU datacube from provided location (filename)
    Returns the 3D unumpy array (with values and corresponding uncertainties)
    '''
    print(f'\nReading in {filename}..')
    data = fits.open(filename)
    try: sci_ext = data[0].header['FLUXEXT']
    except: sci_ext = 'SCI'
    try: err_ext = data[0].header['ERREXT']
    except: err_ext = 'ERR'

    val = data[sci_ext].data
    err = data[err_ext].data
    val_header = data[sci_ext].header

    # -----------getting the wavelength array----------
    if 'wave' in val_header[f'CTYPE{args.wave_axis}'].lower():
        wave_ref_pix = val_header[f'CRPIX{args.wave_axis}']
        wave_ref_val = val_header[f'CRVAL{args.wave_axis}'] # assuming wave is in microns
        wave_delta = val_header[f'CDELT{args.wave_axis}']
        nwave = val_header[f'NAXIS{args.wave_axis}']
        wave_arr = (np.arange(nwave) - (wave_ref_pix - 1)) * wave_delta + wave_ref_val
    else:
        raise TypeError('Axis 3 is not wavlength. Therefore cannot obtain wavelength array. Try a different value(1 or 2) with --wave_axis option.')

    # -------conditonally returning stuff----------
    if return_components:
        primary_header = data[0].header
        err_header = data[err_ext].header
        return primary_header, val, err, val_header, err_header, wave_arr
    else:
        # ---------converting flux into decent units----------
        args.pix_size_arcsec = val_header['CDELT2'] * 3600 # deg to arcsec
        args.pix_scale_kpc = args.pix_size_arcsec / cosmo.arcsec_per_kpc_proper(args.z).value # kpc
        args.distance_kpc = cosmo.luminosity_distance(args.z).to('kpc').value
        steradian = (args.pix_scale_kpc / args.distance_kpc) ** 2 # to convert pixel size in degree to steradian
        args.flux_conversion_factor = val_header['PHOTMJSR'] * 1e6 * steradian * 2.99792458E-05 / (wave_arr * 1e-4) ** 2 # to convert to MJy/sr and then to Jy and then to ergs/s/cm^2/A
        val = (val.T * args.flux_conversion_factor).T # ergs/s/cm^2/A
        err = (err.T * args.flux_conversion_factor).T # ergs/s/cm^2/A

        mask = np.isnan(val) | np.isnan(err) | (val < 0) | (err <= 0)
        err[mask] = 999. # dummy high value
        flux = np.ma.masked_where(mask, val) # ergs/s/cm^2/A    
        err = np.ma.masked_where(mask, err) # ergs/s/cm^2/A    
        return flux, err, wave_arr

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # --------setting up global variables----------
    bkg_region = [[10, 14], [16, 20]]
    redshift = 3.47

    # --------setting up filenames and directory structure----------
    ifu_cube_filename = Path(args.input_file)
    args.output_dir = ifu_cube_filename.parent
    if 'bkg' in ifu_cube_filename.stem:
        ifu_cube_filename_bkg_sub = ifu_cube_filename
        do_bkg = False
    else:
        ifu_cube_filename_bkg_sub = Path(str(ifu_cube_filename).replace('.fits', '_bkg.fits'))
        do_bkg = True

    args.output_fits_dir = args.output_dir / 'fits'
    args.output_fig_dir = args.output_dir / 'plots'
    args.output_fits_dir.mkdir(parents=True, exist_ok=True)
    args.output_fig_dir.mkdir(parents=True, exist_ok=True)

    line_cube_filename = args.output_fits_dir / 'fitted_lines.fits'
    args.z = redshift
    args.distance = cosmo.luminosity_distance(args.z)

    args.arcsec_limit = 1.5 # half of NIRSpec FoV
    offset = 0
    args.extent = (-args.arcsec_limit - offset, args.arcsec_limit - offset, -args.arcsec_limit - offset, args.arcsec_limit - offset)

    # -------------------handling data cube------------------------------

    # ------reading in background subtracted input IFU data cube---------------
    if not os.path.exists(ifu_cube_filename_bkg_sub):# or args.clobber:
        background_subtract(ifu_cube_filename, ifu_cube_filename_bkg_sub, args, bkg_region=bkg_region, smooth=True)
    ifu_cube, ifu_err_cube, wave_arr = read_datacube(ifu_cube_filename_bkg_sub, args)
    rest_wave_arr = wave_arr * 1e4 / (1 + redshift) # microns to Angstroms

    # -------reading the linelist-------------------
    df_linelist = get_linelist()
    df_linelist = df_linelist[df_linelist['restwave'].between(np.min(rest_wave_arr), np.max(rest_wave_arr))].reset_index(drop=True)

    if not os.path.exists(line_cube_filename) or args.clobber:
        # ------fitting IFU data cube------------------
        flux_cube, vel_cube, vdisp_cube = linefit_cube(ifu_cube, ifu_err_cube, rest_wave_arr, df_linelist, contmask_lambda_window=args.contmask_lambda, args=args)

        # -------saving the fitted line cube----------
        save_fitted_cube(flux_cube, vel_cube, vdisp_cube, df_linelist, line_cube_filename, args)
       
    # ------reading in the fitted line cube----------
    print(f'Reading existing fitted cube from {line_cube_filename}')
    line_cube = fits.open(line_cube_filename)
    df_linelist[df_linelist['LineID'].isin(line_cube[0].header['LINES'].split(','))]
    
    # -------plotting emission lines------------------
    nrow, ncol = 2, 3
    fig, axes = plt.subplots(nrow, ncol, figsize=(8, 6), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.07, right=0.95, top=0.95, bottom=0.07, wspace=0.1, hspace=0.1)
    
    args.EB_V = 0. #get_EB_V(line_cube)
    line_label_arr = ['OIII5007', 'HBeta', 'NII6584', 'H6562', 'SII6717', 'SII6730']
    vmin, vmax = -1, 3
    line_map_arr = get_emission_line_maps(line_cube, line_label_arr, EB_V=args.EB_V, silent=args.only_integrated)
    for index, line_map in enumerate(line_map_arr):
        i, j = int(index / ncol), index % ncol
        ax = axes[i][j]
        ax, _ = plot_2D_map(unp.nominal_values(line_map), ax, args, takelog=True, label=line_label_arr[index], cmap='summer', vmin=vmin, vmax=vmax, hide_xaxis=i != nrow - 1, hide_yaxis=j > 0, hide_cbar=False, image_err=unp.std_devs(line_map), fontsize=args.fontsize)

    outfigname = args.output_fig_dir / 'line_maps.png'
    save_fig(fig, outfigname, args)
    
    # ------computing properties-------------
    sfr_map = get_SFR(line_cube, args.distance, EB_V=args.EB_V)
    logOH_map = get_Z_C19(line_cube, args)
    vel_map = get_emission_line_map('H6562', line_cube, EB_V=0, ext='VEL')
    vdisp_map = get_emission_line_map('H6562', line_cube, EB_V=0, ext='VDISP')

    # -------plotting properties------------------
    fig, axes = plt.subplots(1, 3, figsize=(10, 4), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.07, right=0.95, top=0.95, bottom=0.07, wspace=0.1, hspace=0.1)
    
    #axes[0][0], _ = plot_2D_map(args.EB_V, axes[0][0], args)
    #axes[0], _ = plot_2D_map(logOH_map, axes[0], args, takelog=False, label='log O/H + 12', cmap='winter', vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=False, image_err=None, fontsize=args.fontsize)
    axes[0], _ = plot_2D_map(vel_map, axes[0], args, takelog=False, label='Vel (km/s)', cmap='coolwarm', vmin=-500, vmax=500, hide_xaxis=False, hide_yaxis=False, hide_cbar=False, image_err=None, fontsize=args.fontsize)
    axes[1], _ = plot_2D_map(vdisp_map, axes[1], args, takelog=False, label='Vel disp (km/s)', cmap='plasma', vmin=0, vmax=5000, hide_xaxis=False, hide_yaxis=True, hide_cbar=False, image_err=None, fontsize=args.fontsize)
    axes[2], _ = plot_2D_map(sfr_map.data, axes[2], args, takelog=True, label='SFR', cmap='magma', vmin=18, vmax=21, hide_xaxis=False, hide_yaxis=True, hide_cbar=False, image_err=None, fontsize=args.fontsize)

    outfigname = args.output_fig_dir / 'properties.png'
    save_fig(fig, outfigname, args)
    
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
