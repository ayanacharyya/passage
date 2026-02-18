'''
    Filename: read_line_maps.py
    Notes: Reads in grizli line flux maps from any given filename, and computes various quantitites
    Author : Ayan
    Created: 18-02-26
    Example: run read_line_maps.py --Zdiag R23 --use_C25 --Zbranch low
             run read_line_maps.py --Zdiag R23 --use_C25
'''

from header import *
from util import *
from make_diagnostic_maps import trim_image, take_safe_log_ratio, take_safe_log_sum, get_dereddened_flux
from plots_for_zgrad_paper import odr_fit, plot_fitted_line

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_emission_line_map(line, full_hdu, args, dered=True, silent=True):
    '''
    Retrieve the emission map for a given line from the HDU
    Returns the 2D line map and the corresponding 2D uncertainty map
    '''
    # ---------getting the spatially resolved line flux map----------------
    try:
        line_hdu = full_hdu['LINE', line]
        line_wht = full_hdu['LINEWHT', line]
    except KeyError:
        if line == 'OIII':
            line_hdu = full_hdu['LINE', 'OIII-5007']
            line_wht = full_hdu['LINEWHT', 'OIII-5007']

    line_map = line_hdu.data * 1e-17 # in units of ergs/s/cm^2
    line_map_err = 1e-17 / (line_wht.data ** 0.5)  # ERR = 1/sqrt(LINEWHT) = flux uncertainty; in units of ergs/s/cm^2
    line_wave = line_hdu.header['RESTWAVE'] # in Angstrom

    # ----------deblending flux--------------------
    factor = 1.0
    if not args.do_not_correct_flux:
        if line == 'OIII': # special treatment for OIII 5007 line, in order to account for and remove the OIII 4959 component
            ratio_5007_to_4959 = 2.98 # from grizli source code
            factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
            if not silent: print(f'Correcting OIII for 4959 component, by factor of {factor:.3f}')
        elif line == 'Ha': # special treatment for Ha line, in order to account for and remove the NII component
            factor = 0.823 # from James et al. 2023?
            if not silent: print(f'Correcting Ha for NII component, by factor of {factor:.3f}')

    line_map = line_map * factor
    line_map_err = line_map_err * factor

    # ----------getting a smaller cutout around the object center-----------
    line_map = trim_image(line_map, args, skip_re_trim=True)
    line_map_err = trim_image(line_map_err, args, skip_re_trim=True)

    # ----------re-center line map-----------
    line_map = ndimage.shift(line_map, [args.ndelta_xpix, args.ndelta_ypix], order=0, cval=np.nan)
    line_map_err = ndimage.shift(line_map_err, [args.ndelta_xpix, args.ndelta_ypix], order=0, cval=np.nan)

    # -------masking outside segmentation map---------
    line_map = np.ma.masked_where(args.segmentation_map != args.id, line_map)
    line_map_err = np.ma.masked_where(args.segmentation_map != args.id, line_map_err)

    # -----------getting the dereddened flux value-----------------
    if dered and args.EB_V != 0:
        line_map_quant = get_dereddened_flux(unp.uarray(line_map, line_map_err), line_wave, args.EB_V)
        line_map = unp.nominal_values(line_map_quant)
        line_map_err = unp.std_devs(line_map_quant)

    # ---------converting from flux to surface brightness units------------
    pixscale_kpc = args.pix_size_arcsec/ cosmo.arcsec_per_kpc_proper(args.z).value # kpc
    line_map /= (pixscale_kpc ** 2) # this is now in ergs/s/cm^2/kpc^2 (i.e. Surface Brightness)
    line_map_err /= (pixscale_kpc ** 2) # this is now in ergs/s/cm^2/kpc^2 (i.e. Surface Brightness)

    # -----------------------------------------------------
    line_map = unp.uarray(line_map, line_map_err)
    line_int = ufloat(np.nan, np.nan)
    line_sum = ufloat(np.nan, np.nan)
    line_ew = ufloat(np.nan, np.nan)

    return line_map, line_wave, line_int, line_sum, line_ew

# --------------------------------------------------------------------------------------------------------------------
def compute_Z_NB(line_label_array, line_waves_array, line_flux_array):
    '''
    Calculates and returns the NebulaBayes metallicity given a list of observed line fluxes
    '''
    line_flux_array = [np.atleast_1d(item) for item in line_flux_array]
    npixels = len(line_flux_array[0])
    if args.vorbin and npixels > 1: # getting how many unique IDs present, so that NB does not have to run unnecessary repeats
        IDs_array = args.voronoi_bin_IDs.astype(int).flatten()
    else:
        IDs_array = np.arange(npixels).flatten()
    unique_IDs_array = np.unique(IDs_array)
    print(f'\nAbout to start running NB, with {len(line_label_array)} lines: {line_label_array}..\n')
    if len(unique_IDs_array) > 60: print(f'This might take ~{int(len(unique_IDs_array) / 60)} min')

    # -----making a "net" mask array and separating out the line fluxes form the input unumpy arrays---------
    net_mask = np.zeros(np.shape(line_flux_array[0]), dtype=bool)
    obs_flux_array, obs_err_array = [], []

    for index in range(len(line_flux_array)):
        if np.ma.isMaskedArray(line_flux_array[index]):
            net_mask = net_mask | line_flux_array[index].mask
            obs_flux_array.append(unp.nominal_values(line_flux_array[index].data).flatten())
            obs_err_array.append(unp.std_devs(line_flux_array[index].data).flatten())
        else:
            obs_flux_array.append(unp.nominal_values(line_flux_array[index]).flatten())
            obs_err_array.append(unp.std_devs(line_flux_array[index]).flatten())

    obs_flux_array = np.array(obs_flux_array)
    obs_err_array = np.array(obs_err_array)
    net_mask_array = net_mask.flatten()

    # -----loading the NB HII region model grid---------
    NB_Model_HII = NB_Model("HII", line_list=line_label_array)

    out_dir = args.output_dir / args.field / f'{args.id:05d}_NB_results'
    if args.vorbin: out_dir = out_dir / f'{vorbin_text}'
    out_subdirs = [out_dir / 'prior_plots', out_dir / 'likelihood_plots', out_dir / 'posterior_plots', out_dir / 'best_model_catalogs', out_dir / 'param_estimates_catalogs']
    for this_out_subdir in out_subdirs: this_out_subdir.mkdir(exist_ok=True, parents=True)

    # -----looping over each pixel to calculate NB metallicity--------
    logOH_array = []
    logOH_dict_unique_IDs = {}
    counter = 0
    start_time3 = datetime.now()

    for index in range(len(obs_flux_array[0])):
        this_ID = IDs_array[index]
        if net_mask_array[index]: # no need to calculate for those pixels that are already masked
            #print(f'Skipping NB for masked pixel {index + 1} out of {len(obs_flux_array[0])}..')
            logOH = ufloat(np.nan, np.nan)
        elif this_ID in logOH_dict_unique_IDs.keys():
            #print(f'Skipping NB due to existing measurement from unique ID {this_ID} for pixel {index + 1} out of {len(obs_flux_array[0])}..')
            logOH = logOH_dict_unique_IDs[this_ID]
        else:
            start_time4 = datetime.now()

            # ------getting all line fluxes-------------
            obs_fluxes = obs_flux_array[:,index]
            obs_errs = obs_err_array[:, index]

            # ------discarding lines with negative fluxes-------------
            good_obs = obs_fluxes >= 0
            obs_fluxes = obs_fluxes[good_obs]
            obs_errs = obs_errs[good_obs]
            line_labels = list(np.array(line_label_array)[good_obs])
            line_waves = np.array(line_waves_array)[good_obs]

            if len(line_labels) > 1:
                # -------setting up NB parameters----------
                dered = 'Hbeta' in line_labels and 'Halpha' in line_labels and args.dered_in_NB
                norm_line = 'Hbeta' if 'Hbeta' in line_labels else 'OIII5007' if 'OIII5007' in line_labels else 'NII6583_Halpha' if 'NII6583_Halpha' in line_labels else 'Halpha' if 'Halpha' in line_labels else line_labels[0]
                kwargs = {'prior_plot': os.path.join(out_dir, 'prior_plots', f'{this_ID}_HII_prior_plot.pdf'),
                        'likelihood_plot': os.path.join(out_dir, 'likelihood_plots', f'{this_ID}_HII_likelihood_plot.pdf'),
                        'posterior_plot': os.path.join( out_dir, 'posterior_plots', f'{this_ID}_HII_posterior_plot.pdf'),
                        'estimate_table': os.path.join(out_dir, 'best_model_catalogs', f'{this_ID}_HII_param_estimates.csv'),
                        'best_model_table': os.path.join(out_dir, 'param_estimates_catalogs', f'{this_ID}_HII_best_model.csv'),
                        'verbosity': 'ERROR',
                        'norm_line':norm_line,
                        'deredden': dered,
                        'propagate_dered_errors': dered,
                        'obs_wavelengths': line_waves if dered else None
                        }

                # -------running NB--------------
                print(f'Deb1576: binID {this_ID}: nlines={len(obs_fluxes)}, {dict(zip(line_labels, obs_fluxes))}, norm_line = {norm_line}, dereddening on the fly? {dered}') ##
                Result = NB_Model_HII(obs_fluxes, obs_errs, line_labels, **kwargs)

                # -------estimating the resulting logOH, and associated uncertainty-----------
                df_estimates = Result.Posterior.DF_estimates # pandas DataFrame
                logOH_est = df_estimates.loc['12 + log O/H', 'Estimate']
                logOH_low = df_estimates.loc['12 + log O/H', 'CI68_low']
                logOH_high = df_estimates.loc['12 + log O/H', 'CI68_high']
                logOH_err = np.mean([logOH_est - logOH_low, logOH_high - logOH_est])
                logOH = ufloat(logOH_est, logOH_err)
                
                print(f'Ran NB for unique ID {this_ID} ({counter + 1} out of {len(unique_IDs_array)}) with {len(obs_fluxes)} good fluxes in {timedelta(seconds=(datetime.now() - start_time4).seconds)}')

                best_model_dict = Result.Posterior.best_model
                if 'Halpha' in line_labels and 'NII6583' in line_labels:
                    N2Ha_model = best_model_dict['table'].loc['NII6583', 'Model'] / best_model_dict['table'].loc['Halpha', 'Model']
                    print(f'N2Ha = {N2Ha_model}')
            
            else:
                logOH = ufloat(np.nan, np.nan)
                print(f'Could not run NB for unique ID {this_ID} ({counter + 1} out of {len(unique_IDs_array)}) with only {len(obs_fluxes)} good fluxes')

            counter += 1
            logOH_dict_unique_IDs.update({this_ID: logOH}) # updating to unique ID dictionary once logOH has been calculated for this unique ID

            if args.only_integrated: ##
                print('\n') ##
                norm_flux = obs_fluxes[np.where(np.array(line_labels)==kwargs['norm_line'])[0][0]]
                for index, line in enumerate(line_labels): print(f'{args.id}: {line}={obs_fluxes[index]/norm_flux: .1f}+/-{obs_errs[index]/norm_flux: .1f}')  ##
                print(f'{args.id}: log O/H + 12={logOH: .1f}')

        logOH_array.append(logOH)
    print(f'\nRan NB for total {counter} unique pixels out of {len(obs_flux_array[0])}, in {timedelta(seconds=(datetime.now() - start_time3).seconds)}\n')

    # ---------collating all the metallicities computed---------
    log_OH = np.ma.masked_where(net_mask, np.reshape(logOH_array, np.shape(line_flux_array[0])))
    if len(log_OH) == 1: log_OH = log_OH.data[0]

    return log_OH

# --------------------------------------------------------------------------------------------------------------------
def get_Z_NB(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated metallicity from a given HDU, based on Bayesian
    statistics, using NebulaBayes
    '''
    # -----dict for converting line label names to those acceptable to NB---------
    if args.use_original_NB_grid: line_label_dict = {'OII':'OII3726_29', 'Hb':'Hbeta', 'OIII':'OIII5007', 'OIII-4363':'OIII4363', 'OI-6302':'OI6300', \
                       'Ha':'Halpha', 'SII':'SII6716', 'NeIII-3867':'NeIII3869'}
    else: line_label_dict = {'OII':'OII3726_29', 'Hb':'Hbeta', 'OIII':'OIII5007', 'OIII-4363':'OIII4363', 'OI-6302':'OI6300', \
                       'Ha':'Halpha' if args.dered_in_NB else 'NII6583_Halpha', 'SII':'SII6716_31', 'NeIII-3867':'NeIII3869', 'OII-7325':'OII7320', \
                       'HeI-5877':'HeI5876', 'Hg':'Hgamma', 'Hd':'Hdelta', 'NII':'NII6583'}

    line_map_array, line_int_array, line_sum_array, line_label_array, line_waves_array = [], [], [], [], []
    for line in args.available_lines:
        line_map, line_wave, line_int, line_sum, _ = get_emission_line_map(line, full_hdu, args, dered=not args.dered_in_NB, silent=False)
        if not args.do_not_correct_flux:
            factor = 1.
            if args.use_original_NB_grid:
                if line == 'SII':
                    factor = 2.  # from grizli source code
                    print(f'Correcting SII to exclude the other SII component, i.e. dividing by factor {factor:.3f}')
            else:
                if line == 'OIII':
                    ratio_5007_to_4959 = 2.98  # from grizli source code
                    factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
                    print(f'Un-correcting OIII to include the 4959 component back, i.e. dividing by factor {factor:.3f} because in NB, OIII = 5007+4959')
                if line == 'Ha' and not args.dered_in_NB:
                    factor = 0.823  # from grizli source code
                    print(f'Un-correcting Ha to include the NII component back, i.e. dividing by factor {factor:.3f} because NB can use line summation')
            if line_map is not None: line_map = np.ma.masked_where(line_map.mask, line_map.data / factor)
            line_int = line_int / factor
            line_sum = line_sum / factor

        #line_snr = line_int.n / line_int.s
        line_snr = line_sum.n / line_sum.s

        if args.only_integrated:
            print(f'Deb1485: {line}: sum={line_sum:.2e}, sum-snr={line_sum.n / line_sum.s:.1f}, int={line_int:.2e}, int-snr={line_int.n / line_int.s:.1f}')

        #if line_snr > 2 and line in line_label_dict.keys() and line not in args.exclude_lines:
        if line in line_label_dict.keys() and line not in args.exclude_lines:
            line_map_array.append(line_map)
            line_int_array.append(line_int)
            line_sum_array.append(line_sum)
            line_label_array.append(line_label_dict[line])
            line_waves_array.append(rest_wave_dict[line] * 10) # factor of 10 to convert from nm to Angstroms

    # -------also adding NII flux to list if Halpha and NII are being used separately by NB----------
    add_nii = False
    if add_nii and 'Ha' in args.available_lines and args.dered_in_NB:
        print(f'Providing the NII components separately to NB, since Halpha was already corrected for NII')
        line = 'Ha'
        line_map, line_wave, line_int, line_sum, _ = get_emission_line_map(line, full_hdu, args, dered=not args.dered_in_NB, silent=False)
        factor = 0.823  # from grizli source code
        line_map = np.ma.masked_where(line_map.mask, line_map.data * (1 - factor) / factor)
        line_int = line_int  * (1 - factor) / factor
        line_sum = line_sum  * (1 - factor) / factor

        line_map_array.append(line_map)
        line_int_array.append(line_int)
        line_sum_array.append(line_sum)
        line_label_array.append('NII6583')
        line_waves_array.append(6583.454) # factor of 10 to convert from nm to Angstroms

    # ----------calling NB----------------------
    logOH_map = compute_Z_NB(line_label_array, line_waves_array, line_map_array) if not args.only_integrated else None
    # if (np.array(line_int_array) > 0).all(): logOH_int = compute_Z_NB(line_label_array, line_waves_array, line_int_array)
    # else: logOH_int = ufloat(np.nan, np.nan)
    # if (np.array(line_int_array) > 0).all(): logOH_sum = compute_Z_NB(line_label_array, line_waves_array, line_sum_array)
    # else: logOH_sum = ufloat(np.nan, np.nan)
    logOH_int = compute_Z_NB(line_label_array, line_waves_array, line_int_array)
    logOH_sum = compute_Z_NB(line_label_array, line_waves_array, line_sum_array)

    return logOH_map, logOH_int, logOH_sum, line_label_array

# --------------------------------------------------------------------------------------------------------------------
def compute_Z_C19(ratio, coeff, ax=None, branch='high', use_C25=False, silent=False):
    '''
    Calculates and returns the metallicity given observed line fluxes ratio and coefficient, according to Curti+2019
    '''
    # -----handling turnover situations, where measured ratio is beyond model peak ratio---------
    metallicity_offset = 0 if use_C25 else 8.69
    reasonable_Z_limit = [7.0, 8.4] if use_C25 else [7.6, 8.9] # Z limits within which each calibration is valid
    #reasonable_Z_limit = [6.5, 9.1]
    model = np.poly1d(coeff)
    model_diff = np.polyder(model, m=1)
    model_turnovers = np.roots(model_diff) + metallicity_offset
    possible_model_turnovers = model_turnovers[(model_turnovers > reasonable_Z_limit[0]) & (model_turnovers < reasonable_Z_limit[1])]

    if len(possible_model_turnovers) > 0:
        logOH_turnover = np.max(possible_model_turnovers) # based on solving the differential of polynomial with above coefficients, this is the value of Z where the relation peaks
        ratio_turnover = model(logOH_turnover - metallicity_offset)
    else:
        logOH_turnover, ratio_turnover = np.nan, np.nan

    if ax is not None and np.isfinite(logOH_turnover):
        for index, ax1 in enumerate(ax):
            ax1.axvline(logOH_turnover, ls='--', c='k', lw=1, label='Turnover location' if index == 0 else None)
            ax1.axhline(ratio_turnover, ls='--', c='k', lw=1)
            #ax1.fill_betweenx([-5, 5], reasonable_Z_limit[0], reasonable_Z_limit[1], color='cyan', alpha=0.1, lw=0)

    # ---------getting the line ratio limits of this diagnostics-----------
    if not silent:
        print(f'\nRatio at turnover = {ratio_turnover}')
        print(f'Ratio at min valid metallicity = {model(reasonable_Z_limit[0] - metallicity_offset)}')
        print(f'Ratio at max valid metallicity = {model(reasonable_Z_limit[1] - metallicity_offset)}')

     # --------determining data and masks------------
    if np.ma.isMaskedArray(ratio):
        ratio_arr = np.atleast_1d(ratio.data).flatten()
        mask = ratio.mask
    else:
        ratio_arr = np.atleast_1d(ratio).flatten()
        mask = None

    # --------computing the metallicitities------------
    log_OH = []
    no_solution_not_labelled = True

    for index2, this_ratio in enumerate(ratio_arr):
        this_ratio = unp.nominal_values(this_ratio)
        if np.isfinite(ratio_turnover) and this_ratio > ratio_turnover:
            log_OH.append(ufloat(logOH_turnover, 0.))
            if ax is not None:
                ax[0].axhline(this_ratio, lw=0.5, ls='solid', c='sienna', alpha=0.5, label='No solution' if no_solution_not_labelled else None)
                no_solution_not_labelled = False
        else:
            try:
                poly_to_solve = np.hstack([coeff[:-1], [coeff[-1] - this_ratio]])
                roots = np.roots(poly_to_solve)
                real_roots = np.sort(roots[np.isreal(roots)]) + metallicity_offset # see Table 1 caption in Curti+19
                possible_roots = real_roots[(real_roots > reasonable_Z_limit[0]) & (real_roots < reasonable_Z_limit[1])]
                
                impossible_roots = list(set(real_roots) - set(possible_roots))
                if ax is not None:
                    for index, real_root in enumerate(impossible_roots): ax[0].scatter(real_root, this_ratio, lw=0, s=10, c='grey')

                if branch == 'high': # THIS IS WHERE MAKING THE CHOICE TO GO WITH THE HIGHER METALLICITY BRANCH, WHEREVER THE CHOICE NEEDS TO BE MADE
                    this_log_OH = np.max(possible_roots)
                elif branch == 'low':
                    this_log_OH = np.min(possible_roots)
                log_OH.append(ufloat(this_log_OH, 0.))
            except:
                this_log_OH = np.nan
                log_OH.append(ufloat(this_log_OH, 0))
                if ax is not None:
                    ax[0].axhline(this_ratio, ls='solid', lw=0.5, c='sienna', alpha=0.5, label='No solution' if no_solution_not_labelled else None)
                    no_solution_not_labelled = False

            if np.isfinite(this_log_OH) and ax is not None:
                col_arr = ['cornflowerblue', 'limegreen', 'k']
                for index, real_root in enumerate(possible_roots): ax[0].scatter(real_root, this_ratio, ec='k', lw=0.5, s=50, c=col_arr[index] if len(real_roots) > 1 else 'salmon', zorder=100)
                try:
                    this_model_turnover = model_turnovers[-2]
                    ratio_turnover2 = max(model(this_model_turnover - metallicity_offset), model(reasonable_Z_limit[0] - metallicity_offset))
                    ax[1].scatter(this_log_OH, this_ratio, ec='k', lw=0.5, s=50, c='limegreen' if branch == 'low' and this_ratio > ratio_turnover2 else 'cornfloweblue' if branch == 'high' and this_ratio > ratio_turnover2 else 'salmon', zorder=100)
                except:
                    pass
    log_OH = np.reshape(log_OH, np.shape(ratio))
    if mask is not None: log_OH = np.ma.masked_where(mask, log_OH)
    if ax is not None:
        handles, labels = ax[0].get_legend_handles_labels()
        fig = ax[0].figure
        fig.legend(handles, labels, bbox_to_anchor=(0.9, 0.1), loc='lower right', fontsize=args.fontsize)
        #ax[0].legend(loc='lower left')

    return log_OH

# --------------------------------------------------------------------------------------------------------------------
def get_emission_line_maps(full_hdu, line_labels, args, silent=False):
    '''
    Computes and returns the spatially resolved as well as intregrated metallicity based on a given Curti+2019 calibration, from a given HDU
    '''
    if not all([line in args.available_lines for line in line_labels]):
        print(f'All lines in {line_labels} not available in object {args.id}')
        return None, None

    line_map_arr, line_int_arr, line_sum_arr = [], [], []
    for line in line_labels:
        line_map, line_wave, line_int, line_sum, _ = get_emission_line_map(line, full_hdu, args, silent=silent)

        if not args.do_not_correct_flux and args.Zdiag == 'R23':
            factor = 1.
            if line == 'OIII': # special treatment for OIII 5007 line, in order to account for and ADD the OIII 4959 component back
                ratio_5007_to_4959 = 2.98  # from grizli source code
                factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
                print(f'Un-correcting OIII to include the 4959 component, for computing R23 metallicity, by factor of {factor:.3f}')

            if line_map is not None: line_map = np.ma.masked_where(line_map.mask, line_map.data / factor)
            line_int = line_int / factor
            line_sum = line_sum / factor

        line_map_arr.append(line_map)
        line_int_arr.append(line_int)
        line_sum_arr.append(line_sum)


    return line_map_arr, line_int_arr, line_sum_arr

# --------------------------------------------------------------------------------------------------------------------
def get_Z_C19(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated metallicity based on a given Curti+2019 and Cataldi+2025 calibrations, from a given HDU
    '''
    # ------getting appropriate emission lines and calibration coefficients--------------
    if args.Zdiag == 'O3O2':
        line_map_arr, line_int_arr, line_sum_arr = get_emission_line_maps(full_hdu, ['OIII', 'OII'], args, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        ratio_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1])
        try: ratio_int = unp.log10(line_int_arr[0] / line_int_arr[1])
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        try: ratio_sum = unp.log10(line_sum_arr[0] / line_sum_arr[1])
        except ValueError: ratio_sum = ufloat(np.nan, np.nan)
        if args.use_C25: coeff = [0.01124, 0.03072, 0.1251, -0.01470] # from Table 4 of Cataldi+25
        else: coeff = [-0.691, -2.944, -1.308]  # c0-2 parameters from Table 2 of Curti+19 3rd row (O3O2)

    elif args.Zdiag == 'R3':
        line_map_arr, line_int_arr, line_sum_arr = get_emission_line_maps(full_hdu, ['OIII', 'Hb'], args, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        ratio_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1])
        try: ratio_int = unp.log10(line_int_arr[0] / line_int_arr[1])
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        try: ratio_sum = unp.log10(line_sum_arr[0] / line_sum_arr[1])
        except ValueError: ratio_sum = ufloat(np.nan, np.nan)
        if args.use_C25: coeff = [-3.587, -4.475, 1.345, -0.08951] # from Table 4 of Cataldi+25
        else: coeff = [-0.277, -3.549, -3.593, -0.981]  # c0-3 parameters from Table 2 of Curti+19 2nd row (R3)

    elif args.Zdiag == 'R23':
        line_map_arr, line_int_arr, line_sum_arr = get_emission_line_maps(full_hdu, ['OIII', 'OII', 'Hb'], args, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        R1_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[2], skip_log=True)
        R2_map = take_safe_log_ratio(line_map_arr[1], line_map_arr[2], skip_log=True)
        ratio_map = take_safe_log_sum(R1_map, R2_map)
        try: ratio_int = unp.log10((line_int_arr[0] + line_int_arr[1]) / line_int_arr[2])
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        try: ratio_sum = unp.log10((line_sum_arr[0] + line_sum_arr[1]) / line_sum_arr[2])
        except ValueError: ratio_sum = ufloat(np.nan, np.nan)
        if args.use_C25: coeff = [-2.555, -3.192, 0.9630, -0.06351] # from Table 4 of Cataldi+25
        else: coeff = [0.527, -1.569, -1.652, -0.421]  # c0-3 parameters from Table 2 of Curti+19 4th row (R23)

    elif args.Zdiag == 'O3S2':
        line_map_arr, line_int_arr, line_sum_arr = get_emission_line_maps(full_hdu, ['OIII', 'Hb', 'SII', 'Ha'], args, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        R1_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1], skip_log=True)
        R2_map = take_safe_log_ratio(line_map_arr[2], line_map_arr[3], skip_log=True)
        ratio_map = take_safe_log_ratio(R1_map, R2_map)
        try: ratio_int = unp.log10((line_int_arr[0] / line_int_arr[1]) / (line_int_arr[2] / line_int_arr[3]))
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        try: ratio_sum = unp.log10((line_sum_arr[0] / line_sum_arr[1]) / (line_sum_arr[2] / line_sum_arr[3]))
        except ValueError: ratio_sum = ufloat(np.nan, np.nan)
        coeff = [0.191, -4.292, -2.538, 0.053, 0.332]  # c0-4 parameters from Table 2 of Curti+19 last row (O3S2)

    elif args.Zdiag == 'S2':
        line_map_arr, line_int_arr, line_sum_arr = get_emission_line_maps(full_hdu, ['SII', 'Ha'], args, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        ratio_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1])
        try: ratio_int = unp.log10(line_int_arr[0] / line_int_arr[1])
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        try: ratio_sum = unp.log10(line_sum_arr[0] / line_sum_arr[1])
        except ValueError: ratio_sum = ufloat(np.nan, np.nan)
        coeff = [-0.442, -0.360, -6.271, -8.339, -3.559]  # c0-3 parameters from Table 2 of Curti+19 3rd-to-last row (S2)

    elif args.Zdiag == 'RS32':
        line_map_arr, line_int_arr, line_sum_arr = get_emission_line_maps(full_hdu, ['OIII', 'Hb', 'SII', 'Ha'], args, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        R1_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1], skip_log=True)
        R2_map = take_safe_log_ratio(line_map_arr[2], line_map_arr[3], skip_log=True)
        ratio_map = take_safe_log_sum(R1_map, R2_map)
        try: ratio_int = unp.log10((line_int_arr[0] / line_int_arr[1]) + (line_int_arr[2] / line_int_arr[3]))
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        try: ratio_sum = unp.log10((line_sum_arr[0] / line_sum_arr[1]) + (line_sum_arr[2] / line_sum_arr[3]))
        except ValueError: ratio_sum = ufloat(np.nan, np.nan)
        coeff = [-0.054, -2.546, -1.970, 0.082, 0.222]  # c0-3 parameters from Table 2 of Curti+19 2nd-to-last row (RS32)

    elif args.Zdiag == 'R2':
        line_map_arr, line_int_arr, line_sum_arr = get_emission_line_maps(full_hdu, ['OII', 'Hb'], args, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan)
        ratio_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1])
        try: ratio_int = unp.log10(line_int_arr[0] / line_int_arr[1])
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        try: ratio_sum = unp.log10(line_sum_arr[0] / line_sum_arr[1])
        except ValueError: ratio_sum = ufloat(np.nan, np.nan)
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
    logOH_map = compute_Z_C19(ratio_map, coeff, ax=ax, branch=args.Zbranch, use_C25=args.use_C25) if not args.only_integrated else None
    logOH_int = compute_Z_C19(ratio_int, coeff, branch=args.Zbranch, use_C25=args.use_C25)
    logOH_int = ufloat(unp.nominal_values(np.atleast_1d(logOH_int))[0], unp.std_devs(np.atleast_1d(logOH_int))[0])
    logOH_sum = compute_Z_C19(ratio_sum, coeff, branch=args.Zbranch, use_C25=args.use_C25)
    logOH_sum = ufloat(unp.nominal_values(np.atleast_1d(logOH_sum))[0], unp.std_devs(np.atleast_1d(logOH_sum))[0])

    # -------saving the debugging plots---------------
    if ax is not None:
        Zbranch_text = '' if args.Zdiag in ['NB', 'P25', 'Te'] else f'-{args.Zbranch}'
        figname = fig_dir / f'{args.field}_{args.id:05d}_metallicity_debug{snr_text}{vorbin_text}_Zdiag_{args.Zdiag}{Zbranch_text}.png'
        fig.savefig(figname, transparent=args.fortalk, dpi=200)
        print(f'\nSaved figure at {figname}')

    return logOH_map, logOH_int, logOH_sum

# --------------------------------------------------------------------------------------------------------------------
def get_Z(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated metallicity from a given HDU
    Also saves the metallicity map as a fits file
    '''
    # -----------determining output fits file name---------
    NB_text = '_orig_grid' if args.use_original_NB_grid and args.Zdiag == 'NB' else ''
    exclude_text = f'_without_{",".join(args.exclude_lines)}' if len(args.exclude_lines) > 0 and args.Zdiag == 'NB' else ''
    Zbranch_text = '' if args.Zdiag in ['NB', 'P25'] else f'-{args.Zbranch}'
    AGN_diag_text = f'_AGNdiag_{args.AGN_diag}'if args.AGN_diag != 'None' else ''
    extent_text = f'{args.arcsec_limit}arcsec' if args.re_limit is None else f'{args.re_limit}re'
    C25_text = '_wC25' if args.use_C25 and 'NB' not in args.Zdiag else ''
    dered_text = '_dered' if args.dered_in_NB and 'NB' in args.Zdiag else ''
    
    # --------deriving the metallicity map-------------
    if args.Zdiag == 'NB':
        logOH_map, logOH_int, logOH_sum, line_label_array = get_Z_NB(full_hdu, args)
    else:
        logOH_map, logOH_int, logOH_sum = get_Z_C19(full_hdu, args)

    return logOH_map, logOH_int, logOH_sum

# --------------------------------------------------------------------------------------------------------------------
def plot_radial_profile(df, ax, args, ylim=None, xlim=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, skip_annotate=False, quant='log_OH'):
    '''
    Plots and fits the radial profile from a given dataframe in a given axis
    Returns the axis handle and the linefit
    '''
    label_dict = smart_dict({'SFR': r'$\log$ $\Sigma_{\rm SFR}$ (M$_{\odot}$/yr/kpc$^2$)', 'logOH': r'$\log$ (O/H) + 12', 'Z': r'$\log$ (O/H) + 12', 'log_OH': r'$\log$ (O/H) + 12'})

    df = df.sort_values(by='distance').reset_index(drop=True)

    # -------fitting by various methods--------------    
    linefit_odr = odr_fit(df, quant_x='distance', quant_y=quant)

    # -----------plot data---------------
    quant_x, quant_y = 'distance', quant
    col='darkgoldenrod'
    ax.scatter(df[quant_x], df[quant_y], c=col, s=20, lw=1, edgecolor='k', alpha=1)
    if quant_y + '_u' in df: ax.errorbar(df[quant_x], df[quant_y], yerr=df[quant_y + '_u'], c=col, fmt='none', lw=0.5, alpha=0.2)

    # ----------plot fitted line-------------
    xarr = df[quant_x]
    ax = plot_fitted_line(ax, linefit_odr, xarr, col, args, quant=quant, short_label=False, index=index)

    # --------annotating axis--------------
    if xlim is not None: ax.set_xlim(xlim[0], xlim[1]) # kpc
    if ylim is not None: ax.set_ylim(ylim[0], ylim[1])
    if not skip_annotate: ax = annotate_axes(ax, 'Radius (kpc)' if args.re_limit is None else r'Radius (R$_e$)', label_dict[quant], args, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar)
    
    return ax, linefit_odr

# --------------------------------------------------------------------------------------------------------------------
def plot_2D_map(image, ax, label, args, cmap='viridis', clabel='', takelog=True, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, hide_cbar_ticks=False, cticks_integer=True):
    '''
    Plots a given 2D image in a given axis
    Returns the axis handle
    '''
    image = np.ma.masked_where(image.mask, unp.nominal_values(image.data))

    if takelog:
        new_mask = image <= 0
        image[new_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        image = np.ma.masked_where(image.mask | new_mask, np.log10(image.data))

    p = ax.imshow(image, cmap=cmap, origin='lower', extent=args.extent, vmin=vmin, vmax=vmax)
    ax.scatter(0, 0, marker='x', s=10, c='grey')
    ax.set_aspect('equal') 
    
    units = 'arcsec'
    ax = annotate_axes(ax, f'Offset ({units})', f'Offset ({units})', args, label=label, clabel=clabel, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar, p=p, hide_cbar_ticks=hide_cbar_ticks, cticks_integer=cticks_integer)

    return ax

# --------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.re_limit is None: args.re_limit = 2.
    args.fontsize = 10
    args.fontfactor = 1
                        
    line_list = ['OII', 'Hb', 'OIII-4959', 'OIII', 'Ha', 'SII']
    dered = True
    # --------------defining filenames and directories-------------------------------
    input_dir = Path('/Users/acharyya/Work/astro/passage/glass_data/resolved_sfh')
    fig_dir = input_dir / 'figs'
    files = glob.glob(str(input_dir) + '/regions_*arcsec.line.fits')

    for ind, full_filename in enumerate(files):
        print(f'\nDoing ({ind + 1}/{len(files)})..')
        switched_diag = False
        # -----------------------------------------------
        full_filename = Path(full_filename)
        full_hdu = fits.open(full_filename)

        # ----------determining object parameters------------
        args.id = full_hdu[0].header['ID']
        args.available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
        args.available_lines = np.array(['OIII' if item == 'OIII-5007' else item for item in args.available_lines]) # replace 'OIII-5007' with 'OIII'
        args.z = full_hdu[0].header['REDSHIFT']
        args.distance = cosmo.luminosity_distance(args.z)
        args.pix_size_arcsec = full_hdu[5].header['PIXASEC']
        args.pix_size_kpc = args.pix_size_arcsec / cosmo.arcsec_per_kpc_proper(args.z).value
        
        args.EB_V = 0. # until it gets over-written, if both H alpha and H beta lines are present
        offset = args.pix_size_arcsec / 2 # half a pixel offset to make sure cells in 2D plot are aligned with centers and not edges
        args.extent = (-args.arcsec_limit - offset, args.arcsec_limit - offset, -args.arcsec_limit - offset, args.arcsec_limit - offset)
        
        # --------determining true center of object rom direct image---------------------
        args.ndelta_xpix, args.ndelta_ypix = 0, 0

        # ---------------segmentation map---------------
        segmentation_map = full_hdu['SEG'].data
        segmentation_map = trim_image(segmentation_map, args, skip_re_trim=True)
        args.segmentation_map = segmentation_map

        # ------------setup figure-------------
        nrows, ncols = 3, 3
        npanels = len(line_list) + 2
        fig, axes = plt.subplots(nrows, ncols, figsize=(7, 7), layout='constrained')

        # ----------------line maps------------------
        for index, line in enumerate(line_list):
            row = index // ncols
            col = index % ncols
            ax = axes[row][col]
            if line in args.available_lines:
                line_map, _, _, _, _ = get_emission_line_map(line, full_hdu, args, dered=dered, silent=True)
                ax = plot_2D_map(line_map, ax, f'{line}', args, hide_xaxis=row < nrows - 1, hide_yaxis=col > 0)
            else:
                ax.remove()

        # -----------------metallicity maps-------------------
        index += 1
        ax = axes[index // ncols][index % ncols]
        if args.Zdiag == 'R23' and 'OII' not in args.available_lines:
            args.Zdiag = 'R3'
            switched_diag = True
            print(f'No OII for object {ind + 1}, so switching to R3 (instead of R23)..')
        
        logOH_map, _, _ = get_Z(full_hdu, args)
        ax = plot_2D_map(logOH_map, ax, f'Z ({args.Zdiag})', args)
        
        if switched_diag: args.Zdiag = 'R23'

        # -----------------radial profile-------------------
        index += 1
        ax = axes[index // ncols][index % ncols]
        quant = 'log_OH'
        shape = np.shape(logOH_map.data)
        center_xpix, center_ypix = shape[0] / 2., shape[1] / 2.
        distance_map = np.array([[np.sqrt((i - center_xpix)**2 + (j - center_ypix)**2) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size_kpc # Re or kpc

        logOH_df = pd.DataFrame({'distance': distance_map.flatten(), f'{quant}': unp.nominal_values(logOH_map.data).flatten(), f'{quant}_u': unp.std_devs(logOH_map.data).flatten()})
        logOH_df = logOH_df.dropna(subset=[f'{quant}', f'{quant}_u'], axis=0)

        ax, radial_linefit_odr = plot_radial_profile(logOH_df, ax, args, ylim=None, xlim=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, skip_annotate=False, quant=quant)
        
        # ------------deleting remaining axes---------
        for i in range(index + 1, nrows * ncols):
            ax = axes[i // ncols][i % ncols]
            ax.remove()
   
        # -----------save figure---------------
        plt.title(f'{full_filename.stem}')
        figname = f'{full_filename.stem}.png'
        save_fig(fig, fig_dir, figname, args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
