'''
    Filename: plot_stacked_maps.py
    Notes: Reads stacked (in bins of mass and/or SFR) 2D emission line maps for objects in a given field/s, computes metallicity gradient, and plots them
    Author : Ayan
    Created: 18-01-26
    Example: run plot_stacked_maps.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par28
             run plot_stacked_maps.py --field Par28 --Zdiag R23 --use_C25 --plot_line_maps --plot_metallicity --plot_radial_profiles --debug_bin
             run plot_stacked_maps.py --field Par28 --Zdiag R23 --use_C25 --fold_maps --plot_line_and_metallicity --debug_bin
             run plot_stacked_maps.py --field Par28 --Zdiag R23 --use_C25 --fold_maps --plot_line_and_metallicity --debug_bin --skip_deproject --skip_re_scaling
             run plot_stacked_maps.py --field Par28 --Zdiag R23 --use_C25 --fold_maps --plot_metallicity --debug_bin
             run plot_stacked_maps.py --field Par28 --Zdiag R23 --use_C25 --fold_maps --plot_radial_profiles --debug_bin
             run plot_stacked_maps.py --field Par28 --Zdiag R23 --use_C25 --adaptive_bins --fold_maps --plot_radial_profiles --debug_bin
             run plot_stacked_maps.py --field Par28 --Zdiag R23 --use_C25 --plot_radial_profiles --plot_minor_major_profile
             run plot_stacked_maps.py --field Par28 --Zdiag R23 --use_C25 --plot_radial_profiles
             run plot_stacked_maps.py --system ssd --do_all_fields --Zdiag R23 --use_C25 --plot_radial_profiles --adaptive_bins --bin_by_distance --fold_maps
             run plot_stacked_maps.py --system ssd --do_all_fields --Zdiag R23 --use_C25 --plot_line_and_metallicity --plot_radial_profiles --adaptive_bins --bin_by_distance_mass --fold_maps
             run plot_stacked_maps.py --system ssd --do_all_fields --plot_line_and_bpt --AGN_diag Ne3O2 --adaptive_bins --bin_by_distance_mass --fold_maps
'''

from header import *
from util import *
from make_diagnostic_maps import compute_Z_C19, compute_Z_KD02_R23, compute_Z_P25, compute_Z_Te, compute_Te, take_safe_log_ratio, take_safe_log_sum, myimshow, get_AGN_func_methods
from make_sfms_bins import get_binned_df, z_lim
from stack_emission_maps import read_stacked_maps
from plots_for_zgrad_paper import plot_fitted_line, odr_fit, plot_AGN_demarcation_ax, get_distance_map_from_AGN_line, get_ratio_labels, overplot_AGN_line_on_BPT

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def plot_stacked_line_maps(line_dict, args, bin_text='', cmin=None, cmax=None, takelog=True):
    '''
    Makes a nice plot of all emission lines present in a given list of stacked line maps
    Returns figure handle and the list of lines plotted
    '''
    # ----------getting line list from fits file-------------
    line_list = [item for item in list(line_dict.keys()) if '_nobj' not in item and '_id' not in item]
    n_lines = len(line_list)
    if n_lines == 0:
        fig = None
    else:
        # -----------------setup the figure---------------
        ncols = n_lines
        nrows = 2 if args.plot_radial_profiles else 1
        
        fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3.5 * nrows))
        axes = np.atleast_2d(axes)
        fig.subplots_adjust(left=0.07, right=0.96, top=0.92, bottom=0.1, wspace=0., hspace=0.2)
        fig.text(0.05, 0.95, f'{bin_text}', fontsize=args.fontsize, c='k', ha='left', va='top')

        # -----------------plot emission line maps of this bin---------------
        for index3, this_line in enumerate(line_list):
            this_map, _ , _= get_emission_line_map(this_line, line_dict, args)
            nobj = line_dict[f'{this_line}_nobj']
            if type(unp.nominal_values(this_map.data)) == np.ndarray:
                axes[0][index3] = plot_2D_map(this_map, axes[0][index3], f'Stacked {this_line}: {nobj}', args, cmap='cividis', clabel='', takelog=takelog, vmin=cmin, vmax=cmax, hide_xaxis=False, hide_yaxis=index3, hide_cbar=index3 < ncols - 1, hide_cbar_ticks=False, cticks_integer=True)

            # ----------getting distance maps, in case radial profile plotting-------------
            if args.plot_radial_profiles:
                shape = np.shape(this_map)
                if args.fold_maps:
                    center_xpix, center_ypix = (args.npix_side % 2 == 0) * 0.5, (args.npix_side % 2 == 0) * 0.5 # this yields 0.5 (instead of 0) pixel offset for even-sized stacked maps, because the center of the map is in the center (and not the edge) of the first pixel
                else:
                    center_xpix, center_ypix = shape[0] / 2., shape[1] / 2.
                distance_map = np.array([[np.sqrt((i - center_xpix)**2 + (j - center_ypix)**2) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size # Re or kpc
                minor_distance_map = np.array([[np.abs(j - center_ypix) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size # Re or kpc
                major_distance_map = np.array([[np.abs(i - center_xpix) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size # Re or kpc
                this_df = pd.DataFrame({'distance': distance_map.flatten(), 'minor_distance': minor_distance_map.flatten(), 'major_distance': major_distance_map.flatten()})

                this_df[f'{this_line}'] = unp.nominal_values(this_map.data).flatten()
                this_df[f'{this_line}_u'] = unp.std_devs(this_map.data).flatten()
                this_df = this_df.dropna(subset=[f'{this_line}', f'{this_line}_u'], axis=0)
                this_df = this_df[this_df[f'{this_line}'] > 0] # drop negative fluxes before taking log for the radial profile plot

                axes[1][index3], _, _, _ = plot_radial_profile(this_df, axes[1][index3], args, ylim=[cmin, cmax], takelog=takelog, xlim=None, hide_xaxis=False, hide_yaxis=index3, hide_cbar=True, skip_annotate=False, quant=this_line, skip_fitting=True)

    return fig, line_list         

# --------------------------------------------------------------------------------------------------------------------
def get_integrated_map(line_map):
    '''
    Integrate the given emission line map (with uncertainty) to compute a variance-weighted sum
    Inputs a unumpy.uarray object
    Returns integrated flux (along with uncertainty) as a unumpy.ufloat object
    '''
    flat_map = line_map.flatten()
    mask = ~unp.isnan(flat_map)
    clean_map = flat_map[mask]
    
    weights = 1.0 / (unp.std_devs(clean_map)**2)
    weighted_sum = (clean_map * weights).sum() / weights.sum()

    return weighted_sum

# --------------------------------------------------------------------------------------------------------------------
def get_emission_line_map(line, line_dict, args, silent=True):
    '''
    Retrieve the emission map for a given line from the given dictionary of emission lines
    Returns the 2D line map and the corresponding 2D uncertainty map
    '''
    line = line.upper()
    # ---------getting the spatillay resolved flux---------------
    line_map_object = line_dict[line]
    line_map_mask = line_map_object.mask
    line_map = line_map_object.data

    # ---------getting the integrated quanitty---------------
    line_int = get_integrated_map(line_map)
    nobj = line_dict[f'{line}_nobj']

    # ----------deblending flux--------------------
    factor = 1.0
    if not args.do_not_correct_flux:
        if line == 'OIII': # special treatment for OIII 5007 line, in order to account for and remove the OIII 4959 component
            ratio_5007_to_4959 = 2.98 # from grizli source code
            factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
            if not silent: print(f'Correcting OIII for 4959 component, by factor of {factor:.3f}')
        elif line == 'HA': # special treatment for Ha line, in order to account for and remove the NII component
            factor = 0.823 # from James et al. 2023?
            if not silent: print(f'Correcting Ha for NII component, by factor of {factor:.3f}')

    line_map = line_map * factor
    line_int = line_int * factor
    line_map = np.ma.masked_where(line_map_mask, line_map)

    return line_map, line_int, nobj

# --------------------------------------------------------------------------------------------------------------------
def get_Te(line_dict, args):
    '''
    Computes and returns the spatially resolved as well as intregrated Te from a given dictionary of emission lines
    '''
    OIII4363_map, OIII4363_int, _ = get_emission_line_map('OIII-4363', line_dict, args)
    OIII5007_map, OIII5007_int, _ = get_emission_line_map('OIII', line_dict, args)

    Te_map = compute_Te(OIII4363_map, OIII5007_map)
    Te_int = compute_Te(OIII4363_int, OIII5007_int)

    return Te_map, Te_int

# --------------------------------------------------------------------------------------------------------------------
def get_Z_Te(line_dict, args):
    '''
    Computes and returns the spatially resolved as well as intregrated Te metallicity from a given dictionary of emission lines
    '''
    OII3727_map, OII3727_int, OII3727_nobj = get_emission_line_map('OII', line_dict, args)
    OIII5007_map, OIII5007_int, OIII5007_nobj = get_emission_line_map('OIII', line_dict, args)
    Hbeta_map, Hbeta_int, Hbeta_nobj = get_emission_line_map('Hb', line_dict, args)
    nobj = min(OII3727_nobj, OIII5007_nobj, Hbeta_nobj)

    Te_map, Te_int = get_Te(line_dict, args)

    logOH_map = compute_Z_Te(OII3727_map, OIII5007_map, Hbeta_map, Te_map)
    logOH_int = compute_Z_Te(OII3727_int, OIII5007_int, Hbeta_int, Te_int)

    return logOH_map, logOH_int, nobj

# --------------------------------------------------------------------------------------------------------------------
def get_Z_P25(line_dict, args, branch='low'):
    '''
    Computes and returns the spatially resolved as well as intregrated Peluso+2025 metallicity from a given dictionary of emission lines
    '''
    SII6717_map, SII6717_int, SII6717_nobj = get_emission_line_map('SII', line_dict, args)
    OIII5007_map, OIII5007_int, OIII5007_nobj = get_emission_line_map('OIII', line_dict, args)
    Halpha_map, Halpha_int, Halpha_nobj = get_emission_line_map('Ha', line_dict, args)
    nobj = min(SII6717_nobj, OIII5007_nobj, Halpha_nobj)

    if not args.do_not_correct_flux:
        # special treatment for H-alpha line, in order to account for NII 6584 component
        factor = 0.823  # from grizli source code
    else:
        factor = 1.

    NII6584_map = np.ma.masked_where(Halpha_map.mask, Halpha_map.data * (1 - 0.823) / factor)
    NII6584_int = Halpha_int * (1 - 0.823) / factor

    logOH_map = compute_Z_P25(OIII5007_map, NII6584_map, SII6717_map, args.distance_from_AGN_line_map)
    logOH_int = compute_Z_P25(OIII5007_int, NII6584_int, SII6717_int, args.distance_from_AGN_line_int)

    return logOH_map, logOH_int, nobj
 
# --------------------------------------------------------------------------------------------------------------------
def get_Z_KD02_R23(line_dict, args, branch='low'):
    '''
    Computes and returns the spatially resolved as well as intregrated R23 metallicity from a given dictionary of emission lines
    '''
    OII3727_map, OII3727_int, OII3727_nobj = get_emission_line_map('OII', line_dict, args)
    OIII5007_map, OIII5007_int, OIII5007_nobj = get_emission_line_map('OIII', line_dict, args)
    Hbeta_map, Hbeta_int, Hbeta_nobj = get_emission_line_map('Hb', line_dict, args)
    nobj = min(OII3727_nobj, OIII5007_nobj, Hbeta_nobj)

    if not args.do_not_correct_flux:
        # special treatment for OIII 5007 line, in order to account for and ADD the OIII 4959 component back
        ratio_5007_to_4959 = 2.98  # from grizli source code
        factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
        print(f'Un-correcting OIII to include the 4959 component, for computing R23 metallicity, by factor of {factor:.3f}')
        OIII5007_map = np.ma.masked_where(OIII5007_map.mask, OIII5007_map.data / factor)
        OIII5007_int = OIII5007_int / factor

    logOH_map = compute_Z_KD02_R23(OII3727_map, OIII5007_map, Hbeta_map, branch=branch)
    logOH_int = compute_Z_KD02_R23(OII3727_int, OIII5007_int, Hbeta_int, branch=branch)

    return logOH_map, logOH_int, nobj

# --------------------------------------------------------------------------------------------------------------------
def get_emission_line_maps(line_dict, line_labels, args, silent=False):
    '''
    Computes and returns the spatially resolved as well as intregrated metallicity based on a given Curti+2019 calibration, from a given dictionary of emission lines
    '''

    # ----------getting line list from fits file-------------
    line_list = [item for item in list(line_dict.keys()) if '_nobj' not in item and '_id' not in item]
    if not set(list(map(str.upper, line_labels))) <= set(line_list):
        print(f'All lines in {line_labels} not available in this stack')
        return None, None, 0

    line_map_arr, line_int_arr, nobj_arr = [], [], []
    for line in line_labels:
        line_map, line_int, nobj = get_emission_line_map(line, line_dict, args, silent=silent)

        if not args.do_not_correct_flux and args.Zdiag == 'R23':
            factor = 1.
            if line == 'OIII': # special treatment for OIII 5007 line, in order to account for and ADD the OIII 4959 component back
                ratio_5007_to_4959 = 2.98  # from grizli source code
                factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
                print(f'Un-correcting OIII to include the 4959 component, for computing R23 metallicity, by factor of {factor:.3f}')

            if line_map is not None: line_map = np.ma.masked_where(line_map.mask, line_map.data / factor)
            line_int = line_int / factor

        line_map_arr.append(line_map)
        line_int_arr.append(line_int)
        nobj_arr.append(nobj)
    
    nobj = np.min(np.array(nobj_arr))

    return line_map_arr, line_int_arr, nobj

# --------------------------------------------------------------------------------------------------------------------
def get_Z_C19(line_dict, args):
    '''
    Computes and returns the spatially resolved as well as intregrated metallicity based on a given Curti+2019 and Cataldi+2025 calibrations, from a given dictionary of emission lines
    '''
    # ------getting appropriate emission lines and calibration coefficients--------------
    if args.Zdiag == 'O3O2':
        line_map_arr, line_int_arr, nobj = get_emission_line_maps(line_dict, ['OIII', 'OII'], args, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan), 0
        ratio_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1])
        try: ratio_int = unp.log10(line_int_arr[0] / line_int_arr[1])
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        if args.use_C25: coeff = [0.01124, 0.03072, 0.1251, -0.01470] # from Table 4 of Cataldi+25
        else: coeff = [-0.691, -2.944, -1.308]  # c0-2 parameters from Table 2 of Curti+19 3rd row (O3O2)

    elif args.Zdiag == 'R3':
        line_map_arr, line_int_arr, nobj = get_emission_line_maps(line_dict, ['OIII', 'Hb'], args, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan), 0
        ratio_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1])
        try: ratio_int = unp.log10(line_int_arr[0] / line_int_arr[1])
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        if args.use_C25: coeff = [-3.587, -4.475, 1.345, -0.08951] # from Table 4 of Cataldi+25
        else: coeff = [-0.277, -3.549, -3.593, -0.981]  # c0-3 parameters from Table 2 of Curti+19 2nd row (R3)

    elif args.Zdiag == 'R23':
        line_map_arr, line_int_arr, nobj = get_emission_line_maps(line_dict, ['OIII', 'OII', 'Hb'], args, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan), 0
        R1_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[2], skip_log=True)
        R2_map = take_safe_log_ratio(line_map_arr[1], line_map_arr[2], skip_log=True)
        ratio_map = take_safe_log_sum(R1_map, R2_map)
        try: ratio_int = unp.log10((line_int_arr[0] + line_int_arr[1]) / line_int_arr[2])
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        if args.use_C25: coeff = [-2.555, -3.192, 0.9630, -0.06351] # from Table 4 of Cataldi+25
        else: coeff = [0.527, -1.569, -1.652, -0.421]  # c0-3 parameters from Table 2 of Curti+19 4th row (R23)

    elif args.Zdiag == 'O3S2':
        line_map_arr, line_int_arr, nobj = get_emission_line_maps(line_dict, ['OIII', 'Hb', 'SII', 'Ha'], args, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan), 0
        R1_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1], skip_log=True)
        R2_map = take_safe_log_ratio(line_map_arr[2], line_map_arr[3], skip_log=True)
        ratio_map = take_safe_log_ratio(R1_map, R2_map)
        try: ratio_int = unp.log10((line_int_arr[0] / line_int_arr[1]) / (line_int_arr[2] / line_int_arr[3]))
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        coeff = [0.191, -4.292, -2.538, 0.053, 0.332]  # c0-4 parameters from Table 2 of Curti+19 last row (O3S2)

    elif args.Zdiag == 'S2':
        line_map_arr, line_int_arr, nobj = get_emission_line_maps(line_dict, ['SII', 'Ha'], args, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan), 0
        ratio_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1])
        try: ratio_int = unp.log10(line_int_arr[0] / line_int_arr[1])
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        coeff = [-0.442, -0.360, -6.271, -8.339, -3.559]  # c0-3 parameters from Table 2 of Curti+19 3rd-to-last row (S2)

    elif args.Zdiag == 'RS32':
        line_map_arr, line_int_arr, nobj = get_emission_line_maps(line_dict, ['OIII', 'Hb', 'SII', 'Ha'], args, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan), 0
        R1_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1], skip_log=True)
        R2_map = take_safe_log_ratio(line_map_arr[2], line_map_arr[3], skip_log=True)
        ratio_map = take_safe_log_sum(R1_map, R2_map)
        try: ratio_int = unp.log10((line_int_arr[0] / line_int_arr[1]) + (line_int_arr[2] / line_int_arr[3]))
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        coeff = [-0.054, -2.546, -1.970, 0.082, 0.222]  # c0-3 parameters from Table 2 of Curti+19 2nd-to-last row (RS32)

    elif args.Zdiag == 'R2':
        line_map_arr, line_int_arr, nobj = get_emission_line_maps(line_dict, ['OII', 'Hb'], args, silent=args.only_integrated)
        if line_map_arr is None: return None, ufloat(np.nan, np.nan), 0
        ratio_map = take_safe_log_ratio(line_map_arr[0], line_map_arr[1])
        try: ratio_int = unp.log10(line_int_arr[0] / line_int_arr[1])
        except ValueError: ratio_int = ufloat(np.nan, np.nan)
        if args.use_C25: coeff = [-6.481, 0.8163] # from Table 4 of Cataldi+25
        else: coeff = [0.435, -1.362, -5.655, -4.851, -0.478, 0.736]  # c0-3 parameters from Table 2 of Curti+19 1st row (R2)

    else:
        print(f'Could not apply any of the metallicity diagnostics, so returning NaN metallicities')
        return None, ufloat(np.nan, np.nan), 0

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
    logOH_map = compute_Z_C19(ratio_map, coeff, ax=ax, branch=args.Zbranch, use_C25=args.use_C25, silent=True)
    logOH_int = compute_Z_C19(ratio_int, coeff, branch=args.Zbranch, use_C25=args.use_C25, silent=True)
    logOH_int = ufloat(unp.nominal_values(np.atleast_1d(logOH_int))[0], unp.std_devs(np.atleast_1d(logOH_int))[0])

    # -------saving the debugging plots---------------
    if ax is not None:
        Zbranch_text = '' if args.Zdiag in ['NB', 'P25', 'Te'] else f'-{args.Zbranch}'
        figname = args.fig_dir / f'stacked_metallicity_debug_Zdiag_{args.Zdiag}{Zbranch_text}.png'
        fig.savefig(figname, transparent=args.fortalk, dpi=200)
        print(f'\nSaved figure at {figname}')

    return logOH_map, logOH_int, nobj

# --------------------------------------------------------------------------------------------------------------------
def compute_Z_NB(line_label_array, line_waves_array, line_flux_array):
    '''
    Calculates and returns the NebulaBayes metallicity given a list of observed line fluxes
    '''
    line_flux_array = [np.atleast_1d(item) for item in line_flux_array]
    
    map_shape = np.shape(line_flux_array[0])
    if len(map_shape) == 1: npixels = map_shape[0]
    else: npixels = map_shape[0] * map_shape[1]

    IDs_array = np.arange(npixels).flatten()
    unique_IDs_array = np.unique(IDs_array)  

    if len(line_label_array) <= 1:
        if hasattr(line_flux_array[0][0], "__len__"): return np.ma.masked_where(True, unp.uarray(np.ones(np.shape(line_flux_array[0])) * np.nan, np.ones(np.shape(line_flux_array[0])) * np.nan))
        else: return ufloat(np.nan, np.nan)

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

    out_dir = args.fits_dir / f'{bin_text}_NB_results'
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
            good_obs = (obs_fluxes >= 0) & (obs_errs > 0)
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
                #print(f'Deb1576: binID {this_ID}: nlines={len(obs_fluxes)}, {dict(zip(line_labels, obs_fluxes))}, norm_line = {norm_line}, dereddening on the fly? {dered}') ##
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
            
                counter += 1
            else:
                logOH = ufloat(np.nan, np.nan)
                print(f'Could not run NB for unique ID {this_ID} ({counter + 1} out of {len(unique_IDs_array)}) with only {len(obs_fluxes)} good fluxes')

            logOH_dict_unique_IDs.update({this_ID: logOH}) # updating to unique ID dictionary once logOH has been calculated for this unique ID

        logOH_array.append(logOH)
    print(f'\nRan NB for total {counter} unique pixels out of {len(obs_flux_array[0])}, in {timedelta(seconds=(datetime.now() - start_time3).seconds)}\n')

    # ---------collating all the metallicities computed---------
    log_OH = np.ma.masked_where(net_mask, np.reshape(logOH_array, np.shape(line_flux_array[0])))
    if len(log_OH) == 1: log_OH = log_OH.data[0]

    return log_OH

# --------------------------------------------------------------------------------------------------------------------
def get_Z_NB(line_dict, args):
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

    line_list = [item for item in list(line_dict.keys()) if '_nobj' not in item and '_id' not in item]
    line_map_array, line_int_array, line_label_array, line_waves_array, nobj_array = [], [], [], [], []
    
    for line in line_list:
        line_map, line_int, nobj = get_emission_line_map(line, line_dict, args, silent=False)
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

        if line in line_label_dict.keys() and line not in args.exclude_lines:
            nobj_array.append(nobj)
            line_map_array.append(line_map)
            line_int_array.append(line_int)
            line_label_array.append(line_label_dict[line])
            line_waves_array.append(rest_wave_dict[line] * 10) # factor of 10 to convert from nm to Angstroms

    # ----------calling NB----------------------
    logOH_map = compute_Z_NB(line_label_array, line_waves_array, line_map_array)
    logOH_int = compute_Z_NB(line_label_array, line_waves_array, line_int_array)

    return logOH_map, logOH_int, line_label_array, np.min(nobj_array)

# --------------------------------------------------------------------------------------------------------------------
def get_metallicity_map(line_dict, args):
    '''
    Computes 2D metallicity and integrated metallicity from a given dictionary of several stacked lines 
    Returns metallicity map (2D array) along with uncertainty
    '''
    # ----------getting line list from fits file-------------
    line_list = [item for item in list(line_dict.keys()) if '_nobj' not in item and '_id' not in item]

    # --------deriving the metallicity map-------------
    if args.Zdiag == 'NB':
        logOH_map, logOH_int, line_label_array, nobj = get_Z_NB(line_dict, args)
    elif args.Zdiag == 'KD02_R23' and set(['OIII', 'OII', 'HB']) <= set(line_list):
        logOH_map, logOH_int, nobj = get_Z_KD02_R23(line_dict, args, branch=args.Zbranch)
    elif args.Zdiag == 'Te' and set(['OIII', 'OIII-4363', 'OII', 'HB']) <= set(line_list):
        logOH_map, logOH_int, nobj = get_Z_Te(line_dict, args)
    elif args.Zdiag == 'P25' and set(['OIII', 'HA', 'SII']) <= set(line_list):
        logOH_map, logOH_int, nobj = get_Z_P25(line_dict, args)
    else:
        logOH_map, logOH_int, nobj = get_Z_C19(line_dict, args)

    return logOH_map, logOH_int, nobj

# --------------------------------------------------------------------------------------------------------------------
def write_metallicity_map(logOH_map, logOH_int, outfilename, args, nobj=None):
    '''
    Writes the 2D metallicity map and the integrated metallicity (and uncertainties) as a fits file
    '''
    logOH_map_val = np.where(logOH_map.mask, np.nan, unp.nominal_values(logOH_map.data))
    logOH_map_err = np.where(logOH_map.mask, np.nan, unp.std_devs(logOH_map.data))

    # --------processing the dataframe in case voronoi binning has been performed and there are duplicate data values------
    hdr1, hdr2 = fits.Header(), fits.Header()
    primary_hdu = fits.PrimaryHDU(header=hdr1)

    hdr2['z_diag'] = args.Zdiag
    hdr2['zbranch'] = args.Zbranch
    hdr2['log_oh_int'] = None if ~np.isfinite(logOH_int.n) else logOH_int.n
    hdr2['log_oh_int_err'] = None if ~np.isfinite(logOH_int.s) else logOH_int.s
    if 'NB' in args.Zdiag:
        hdr2['nb_old_grid'] = True if args.use_original_NB_grid and args.Zdiag == 'NB' else False
    if nobj is not None:
        hdr2['nobj'] = f'{nobj}'


    logOH_val_hdu = fits.ImageHDU(data=logOH_map_val, name='log_OH', header=hdr2)
    logOH_err_hdu = fits.ImageHDU(data=logOH_map_err, name='log_OH_u')
    hdul = fits.HDUList([primary_hdu, logOH_val_hdu, logOH_err_hdu])
    hdul.writeto(outfilename, overwrite=True)

    print(f'Saved metallicity maps in {outfilename}')

# --------------------------------------------------------------------------------------------------------------------
def read_metallicity_map(infilename):
    '''
    Reads the 2D metallicity map and the integrated metallicity (and uncertainties) from a fits file
    '''
    print(f'\nReading from existing file {infilename}')

    hdul = fits.open(infilename)
    
    logOH_map_val = hdul['log_OH'].data
    logOH_map_err = hdul['log_OH_u'].data
    mask = np.isnan(logOH_map_val)
    logOH_map = np.ma.masked_where(mask, unp.uarray(logOH_map_val, logOH_map_err))

    header = hdul['log_OH'].header
    if 'log_oh_int' in header:
        val, err = header['log_oh_int'], header['log_oh_int_err']
        if val is not None and err is not None: logOH_int = ufloat(float(val), float(err))
        else: logOH_int = ufloat(np.nan, np.nan)
    else:
        logOH_int = ufloat(np.nan, np.nan)
    if 'nobj' in header: nobj = header['nobj']
    else: nobj = np.nan

    return logOH_map, logOH_int, nobj

# --------------------------------------------------------------------------------------------------------------------
def plot_profile(df, ax, linefit_odr, quant_x, quant_y, col='grey', index=0, skip_fitting=False):
    '''
    Plots the radial profile from a given dataframe in a given axis
    Returns the axis handle
    '''
    # -------plotting the data and the fits--------
    ax.scatter(df[quant_x], df[quant_y], c=col, s=20, lw=1, edgecolor='k', alpha=1)
    if quant_y + '_u' in df: ax.errorbar(df[quant_x], df[quant_y], yerr=df[quant_y + '_u'], c=col, fmt='none', lw=0.5, alpha=0.2)

    if not skip_fitting:
        xarr = df[quant_x]
        ax = plot_fitted_line(ax, linefit_odr, xarr, col, args, quant=quant_y, short_label=False, index=index)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_radial_profile(df, ax, args, ylim=None, xlim=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, skip_annotate=False, takelog=False, quant='log_OH', skip_fitting=False, yaxis_on_right=False):
    '''
    Plots and fits the radial profile from a given dataframe in a given axis
    Returns the axis handle and the linefit
    '''
    label_dict = smart_dict({'SFR': r'$\log$ $\Sigma_{\rm SFR}$ (M$_{\odot}$/yr/kpc$^2$)', 'logOH': r'$\log$ (O/H) + 12', 'Z': r'$\log$ (O/H) + 12', 'log_OH': r'$\log$ (O/H) + 12'})

    df = df.sort_values(by='distance').reset_index(drop=True)
    if takelog:
        if f'{quant}_u' in df:
            value_with_unc = unp.uarray(df[quant], df[f'{quant}_u'])
            log_value_with_unc = unp.log10(value_with_unc)
            df[quant] = unp.nominal_values(log_value_with_unc)
            df[f'{quant}_u'] = unp.std_devs(log_value_with_unc)
        else:
            df[quant] = np.log10(df[quant])

    # -------fitting by various methods--------------
    width = 1 if args.re_limit is None else 0.2 # 1 kpc or 0.2 Re
    df_minor = df[df['major_distance'] <= width]
    minor_linefit_odr = odr_fit(df_minor, quant_x='minor_distance', quant_y=quant)
    
    df_major = df[df['minor_distance'] <= width]
    major_linefit_odr = odr_fit(df_major, quant_x='major_distance', quant_y=quant)
    
    radial_linefit_odr = odr_fit(df, quant_x='distance', quant_y=quant)

    if args.plot_minor_major_profile:
        ax = plot_profile(df_minor, ax, minor_linefit_odr, 'minor_distance', quant, col='salmon', index=0, skip_fitting=skip_fitting)
        ax = plot_profile(df_major, ax, major_linefit_odr, 'major_distance', quant, col='cornflowerblue', index=1, skip_fitting=skip_fitting)
    else:
        ax = plot_profile(df, ax, radial_linefit_odr, 'distance', quant, col='darkgoldenrod' if args.fortalk else 'grey', skip_fitting=skip_fitting)

    ax.set_aspect('auto') 

    # --------annotating axis--------------
    if not skip_annotate: ax = annotate_axes(ax, 'Radius (kpc)' if args.re_limit is None else r'Radius (R$_e$)', label_dict[quant], args=args, xlim=xlim, ylim=ylim, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar, yaxis_on_right=yaxis_on_right)
    
    return ax, minor_linefit_odr, major_linefit_odr, radial_linefit_odr

# --------------------------------------------------------------------------------------------------------------------
def plot_2D_map(image, ax, label, args, cmap='cividis', clabel='', takelog=True, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, hide_cbar_ticks=False, cticks_integer=True):
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
    
    units = 'kpc' if args.skip_re_scaling else r'R$_e$' 
    ax = annotate_axes(ax, f'Offset ({units})', f'Offset ({units})', args=args, label=label, labelx=0.05, labely=0.92, clabel=clabel, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar, p=p, hide_cbar_ticks=hide_cbar_ticks, cticks_integer=cticks_integer)

    if args.plot_radial_profiles and args.plot_minor_major_profile:
        width = 1 if args.re_limit is None else 0.2 # 1 kpc or 0.2 Re
        if args.fold_maps:
            ax.add_patch(plt.Rectangle((0, 0), args.extent[1], width, lw=2, color='salmon', fill=False))
            ax.add_patch(plt.Rectangle((0, 0), width, args.extent[3], lw=2, color='cornflowerblue', fill=False))
        else:
            ax.add_patch(plt.Rectangle((args.extent[0], -width/2), args.extent[1] - args.extent[0], width, lw=2, color='salmon', fill=False))
            ax.add_patch(plt.Rectangle((-width/2, args.extent[2]), width, args.extent[3] - args.extent[2], lw=2, color='cornflowerblue', fill=False))

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_map(logOH_map, args, bin_text='', Zmin=None, Zmax=None):
    '''
    Makes a nice plot of a given metallicity map
    Returns figure handle
    '''
    # -----------------setup the figure---------------
    fig, axes = plt.subplots(1, 2 if args.plot_radial_profiles else 1, figsize=(10, 4) if args.plot_radial_profiles else (8, 6))
    fig.subplots_adjust(left=0.1, right=0.98 if args.plot_radial_profiles else 0.85, top=0.98, bottom=0.14, wspace=0.6, hspace=0.)
    axes = np.atleast_1d(axes)
    
    # -----------------plot metallicity maps of this bin---------------
    axes[0] = plot_2D_map(logOH_map, axes[0], f'Stacked logOH {bin_text}', args, cmap='plasma', clabel=r'$\log$(O/H) + 12', takelog=False, vmin=Zmin, vmax=Zmax, hide_cbar=False, cticks_integer=True)
        
    if args.plot_radial_profiles:
        quant = 'log_OH'
        shape = np.shape(logOH_map)
        if args.fold_maps:
            center_xpix, center_ypix = (args.npix_side % 2 == 0) * 0.5, (args.npix_side % 2 == 0) * 0.5 # this yields 0.5 (instead of 0) pixel offset for even-sized stacked maps, because the center of the map is in the center (and not the edge) of the first pixel
        else:
            center_xpix, center_ypix = shape[0] / 2., shape[1] / 2.
        
        distance_map = np.array([[np.sqrt((i - center_xpix)**2 + (j - center_ypix)**2) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size # Re or kpc
        minor_distance_map = np.array([[np.abs(j - center_ypix) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size # Re or kpc
        major_distance_map = np.array([[np.abs(i - center_xpix) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size # Re or kpc
        df = pd.DataFrame({'distance': distance_map.flatten(), \
                           'minor_distance': minor_distance_map.flatten(), \
                           'major_distance': major_distance_map.flatten(), \
                           f'{quant}': unp.nominal_values(logOH_map).flatten(), \
                           f'{quant}_u': unp.std_devs(logOH_map).flatten(), \
                           })
        df = df.dropna(subset=[f'{quant}', f'{quant}_u'], axis=0)
        axes[1], minor_linefit_odr, major_linefit_odr, radial_linefit_odr = plot_radial_profile(df, axes[1], args, ylim=[Zmin, Zmax], xlim=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, skip_annotate=False, quant=quant, skip_fitting=False)
    else:
        minor_linefit_odr, major_linefit_odr, radial_linefit_odr = None, None, None

    return fig, minor_linefit_odr, major_linefit_odr, radial_linefit_odr      

# --------------------------------------------------------------------------------------------------------------------
def get_interval_from_filename(filename):
    '''
    Accepts a filename of the pattern stacked_maps_logmassbin_X-Y_logsfrbin_A-B.fits and extracts the mass and SFR interval from the filename
    Returns object: pair of Intervals
    '''
    filename = Path(filename).stem
    pattern1 = r'logmassbin_([-?\d.]+)-([-?\d.]+)_logsfrbin_([-?\d.]+)-([-?\d.]+)'
    match1 = re.search(pattern1, filename)
    pattern2 = r'delta_sfms_bin_([-]?\d+(?:\.\d+)?)-([-]?\d+(?:\.\d+)?)'
    match2 = re.search(pattern2, filename)

    if match1:
        m_low, m_high, s_low, s_high = map(float, match1.groups())        
        mass_int = pd.Interval(m_low, m_high, closed='left')
        sfr_int = pd.Interval(s_low, s_high, closed='left')
        
        combined_interval = (mass_int, sfr_int)
    elif match2:
        delta_sfms_low, delta_sfms_high = map(float, match2.groups())
        delta_sfms_interval = pd.Interval(delta_sfms_low, delta_sfms_high, closed='left')
        
        combined_interval = delta_sfms_interval
    else:
        raise ValueError(f'The file {filename} does not have the expected pattern. So could not extract intervals from it.')
    
    return combined_interval

# --------------------------------------------------------------------------------------------------------------------
def fold_line_map(line_map, args, line=''):
    '''
    Fold a given emission line map along vertical and horizontal axis (along with the uncertainty map), by calling fold_image_to_quadrant()
    Returns folded 2D array (along with uncertainties) which is a 'quarter' map
    '''
    folded_data = fold_image_to_quadrant(unp.nominal_values(line_map.data), args, line=f'{line}')
    folded_err = fold_image_to_quadrant(unp.std_devs(line_map.data), args, line=f'{line}_err')
    if np.ndim(line_map.mask) > 0:
        folded_mask = fold_image_to_quadrant(line_map.mask, args, line=f'{line}_mask')
    else:
        folded_mask = line_map.mask
    
    folded_line_map = np.ma.masked_where(folded_mask, unp.uarray(folded_data, folded_err))

    return folded_line_map

# --------------------------------------------------------------------------------------------------------------------
def fold_image_to_quadrant(data_map, args, line=''):
    '''
    Fold a given image along vertical and horizontal axis (this assumes that the major and minor axes align with the vertical and horizontal axes, respectively)
    Returns 2D array which is a 'quarter' map
    '''
    mid_x, mid_y = data_map.shape[0] // 2, data_map.shape[1] // 2
    
    # ----extract the four quadrants: Top-Right (TR), Top-Left (TL), Bottom-Left (BL), Bottom-Right (BR)---
    tr = data_map[mid_y:, mid_x:]            # Already oriented correctly
    tl = np.fliplr(data_map[mid_y:, :mid_x]) # Flip horizontally
    bl = np.flipud(np.fliplr(data_map[:mid_y, :mid_x])) # Flip both
    br = np.flipud(data_map[:mid_y, mid_x:]) # Flip vertically
    
    # ---------taking their median-----------
    folded_quadrant = np.median([tr, tl, bl, br], axis=0)

    # ----------for diagnostics------------
    if args.debug_folding:
        print(f'Deb467: shapes: {line} map = {np.shape(data_map)}, tr = {np.shape(tr)}, tl = {np.shape(tl)}, bl = {np.shape(bl)}, br = {np.shape(br)}, folded = {np.shape(folded_quadrant)}') ##
        fig, axes = plt.subplots(1, 6, figsize=(13, 3))
        fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.1, wspace=0.3, hspace=0.1)
        cmap = 'viridis'
        fig.text(0.05, 0.98, f'{args.field}: ID {args.id}: {line} map: Folding diagnostics', fontsize=args.fontsize, c='k', ha='left', va='top')
        
        axes[0] = myimshow(data_map, axes[0], label='Original', cmap=cmap, col='k')
        axes[1] = myimshow(tr, axes[1], label='TR', cmap=cmap, col='k')
        axes[2] = myimshow(tl, axes[2], label='TL', cmap=cmap, col='k')
        axes[3] = myimshow(br, axes[3], label='BR', cmap=cmap, col='k')
        axes[4] = myimshow(bl, axes[4], label='BL', cmap=cmap, col='k')
        axes[5] = myimshow(folded_quadrant, axes[5], label='Folded', cmap=cmap, col='k')

        plt.show(block=False)
        sys.exit(f'Exiting here because of --debug_folding mode; if you want to run the full code as usual then remove the --debug_folding option and re-run')
    
    return folded_quadrant

# --------------------------------------------------------------------------------------------------------------------
def plot_line_and_metallicity_maps(line_dict, logOH_map, args, bin_text='', cmin=None, cmax=None, Zmin=None, Zmax=None, takelog=True):
    '''
    Makes a nice plot of all emission lines present in a given list of stacked line maps along with the metallicity map
    Returns figure handle and the list of lines plotted
    '''
    # ----------getting line list from fits file-------------
    line_list = [item for item in list(line_dict.keys()) if '_nobj' not in item and '_id' not in item]
    n_lines = len(line_list)
    if n_lines == 0:
        fig = None
    else:
        # ----------getting distance maps, in case radial profile plotting-------------
        if args.plot_radial_profiles:
            shape = np.shape(logOH_map)
            if args.fold_maps:
                center_xpix, center_ypix = (args.npix_side % 2 == 0) * 0.5, (args.npix_side % 2 == 0) * 0.5 # this yields 0.5 (instead of 0) pixel offset for even-sized stacked maps, because the center of the map is in the center (and not the edge) of the first pixel
            else:
                center_xpix, center_ypix = shape[0] / 2., shape[1] / 2.
            distance_map = np.array([[np.sqrt((i - center_xpix)**2 + (j - center_ypix)**2) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size # Re or kpc
            minor_distance_map = np.array([[np.abs(j - center_ypix) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size # Re or kpc
            major_distance_map = np.array([[np.abs(i - center_xpix) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size # Re or kpc
            df = pd.DataFrame({'distance': distance_map.flatten(), 'minor_distance': minor_distance_map.flatten(), 'major_distance': major_distance_map.flatten()})

        # -----------------setup the figure---------------
        ncols = n_lines + 1
        nrows = 2 if args.plot_radial_profiles else 1
        
        fig, axes = plt.subplots(nrows, ncols, figsize=(2.5 * ncols, 3.5 * nrows))
        axes = np.atleast_2d(axes)
        fig.subplots_adjust(left=0.06, right=0.94, top=0.94, bottom=0.08, wspace=0., hspace=0.2)
        fig.text(0.05, 0.95, f'{bin_text}', fontsize=args.fontsize, c='k', ha='left', va='top')
        
        # -----------------plot emission line maps of this bin---------------
        for index, this_line in enumerate(line_list):
            this_map, _ , _= get_emission_line_map(this_line, line_dict, args)
            nobj = line_dict[f'{this_line}_nobj']
            if type(unp.nominal_values(this_map.data)) == np.ndarray:
                axes[0][index] = plot_2D_map(this_map, axes[0][index], f'{this_line}: {nobj}', args, cmap='cividis', clabel='', takelog=takelog, vmin=cmin, vmax=cmax, hide_xaxis=False, hide_yaxis=index, hide_cbar=True, hide_cbar_ticks=False, cticks_integer=True)

            if args.plot_radial_profiles:
                this_df = df.copy()
                this_df[f'{this_line}'] = unp.nominal_values(this_map.data).flatten()
                this_df[f'{this_line}_u'] = unp.std_devs(this_map.data).flatten()
                this_df = this_df.dropna(subset=[f'{this_line}', f'{this_line}_u'], axis=0)
                this_df = this_df[this_df[f'{this_line}'] > 0] # drop negative fluxes before taking log for the radial profile plot
                
                axes[1][index], _, _, _ = plot_radial_profile(this_df, axes[1][index], args, ylim=[cmin, cmax], takelog=takelog, xlim=None, hide_xaxis=False, hide_yaxis=index, hide_cbar=True, skip_annotate=False, quant=this_line, skip_fitting=True)

        # -----------------plot metallicity maps of this bin---------------
        axes[0][index + 1] = plot_2D_map(logOH_map, axes[0][index + 1], f'logOH {args.Zdiag}', args, cmap='plasma', clabel=r'$\log$(O/H) + 12', takelog=False, vmin=Zmin, vmax=Zmax, hide_yaxis=True, hide_cbar=False, cticks_integer=True)
        
        if args.plot_radial_profiles:
            quant = 'log_OH'
            logOH_df = df.copy()
            logOH_df[f'{quant}'] = unp.nominal_values(logOH_map).flatten()
            logOH_df[f'{quant}_u'] = unp.std_devs(logOH_map).flatten()
            logOH_df = logOH_df.dropna(subset=[f'{quant}', f'{quant}_u'], axis=0)

            axes[1][index + 1], minor_linefit_odr, major_linefit_odr, radial_linefit_odr = plot_radial_profile(logOH_df, axes[1][index + 1], args, ylim=[Zmin, Zmax], xlim=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, skip_annotate=False, quant=quant, skip_fitting=False, yaxis_on_right=True)
        else:
            minor_linefit_odr, major_linefit_odr, radial_linefit_odr = None, None, None

    return fig, line_list, minor_linefit_odr, major_linefit_odr, radial_linefit_odr

# --------------------------------------------------------------------------------------------------------------------
def plot_AGN_demarcation_ax(x_num, x_den, y_num, y_den, ax, args, color=None, marker='o', size=20, lw=0.5):
    '''
    Plots line ratio vs line ratio on a given axis
    Returns axis handle and the scatter plot handle
    '''
    shape = np.shape(y_num)
    if args.fold_maps:
        center_xpix, center_ypix = (args.npix_side % 2 == 0) * 0.5, (args.npix_side % 2 == 0) * 0.5 # this yields 0.5 (instead of 0) pixel offset for even-sized stacked maps, because the center of the map is in the center (and not the edge) of the first pixel
    else:
        center_xpix, center_ypix = shape[0] / 2., shape[1] / 2.
    distance_map = np.array([[np.sqrt((i - center_xpix)**2 + (j - center_ypix)**2) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size # Re or kpc

    df = pd.DataFrame({'xnum': unp.nominal_values(np.atleast_1d(x_num)).flatten(), \
                       'xnum_u': unp.std_devs(np.atleast_1d(x_num)).flatten(), \
                       'xden': unp.nominal_values(np.atleast_1d(x_den)).flatten(), \
                       'xden_u': unp.std_devs(np.atleast_1d(x_den)).flatten(), \
                       'ynum': unp.nominal_values(np.atleast_1d(y_num)).flatten(), \
                       'ynum_u': unp.std_devs(np.atleast_1d(y_num)).flatten(), \
                       'yden': unp.nominal_values(np.atleast_1d(y_den)).flatten(), \
                       'yden_u': unp.std_devs(np.atleast_1d(y_den)).flatten(), \
                       'distance': np.atleast_1d(distance_map).flatten()
                       })
    df = df.drop_duplicates().reset_index(drop=True)
    df = df[(df['xnum'] > 0) & (df['xden'] > 0) & (df['ynum'] > 0) & (df['yden'] > 0)]

    y_ratio = unp.log10(unp.uarray(df['ynum'], df['ynum_u']) / unp.uarray(df['yden'], df['yden_u']))
    x_ratio = unp.log10(unp.uarray(df['xnum'], df['xnum_u']) / unp.uarray(df['xden'], df['xden_u']))

    if color is None:
        color = get_distance_map_from_AGN_line(df['xnum'], df['xden'], df['ynum'], df['yden'], args).data
        cmap, vmin, vmax = args.diverging_cmap, -1, 1
    elif color == 'distance':
        color = df['distance']
        cmap, vmin, vmax = 'viridis', None, None

    p = ax.scatter(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), c=color, cmap=cmap, vmin=vmin, vmax=vmax, marker=marker, s=size, lw=lw, edgecolor='w' if args.fortalk else 'k', zorder=10)
    ax.errorbar(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), xerr=unp.std_devs(x_ratio), yerr=unp.std_devs(y_ratio), c='gray', fmt='none', lw=lw, alpha=0.5)

    return ax, p

# --------------------------------------------------------------------------------------------------------------------
def plot_AGN_demarcation_object(line_dict, args, ax, marker='o', size=50, ax_inset=None, hide_cbar=True):
    '''
    Plots the spatially resolved AGN demarcation for a given object on a given axis
    Returns axis handle
    '''    
    theoretical_lines, line_labels = get_AGN_func_methods(args)

    # -----------getting the fluxes------------------
    try:
        ynum_map, _, _ = get_emission_line_map(args.ynum_line, line_dict, args, silent=True)
        yden_map, _, _ = get_emission_line_map(args.yden_line, line_dict, args, silent=True)

        xnum_map, _, _ = get_emission_line_map(args.xnum_line, line_dict, args, silent=True)
        xden_map, _, _ = get_emission_line_map(args.xden_line, line_dict, args, silent=True)
    except:
        print(f'Required emission lines not available for {args.id} with {args.AGN_diag} AGN diagnostic. So skipping this object')
        return ax, None, None
    
    if not args.do_not_correct_flux and args.AGN_diag in ['H21', 'B22'] and args.xden_line == 'Ha': # special treatment for H-alpha line, in order to add the NII 6584 component back
        factor = 0.823  # from grizli source code
        print(f'Adding the NII component back to Ha, i.e. dividing by factor {factor} because using Henry+21 for AGN-SF separation')
        xden_map = np.ma.masked_where(xden_map.mask, xden_map.data / factor)

   # -----------spatially_resolved-----------------------
    ax, scatter_plot_handle = plot_AGN_demarcation_ax(xnum_map, xden_map, ynum_map, yden_map, ax, args, marker=marker, size=size, lw=0.5, color=args.AGN_colorby if 'AGN_colorby' in args else None)
    ax.set_aspect('auto')

    # -----------2D map inset-----------------------
    if ax_inset is not None:
        distance_from_AGN_line_map = get_distance_map_from_AGN_line(xnum_map, xden_map, ynum_map, yden_map, args)
        ax_inset = plot_2D_map(distance_from_AGN_line_map, ax_inset, '', args, takelog=False, cmap=args.diverging_cmap, clabel='Distance from\nNB line', vmin=-1, vmax=1, hide_xaxis=True, hide_yaxis=True, hide_cbar=True)

    # -----------annotating axes-----------------------
    if not hide_cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='2%', pad=0.05)
        cbar = plt.colorbar(scatter_plot_handle, cax=cax)
        unit_text = r'kpc' if args.re_limit is None else r'R_e'
        cbar.set_label(f'Distance from center ({unit_text})' if 'AGN_colorby' in args and args.AGN_colorby == 'distance' else f'Distance from {theoretical_lines[0]} line', fontsize=args.fontsize)
        cbar.ax.tick_params(labelsize=args.fontsize)

    ax.set_xlim(-1, 0.2)
    ax.set_ylim(-1., 1.)
    ax.set_xlabel(f'Log {get_ratio_labels("NeIII-3867/OII")}' if args.AGN_diag == 'Ne3O2' else f'Log {get_ratio_labels("SII/NII,Ha")}' if args.AGN_diag == 'H21' else f'Log {get_ratio_labels(f"{args.xnum_line}/{args.xden_line}")}', fontsize=args.fontsize)
    ax.set_ylabel(f'Log {get_ratio_labels("OIII/Hb")}', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    # ---------adding literature AGN demarcation lines----------
    color_arr = ['brown', 'darkgreen', 'dodgerblue', 'cyan', 'sienna']
    for index, (theoretical_line, line_label) in enumerate(zip(theoretical_lines, line_labels)):
        if line_label == 'This work': line_label = 'NB'
        overplot_AGN_line_on_BPT(ax, theoretical_line, label=line_label, label_loc='upper left', color=color_arr[index], fontsize=args.fontsize / args.fontfactor, lw=0.5 if index else 1, ls='solid')

    return ax, scatter_plot_handle, ax_inset

# --------------------------------------------------------------------------------------------------------------------
def plot_line_and_bpt_maps(line_dict, args, bin_text='', cmin=None, cmax=None, takelog=True):
    '''
    Makes a nice plot of all emission lines present in a given list of stacked line maps along with the metallicity map
    Returns figure handle and the list of lines plotted
    '''
    # ----------getting line list from fits file-------------
    line_list = ['OIII', 'HB', 'OII', 'NEIII-3867']
    absent_lines = [item for item in line_list if item.upper() not in line_dict.keys()]
    if len(absent_lines) > 0:
        fig = None
        print(f'\n{absent_lines} line/s not availble in the stack, so cannot make OHNO BPT\n')
    else:
        # -------setting up the figure--------------------
        fig = plt.figure(figsize=(10, 5))
        fig.subplots_adjust(left=0.07, right=0.9, bottom=0.12, top=0.9, wspace=0., hspace=0.)
        outer_gs = gridspec.GridSpec(1, 2, width_ratios=[1.2, 1], figure=fig, wspace=0.3, hspace=0.)

        # --------setting up the sub-figure------------
        map_gs = outer_gs[0].subgridspec(2, 2, wspace=0., hspace=0.)
        axes = [fig.add_subplot(map_gs[0, 0]), fig.add_subplot(map_gs[0, 1]), fig.add_subplot(map_gs[1, 0]), fig.add_subplot(map_gs[1, 1])]
        bpt_gs = outer_gs[1].subgridspec(1, 1, wspace=0., hspace=0.)
        axes += [fig.add_subplot(bpt_gs[0, 0])]
        
        fig.text(0.05, 0.97, f'{bin_text}', fontsize=args.fontsize, c='k', ha='left', va='top')
        
        # -----------------plot emission line maps of this bin---------------
        for index, this_line in enumerate(line_list):
            this_map, _ , _= get_emission_line_map(this_line, line_dict, args)
            nobj = line_dict[f'{this_line}_nobj']
            if type(unp.nominal_values(this_map.data)) == np.ndarray:
                axes[index] = plot_2D_map(this_map, axes[index], f'{this_line}: {nobj}', args, cmap='cividis', clabel='', takelog=takelog, vmin=cmin, vmax=cmax, hide_xaxis=index < 2, hide_yaxis=index%2 == 1, hide_cbar=True, hide_cbar_ticks=False, cticks_integer=True)

        # -----------------plot BPT diagram of this bin---------------
        ax = axes[-1]
        ax_inset = ax.inset_axes([0.65, 0.50, 0.3, 0.3])
        ax, scatter_plot_handle, ax_inset = plot_AGN_demarcation_object(line_dict, args, ax, marker='o', size=30, ax_inset=ax_inset, hide_cbar=False)
        
    return fig

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.re_limit is None: args.re_limit = 2.
    args.fontfactor = 1.5

    # ------------reading and binning dataframe-------------
    df, bin_list, args = get_binned_df(args, z_lim=z_lim)
    
    # ------------setting up master dataframe----------------------
    common_cols = ['nobj', 'logOH_int', 'logOH_int_u', 'minor_logOH_grad', 'minor_logOH_grad_u', 'major_logOH_grad', 'major_logOH_grad_u', 'radial_logOH_grad', 'radial_logOH_grad_u']
    if args.bin_by_distance_mass:
        bin_cols = ['delta_sfms_bin', 'log_mass_bin']
    elif args.bin_by_distance:
        bin_cols = ['delta_sfms_bin']
    else:
        bin_cols = ['log_mass_bin', 'log_sfr_bin']
    df_grad = pd.DataFrame(columns=common_cols + bin_cols)
    nbin_good = 0
    
    # ------------looping over each bin-----------------------
    for index2, this_bin in enumerate(bin_list):
        if args.debug_bin and nbin_good > 0: continue
        start_time3 = datetime.now()
        if args.bin_by_distance:
            bin_text = f'delta_sfms_bin_{this_bin.left}-{this_bin.right}'
        elif args.bin_by_distance_mass:
            this_delta_sfms_bin = this_bin[0]
            this_mass_bin = this_bin[1]
            bin_text = f'delta_sfms_bin_{this_delta_sfms_bin.left}-{this_delta_sfms_bin.right}_mass_bin_{this_mass_bin.left}-{this_mass_bin.right}'
        else:
            bin_text = f'logmassbin_{this_bin[0].left}-{this_bin[0].right}_logsfrbin_{this_bin[1].left}-{this_bin[1].right}'
        print(f'\n\tStarting ({index2 + 1}/{len(bin_list)}) {bin_text}..', end=' ')
        
        # -------reading previously saved stacked fits file------------
        stacked_filename = args.fits_dir / f'stacked_maps{args.deproject_text}{args.rescale_text}_{bin_text}.fits'
        if not stacked_filename.exists():
            print(f'No stacked fits file found for {bin_text}, so skipping this bin.')
            continue
        line_dict = read_stacked_maps(stacked_filename, args)
        if not line_dict:
            print(f'No lines found for {bin_text}. So Skipping.')
            continue
        nbin_good += 1

        # ---------fold stacked maps along major and minor axis--------------------
        if args.fold_maps:
            line_list = [item for item in list(line_dict.keys()) if '_nobj' not in item and '_id' not in item]
            print(f'Folding {len(line_list)} stacked emission maps for this mass-sfr bin..')
            for line in line_list:
                line_dict[line] = fold_line_map(line_dict[line], args, line=line)

        # ---------------plot emission line maps of this bin---------------------
        if args.plot_line_maps:
            fig_em, line_list = plot_stacked_line_maps(line_dict, args, bin_text=bin_text, takelog=True, cmin=-3, cmax=-2)
            save_fig(fig_em, args.fig_dir, f'stacked{args.fold_text}_line_maps{args.deproject_text}{args.rescale_text}_{bin_text}.png', args) # saving the figure

        # -----------------computing metallicity maps of this bin---------------
        if args.plot_metallicity or args.plot_line_and_metallicity:
            metallicity_map_fits_file = args.fits_dir / f'stacked{args.fold_text}_metallicity_map_{args.Zdiag}{args.C25_text}{args.deproject_text}{args.rescale_text}_{bin_text}.fits'
            if not os.path.exists(metallicity_map_fits_file) or args.clobber:
                logOH_map, logOH_int, nobj = get_metallicity_map(line_dict, args)
                if logOH_map is None:
                    print(f'Unable to compute {args.Zdiag} metallicity for {bin_text}. So Skipping.')
                    continue
                else:
                    write_metallicity_map(logOH_map, logOH_int, metallicity_map_fits_file, args) # saving the metallicity maps as fits files
            else:
                logOH_map, logOH_int, nobj = read_metallicity_map(metallicity_map_fits_file)
        
        # -----------------plot metallicity maps of this bin---------------
        if args.plot_metallicity:
            fig_met, minor_linefit_odr, major_linefit_odr, radial_linefit_odr = plot_metallicity_map(logOH_map, args, bin_text=bin_text, Zmin=None, Zmax=None)
            save_fig(fig_met, args.fig_dir, f'stacked{args.fold_text}_metallicity_map_{args.Zdiag}{args.C25_text}{args.deproject_text}{args.rescale_text}_{bin_text}.png', args) # saving the figure

        # -----------------plot line maps and metallicity maps of this bin---------------
        if args.plot_line_and_metallicity:
            fig_met, line_list, minor_linefit_odr, major_linefit_odr, radial_linefit_odr = plot_line_and_metallicity_maps(line_dict, logOH_map, args, bin_text=bin_text, takelog=True, cmin=-3, cmax=-2, Zmin=None, Zmax=None)
            save_fig(fig_met, args.fig_dir, f'stacked{args.fold_text}_line_and_metallicity_map_{args.Zdiag}{args.C25_text}{args.deproject_text}{args.rescale_text}_{bin_text}.png', args) # saving the figure

        # -------------save fit results to dataframe-----------------------
        if args.plot_radial_profiles and (args.plot_metallicity or args.plot_line_and_metallicity):
            if args.bin_by_distance_mass:
                thisrow = {'delta_sfms_bin':this_delta_sfms_bin, 'log_mass_bin':this_mass_bin}
            elif args.bin_by_distance:
                thisrow = {'delta_sfms_bin':this_bin} 
            else:
                thisrow = {'log_mass_bin':this_bin[0], 'log_sfr_bin':this_bin[1]}
            
            thisrow.update({'nobj': nobj, \
                            'logOH_int': logOH_int.n, \
                            'logOH_int_u': logOH_int.s, \
                            'minor_logOH_grad': minor_linefit_odr[0].n, \
                            'minor_logOH_grad_u': minor_linefit_odr[0].s, \
                            'major_logOH_grad': major_linefit_odr[0].n, \
                            'major_logOH_grad_u': major_linefit_odr[0].s, \
                            'radial_logOH_grad': radial_linefit_odr[0].n, \
                            'radial_logOH_grad_u': radial_linefit_odr[0].s, \
                            })
            df_grad.loc[len(df_grad)] = thisrow

        # -----------------plot line maps and metallicity maps of this bin---------------
        if args.plot_line_and_bpt:
            fig_bpt = plot_line_and_bpt_maps(line_dict, args, bin_text=bin_text, takelog=True, cmin=-3.5, cmax=-2.5)
            if fig_bpt is not None:
                save_fig(fig_bpt, args.fig_dir, f'stacked{args.fold_text}_line_and_bpt_map_{args.AGN_diag}{args.deproject_text}{args.rescale_text}_{bin_text}.png', args) # saving the figure

        if len(bin_list) > 5: plt.close('all')
        print(f'\nCompleted bin {bin_text} in {timedelta(seconds=(datetime.now() - start_time3).seconds)}, {len(bin_list) - index2 - 1} to go!')
    
    # --------------------save master dataframe-----------------------------------
    if nbin_good > 1 and args.plot_radial_profiles:
        for thiscol in ['log_mass_bin', 'log_sfr_bin', 'delta_sfms_bin']:
            if thiscol in df_grad: df_grad[thiscol] = df_grad[thiscol].astype(str) # otherwise FITS cannot save 'Interval' datatype

        Table.from_pandas(df_grad).write(args.grad_filename, format='fits', overwrite=True)
        print(f'Saved fitting data from all bins in {args.grad_filename}')
        

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

