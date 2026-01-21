'''
    Filename: plot_stacked_maps.py
    Notes: Reads stacked (in bins of mass and/or SFR) 2D emission line maps for objects in a given field/s, computes metallicity gradient, and plots them
    Author : Ayan
    Created: 18-01-26
    Example: run plot_stacked_maps.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par28
             run plot_stacked_maps.py --field Par28 --Zdiag R23 --use_C25 --plot_line_maps --plot_metallicity --plot_radial_profiles --debug_bin
             run plot_stacked_maps.py --field Par28 --Zdiag R23 --use_C25 --fold_maps --plot_metallicity --debug_bin
             run plot_stacked_maps.py --field Par28 --Zdiag R23 --use_C25 --fold_maps --plot_radial_profiles --debug_bin
             run plot_stacked_maps.py --field Par28 --Zdiag R23 --use_C25 --adaptive_bins --fold_maps --plot_radial_profiles --debug_bin
             run plot_stacked_maps.py --field Par28 --Zdiag R23 --use_C25 --plot_radial_profiles
'''

from header import *
from util import *
from make_diagnostic_maps import compute_Z_C19, compute_Z_KD02_R23, compute_Z_P25, compute_Z_Te, compute_Te, take_safe_log_ratio, take_safe_log_sum, myimshow
from stack_emission_maps import read_stacked_maps
from plots_for_zgrad_paper import plot_2D_map, plot_radial_profile

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def save_fig(fig, fig_dir, figname, args):
    '''
    Saves a given figure handle as a given output filename
    '''

    if args.fortalk:
        #mplcyberpunk.add_glow_effects()
        #try: mplcyberpunk.make_lines_glow()
        #except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    fig_dir.mkdir(exist_ok=True, parents=True)
    figname = fig_dir / figname
    fig.savefig(figname, transparent=args.fortalk)
    print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    return

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
        ncols = min(n_lines, max_ncols)
        nrows = int(np.ceil(n_lines / ncols))
        
        fig, axes2d = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows), sharex=True, sharey=True)
        fig.subplots_adjust(left=0.1, right=0.93, top=0.9, bottom=0.1, wspace=0., hspace=0.)
        fig.text(0.05, 0.95, f'{bin_text}', fontsize=args.fontsize, c='k', ha='left', va='top')
        
        if nrows == 1: axes2d = axes2d[np.newaxis, :]
        if ncols == 1: axes2d = axes2d[:, np.newaxis]
        axes = np.atleast_1d(axes2d).flatten()

        # -----------------plot emission line maps of this bin---------------
        for index3, this_line in enumerate(line_list):
            row = index3 // ncols
            col = index3 % ncols
            ax = axes2d[row, col]

            this_map, _ , _= get_emission_line_map(this_line, line_dict, args)
            nobj = line_dict[f'{this_line}_nobj']
            if type(unp.nominal_values(this_map.data)) == np.ndarray:
                ax = plot_2D_map(this_map, ax, f'Stacked {this_line}: {nobj}', args, cmap='cividis', clabel='', takelog=takelog, vmin=cmin, vmax=cmax, hide_xaxis=row < nrows - 1, hide_yaxis=col, hide_cbar=col < ncols - 1, skip_annotate=False, hide_cbar_ticks=False, cticks_integer=True)

        for j in range(index3 + 1, len(axes)): axes[j].axis('off') # hide unused subplots

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
        figname = fig_dir / f'stacked_metallicity_debug_Zdiag_{args.Zdiag}{Zbranch_text}.png'
        fig.savefig(figname, transparent=args.fortalk, dpi=200)
        print(f'\nSaved figure at {figname}')

    return logOH_map, logOH_int, nobj

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
def write_metallicity_map(logOH_map, logOH_int, outfilename, args):
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

    logOH_val_hdu = fits.ImageHDU(data=logOH_map_val, name='log_OH', header=hdr2)
    logOH_err_hdu = fits.ImageHDU(data=logOH_map_err, name='log_OH_u')
    hdul = fits.HDUList([primary_hdu, logOH_val_hdu, logOH_err_hdu])
    hdul.writeto(outfilename, overwrite=True)

    print(f'Saved metallicity maps in {outfilename}')

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_map(logOH_map, args, bin_text='', Zmin=None, Zmax=None):
    '''
    Makes a nice plot of a given metallicity map
    Returns figure handle
    '''
    # -----------------setup the figure---------------
    fig, axes = plt.subplots(1, 2 if args.plot_radial_profiles else 1, figsize=(10, 4) if args.plot_radial_profiles else (8, 6))
    fig.subplots_adjust(left=0.12, right=0.95 if args.plot_radial_profiles else 0.85, top=0.9, bottom=0.13, wspace=0.6, hspace=0.)
    axes = np.atleast_1d(axes)
    
    # -----------------plot metallicity maps of this bin---------------
    axes[0] = plot_2D_map(logOH_map, axes[0], f'Stacked logOH {bin_text}', args, cmap='plasma', clabel=r'$\log$(O/H) + 12', takelog=False, vmin=Zmin, vmax=Zmax, hide_cbar=False, cticks_integer=True)

    return fig         

# --------------------------------------------------------------------------------------------------------------------
def get_interval_from_filename(filename):
    '''
    Accepts a filename of the pattern stacked_maps_logmassbin_X-Y_logsfrbin_A-B.fits and extracts the mass and SFR interval from the filename
    Returns object: pair of Intervals
    '''
    filename = Path(filename).stem
    pattern = r'logmassbin_([-?\d.]+)-([-?\d.]+)_logsfrbin_([-?\d.]+)-([-?\d.]+)'
    match = re.search(pattern, filename)

    if match:
        m_low, m_high, s_low, s_high = map(float, match.groups())        
        mass_int = pd.Interval(m_low, m_high, closed='left')
        sfr_int = pd.Interval(s_low, s_high, closed='left')
        
        combined_interval = (mass_int, sfr_int)
    else:
        raise ValueError(f'This file {filename} does not have the expected pattern. So could not extract intervals from it.')
    
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
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.re_limit is None: args.re_limit = 2.
    if args.plot_radial_profiles: args.plot_metallicity = True

    # -------------for figure annotations--------------------
    max_ncols = 4
    args.fontfactor = 1.5
    args.pix_size_re = (2 * args.re_limit) / args.npix_side
    offset = args.pix_size_re / 2 # half a pixel offset to make sure cells in 2D plot are aligned with centers and not edges
    args.extent = (-args.re_limit - offset, args.re_limit - offset, -args.re_limit - offset, args.re_limit - offset)
    fold_text = '_folded' if args.fold_maps else ''

    # ---------determining list of fields----------------
    if args.do_all_fields:
        field_list = [os.path.split(item[:-1])[1] for item in glob.glob(str(args.input_dir / 'Par*') + '/')]
        field_list.sort(key=natural_keys)
    else:
        field_list = args.field_arr
    
    # --------loop over all fields------------------
    for index, field in enumerate(field_list):
        start_time2 = datetime.now()
        args.field = f'Par{int(field[3:]):03}' if len(field) < 6 else field
        print(f'\nStarting field {args.field} which is {index + 1} of {len(field_list)}..')

        # ------determining field-specific paths, etc-----------
        product_dir = args.input_dir / args.field / 'Products'
        output_dir = args.output_dir / args.field / 'stacking'
        if args.adaptive_bins: output_dir = Path(str(output_dir).replace('stacking', 'stacking_adaptive'))
        output_dir.mkdir(parents=True, exist_ok=True)
        fig_dir = output_dir / 'plots'
        fig_dir.mkdir(parents=True, exist_ok=True)
        fits_dir = output_dir / 'maps'
        fits_dir.mkdir(parents=True, exist_ok=True)

        C25_text = '_wC25' if args.use_C25 and 'NB' not in args.Zdiag else ''
        output_filename = fits_dir / f'stacked{fold_text}_fits_allbins_Zdiag_{args.Zdiag}{C25_text}.fits'

        # --------getting the list of bins from available fits files---------------
        stacked_files = glob.glob(str(fits_dir) + '/stacked_maps_*.fits')
        bin_list = [get_interval_from_filename(item) for item in stacked_files]
        bin_list.sort(key=lambda x: (x[0].left, x[1].left))

        # ------------setting up master dataframe-----------------------
        df_grad = pd.DataFrame(columns=['log_mass_bin', 'log_sfr_bin', 'nobj', 'logOH_int', 'logOH_int_u', 'logOH_grad', 'logOH_grad_u'])
        nbin_good = 0
        
        # ------------looping over each bin-----------------------
        for index2, this_mass_sfr_bin in enumerate(bin_list):
            if args.debug_bin and nbin_good > 0: continue
            start_time3 = datetime.now()
            bin_text = f'logmassbin_{this_mass_sfr_bin[0].left}-{this_mass_sfr_bin[0].right}_logsfrbin_{this_mass_sfr_bin[1].left}-{this_mass_sfr_bin[1].right}'
            print(f'\tStarting ({index2 + 1}/{len(bin_list)}) {bin_text}..', end=' ')
            
            # -------reading previously saved stacked fits file------------
            stacked_filename = fits_dir / f'stacked_maps_{bin_text}.fits'
            if not stacked_filename.exists():
                print(f'No stacked fits file found for {bin_text}, so skipping this bin.')
                continue
            line_dict = read_stacked_maps(stacked_filename, args)
            nbin_good += 1

            # ---------fold stacked maps along major and minor axis--------------------
            if args.fold_maps:
                line_list = [item for item in list(line_dict.keys()) if '_nobj' not in item and '_id' not in item]
                print(f'Folding {len(line_list)} stacked emission maps for this mass-sfr bin..')
                for line in line_list:
                    line_dict[line] = fold_line_map(line_dict[line], args, line=line)

            # ---------------plot emission line maps of this bin---------------------
            if args.plot_line_maps:
                fig_em, line_list = plot_stacked_line_maps(line_dict, args, bin_text=bin_text, takelog=True, cmin=-20.5, cmax=-17.8)
                if len(line_list) == 0:
                    print(f'No lines found for {bin_text}. So Skipping.')
                    continue
                save_fig(fig_em, fig_dir, f'stacked{fold_text}_line_maps_{bin_text}.png', args) # saving the figure

            # -----------------computing metallicity maps of this bin---------------
            logOH_map, logOH_int, nobj = get_metallicity_map(line_dict, args)
            if logOH_map is None:
                print(f'Unable to compute {args.Zdiag} metallicity for {bin_text}. So Skipping.')
                continue
            else:
                write_metallicity_map(logOH_map, logOH_int, fits_dir / f'stacked{fold_text}_metallicity_map_{bin_text}.fits', args) # saving the metallicity maps as fits files
            
            # -----------------plot metallicity maps of this bin---------------
            if args.plot_metallicity:
                fig_met = plot_metallicity_map(logOH_map, args, bin_text=bin_text, Zmin=None, Zmax=None)
                save_fig(fig_met, fig_dir, f'stacked{fold_text}_metallicity_map_{bin_text}.png', args) # saving the figure

            # -----------------computing metallicity gradient of this bin---------------
            if args.plot_radial_profiles:
                ax = fig_met.axes[1]
                quant = 'log_OH'
                shape = np.shape(logOH_map)

                if args.fold_maps:
                    center_xpix, center_ypix = (args.npix_side % 2 == 0) * 0.5, (args.npix_side % 2 == 0) * 0.5 # this yields 0.5 (instead of 0) pixel offset for even-sized stacked maps, because the center of the map is in the center (and not the edge) of the first pixel
                else:
                    center_xpix, center_ypix = shape[0] / 2., shape[1] / 2.
                distance_map = np.array([[np.sqrt((i - center_xpix)**2 + (j - center_ypix)**2) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size_re # Re
                
                logOH_df = pd.DataFrame({'distance': distance_map.flatten(), f'{quant}': unp.nominal_values(logOH_map).flatten(), f'{quant}_u': unp.std_devs(logOH_map).flatten()})
                logOH_df = logOH_df.dropna(subset=[f'{quant}', f'{quant}_u'], axis=0)

                ax, [linefit_original, linefit_odr, linefit_lenstronomy, params_llim, params_median, params_ulim] = plot_radial_profile(logOH_df, ax, args, ylim=None, xlim=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=True, skip_annotate=False, short_label=False, quant=quant, Zdiag=args.Zdiag, do_mcmc=False, already_in_re=True)
                save_fig(fig_met, fig_dir, f'stacked{fold_text}_metallicity_map_{bin_text}.png', args) # saving the figure

                # -------------save fit results to dataframe-----------------------
                thisrow = [this_mass_sfr_bin[0], this_mass_sfr_bin[1], nobj, logOH_int.n, logOH_int.s, linefit_odr[0].n, linefit_odr[0].s]
                df_grad.loc[len(df_grad)] = thisrow

            print(f'\nCompleted bin mass={this_mass_sfr_bin[0]}, sfr={this_mass_sfr_bin[1]} in {timedelta(seconds=(datetime.now() - start_time3).seconds)}, {len(bin_list) - index2 - 1} to go!')
        
        # --------------------save master dataframe-----------------------------------
        if nbin_good > 1 and args.plot_radial_profiles:
            for thiscol in ['log_mass_bin', 'log_sfr_bin']: df_grad[thiscol] = df_grad[thiscol].astype(str) # otherwise FITS cannot save 'Interval' datatype

            Table.from_pandas(df_grad).write(output_filename, format='fits', overwrite=True)
            print(f'Saved fitting data from all bins in {output_filename}')
        
        print(f'Completed field {field} ({nbin_good} / {len(bin_list)} bins) in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(field_list) - index - 1} to go!')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

