'''
    Filename: plot_stacked_dig.py
    Notes: Reads stacked (in bins of mass and/or SFR) 2D emission line maps for objects in a given field/s, computes resolve DIG and SFR maps, and plots them
    Author : Ayan
    Created: 18-01-26
    Example: run plot_stacked_dig.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par28
             run plot_stacked_dig.py --system ssd --do_all_fields --plot_radial_profiles --cut_z_flag 4 --adaptive_bins --bin_by_distance_mass --fold_maps --skip_deproject --plot_quant cdig,sfr,logOH
'''

from header import *
from util import *
setup_plot_style()
from make_sfms_bins import get_binned_df, required_lines, sfms
from stack_emission_maps import read_stacked_maps
from plot_stacked_maps import get_metallicity_map, write_metallicity_map, read_metallicity_map, get_emission_line_map, plot_2D_map, fold_line_map
from make_diagnostic_maps import compute_SFR, take_safe_log_ratio

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def invert_Cataldi25_N2(logOH):
    '''
    Inverts the Cataldi et al. (2025) calibration Z(N2) = c0 + c1*N2 + c2*N2^2 + ...
    to find N2 = log10([NII]6583 / Halpha) for a given metallicity value.
    '''
    if np.isnan(logOH.n) or logOH <= 0:
        return ufloat(np.nan, np.nan)

    coeff_N2 = [-5.2414, 0.4951] # from Table 4 second-to-last row of C25
    N2 = np.poly1d(coeff_N2[::-1])(logOH)

    return N2

# --------------------------------------------------------------------------------------------------------------------
def deblend_nii_ha_maps(nii_plus_ha_map, logOH_map_R3, args):
    '''
    De-blends an unresolved Halpha + [NII] map into uncontaminated Halpha 
    and [NII]6583 maps.
    
    Parameters:
    -----------
    nii_plus_ha_map : numpy.ndarray
        Unresolved line flux map (Halpha + [NII]6548 + [NII]6583).
    logOH_map_R3 : numpy.ndarray
        Metallicity map (12+log(O/H)) derived from R3 ratio.
    args : argparse.Namespace
        Command-line arguments.
        
    Returns:
    --------
    ha_map : numpy.ndarray
        Uncontaminated H-alpha map.
    nii_6583_map : numpy.ndarray
        Uncontaminated [NII] 6583 map.
    '''    
    print(f'\n\tDeblending NII and Halpha..')

    # Invert Cataldi+25 relation pixel-by-pixel to retrieve N2 index
    log_nii_over_ha_map = np.full_like(logOH_map_R3, np.nan)
    flat_logOH = logOH_map_R3.data.flatten()

    flat_n2 = np.array([invert_Cataldi25_N2(logOH) for logOH in flat_logOH])
    log_nii_over_ha_map = flat_n2.reshape(logOH_map_R3.shape)
    
    # De-blend using theoretical line ratios
    nii6583_over_ha_map = 10**log_nii_over_ha_map
    nii_doublet_factor = 1.337  # 1 + 1/3 for [NII]6548 + [NII]6583
    
    ha_map = nii_plus_ha_map / (1.0 + nii_doublet_factor * nii6583_over_ha_map)
    nii_6583_map = ha_map * nii6583_over_ha_map

    ha_map = np.ma.masked_where(logOH_map_R3.mask, ha_map)
    nii_6583_map = np.ma.masked_where(logOH_map_R3.mask, nii_6583_map)

    return nii_6583_map, ha_map

# --------------------------------------------------------------------------------------------------------------------
def plot_stacked_line_maps(line_dict, args, cmin=None, cmax=None, takelog=True, axes=None):
    '''
    Makes a nice plot of all emission lines present in a given list of stacked line maps on a given axis array
    Returns axes handle
    '''
    print(f'\n\tPlotting stacked line maps..')

    line_list = [item for item in list(line_dict.keys()) if '_nobj' not in item and '_id' not in item]
    n_lines = len(line_list)
    ncols = n_lines
    nrows = 1

    if axes is None:    
        fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3.5 * nrows))
        axes = np.atleast_2d(axes)
        fig.subplots_adjust(left=0.07, right=0.96, top=0.92, bottom=0.1, wspace=0., hspace=0.2)
        fig.text(0.05, 0.95, f'{bin_text}', fontsize=args.fontsize, c='k', ha='left', va='top')

    # -----------------plot emission line maps of this bin---------------
    for index, this_line in enumerate(line_list):
        this_map, _ , _= get_emission_line_map(this_line, line_dict, args)
        nobj = line_dict[f'{this_line}_nobj']
        if type(unp.nominal_values(this_map.data)) == np.ndarray:
            axes[index] = plot_2D_map(this_map, axes[index], f'Stacked {this_line}: {nobj}', args, cmap='cividis', clabel='', takelog=takelog, vmin=cmin, vmax=cmax, hide_xaxis=True, hide_yaxis=index, hide_cbar=index < ncols - 1, hide_cbar_ticks=False, cticks_integer=True)

    return axes        

# --------------------------------------------------------------------------------------------------------------------
def compute_dig_contribution(line_dict, logOH_map, args):
    '''
    Plots spatially resolved DIG contribution based on Tomicic+21, on an existing axis
    Returns axis handle and the handle of the spatially resolved scatter plot
    '''
    print(f'\n\tComputing DIG contribution..')

    # ---------getting the Ha and SII maps-------------
    ha_map, _ , _= get_emission_line_map('HA', line_dict, args)
    sii_map, _ , _= get_emission_line_map('SII', line_dict, args)

    # ---------getting the S2Ha ratio map-------------
    S2Ha_map = take_safe_log_ratio(sii_map, ha_map, skip_log=True)

    # -----correcting S2Ha by metallicity, following Tomicic+21----------------
    Z_MW = 8.69 # T21 (and Franchetto+20) assumed MW abundance = solar neighbourhood abundance
    Z_ratio = 10 ** (logOH_map.data - Z_MW)
    S2Ha_map_corr = np.ma.masked_where(S2Ha_map.mask | logOH_map.mask, (S2Ha_map.data / Z_ratio))

    # ------converting units of Ha map from ergs/s/cm^2 to ergs/s/kpc^2---------
    quant = (unp.nominal_values(ha_map.data) * u.erg / u.second / u.cm ** 2).to(u.erg / u.second / (u.kpc ** 2)).value
    quant_u = (unp.std_devs(ha_map.data) * u.erg / u.second / u.cm ** 2).to(u.erg / u.second / (u.kpc ** 2)).value
    ha_map = np.ma.masked_where(ha_map.mask, unp.uarray(quant, quant_u))

    # ---------plotting S2Ha_corr vs Ha SB-------------
    bad_mask = S2Ha_map_corr.mask | ha_map.mask
    S2Ha_corr_array = np.ma.compressed(np.ma.masked_where(bad_mask, S2Ha_map_corr))
    ha_array = np.ma.compressed(np.ma.masked_where(bad_mask, ha_map))
    log_ha_array = unp.log10(ha_array)

    # ---------computing spatially resolved C_DIG, following Tomicic+21 S3.1-------------
    df = pd.DataFrame({'S2Ha_corr': unp.nominal_values(S2Ha_corr_array), 'S2Ha_corr_u': unp.std_devs(S2Ha_corr_array), \
                       'log_Ha_SB':unp.nominal_values(log_ha_array), 'log_Ha_SB_u':unp.std_devs(log_ha_array), \
                       })
    df['pct_Ha'] = df['log_Ha_SB'].rank(pct=True)
    #S2Ha_corr_DIG = np.median(df[df['pct_Ha'] <= 0.05]['S2Ha_corr']) # percentile valyes from Tomicic+21
    #S2Ha_corr_dense = np.median(df[df['pct_Ha'] >= 0.95]['S2Ha_corr'])
    S2Ha_corr_DIG = np.median(df[df['pct_Ha'] <= 0.1]['S2Ha_corr']) # using different percentile values because our data did not seem to have below 5% etc
    S2Ha_corr_dense = np.median(df[df['pct_Ha'] >= 0.9]['S2Ha_corr'])
    c_dig_array = (S2Ha_corr_dense - S2Ha_corr_array) / (S2Ha_corr_dense - S2Ha_corr_DIG)

    # ---------fitting C_DIG vs Ha SB------------------
    def cdig_func(x, *popt): return np.piecewise(x, [x <= popt[0], x > popt[0]], [1, lambda x: (popt[0]/x) ** popt[1]])

    try:
        popt, pcov = curve_fit(cdig_func, unp.nominal_values(ha_array), unp.nominal_values(c_dig_array), p0=[1e24, 0.7], sigma=unp.std_devs(c_dig_array), absolute_sigma=True)

        # ---------deriving the fitted C_DIG map and uncertainty------------------
        c_dig_map = cdig_func(unp.nominal_values(ha_map), *popt)
        c_dig_map_lowlim = cdig_func(unp.nominal_values(ha_map) + unp.std_devs(ha_map), *popt)
        c_dig_map_uplim = cdig_func(unp.nominal_values(ha_map) - unp.std_devs(ha_map), *popt)
        c_dig_map_u = np.mean(np.array([c_dig_map_lowlim, c_dig_map_uplim]), axis=0)
        c_dig_map = np.ma.masked_where(S2Ha_map_corr.mask, unp.uarray(c_dig_map, c_dig_map_u))
    except ValueError:
        print(f'\tFailed to compute DIG contribution because min and max Ha SB percentiles available were {df["pct_Ha"].min()} and {df["pct_Ha"].max()}, respectively.')
        c_dig_map = np.ma.masked_where(S2Ha_map_corr.mask, unp.uarray(np.full_like(logOH_map.data, np.nan), np.full_like(logOH_map.data, np.nan)))

    return c_dig_map

# --------------------------------------------------------------------------------------------------------------------
def plot_stacked_rad_profiles(df_list, args, quant_arr=['dig']):
    '''
    Plots the stacked radial profiles of a given quantity from a list of dataframes containing radial profiles for each bin
    Returns figure handle
    '''
    xcol_arr = ['distance', 'minor_distance', 'major_distance']
    colorcol_arr = ['mass_intervals', 'bin_intervals']
    cmap_arr = ['viridis', 'RdBu']
    clim_dict = {'mass_intervals': [7, 11], 'bin_intervals': [-1, 1]}
    clabel_dict = {'mass_intervals': r'$\log{(M_*/M_{\odot})}$', 'bin_intervals': r'$\delta_{SFMS}$', 'sfr': r'SFR [M$_\odot$/yr]', 'cdig': 'DIG contribution', 'logOH': r'$\log{(O/H)}$ + 12'}

    # -------looping over quants-------------
    for quant in quant_arr:
        print(f'\nPlotting radial {quant} profiles..')

        # ----------setup figure---------
        fig = plt.figure(figsize=(8, 7.5))

        gs = fig.add_gridspec(
            nrows=len(xcol_arr) + 1, 
            ncols=len(colorcol_arr), 
            height_ratios=[0.035, 1, 1, 1],  # Small top row for colorbars
            wspace=0.02,                      # Near-zero horizontal gap between columns
            hspace=0.05                       # Minimal vertical gap between panels
        )

        # --------looping over color coding----------
        axes = np.empty((len(xcol_arr), len(colorcol_arr)), dtype=object)
        for c, colorcol in enumerate(colorcol_arr):
            cmap = mplcolormaps[cmap_arr[c]]
            cmin = clim_dict[colorcol][0]
            cmax = clim_dict[colorcol][1]
            norm = mplcolors.Normalize(vmin=cmin, vmax=cmax)

            # -------looping over dist, maj and min dist-------
            for r, xcol in enumerate(xcol_arr):
                sharex_ax = axes[0, c] if r > 0 else None
                sharey_ax = axes[0, 0] if (r > 0 or c > 0) else None
                axes[r, c] = fig.add_subplot(gs[r + 1, c], sharex=sharex_ax, sharey=sharey_ax)
                    
                # -------looping over bins---------
                for df in df_list:
                    df[colorcol] = df[colorcol].astype('str')
                    df[colorcol] = pd.IntervalIndex.from_tuples(df[colorcol].apply(lambda x: pd.to_numeric(x.strip('()[]').split(', '))).map(tuple), closed='right')

                    color_value = (df[colorcol].values[0].left + df[colorcol].values[0].right) / 2.
                    color = cmap(norm(color_value))            

                    # # ------plotting individual radial profile--------
                    # axes[r][c].scatter(df[xcol], df[f'{quant}'], c=color, s=5, lw=1, ec='k')
                    # axes[r][c].errorbar(df[xcol], df[f'{quant}'], yerr=df[f'{quant}_u'], c=color, lw=0.5)

                    # ------radially binning df---------
                    df['xcol_bin'] = pd.cut(df[xcol], bins=10)
                    df['w'] = 1.0 / (df[f'{quant}_u']**2)
                    df['wy'] = df[f'{quant}'] * df['w']

                    # Aggregate using a single groupby pass
                    binned_df = df.groupby('xcol_bin', observed=False).agg(
                        x_mid=(xcol, 'mean'),   # Bin center (or use .median())
                        sum_w=('w', 'sum'),
                        sum_wy=('wy', 'sum'),
                        n_points=(f'{quant}', 'count')
                    ).reset_index()

                    # Compute inverse-variance weighted mean and error
                    binned_df['quant_weighted'] = binned_df['sum_wy'] / binned_df['sum_w']
                    binned_df['quant_weighted_u'] = 1.0 / np.sqrt(binned_df['sum_w'])
                    
                    # ------plotting binned radial profile--------
                    line, = axes[r][c].plot(binned_df['x_mid'], binned_df['quant_weighted'], c=color, lw=2)
                    axes[r][c].errorbar(binned_df['x_mid'], binned_df['quant_weighted'], yerr=binned_df['quant_weighted_u'], c=color, lw=0.5)
                    line.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

                    #if quant == 'sfr': axes[r][c].set_yscale('log')
                    axes[r][c] = annotate_axes(axes[r][c], '', '', label= f'{xcol}' if c == 0 else '', labelx=0.6, xlim=[0, 2], ylim=None, args=args, clabel='', hide_xaxis=r < len(xcol_arr) - 1, hide_yaxis=c > 0, hide_cbar=True)

            # ---------annotating colorbar------------
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            cax = fig.add_subplot(gs[0, c])
            cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
            cax.xaxis.set_ticks_position('top')
            cax.xaxis.set_label_position('top')
            cbar.set_label(clabel_dict[colorcol], fontsize=args.fontsize, labelpad=6)
            cbar.ax.tick_params(labelsize=args.fontsize)
            cbar.locator = ticker.MaxNLocator(integer=False, nbins=4)#, prune='both')
            cbar.update_ticks()

        fig.supylabel(clabel_dict[quant], fontsize=args.fontsize, x=0.02)
        fig.supxlabel(r'Radius (R$_e$)', fontsize=args.fontsize, y=0.01)
        plt.subplots_adjust(left=0.10, right=0.98, top=0.91, bottom=0.07)

        save_fig(fig, args.fig_dir, f'stacked{args.fold_text}_{quant}_profiles{args.deproject_text}{args.rescale_text}.png', args)

    return fig

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.re_limit is None: args.re_limit = 2.
    args.fontfactor = 1.5
    args.Zdiag = 'R3'
    args.use_C25 = True
    args.dered_in_NB = True

    # ------------reading and binning dataframe-------------
    df, bin_list, args = get_binned_df(args, required_lines=required_lines, sfms=sfms)
    
    # ------------setting up master dataframe----------------------
    nbin_good = 0
    rad_profiles = []
    
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
        elif args.bin_by_sfh_mass:
            this_delta_sfms_bin = this_bin[0]
            this_mass_bin = this_bin[1]
            bin_text = f'tform_ratio_bin_{this_delta_sfms_bin.left}-{this_delta_sfms_bin.right}_mass_bin_{this_mass_bin.left}-{this_mass_bin.right}'
        else:
            bin_text = f'logmassbin_{this_bin[0].left}-{this_bin[0].right}_logsfrbin_{this_bin[1].left}-{this_bin[1].right}'
        print(f'\n\tStarting ({index2 + 1}/{len(bin_list)}) {bin_text}..', end=' ')

        # -------reading previously saved radial profile dataframes-----------
        nbin_good += 1
        stacked_dfname = args.fits_dir / f'stacked_radial_df{args.deproject_text}{args.rescale_text}_{bin_text}.csv'
        if not os.path.exists(stacked_dfname) or args.clobber:
            print(f'\n\tStacked df does not exist, making a new one..')

            # -------reading previously saved stacked fits file------------
            stacked_filename = args.fits_dir / f'stacked_maps{args.deproject_text}{args.rescale_text}_{bin_text}.fits'
            if not stacked_filename.exists():
                print(f'No stacked fits file found for {bin_text}, so skipping this bin.')
                continue
            line_dict = read_stacked_maps(stacked_filename, args)
            if not line_dict:
                print(f'No lines found for {bin_text}. So Skipping.')
                continue

            # ---------fold stacked maps along major and minor axis--------------------
            if args.fold_maps:
                line_list = [item for item in list(line_dict.keys()) if '_nobj' not in item and '_id' not in item]
                print(f'Folding {len(line_list)} stacked emission maps for this mass-sfr bin..')
                for line in line_list:
                    line_dict[line] = fold_line_map(line_dict[line], args, line=line)

            # -----------------setup the figure---------------
            line_list = [item for item in list(line_dict.keys()) if '_nobj' not in item and '_id' not in item]
            n_lines = len(line_list)
            ncols = max(n_lines, 5)
            nrows = 2
            
            fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3.5 * nrows))
            axes = np.atleast_2d(axes)
            fig.subplots_adjust(left=0.07, right=0.96, top=0.92, bottom=0.1, wspace=0., hspace=0.2)
            fig.text(0.05, 0.95, f'{bin_text}', fontsize=args.fontsize, c='k', ha='left', va='top')

            #----------------plot line maps of this bin---------------------
            axes0 = plot_stacked_line_maps(line_dict, args, cmin=-3, cmax=-2, takelog=True, axes=axes[0])
            fig.delaxes(axes[0][4])

            # ---------------deriving and plotting R3 metallicity---------------------
            metallicity_map_fits_file = args.fits_dir / f'stacked{args.fold_text}_metallicity_map_R3_wC25{args.deproject_text}{args.rescale_text}_{bin_text}.fits'
            if not os.path.exists(metallicity_map_fits_file) or args.clobber:
                logOH_map, logOH_int, nobj = get_metallicity_map(line_dict, args)
                if logOH_map is None:
                    print(f'Unable to compute R3 metallicity for {bin_text}. So Skipping.')
                    continue
                else:
                    write_metallicity_map(logOH_map, logOH_int, metallicity_map_fits_file, args, nobj=nobj) # saving the metallicity maps as fits files
            
            logOH_map_R3, logOH_int_R3, nobj = read_metallicity_map(metallicity_map_fits_file, args)
            Zmax, Zmin = None, None
            axes[1][0] = plot_2D_map(logOH_map_R3, axes[1][0], rf'$\log$(O/H) + 12: {args.Zdiag}', args, cmap='plasma', takelog=False, vmin=Zmin, vmax=Zmax, hide_cbar=True, cticks_integer=True)

            # ---------------correct for NII/Ha based on R3 metallicity---------------------
            corrected_ha_map = line_dict['HA'].data
            factor = 0.823 # from James et al. 2023
            nii_plus_ha_map = corrected_ha_map / factor # undoing the correction applied in stack_emission_maps.py

            nii_map, ha_map = deblend_nii_ha_maps(nii_plus_ha_map, logOH_map_R3, args)
            line_dict['NII'] = np.ma.masked_array(nii_map, mask=line_dict['HA'].mask)
            line_dict['HA'] = np.ma.masked_array(ha_map, mask=line_dict['HA'].mask)

            # ---------------plotting derived NII6548 and Halpha maps---------------------
            axes[1][1] = plot_2D_map(line_dict['NII'], axes[1][1], f'NII6548', args, cmap='cividis', clabel='', takelog=True, vmin=-3, vmax=-2, hide_xaxis=False, hide_yaxis=True, hide_cbar=True)
            axes[1][2] = plot_2D_map(line_dict['HA'], axes[1][2], f'Halpha', args, cmap='cividis', clabel='', takelog=True, vmin=-3, vmax=-2, hide_xaxis=False, hide_yaxis=True, hide_cbar=True)

            # ---------------making radial profiles---------------------
            shape = np.shape(logOH_map_R3)
            if args.fold_maps:
                center_xpix, center_ypix = (args.npix_side % 2 == 0) * 0.5, (args.npix_side % 2 == 0) * 0.5 # this yields 0.5 (instead of 0) pixel offset for even-sized stacked maps, because the center of the map is in the center (and not the edge) of the first pixel
            else:
                center_xpix, center_ypix = shape[0] / 2., shape[1] / 2.
            distance_map = np.array([[np.sqrt((i - center_xpix)**2 + (j - center_ypix)**2) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size # Re or kpc
            minor_distance_map = np.array([[np.abs(j - center_ypix) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size # Re or kpc
            major_distance_map = np.array([[np.abs(i - center_xpix) for j in range(shape[1])] for i in range(shape[0])]) * args.pix_size # Re or kpc
            this_df = pd.DataFrame({'distance': distance_map.flatten(), 'minor_distance': minor_distance_map.flatten(), 'major_distance': major_distance_map.flatten()})
            this_df['bin_intervals'] = this_bin[0]
            this_df['mass_intervals'] = this_bin[1]

            # ---------------plot SFR maps of this bin---------------------
            df_thisbin = df[(df['bin_intervals'] == this_bin[0]) & (df['mass_intervals'] == this_bin[1])]
            mean_redshift = df_thisbin.redshift.mean()
            distance = cosmo.luminosity_distance(mean_redshift)
            sfr_map = compute_SFR(line_dict['HA'], distance)

            axes[1][3] = plot_2D_map(sfr_map, axes[1][3], f'SFR', args, cmap='Blues', takelog=False, vmin=None, vmax=None, hide_cbar=False, hide_yaxis=True)

            this_df[f'sfr'] = unp.nominal_values(sfr_map.data).flatten()
            this_df[f'sfr_u'] = unp.std_devs(sfr_map.data).flatten()

            # ---------------plot DIG maps of this bin---------------------
            dig_map = compute_dig_contribution(line_dict, logOH_map_R3, args)
            dig_cmap = get_combined_cmap([0., 0.3, 0.7, 1.], ['binary', 'spring', 'hot']) # making C_DIG colormap, broken at 0.3 and 0.7
            axes[1][4] = plot_2D_map(dig_map, axes[1][4], r'C$_{\mathrm{DIG}}$', args, cmap=dig_cmap, takelog=False, vmin=0, vmax=1, hide_cbar=False, hide_yaxis=True)
            save_fig(fig, args.fig_dir, f'stacked{args.fold_text}_sfr_dig_maps{args.deproject_text}{args.rescale_text}_{bin_text}.png', args) # saving the figure

            this_df[f'cdig'] = unp.nominal_values(dig_map.data).flatten()
            this_df[f'cdig_u'] = unp.std_devs(dig_map.data).flatten()

            # ---------------add NB metallicity radial profile of this bin---------------------
            metallicity_map_fits_file = args.fits_dir / f'stacked{args.fold_text}_metallicity_map_NB{args.deproject_text}{args.rescale_text}_{bin_text}.fits' 
            logOH_map_NB, _, _ = read_metallicity_map(metallicity_map_fits_file, args)

            this_df[f'logOH'] = unp.nominal_values(logOH_map_NB.data).flatten()
            this_df[f'logOH_u'] = unp.std_devs(logOH_map_NB.data).flatten()

            # ----------writing df to file-------------
            this_df.to_csv(stacked_dfname, index=None)
        else:
            print(f'\n\tReading from existing dataframe {stacked_dfname}..')
            this_df = pd.read_csv(stacked_dfname, comment='#')

        rad_profiles.append(this_df)

        if len(bin_list) > 5: plt.close('all')
        print(f'\nCompleted bin {bin_text} in {timedelta(seconds=(datetime.now() - start_time3).seconds)}, {len(bin_list) - index2 - 1} to go!')
    
    # --------------------plotting all profiles in one plot-----------------------------------
    if nbin_good > 1:
        fig_radprof = plot_stacked_rad_profiles(rad_profiles, args, quant_arr=args.plot_quant.split(','))
        
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
