'''
    Filename: make_diagnostic_maps.py
    Notes: Plots integrated 1D spectra and 2D emission line maps (from existing .full.fits file), for a given object/s in a given field, while accounting for ALL uncertainties via error propagation
    Author : Ayan
    Created: 17-07-24
    Example: run make_diagnostic_maps.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par50 --id 3667
             run make_diagnostic_maps.py --field Par50 --id 823 --plot_radial_profiles
             run make_diagnostic_maps.py --field Par50 --id 823 --pixscale 0.2 --vorbin --voronoi_line Ha --voronoi_snr 3 --plot_radial_profiles
             run make_diagnostic_maps.py --field Par50 --id 823 --plot_radial_profiles --only_seg --snr_cut 3 --plot_mappings
             run make_diagnostic_maps.py --field Par51 --do_all_obj --plot_radial_profiles --only_seg --snr_cut 3 --write_file
             run make_diagnostic_maps.py --field Par51 --re_extract --do_all_obj --plot_radial_profiles --only_seg --snr_cut 3 --write_file
             run make_diagnostic_maps.py --field Par28 --id 58,1457,1585,1588 --plot_radial_profiles --only_seg --plot_mappings
             run make_diagnostic_maps.py --field Par28 --id 58,1457,1585,1588 --plot_radial_profiles --only_seg --plot_mappings --vorbin --voronoi_snr 5
             run make_diagnostic_maps.py --field Par28 --id 58,1457,1585,1588 --plot_starburst --vorbin --voronoi_snr 5 --plot_radial_profile --only_seg
             run make_diagnostic_maps.py --field Par28 --id 58,1457,1588 --plot_metallicity --vorbin --voronoi_snr 5 --plot_radial_profile --only_seg
    Afterwards, to make the animation: run /Users/acharyya/Work/astro/ayan_codes/animate_png.py --inpath /Volumes/Elements/acharyya_backup/Work/astro/passage/passage_output/Par028/all_diag_plots_wradprof_snr3.0_onlyseg/ --rootname Par028_*_all_diag_plots_wradprof_snr3.0_onlyseg.png --delay 0.1
'''

from header import *
from util import *
from matplotlib import cm as mpl_cm
import imageio

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def annotate_PAs(pa_arr, ax, fontsize=10):
    '''
    Annotates a given plot with the PAs of all available filters in the given axis
    Returns the axis handle
    '''
    color = 'grey'
    x_cen = ax.get_xlim()[1] - 0.1 * np.diff(ax.get_xlim())[0]
    y_cen = ax.get_ylim()[1] - 0.15 * np.diff(ax.get_ylim())[0]
    len = np.diff(ax.get_xlim())[0] * 0.1

    for pa in pa_arr:
        x_comp = len * np.sin(pa * np.pi / 180)
        y_comp = len * np.cos(pa * np.pi / 180)
        ax.plot([x_cen, x_cen - x_comp], [y_cen, y_cen + y_comp], lw=1, c=color)
        ax.text(x_cen - x_comp - 0.02, y_cen + y_comp + 0.02, r'%d$^\circ$' % pa, fontsize=fontsize, ha='center', va='center', rotation=pa)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_direct_image(full_hdu, ax, args, hide_xaxis=False, hide_yaxis=False):
    '''
    Plots the combined direct image of all available filters in the given axis
    Returns the axis handle
    '''
    cmap_arr = ['Reds', 'Blues', 'Greens']

    # -------plot direct image for each filter---------
    for index in range(args.ndfilt):
        filt = full_hdu[0].header[f'DFILT{index + 1:02d}']
        print(f'Plotting direct image for filter {filt} which is {index+1} of {args.ndfilt}..')

        ext = 5 + index * 2
        image = full_hdu[ext].data
        image = trim_image(image, args)

        p = ax.imshow(image, cmap=cmap_arr[index], origin='lower', extent=args.extent, alpha=1)#, vmin=0, vmax=0.03)

        ax.set_xlim(-args.arcsec_limit, args.arcsec_limit)  # arcsec
        ax.set_ylim(-args.arcsec_limit, args.arcsec_limit)  # arcsec

        textcolor = mpl_cm.get_cmap(cmap_arr[index])(0.9)
        ax.text(ax.get_xlim()[0] * 0.9, ax.get_ylim()[1] * 0.7 - index * 0.1, filt, c=textcolor, fontsize=args.fontsize / 1.5, ha='left', va='top')

    ax.text(ax.get_xlim()[0] * 0.9, ax.get_ylim()[1] * 0.95, f'z={args.z:.2f}', c='k', fontsize=args.fontsize, ha='left', va='top')
    ax.text(ax.get_xlim()[1] * 0.95, ax.get_ylim()[0] * 0.95, f'Mag={args.mag:.1f}', c='k', fontsize=args.fontsize, ha='right', va='bottom')
    ax.scatter(0, 0, marker='x', s=10, c='grey')
    ax = annotate_PAs(args.pa_arr, ax, fontsize=args.fontsize/1.5)
    #cbar = plt.colorbar(p)

    if args.only_seg:
        ax.contour(args.segmentation_map != args.id, levels=0, colors='k', extent=args.extent, linewidths=0.5)

    if hide_xaxis:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel('RA (")', fontsize=args.fontsize)
        ax.tick_params(axis='x', which='major', labelsize=args.fontsize)

    if hide_yaxis:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel('Dec (")', fontsize=args.fontsize)
        ax.tick_params(axis='y', which='major', labelsize=args.fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def get_MAPPINGS_linelist(wave_lim=None):
    '''
    Reads in a MAPPINGS model file to grab the emission line list
    Returns list of lines that are within the given wavelength limits and above some threshold strength, as a pandas dataframe
    '''
    line_list_file = HOME / 'Desktop/Lisa_UV_diag/P_spherical/sp_P70_a05modelfiles/Q700/spec0003.csv'

    lines_df = pd.read_csv(line_list_file, skiprows=55, names=['wave', 'ev', 'flux', 'species', 'kind', 'acc'])
    if wave_lim is not None: lines_df = lines_df[lines_df['wave'].between(wave_lim[0], wave_lim[1])]
    lines_df = lines_df[lines_df['flux'] > 0.003]
    lines_df = lines_df[lines_df['kind'].str.strip().isin(['CM', 'RCAB'])]
    lines_df = lines_df[~lines_df['species'].str.contains('Fe')] # not interested in the numerous Fe lines

    print(f'Found {len(lines_df)} lines in this wavelength regime from {line_list_file}; over-plotting them now..')

    return lines_df

# --------------------------------------------------------------------------------------------------------------------
def plot_MAPPINGS_lines(ax):
    '''
    Plots a list of MAPPINGS emission line wavelengths on the given axis
    Returns axis handle
    '''
    lines_df = get_MAPPINGS_linelist(wave_lim=ax.get_xlim())
    max_flux = lines_df['flux'].max() * 1.001
    min_flux = lines_df['flux'].min() * 0.999

    for index in range(len(lines_df)):
        ax.axvline(lines_df.iloc[index]['wave'], c='cornflowerblue', lw=1, alpha= 0.3 + 0.7 * (lines_df.iloc[index]['flux'] - min_flux) / (max_flux - min_flux))
        xpos = lines_df.iloc[index]['wave'] + np.diff(ax.get_xlim())[0] * 0.01
        ypos = ax.get_ylim()[1] * 0.98 if index % 2 else 0.05 + ax.get_ylim()[0] * 1.02
        ax.text(xpos, ypos, lines_df.iloc[index]['species'].strip(), rotation=90, va='top' if index % 2 else 'bottom', ha='left', fontsize=args.fontsize)

    return ax


# --------------------------------------------------------------------------------------------------------------------
def get_linelist(wave_lim=None):
    '''
    Reads in an emission line list
    Returns list of lines that are within the given wavelength limits, as a pandas dataframe
    '''
    line_list_file = HOME / 'Desktop/mage_plot/labframe.shortlinelist'

    lines_df = pd.read_table(line_list_file, comment='#', delim_whitespace=True)
    lines_df = lines_df[~lines_df['LineID'].str.contains('Fe')] # not interested in the numerous Fe lines
    if wave_lim is not None: lines_df = lines_df[lines_df['restwave'].between(wave_lim[0], wave_lim[1])]

    print(f'Found {len(lines_df)} lines in this wavelength regime from {line_list_file}; over-plotting them now..')

    return lines_df


# --------------------------------------------------------------------------------------------------------------------
def plot_linelist(ax):
    '''
    Plots a list of emission line wavelengths on the given axis
    Returns axis handle
    '''
    lines_df = get_linelist(wave_lim=ax.get_xlim())

    for index in range(len(lines_df)):
        ax.axvline(lines_df.iloc[index]['restwave'], c='cornflowerblue', lw=1)
        xpos = lines_df.iloc[index]['restwave'] + np.diff(ax.get_xlim())[0] * 0.01
        ypos = ax.get_ylim()[1] * 0.98 if index % 2 else 0.02 + ax.get_ylim()[0] * 1.02
        ax.text(xpos, ypos, lines_df.iloc[index]['LineID'].strip(), rotation=90, va='top' if index % 2 else 'bottom', ha='left', fontsize=args.fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_1d_spectra(od_hdu, ax, args):
    '''
    Plots the 1D spectra in the given axis
    Returns the axis handle
    '''
    nfilters = sum(['GRISM' in item for item in list(od_hdu[0].header.keys())])
    filters = [od_hdu[0].header[f'GRISM{item + 1:03d}'] for item in range(nfilters)]

    col_arr = ['palegreen', 'gold', 'coral'] if args.fortalk else ['orange', 'orangered', 'firebrick'] # colors in order: for measured flux, fitted continuum, fitted continuum + line flux
    norm_factor = 1e-19

    # -------plot 1D spectra for each filter-----------
    for index, filter in enumerate(filters):
        print(f'Plotting 1D spectra for filter {filter} which is {index+1} of {nfilters}..')
        table = Table(od_hdu[filter].data).to_pandas()
        table = table[table['wave'].between(table['wave'].min() * (1 + args.trim_filter_by_wavelength_factor), table['wave'].max() * (1 - args.trim_filter_by_wavelength_factor))]
        table['rest_wave'] = table['wave'] / (1 + args.z)
        ax.plot(table['rest_wave'], table['flux'] / table['flat'] / norm_factor, lw=0.5, c=col_arr[index], alpha=0.8) # need to divide all columns with 'flat' to get the right units (ergs/s/cm^2/A)
        ax.plot(table['rest_wave'], table['cont'] / table['flat'] / norm_factor, lw=0.5, c='lightgray' if args.fortalk else 'grey')
        ax.plot(table['rest_wave'], table['line'] / table['flat'] / norm_factor, lw=1, c='cyan' if args.fortalk else 'indigo')
    for index, filter in enumerate(filters): # so that the filter labels come out at right y coordinate
        ax.text(float(filters[0][1:-1]) * 1e2 * 0.85 / (1 + args.z), ax.get_ylim()[1] * 0.98 - (index + 1) * 0.1 * np.diff(ax.get_ylim())[0], filter, c=col_arr[index], fontsize=args.fontsize, ha='left', va='top')

    ax.set_xlabel(r'Rest-frame wavelength ($\AA$)', fontsize=args.fontsize)
    ax.set_ylabel(r'f$_{\lambda}$ ' + '(%.0e ' % norm_factor + r'ergs/s/cm$^2$/A)', fontsize=args.fontsize/1.2)
    ax.set_ylim(0, args.flam_max) # flam_max should be in units of 1e-19 ergs/s/cm^2/A

    # ---observed wavelength axis-------
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(ax.get_xticks())
    ax2.set_xticklabels(['%.2F' % (item * (1 + args.z) / 1e4) for item in ax2.get_xticks()], fontsize=args.fontsize)
    ax2.set_xlabel(r'Observed wavelength ($\mu$)', fontsize=args.fontsize)

    # ---vertical lines for emission line wavelengths------
    if args.plot_mappings: ax = plot_MAPPINGS_lines(ax)
    else: ax = plot_linelist(ax)

    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    # --------for talk plots--------------
    if args.fortalk:
        mplcyberpunk.add_glow_effects()
        try: mplcyberpunk.make_lines_glow()
        except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    return ax

# ---------------------------------------------------------------------------------
def plot_binned_profile(df, ax, color='darkorange', yerr=None, xcol='radius', ycol='data'):
    '''
    Function to overplot binned data on existing plot in a given axis
    Returns axis handle
    '''
    bin_edges = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 10)
    df['binned_cat'] = pd.cut(df[xcol], bin_edges)

    if yerr is not None:
        df['weightcol'] = 1 / yerr.flatten() ** 2
        agg_func = lambda x: np.sum(x * df.loc[x.index, 'weightcol']) / np.sum(df.loc[x.index, 'weightcol']) # function to get weighted mean
        agg_u_func = lambda x: np.sqrt(((np.sum(df.loc[x.index, 'weightcol'] * x**2) / np.sum(df.loc[x.index, 'weightcol'])) - (np.sum(x * df.loc[x.index, 'weightcol']) / np.sum(df.loc[x.index, 'weightcol']))**2) * (np.sum(df.loc[x.index, 'weightcol']**2)) / (np.sum(df.loc[x.index, 'weightcol'])**2 - np.sum(df.loc[x.index, 'weightcol']**2))) # eq 6 of http://seismo.berkeley.edu/~kirchner/Toolkits/Toolkit_12.pdf
    else:
        agg_func, agg_u_func = np.mean, np.std

    y_binned = df.groupby('binned_cat', as_index=False).agg([(ycol, agg_func)])[ycol].values.flatten()
    y_u_binned = df.groupby('binned_cat', as_index=False).agg([(ycol, agg_u_func)])[ycol].values.flatten()
    x_bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2

    # ------getting rid of potential nan values---------
    indices = np.array(np.logical_not(np.logical_or(np.isnan(x_bin_centers), np.isnan(y_binned))))
    x_bin_centers = x_bin_centers[indices]
    y_binned = y_binned[indices]
    y_u_binned = y_u_binned[indices]

    # ----------to fit and plot the binned profile--------------
    try:
        linefit, linecov = np.polyfit(x_bin_centers, y_binned, 1, cov=True, w=1. / (y_u_binned) ** 2)
        y_fitted = np.poly1d(linefit)(x_bin_centers) # in logspace
        ax.plot(x_bin_centers, y_fitted, color=color, lw=1, ls='dashed')
    except Exception:
        print(f'Could not fit radial profile in this case..')
        linefit, linecov = [np.nan, np.nan], None

    # ----------to plot mean binned y vs x profile--------------
    ax.errorbar(x_bin_centers, y_binned, c=color, yerr=y_u_binned, lw=1, ls='none', zorder=1)
    ax.scatter(x_bin_centers, y_binned, c=color, s=20, lw=0.2, ec='black', zorder=10)

    return ax, linefit

# --------------------------------------------------------------------------------------------------------------
def get_distance_map(image_shape, args):
    '''
    Get map of distances from the center, in target rest-frame, on a given 2D grid
    Returns 2D distance map
    '''
    pixscale_kpc = args.pix_size_arcsec/ cosmo.arcsec_per_kpc_proper(args.z).value # kpc
    center_pix = image_shape[0] / 2.
    distance_map = np.array([[np.sqrt((i - center_pix)**2 + (j - center_pix)**2) for j in range(image_shape[1])] for i in range(image_shape[0])]) * pixscale_kpc # kpc

    return distance_map

# --------------------------------------------------------------------------------------------------------------
def trim_image(image, args):
    '''
    Trim a given 2D image to a given arcsecond dimension
    Returns 2D map
    '''
    image_shape = np.shape(image)
    center_pix = int(image_shape[0] / 2.)
    farthest_pix = int(args.arcsec_limit / args.pix_size_arcsec) # both quantities in arcsec

    image = image[center_pix - farthest_pix : center_pix + farthest_pix, center_pix - farthest_pix : center_pix + farthest_pix]
    # print(f'Trimming image of original shape {image_shape} to {args.arcsec_limit} arcseconds, which is from pixels {center_pix - farthest_pix} to {center_pix + farthest_pix}, so new shape is {np.shape(image)}')

    return image

# --------------------------------------------------------------------------------------------------------------------
def plot_radial_profile(image, ax, args, label=None, ymin=None, ymax=None, hide_xaxis=False, hide_yaxis=False):
    '''
    Plots the average radial profile for a given 2D map in the given axis
    Returns the axis handle
    '''
    print(f'Plotting radial profile of {label}..')

    distance_map = get_distance_map(np.shape(image), args)
    try: distance_map = np.ma.masked_where(image.mask, distance_map)
    except AttributeError: distance_map = np.ma.masked_where(False, distance_map)

    # ----making the dataframe before radial profile plot--------------
    xcol, ycol = 'radius', 'data'
    df = pd.DataFrame({xcol: np.ma.compressed(distance_map), ycol: np.ma.compressed(image)})
    df = df[df[xcol] <= args.radius_max]
    df = df.sort_values(by=xcol)

    # --------processing the dataframe in case voronoi binning has been performed and there are duplicate data values------
    df_vorbinned = pd.DataFrame()
    df_vorbinned[xcol] = df.groupby([ycol], as_index=False).agg([(np.mean)])[xcol]['mean']
    df_vorbinned[ycol] = df.groupby([ycol], as_index=False).agg([(np.mean)])[ycol]
    df = df_vorbinned

    # -------proceeding with plotting--------
    ax.scatter(df[xcol], df[ycol], c='grey', s=1, alpha=0.2)

    ax.set_xlim(0, args.radius_max) # kpc
    ax.set_ylim(ymin, ymax)
    ax.set_box_aspect(1)

    ax, linefit = plot_binned_profile(df, ax, xcol=xcol, ycol=ycol)

    if hide_xaxis:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel('Distance (kpc)', fontsize=args.fontsize)
        #ax.set_xticklabels(['%d' % item / cosmo.arcsec_per_kpc_proper(args.z).value for item in ax.get_xticks()], fontsize=args.fontsize)
        ax.tick_params(axis='x', which='major', labelsize=args.fontsize)

    if hide_yaxis:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel(label, fontsize=args.fontsize)
        ax.tick_params(axis='y', which='major', labelsize=args.fontsize)

    return ax, linefit

# --------------------------------------------------------------------------------------------------------------------
def plot_2D_map(image, ax, args, takelog=True, label=None, cmap=None, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=False, radprof_ax=None, vorbin_ax=None, snr_ax=None, image_err=None):
    '''
    Plots the emission map for a given line in the given axis
    Returns the axis handle
    '''
    if image.dtype == 'O': # if it is an uncertainty variable
        if np.ma.isMaskedArray(image): image = np.ma.masked_where(image.mask, unp.nominal_values(image.data))
        else: image = unp.nominal_values(image.data)

    if cmap is None: cmap = 'cividis'
    print(f'Plotting 2D map of {label}..')

    if takelog: image = np.log10(image)

    p = ax.imshow(image, cmap=cmap, origin='lower', extent=args.extent, vmin=vmin, vmax=vmax)

    ax.text(ax.get_xlim()[0] * 0.88, ax.get_ylim()[1] * 0.88, label, c='k', fontsize=args.fontsize, ha='left', va='top', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
    ax.scatter(0, 0, marker='x', s=10, c='grey')

    if hide_xaxis:
        ax.set_xticklabels([])
    else:
        if args.plot_target_frame:
            ax.set_xlabel('Offset (kpc)', fontsize=args.fontsize)
            ax.set_xticklabels(['%d' % (item / cosmo.arcsec_per_kpc_proper(args.z).value) for item in ax.get_xticks()], fontsize=args.fontsize)
        else:
            ax.set_xlabel('RA (")', fontsize=args.fontsize)
        ax.tick_params(axis='x', which='major', labelsize=args.fontsize)

    if hide_yaxis:
        ax.set_yticklabels([])
    else:
        if args.plot_target_frame:
            ax.set_ylabel('Offset (kpc)', fontsize=args.fontsize)
            ax.set_yticklabels(['%d' % (item / cosmo.arcsec_per_kpc_proper(args.z).value) for item in ax.get_xticks()], fontsize=args.fontsize)
        else:
            ax.set_ylabel('Dec (")', fontsize=args.fontsize)
        ax.tick_params(axis='y', which='major', labelsize=args.fontsize)

    if not hide_cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad='2%')
        cbar = plt.colorbar(p, cax=cax, orientation='vertical')
        cbar.ax.tick_params(labelsize=args.fontsize)

    if args.plot_radial_profiles and radprof_ax is not None:
        radius_pix = args.radius_max * cosmo.arcsec_per_kpc_proper(args.z).value # arcsec
        circle = plt.Circle((0, 0), radius_pix, color='k', fill=False, lw=0.5)
        ax.add_patch(circle)
        radprof_ax, radprof_fit = plot_radial_profile(image, radprof_ax, args, label=label.split(r'$_{\rm int}')[0], ymin=vmin, ymax=vmax)
    else:
        radprof_fit = [np.nan, np.nan] # dummy values for when the fit was not performed

    if args.vorbin and args.plot_vorbin and vorbin_ax is not None:
        vorbin_IDs = args.voronoi_bin_IDs
        vorbin_IDs = np.ma.masked_where(image.mask, vorbin_IDs)
        _, _ = plot_2D_map(vorbin_IDs, vorbin_ax, args, takelog=False, label=label + ' vorbin', cmap='rainbow')

    if args.plot_snr and image_err is not None and snr_ax is not None:
        if takelog: image = 10 ** image # undoing the log step that was done initially
        snr_map = image / image_err
        _, _ = plot_2D_map(snr_map, snr_ax, args, takelog=False, label=label + ' SNR', cmap='rainbow')

    return ax, radprof_fit

# --------------------------------------------------------------------------------------------------------------------
def bin_2D(map, bin_IDs, map_err=None):
    '''
    Bin a given 2D map by given bin_IDs
    Returns the binned 2D map (of same shape as input map)
    '''
    binned_map = np.zeros(np.shape(map))
    for id in np.unique(bin_IDs):
        binned_map[bin_IDs == id] = map[bin_IDs == id].mean()

    try: binned_map = np.ma.masked_where(map.mask, binned_map) # propagating the masks from the original 2D image
    except AttributeError: binned_map = np.ma.masked_where(False, binned_map)

    if map_err is not None:
        binned_map_err = np.zeros(np.shape(map_err))
        for id in np.unique(bin_IDs):
            candidates = map_err[bin_IDs == id]
            binned_map_err[bin_IDs == id] = np.sqrt(np.sum(candidates ** 2)) / len(candidates) # this is the apropriate error propagation for mean() operation (which the flux is undergoing above)

        try: binned_map_err = np.ma.masked_where(map.mask, binned_map_err) # propagating the masks from the original 2D image
        except AttributeError: binned_map_err = np.ma.masked_where(False, binned_map_err)

    if map_err is None: return binned_map
    else: return binned_map, binned_map_err

# --------------------------------------------------------------------------------------------------------------------
def get_voronoi_bin_IDs(map, snr_thresh, plot=False, quiet=True):
    '''
    Compute the Voronoi bin IDs a given 2D map and corresponding uncertainty and SNR threshold
    Returns the 2D map (of same shape as input map) with just the IDs
    '''
    x_size, y_size = np.shape(map)
    map_err = unp.std_devs(map)
    map = unp.nominal_values(map)

    x_coords = np.repeat(np.arange(x_size), y_size)
    y_coords = np.tile(np.arange(y_size), x_size)

    map = np.ma.masked_where(~np.isfinite(map_err), map)
    map_err = np.ma.masked_where(~np.isfinite(map_err), map_err)

    map = np.ma.masked_where(map / map_err < 1, map)
    map_err = np.ma.masked_where(map / map_err < 1, map_err)

    binIDs, _, _, _, _, _, _, _ = voronoi_2d_binning(x_coords, y_coords, map.flatten(), map_err.flatten(), snr_thresh, plot=plot, quiet=quiet, cvt=False, pixelsize=1)
    binID_map = binIDs.reshape(np.shape(map))

    return binID_map

# --------------------------------------------------------------------------------------------------------------------
def cut_by_segment(map, args):
    '''
    Mask a given 2D map according to the segmentation map from the HDU
    Returns the masked map
    '''
    cut_map = np.ma.masked_where(args.segmentation_map != args.id, map)

    return cut_map

# --------------------------------------------------------------------------------------------------------------------
def get_emission_line_map(line, full_hdu, args, dered=True, for_vorbin=False):
    '''
    Retrieve the emission map for a given line from the HDU
    Returns the 2D line image
    '''

    line_index = np.where(args.available_lines == line)[0][0]
    ext = 5 + 2 * args.ndfilt + 4 * line_index
    line_map = full_hdu[ext].data * 1e-17 # in units of ergs/s/cm^2
    line_wave = full_hdu[ext].header['RESTWAVE'] # in Angstrom
    line_map_err = 1e-17 / (full_hdu[ext + 3].data ** 0.5)  # 1/sqrt(LINEWHT) = flux uncertainty; in units of ergs/s/cm^2

    line_map = trim_image(line_map, args)
    line_map_err = trim_image(line_map_err, args)

    if args.only_seg:
        line_map = cut_by_segment(line_map, args)
        line_map_err = cut_by_segment(line_map_err, args)

    if args.snr_cut is not None:
        snr_map = line_map / line_map_err
        line_map = np.ma.masked_where(~np.isfinite(snr_map), line_map)
        line_map = np.ma.masked_where(snr_map < args.snr_cut, line_map)

    if args.vorbin and not for_vorbin:
        if args.voronoi_line is None: # No reference emission line specified, so Voronoi IDs need to be computed now
            bin_IDs = get_voronoi_bin_IDs(line_map, line_map_err, args.voronoi_snr)
        else: # Reference emission line specified for Voronoi binning, so bin IDs have been pre-computed
            bin_IDs = args.voronoi_bin_IDs
        line_map, line_map_err = bin_2D(line_map, bin_IDs, map_err=line_map_err)

    line_map = unp.uarray(line_map, line_map_err)

    # -----------getting the integrated flux value-----------------
    line_int = full_hdu[0].header[f'FLUX{line_index + 1:03d}'] # ergs/s/cm^2
    line_int_err = full_hdu[0].header[f'ERR{line_index + 1:03d}'] # ergs/s/cm^2
    line_int = ufloat(line_int, line_int_err)

    # -----------getting the integrated EW value-----------------
    line_index_in_cov = int([item for item in list(full_hdu[2].header.keys()) if full_hdu[0].header[f'FLUX{line_index + 1:03d}'] == full_hdu[2].header[item]][0][5:])
    line_ew = full_hdu[2].header[f'EW50_{line_index_in_cov:03d}']
    line_ew_err = full_hdu[2].header[f'EWHW_{line_index_in_cov:03d}']
    line_ew = ufloat(line_ew, line_ew_err)

    # -----------getting the dereddened flux value-----------------
    if dered:
        line_map = get_dereddened_flux(line_map, line_wave, args.EB_V)
        line_int = get_dereddened_flux(line_int, line_wave, args.EB_V)

    return line_map, line_wave, line_int, line_ew

# --------------------------------------------------------------------------------------------------------------------
def plot_emission_line_map(line, full_hdu, ax, args, cmap='cividis', EB_V=None, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=False):
    '''
    Plots the emission map for a given line in the given axis
    Returns the axes handle
    '''

    line_map, line_wave, line_int, line_ew = get_emission_line_map(line, full_hdu, args, dered=False)
    ax, _ = plot_2D_map(line_map, ax, args, label=r'%s$_{\rm int}$ = %.1e' % (line, line_int.n), cmap=cmap, vmin=vmin, vmax=vmax, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar)
    ax.text(ax.get_xlim()[0] * 0.88, ax.get_ylim()[0] * 0.88, f'EW = {line_ew:.1e}' if line_ew < 1e-3 or line_ew > 1e3 else f'EW = {line_ew:.1f}', c='k', fontsize=args.fontsize, ha='left', va='bottom', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))

    return ax

# -------------------------------------------------------------------------------
def get_kappa(x, i):
    '''
    To calculate kappa according to Clayton, Cardelli & Mathis 1989 dust law
    From ayan_codes/mage_project/ayan/mage.py
    '''
    Rv = 3.1  # Clayton Cardelli Mathis 1989
    x = np.array(x)
    if i == 1:
        a = 0.574 * x ** 1.61
        b = -0.527 * x ** 1.61
    elif i == 2:
        y = x - 1.82
        a = 1 + 0.17699 * y - 0.50447 * y ** 2 - 0.02427 * y ** 3 + 0.72085 * y ** 4 + 0.01979 * y ** 5 - 0.77530 * y ** 6 + 0.32999 * y ** 7
        b = 1.41338 * y + 2.28305 * y ** 2 + 1.07233 * y ** 3 - 5.38434 * y ** 4 - 0.62251 * y ** 5 + 5.30260 * y ** 6 - 2.09002 * y ** 7
    elif i == 3:
        a = 1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341)
        b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263)
    elif i == 4:
        a = 1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341) - 0.04473 * (x - 5.9) ** 2 - 0.009779 * (x - 5.9) ** 3
        b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263) + 0.2130 * (x - 5.9) ** 2 - 0.1207 * (x - 5.9) ** 3
    elif i == 5:
        a = -1.073 - 0.628 * (x - 8) + 0.137 * (x - 8) ** 2 - 0.070 * (x - 8) ** 3
        b = 13.670 + 4.257 * (x - 8) - 0.420 * (x - 8) ** 2 + 0.374 * (x - 8) ** 3
    return a * Rv + b


# -------------------------------------------------------------------------
def get_full_kappa(wave, inAngstrom=True):
    '''
    To calculate kappa for a rang eof wavelengths, according to CCM89 dust law
    From ayan_codes/mage_project/ayan/mage.py
    '''
    flag = 0
    if type(wave) in [float, int, np.float64]:
        wave = [float(wave)]
        flag = 1
    wave = np.array(wave)
    if inAngstrom: wave /= 1e4  # to convert to micron
    x = 1. / wave
    k = np.zeros(len(x))
    k += get_kappa(x, 1) * ((x >= 0.3) & (x <= 1.1))
    k += get_kappa(x, 2) * ((x > 1.1) & (x <= 3.3))
    k += get_kappa(x, 3) * ((x > 3.3) & (x < 5.9))
    k += get_kappa(x, 4) * ((x >= 5.9) & (x <= 8.))
    k += get_kappa(x, 5) * ((x >= 8.) & (x <= 10.))
    if flag: k = k[0]  # if a single wavelength was input as a float, output a float (not array)
    return k

# ------------------Function to calculate extinction and de-redden fluxes-------------
def get_dereddened_flux(obs_flux_map, wavelength, EB_V, inAngstrom=True):
    '''
    Calculates and returns dereddened fluxes
    Does not deal with uncertainties, for now
    From ayan_codes/mage_project/ayan/mage.py
    '''
    if EB_V is None: EB_V = np.zeros(np.shape(obs_flux_map))
    kappa = get_full_kappa(wavelength, inAngstrom=inAngstrom)
    A_map = kappa * EB_V
    flux_corrected_map = obs_flux_map * 10 ** (0.4 * A_map)

    return flux_corrected_map

# ---------------------------------------------------------------------------------------
def compute_EB_V(Ha_flux, Hb_flux,verbose=False):
    '''
    Calculates and returns the color excess given observed H alpha and H beta fluxes
    Based on Eqn 4 of Dominguez+2013 (https://iopscience.iop.org/article/10.1088/0004-637X/763/2/145/pdf)
    '''
    # -----handling masks separately because uncertainty package cannot handle masks---------
    if np.ma.isMaskedArray(Ha_flux):
        net_mask = Ha_flux.mask | Hb_flux.mask
        Ha_flux = Ha_flux.data
        Hb_flux = Hb_flux.data
    else:
        net_mask = False

    theoretical_ratio = 2.86

    if hasattr(Hb_flux, "__len__"): # if it is an array
        # --------computing the ratio and appropriate errors------------
        new_mask = Hb_flux == 0
        Hb_flux[new_mask] = -1 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        obs_ratio = Ha_flux / Hb_flux
        obs_ratio = np.ma.masked_where(new_mask, obs_ratio)

        # --------computing the log of the ratio and appropriate errors------------
        new_mask = obs_ratio <= 0
        obs_ratio[new_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        EB_V = 1.97 * unp.log10(obs_ratio.data / theoretical_ratio)
        EB_V[obs_ratio.data < theoretical_ratio] = 0
        EB_V = np.ma.masked_where(new_mask | obs_ratio.mask | net_mask, EB_V)

    else: # if it is scalar
        try:
            obs_ratio = Ha_flux / Hb_flux
            if obs_ratio < theoretical_ratio:
                EB_V = ufloat(0, 0)
                if verbose: print(f'Based on integrated fluxes Ha = {Ha_flux:.2e}, Hb = {Hb_flux:.2e}, the observed ratio is LOWER than theoretical ratio, so E(B-V) would be unphysical, so just assuming {EB_V}')
            else:
                EB_V = 1.97 * unp.log10(obs_ratio / theoretical_ratio)
                if verbose: print(f'Based on integrated fluxes Ha = {Ha_flux:.2e}, Hb = {Hb_flux:.2e}, determined E(B-V) =', EB_V)
        except:
            EB_V = ufloat(np.nan, np.nan)

    return EB_V

# --------------------------------------------------------------------------------------------------------------------
def get_EB_V(full_hdu, args, verbose=False):
    '''
    Computes and returns the spatially resolved as well as integrated dust extinction map from a given HDU
    Based on Eqn 4 of Dominguez+2013 (https://iopscience.iop.org/article/10.1088/0004-637X/763/2/145/pdf)
    '''

    Ha_map, Ha_wave, Ha_int, _ = get_emission_line_map('Ha', full_hdu, args, dered=False) # do not need to deredden the lines when we are fetching the flux in order to compute reddening
    Hb_map, Hb_wave, Hb_int, _ = get_emission_line_map('Hb', full_hdu, args, dered=False)

    EB_V_map = compute_EB_V(Ha_map, Hb_map)
    EB_V_int = compute_EB_V(Ha_int, Hb_int, verbose=verbose)

    return EB_V_map, EB_V_int

# --------------------------------------------------------------------------------------------------------------------
def plot_EB_V_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the dust extinction map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D E(B-V) map just produced
    '''
    lim, label = [0, 1], 'E(B-V)'

    EB_V_map, EB_V_int = get_EB_V(full_hdu, args)
    ax, EB_V_radfit = plot_2D_map(EB_V_map, ax, args, takelog=False, label=r'%s$_{\rm int}$ = %.1f' % (label, EB_V_int.n), cmap='YlOrBr', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

    return ax, EB_V_map, EB_V_radfit, EB_V_int

# --------------------------------------------------------------------------------------------------------------------
def compute_SFR(Ha_flux, distance):
    '''
    Calculates and returns the SFR given observed H alpha fluxes
    Conversion factor is from Kennicutt 1998 (Eq 2 of https://ned.ipac.caltech.edu/level5/Sept01/Rosa/Rosa3.html)
    '''
    # -----handling masks separately because uncertainty package cannot handle masks---------
    if np.ma.isMaskedArray(Ha_flux):
        net_mask = Ha_flux.mask
        Ha_flux = Ha_flux.data
    else:
        net_mask = False

    Ha_flux = Ha_flux * 4 * np.pi * (distance.to('cm').value) ** 2 # converting to ergs/s
    sfr = Ha_flux * 7.9e-42 # line_map in args/s; SFR in Msun/yr

    if hasattr(Ha_flux, "__len__"): # if it is an array
        sfr = np.ma.masked_where(net_mask, sfr)

    return sfr

# --------------------------------------------------------------------------------------------------------------------
def get_SFR(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated SFR from a given HDU
    '''
    Ha_map, Ha_wave, Ha_int, _ = get_emission_line_map('Ha', full_hdu, args)

    SFR_map = compute_SFR(Ha_map, args.distance)
    SFR_int = compute_SFR(Ha_int, args.distance)

    return SFR_map, SFR_int

# --------------------------------------------------------------------------------------------------------------------
def plot_SFR_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the SFR map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D SFR density map just produced
    '''
    lim, label = [-4, -2], 'SFR'
    SFR_map, SFR_int = get_SFR(full_hdu, args)
    ax, SFR_radfit = plot_2D_map(SFR_map, ax, args, label=r'%s$_{\rm int}$ = %.1f' % (label, SFR_int.n), cmap='Blues', radprof_ax=radprof_ax, vmin=lim[0], vmax=lim[1])

    return ax, SFR_map, SFR_radfit, SFR_int

# --------------------------------------------------------------------------------------------------------------------
def compute_Te(OIII4363_flux, OIII5007_flux):
    '''
    Calculates and returns the Te given observed line fluxes
    Conversion factor is from Nicholls+2017
    '''
    # -----handling masks separately because uncertainty package cannot handle masks---------
    if np.ma.isMaskedArray(OIII4363_flux):
        net_mask = OIII4363_flux.mask | OIII5007_flux.mask
        OIII4363_flux = OIII4363_flux.data
        OIII5007_flux = OIII5007_flux.data
    else:
        net_mask = False

    if hasattr(OIII5007_flux, "__len__"): # if it is an array
        # --------computing the ratio and appropriate errors------------
        new_mask = OIII5007_flux == 0
        OIII5007_flux[new_mask] = -1 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        ratio = OIII4363_flux / OIII5007_flux
        ratio = np.ma.masked_where(new_mask, ratio)

        # --------computing the log of the ratio and appropriate errors------------
        new_mask = ratio <= 0
        ratio[new_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        log_ratio = unp.log10(ratio.data)
        logTe = np.poly1d([0., 9.18962, 3.30355])(log_ratio) / np.poly1d([-0.00935, -0.15034, 2.09136, 1.000])(log_ratio)

        new_mask2 = logTe > 10.
        new_mask3 = logTe < -8.
        logTe[new_mask2 | new_mask3] = 0 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        Te = 10 ** logTe
        Te = np.ma.masked_where(new_mask | new_mask2 | new_mask3 | ratio.mask | net_mask, Te)

    else: # if it is scalar
        try:
            ratio = OIII4363_flux / OIII5007_flux  # in case 'OIII5007_flux' happens to be 0
            log_ratio = unp.log10(ratio.data)
            logTe = np.poly1d([0., 9.18962, 3.30355])(log_ratio) / np.poly1d([-0.00935, -0.15034, 2.09136, 1.000])(log_ratio)
            Te = 10 ** logTe
        except:
            Te = ufloat(np.nan, np.nan)

    return Te

# --------------------------------------------------------------------------------------------------------------------
def get_Te(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated Te from a given HDU
    '''
    OIII4363_map, OIII4363_wave, OIII4363_int, _ = get_emission_line_map('OIII-4363', full_hdu, args)
    OIII5007_map, OIII5007_wave, OIII5007_int, _ = get_emission_line_map('OIII', full_hdu, args)

    Te_map = compute_Te(OIII4363_map, OIII5007_map)
    Te_int = compute_Te(OIII4363_int, OIII5007_int)

    return Te_map, Te_int

# --------------------------------------------------------------------------------------------------------------------
def plot_Te_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the T_e map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D T_e map just produced
    '''
    lim, label = [1, 7], r'T$_e$'
    Te_map, Te_int = get_Te(full_hdu, args)
    ax, Te_radfit = plot_2D_map(Te_map, ax, args, label=r'%s$_{\rm int}$ = %.1e' % (label, Te_int.n), cmap='OrRd_r', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

    return ax, Te_map, Te_radfit, Te_int

# --------------------------------------------------------------------------------------------------------------------
def compute_Z_Te(OII3727_flux, OIII5007_flux, Hbeta_flux, Te, ne=1e3):
    '''
    Calculates and returns the Te metallicity given observed line fluxes
    Conversion factor is from Nicholls+2017
    '''
    # -----handling masks separately because uncertainty package cannot handle masks---------
    if np.ma.isMaskedArray(OII3727_flux):
        net_mask = OII3727_flux.mask | OIII5007_flux.mask | Hbeta_flux.mask | Te.mask
        OII3727_flux = OII3727_flux.data
        OIII5007_flux = OIII5007_flux.data
        Hbeta_flux = Hbeta_flux.data
        Te = Te.data
    else:
        net_mask = False

    # -----------------------------------------------------------
    def poly(R, t, x, a, b, c, d, e):
        return unp.log10(R) + a + b / t - c * unp.log10(t) - d * t + unp.log10(1 + e * x)  # eqn 3 I06 pattern

    t = Te * 1e-4
    x = 1e-4 * ne * unp.sqrt(t)

    if hasattr(Hbeta_flux, "__len__"): # if it is an array
        # --------computing the ratio and appropriate errors------------
        new_mask = Hbeta_flux == 0
        Hbeta_flux[new_mask] = -1 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        ratio1 = OII3727_flux / Hbeta_flux
        ratio2 = OIII5007_flux / Hbeta_flux
        ratio1 = np.ma.masked_where(new_mask, ratio1)
        ratio2 = np.ma.masked_where(new_mask, ratio2)

        # --------computing the log of the ratio and polynomial and appropriate errors------------
        new_mask1 = ratio1 <= 0
        ratio1[new_mask1] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        log_O2H2 = poly(ratio1, t, x, 5.961, 1.676, 0.4, 0.034, 1.35) - 12  # coefficients from eqn 3 I06
        new_mask3 = log_O2H2 > 2.
        log_O2H2[new_mask3] = 0 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        log_O2H2 = np.ma.masked_where(new_mask1 | new_mask3 | ratio1.mask, log_O2H2)

        new_mask2 = ratio2 <= 0
        ratio2[new_mask2] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        log_O3H2 = poly(ratio2, t, x, 6.200, 1.251, -0.55, -0.014, 0.0) - 12  # coefficients from eqn 5 I06
        new_mask4 = log_O3H2 > 2.
        log_O3H2[new_mask4] = 0 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        log_O3H2 = np.ma.masked_where(new_mask2 | new_mask4 | ratio2.mask, log_O3H2)

        # --------computing the combined metallicity and masks, etc.------------
        log_OH = unp.log10(10 ** log_O2H2.data + 10 ** log_O3H2.data) + 12
        log_OH = np.ma.masked_where(log_O3H2.mask | log_O2H2.mask | net_mask, log_OH)

    else: # if it is scalar
        try:
            ratio1 = OII3727_flux / Hbeta_flux # in case 'Hbeta_flux' happens to be 0
            ratio2 = OIII5007_flux / Hbeta_flux

            log_O2H2 = poly(ratio1, t, x, 5.961, 1.676, 0.4, 0.034, 1.35) - 12  # coefficients from eqn 3 I06 # in case 'ratio1' happens to be < 0
            log_O3H2 = poly(ratio2, t, x, 6.200, 1.251, -0.55, -0.014, 0.0) - 12  # coefficients from eqn 5 I06 # in case 'ratio1' happens to be < 0

            log_OH = unp.log10(10 ** log_O2H2 + 10 ** log_O3H2) + 12
        except:
            log_OH = ufloat(np.nan, np.nan)

    return log_OH

# --------------------------------------------------------------------------------------------------------------------
def get_Z_Te(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated Te metallicity from a given HDU
    '''
    OII3727_map, line_wave, OII3727_int, _ = get_emission_line_map('OII', full_hdu, args)
    OIII5007_map, line_wave, OIII5007_int, _ = get_emission_line_map('OIII', full_hdu, args)
    Hbeta_map, line_wave, Hbeta_int, _ = get_emission_line_map('Hb', full_hdu, args)

    Te_map, Te_int = get_Te(full_hdu, args)

    logOH_map = compute_Z_Te(OII3727_map, OIII5007_map, Hbeta_map, Te_map)
    logOH_int = compute_Z_Te(OII3727_int, OIII5007_int, Hbeta_int, Te_int)

    return logOH_map, logOH_int

# --------------------------------------------------------------------------------------------------------------------
def plot_Z_Te_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the T_e-based metallicity map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D metallicity map just produced
    '''
    lim, label = [6, 9], 'Z (Te)'
    logOH_map, logOH_int = get_Z_Te(full_hdu, args)
    ax, logOH_radfit = plot_2D_map(logOH_map, ax, args, takelog=False, label=r'%s$_{\rm int}$ = %.1f' % (label, logOH_int.n), cmap='Greens', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

    return ax, logOH_map, logOH_radfit, logOH_int

# --------------------------------------------------------------------------------------------------------------------
def compute_Z_R23(OII3727_flux, OIII5007_flux, Hbeta_flux):
    '''
    Calculates and returns the R23 metallicity given observed line fluxes
    Conversion factor is from Kewley+2002
    '''
    # -----handling masks separately because uncertainty package cannot handle masks---------
    if np.ma.isMaskedArray(OII3727_flux):
        net_mask = OII3727_flux.mask | OIII5007_flux.mask | Hbeta_flux.mask
        OII3727_flux = OII3727_flux.data
        OIII5007_flux = OIII5007_flux.data
        Hbeta_flux = Hbeta_flux.data
    else:
        net_mask = False

    k = [-44.7026, 10.8052, -0.640113]  # k0-2 parameters for q=8e7 from Table 3 of KD02 last row for q=8e7

    if hasattr(Hbeta_flux, "__len__"): # if it is an array
        # --------computing the ratio and appropriate errors------------
        new_mask = Hbeta_flux == 0
        Hbeta_flux[new_mask] = -1 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        ratio = (OII3727_flux + OIII5007_flux) / Hbeta_flux
        ratio = np.ma.masked_where(new_mask, ratio)

        # --------computing the log of the ratio and appropriate errors------------
        new_mask = ratio <= 0
        ratio[new_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        R23 = unp.log10(ratio.data)
        R23 = np.ma.masked_where(new_mask | ratio.mask, R23)

        # --------computing the polynomial and appropriate errors------------
        p = k[1] ** 2 - 4 * k[2] * (k[0] - R23)
        mask = p < 0
        p[mask] = 0 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        log_OH = (-k[1] + unp.sqrt(p.data))/(2 * k[2])
        log_OH = np.ma.masked_where(mask | R23.mask | net_mask, log_OH)

    else: # if it is scalar
        try:
            ratio =(OII3727_flux + OIII5007_flux) / Hbeta_flux # in case 'Hbeta_flux' happens to be 0
            R23 = unp.log10(ratio)
            p = k[1] ** 2 - 4 * k[2] * (k[0] - R23)
            log_OH = (-k[1] + unp.sqrt(p)) / (2 * k[2])
        except:
            log_OH = ufloat(np.nan, np.nan)

    return log_OH

# --------------------------------------------------------------------------------------------------------------------
def get_Z_R23(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated R23 metallicity from a given HDU
    '''
    OII3727_map, line_wave, OII3727_int, _ = get_emission_line_map('OII', full_hdu, args)
    OIII5007_map, line_wave, OIII5007_int, _ = get_emission_line_map('OIII', full_hdu, args)
    Hbeta_map, line_wave, Hbeta_int, _ = get_emission_line_map('Hb', full_hdu, args)

    logOH_map = compute_Z_R23(OII3727_map, OIII5007_map, Hbeta_map)
    logOH_int = compute_Z_R23(OII3727_int, OIII5007_int, Hbeta_int)

    return logOH_map, logOH_int

# --------------------------------------------------------------------------------------------------------------------
def plot_Z_R23_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the R23 metallicity map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D metallicity map just produced
    '''
    lim, label = [6, 9], 'Z (R23)'
    logOH_map, logOH_int = get_Z_R23(full_hdu, args)
    ax, logOH_radfit = plot_2D_map(logOH_map, ax, args, takelog=False, label=r'%s$_{\rm int}$ = %.1f' % (label, logOH_int.n), cmap='Greens', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

    return ax, logOH_map, logOH_radfit, logOH_int

# --------------------------------------------------------------------------------------------------------------------
def plot_starburst_map(full_hdu, args):
    '''
    Plots the Ha map, direct F115W map and their ratio (starbursty-ness map) in a new figure
    Returns the figure handle and the ratio map just produced
    '''
    nrows = 1
    if args.plot_vorbin: nrows += 1
    if args.plot_snr: nrows += 1
    if args.plot_radial_profiles: nrows += 1

    fig_size_dict = {1: [14, 4, 0.05, 0.97, 0.15, 0.9, 0.2, 0.], 2: [14, 7, 0.02, 0.98, 0.07, 0.95, 0.05, 0.2], 3: [9, 7, 0.02, 0.98, 0.07, 0.95, 0.05, 0.3], 4: [9, 7, 0.02, 0.98, 0.07, 0.95, 0.05, 0.3]} # figsize_w, figsize_h, l, r, b, t, ws, hs

    fig, axes = plt.subplots(nrows, 3, figsize=(fig_size_dict[nrows][0], fig_size_dict[nrows][1]))
    if nrows > 1:
        extra_axes = axes[1:,:]
        axes = axes[0,:]
    if args.plot_vorbin:
        vorbin_axes = extra_axes[0,:].flatten()
        extra_axes = extra_axes[1:,:]
    if args.plot_snr:
        snr_axes = extra_axes[0,:].flatten()
        extra_axes = extra_axes[1:,:]
    if args.plot_radial_profiles:
        radprof_axes = extra_axes[0:, :].flatten()
        extra_axes = extra_axes[1:, :]
    fig.subplots_adjust(left=fig_size_dict[nrows][2], right=fig_size_dict[nrows][3], bottom=fig_size_dict[nrows][4], top=fig_size_dict[nrows][5], wspace=fig_size_dict[nrows][6], hspace=fig_size_dict[nrows][7])

    # ---------getting the Ha map-------------
    ha_map, _, _, _ = get_emission_line_map('Ha', full_hdu, args, dered=True)

    # ---------getting the direct image-------------
    direct_ext = 5
    direct_map = full_hdu[direct_ext].data
    direct_map_err = 1. / (full_hdu[direct_ext + 1].data ** 0.5)

    direct_map = trim_image(direct_map, args)
    direct_map = np.ma.masked_where(~np.isfinite(direct_map), direct_map)

    direct_map_err = trim_image(direct_map_err, args)
    direct_map_err = np.ma.masked_where(~np.isfinite(direct_map), direct_map_err)

    if args.only_seg:
        direct_map = cut_by_segment(direct_map, args)

    if args.vorbin:
        direct_map = bin_2D(direct_map, args.voronoi_bin_IDs)
        direct_map_err = bin_2D(direct_map_err, args.voronoi_bin_IDs)

    direct_map = unp.uarray(direct_map, direct_map_err)

    # ---------getting the ratio and seg maps-------------
    new_mask = unp.nominal_values(direct_map.data) == 0
    direct_map[new_mask] = 1 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
    direct_map = np.ma.masked_where(new_mask, direct_map)

    ratio_map = ha_map.data / direct_map.data
    ratio_map = np.ma.masked_where(ha_map.mask | direct_map.mask | new_mask, ratio_map)

    # --------making arrays for subplots-------------
    maps = [direct_map, ha_map, ratio_map]
    labels = ['Direct', r'H$\alpha$', r'H$\alpha$/Direct']
    lims = [[-4, -1], [-21, -18], [-18, -15]]

    # ---------plotting-------------
    for index, ax in enumerate(axes):
        map_err = np.ma.masked_where(maps[index].mask, unp.std_devs(maps[index].data))
        ax, _ = plot_2D_map(maps[index], ax, args, label=labels[index], cmap='viridis', vmin=lims[index][0], vmax=lims[index][1], radprof_ax=radprof_axes[index], vorbin_ax=vorbin_axes[index] if args.plot_vorbin else None, snr_ax=snr_axes[index] if args.plot_snr else None, image_err=map_err if args.plot_snr else None)
        ax.contour(args.segmentation_map != args.id, levels=0, colors='k', extent=args.extent, linewidths=2)
        if args.plot_vorbin: vorbin_axes[index].contour(args.segmentation_map != args.id, levels=0, colors='k', extent=args.extent, linewidths=2)
        if args.plot_snr: radprof_axes[index].contour(args.segmentation_map != args.id, levels=0, colors='k', extent=args.extent, linewidths=2)

    return fig, ratio_map

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_map(full_hdu, args):
    '''
    Plots the metallicity map, and optionally radial profile, in a new figure
    Returns the figure handle and the metallicity map just produced
    '''
    ncols = 2
    if args.plot_radial_profiles: ncols += 1

    fig_size_dict = {2: [10, 4, 0.05, 0.97, 0.15, 0.9, 0.2, 0.], 3: [14, 4, 0.05, 0.97, 0.15, 0.9, 0.2, 0.]} # figsize_w, figsize_h, l, r, b, t, ws, hs

    fig, axes = plt.subplots(1, ncols, figsize=(fig_size_dict[ncols][0], fig_size_dict[ncols][1]))
    if ncols > 2:
        radprof_ax = axes[2]
        axes = axes[:2]
    else:
        radprof_ax = None
    fig.subplots_adjust(left=fig_size_dict[ncols][2], right=fig_size_dict[ncols][3], bottom=fig_size_dict[ncols][4], top=fig_size_dict[ncols][5], wspace=fig_size_dict[ncols][6], hspace=fig_size_dict[ncols][7])

    if all([line in args.available_lines for line in ['OIII', 'OII', 'Hb']]):
        # --------deriving the metallicity map-------------
        logOH_map, logOH_int = get_Z_R23(full_hdu, args)

        # ---------plotting-------------
        lim, label = [6.5, 8.5], 'log(O/H) (R23)'
        axes[0], logOH_radfit = plot_2D_map(logOH_map, axes[0], args, takelog=False, label=r'%s$_{\rm int}$ = %.1f $\pm$ %.1f' % (label, logOH_int.n, logOH_int.s), cmap='viridis', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])
        axes[0].contour(args.segmentation_map != args.id, levels=0, colors='k', extent=args.extent, linewidths=2)

        logOH_map_err = np.ma.masked_where(logOH_map.mask, unp.std_devs(logOH_map.data))
        axes[1], _ = plot_2D_map(logOH_map_err, axes[1], args, takelog=False, label=r'%s uncertainty' % (label), cmap='cividis')
        axes[1].contour(args.segmentation_map != args.id, levels=0, colors='k', extent=args.extent, linewidths=2)

    else:
        print(f'Not all lines out of OIII, OII and Hb are available, so cannot compute R23 metallicity')
        logOH_map = None

    return fig, logOH_map

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ---------determining filename suffixes-------------------------------
    product_dir = args.input_dir / args.field / 'Products'

    radial_plot_text = '_wradprof' if args.plot_radial_profiles else ''
    snr_text = f'_snr{args.snr_cut}' if args.snr_cut is not None else ''
    only_seg_text = '_onlyseg' if args.only_seg else ''
    pixscale_text = '' if args.pixscale == 0.04 else f'_{args.pixscale}arcsec_pix'
    vorbin_text = '' if not args.vorbin else f'_vorbin_at_{args.voronoi_line}_SNR_{args.voronoi_snr}'
    description_text = f'all_diag_plots{radial_plot_text}{snr_text}{only_seg_text}'

    output_dir = args.output_dir / args.field
    if args.re_extract: output_dir = output_dir / 're_extracted'
    output_dir.mkdir(parents=True, exist_ok=True)
    outfilename = output_dir / f'{args.field}_all_diag_results.txt'

    # ---------prep the catalog file--------------------
    catalog_file = product_dir / f'{args.field}_photcat.fits'
    catalog = GTable.read(catalog_file)

    # --------determine which objects to loop over----------
    if args.do_all_obj:
        if args.re_extract: args.id_arr = ids_to_re_extract_dict[args.field]
        else: args.id_arr = catalog['id']
    else:
        args.id_arr = args.id

    if args.start_id: args.id_arr = args.id_arr[args.start_id - 1:]
    fig_dir = output_dir / f'{description_text}'
    fig_dir.mkdir(parents=True, exist_ok=True)

    # ---------for diplay and amimations----------------
    if len(args.id_arr) > 20: args.hide = True # if too many plots, do not display them, just save them
    if len(args.id_arr) > 20: args.make_anim = False #True
    else: args.make_anim = False

    if args.make_anim:
        outputfile = output_dir / f'{args.field}_{description_text}.mp4'
        duration_per_frame = 0.1 #sec
        writer = imageio.get_writer(outputfile, mode='I', fps=int(1. / duration_per_frame))

    # ----------------------initiliasing dataframe-------------------------------
    all_lines_to_consider = ['Ha', 'Hb', 'OII', 'OIII-4363', 'OIII'] # args.line_list
    measured_quantities_to_plot = ['EB_V', 'SFR', 'Te', 'logOH_Te', 'logOH_R23']

    if args.write_file:
        basic_cols = ['field', 'objid', 'ra', 'dec', 'redshift']
        flag_cols = ['radfit_extent_kpc', 'snr_cut', 'flag_only_seg', 'flag_vorbin', 'vor_snr', 'vor_line']
        cols_in_df = np.hstack([basic_cols, flag_cols, np.hstack([[item + '_int', item + '_EW'] for item in all_lines_to_consider]), np.hstack([[item + '_int', item + '_cen', item + '_slope'] for item in measured_quantities_to_plot])])
        df = pd.DataFrame(columns=cols_in_df)

        # -------checking if about to write the same columns-----------
        if os.path.exists(outfilename) and not args.clobber:
            existing_df = pd.read_table(outfilename, delim_whitespace=True)
            existing_cols = existing_df.columns
            if set(existing_cols) != set(cols_in_df):
                new_cols = set(cols_in_df) - set(existing_cols)
                print(
                f'Existing dataframe at {outfilename} has a different set of columns (the difference being {new_cols}) than currently trying to '
                f'write it in. Either delete/double check the existing file, OR use --clobber to overwrite the '
                f'existing file. Aborting.')
                sys.exit()

    # ------------looping over the provided object IDs-----------------------
    for index, args.id in enumerate(args.id_arr):
        start_time2 = datetime.now()
        print(f'\nCommencing ID {args.id} which is {index+1} of {len(args.id_arr)}..')

        # ------determining directories---------
        output_subdir = output_dir / f'{args.id:05d}{pixscale_text}'
        full_fits_file = output_subdir / f'{args.field}_{args.id:05d}.full.fits'
        maps_fits_file = product_dir / 'maps' / f'{args.field}_{args.id:05d}.maps.fits'

        if os.path.exists(full_fits_file): # if the fits files are in sub-directories for individual objects
            full_filename = full_fits_file
            od_filename = output_subdir/ f'{args.field}_{args.id:05d}.1D.fits'

        elif os.path.exists(maps_fits_file): # if the fits files are in Products/
            full_filename = maps_fits_file
            od_filename = product_dir / 'spec1D' / f'{args.field}_{args.id:05d}.1D.fits'

        else:
            print(f'Could not find {full_fits_file} or {maps_fits_file} for ID {args.id}, so skipping it.')
            continue

        if not os.path.exists(od_filename): od_filename = Path(str(od_filename).replace('.1D.', '.spec1D.'))

        # ------------read in fits files--------------------------------
        if os.path.exists(full_filename):
            full_hdu = fits.open(full_filename)
        else:
            print('Full fits file does not exists, cannot proceed, so skipping..')
            continue

        if os.path.exists(od_filename):
            od_hdu = fits.open(od_filename)
        else:
            print('1D fits file does not exists, cannot proceed, so skipping..')
            continue

        # ----------determining global parameters------------
        args.available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
        args.z = full_hdu[0].header['REDSHIFT']
        args.ndfilt = full_hdu[0].header['NDFILT']
        args.nlines = full_hdu[0].header['NUMLINES']
        args.pix_arcsec = full_hdu[5].header['PIXASEC']
        args.distance = cosmo.comoving_distance(args.z)
        args.pa_arr = np.unique([full_hdu[0].header[item] for item in list(full_hdu[0].header.keys()) if 'PA00' in item])
        try: args.mag = catalog[catalog['id'] == args.id]['mag_auto'].data.data[0]
        except: args.mag = np.nan

        line_wcs = pywcs.WCS(full_hdu['DSCI'].header)
        args.pix_size_arcsec = utils.get_wcs_pscale(line_wcs)
        imsize_arcsec = full_hdu['DSCI'].data.shape[0] * args.pix_size_arcsec
        args.extent = (-args.arcsec_limit, args.arcsec_limit, -args.arcsec_limit, args.arcsec_limit)
        args.EB_V = 0. # until gets over-written, if both H alpha and H beta lines are present

        # ---------------segmentation map---------------
        segmentation_map = full_hdu['SEG'].data
        args.segmentation_map = trim_image(segmentation_map, args)

        # ---------------voronoi binning stuff---------------
        if args.vorbin and args.voronoi_line is not None:
            line_map, _, _, _ = get_emission_line_map(args.voronoi_line, full_hdu, args, for_vorbin=True)
            args.voronoi_bin_IDs = get_voronoi_bin_IDs(line_map, args.voronoi_snr)#, plot=True, quiet=False)

        # ---------------radial profile stuff---------------
        if args.plot_radial_profiles:
            seg_map = full_hdu['SEG'].data
            distance_map = get_distance_map(np.shape(seg_map), args)
            distance_map = np.ma.compressed(np.ma.masked_where(seg_map != args.id, distance_map))
            args.radius_max = np.max(distance_map)
        else:
            args.radius_max = np.nan

        # ---------------dust value---------------
        if all([line in args.available_lines for line in ['Ha', 'Hb']]):
            _, args.EB_V = get_EB_V(full_hdu, args, verbose=True)

        # ---------initialising the starburst figure------------------------------
        if args.plot_starburst:
            fig, ratio_map = plot_starburst_map(full_hdu, args)

            # ---------decorating and saving the figure------------------------------
            fig.text(0.05, 0.98, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
            figname = fig_dir / f'{args.field}_{args.id:05d}_starburst_maps{radial_plot_text}{snr_text}{only_seg_text}{vorbin_text}.png'

        # ---------initialising the metallicity figure------------------------------
        elif args.plot_metallicity:
            fig, logOH_map = plot_metallicity_map(full_hdu, args)

            # ---------decorating and saving the figure------------------------------
            fig.text(0.05, 0.98, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
            figname = fig_dir / f'{args.field}_{args.id:05d}_metallicity_maps{radial_plot_text}{snr_text}{only_seg_text}{vorbin_text}.png'

        # ---------initialising the full figure------------------------------
        else:
            nrow, ncol = 4 if args.plot_radial_profiles else 3, len(all_lines_to_consider)
            fig = plt.figure(figsize=(13/1., 9/1.) if args.plot_radial_profiles else (13, 6), layout='constrained')

            axis_dirimg = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 0), colspan=1)
            axis_1dspec = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 1), colspan=ncol - 1)
            ax_em_lines = [plt.subplot2grid(shape=(nrow, ncol), loc=(1, item), colspan=1) for item in np.arange(ncol)]  # H alpha, H beta, OII, OIII-4363, OIII
            [ax_SFR, ax_EB_V, ax_Te, ax_Z_Te, ax_Z_R23] = [plt.subplot2grid(shape=(nrow, ncol), loc=(2, item), colspan=1) for item in np.arange(ncol)]  # SFR, E(B-V), Te, Z (Te), Z (R23)
            if args.plot_radial_profiles:
                [rax_SFR, rax_EB_V, rax_Te, rax_Z_Te, rax_Z_R23] = [plt.subplot2grid(shape=(nrow, ncol), loc=(3, item), colspan=1) for item in np.arange(ncol)]  # SFR, E(B-V), Te, Z (Te), Z (R23)
            else:
                [rax_SFR, rax_EB_V, rax_Te, rax_Z_Te, rax_Z_R23] = np.tile(None, ncol)
            # ---------direct imaging------------------------------
            axis_dirimg = plot_direct_image(full_hdu, axis_dirimg, args)

            # ---------1D spectra------------------------------
            axis_1dspec = plot_1d_spectra(od_hdu, axis_1dspec, args)

            # -----------------emission line maps---------------
            for ind, line in enumerate(all_lines_to_consider):
                if line in args.available_lines: ax_em_lines[ind] = plot_emission_line_map(line, full_hdu, ax_em_lines[ind], args, cmap='BuPu', vmin=-20, vmax=-18, hide_xaxis=True, hide_yaxis=ind > 0, hide_cbar=False) #ind != len(all_lines_to_consider) - 1) # line_map in ergs/s/cm^2
                else: fig.delaxes(ax_em_lines[ind])

            # ---------------dust map---------------
            if all([line in args.available_lines for line in ['Ha', 'Hb']]):
                ax_EB_V, EB_V_map, EB_V_radfit, EB_V_int = plot_EB_V_map(full_hdu, ax_EB_V, args, radprof_ax=rax_EB_V)
            else:
                fig.delaxes(ax_EB_V)
                if args.plot_radial_profiles: fig.delaxes(rax_EB_V)
                EB_V_map, EB_V_radfit, EB_V_int = np.nan, [np.nan, np.nan], np.nan

            # ---------------SFR map------------------
            if 'Ha' in args.available_lines:
                ax_SFR, SFR_map, SFR_radfit, SFR_int = plot_SFR_map(full_hdu, ax_SFR, args, radprof_ax=rax_SFR)
            else:
                fig.delaxes(ax_SFR)
                if args.plot_radial_profiles: fig.delaxes(rax_SFR)

            # ---------------electron temperature map---------------
            if all([line in args.available_lines for line in ['OIII-4363', 'OIII']]):
                ax_Te, Te_map, Te_radfit, Te_int = plot_Te_map(full_hdu, ax_Te, args, radprof_ax=rax_Te)
            else:
                fig.delaxes(ax_Te)
                if args.plot_radial_profiles: fig.delaxes(rax_Te)

            # ---------------metallicity maps---------------
            if all([line in args.available_lines for line in ['OIII', 'OII', 'Hb']]):
                if 'OIII-4363' in args.available_lines: ax_Z_Te, logOH_Te_map, logOH_Te_radfit, logOH_Te_int = plot_Z_Te_map(full_hdu, ax_Z_Te, args, radprof_ax=rax_Z_Te)
                ax_Z_R23, logOH_R23_map, logOH_R23_radfit, logOH_R23_int = plot_Z_R23_map(full_hdu, ax_Z_R23, args, radprof_ax=rax_Z_R23)
            else:
                fig.delaxes(ax_Z_Te)
                fig.delaxes(ax_Z_R23)
                if args.plot_radial_profiles:
                    fig.delaxes(rax_Z_Te)
                    fig.delaxes(rax_Z_R23)

            # ---------decorating and saving the figure------------------------------
            fig.text(0.05, 0.98, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
            figname = fig_dir / f'{args.field}_{args.id:05d}_{description_text}{vorbin_text}.png'

        # --------for talk plots--------------
        if args.fortalk:
            mplcyberpunk.add_glow_effects()
            try: mplcyberpunk.make_lines_glow()
            except: pass
            try: mplcyberpunk.make_scatter_glow()
            except: pass

        fig.savefig(figname, transparent=args.fortalk)
        print(f'Saved figure at {figname}')
        if args.hide: plt.close('all')
        else: plt.show(block=False)

        # ------------------making animation--------------------------
        if args.make_anim:
            print(f'Appending file {figname} to animation..')  #
            try: writer.append_data(imageio.imread(figname))
            except (ValueError, IOError) as e: print(f'Skipping snapshot due to ' + str(e))

        # ----------appending and writing to catalog file-----------------
        if args.write_file:

            # -------collating all the integrated line fluxes from the HDU header----------
            line_properties = []
            for line in all_lines_to_consider:
                try:
                    line_index = np.where(args.available_lines == line)[0][0]
                    flux = full_hdu[0].header[f'FLUX{line_index + 1:03d}']
                    line_properties += [flux]
                except IndexError:
                    line_properties += [np.nan]
                try:
                    line_index = np.where(args.available_lines == line)[0][0]
                    line_index_in_cov = int([item for item in list(full_hdu[2].header.keys()) if full_hdu[0].header[f'FLUX{line_index + 1:03d}'] == full_hdu[2].header[item]][0][5:])
                    ew = full_hdu[2].header[f'EW50_{line_index_in_cov:03d}']
                    line_properties += [ew]
                except IndexError:
                    line_properties += [np.nan]

            # -------collating all the measured quantities----------
            measured_quants = []
            for quantity in measured_quantities_to_plot:
                if quantity + '_int' in locals():
                    measured_quants += [locals()[quantity + '_int']]
                else:
                    measured_quants += [np.nan]
                if quantity + '_radfit' in locals():
                    measured_quants += [locals()[quantity + '_radfit'][0], locals()[quantity + '_radfit'][1]]
                else:
                    measured_quants += [np.nan, np.nan]

            basic_data = [args.field, f'{args.id:05d}{pixscale_text}', full_hdu[0].header['RA'], full_hdu[0].header['DEC'], args.z]
            flag_data = [args.radius_max, args.snr_cut if args.snr_cut is not None else np.nan, args.only_seg, args.vorbin, args.voronoi_snr if args.vorbin else np.nan, args.voronoi_line if args.vorbin else np.nan]
            this_row = np.hstack([basic_data, flag_data, line_properties, measured_quants])
            this_df = pd.DataFrame(dict(map(lambda i, j: (i, [j]), cols_in_df, this_row)))
            df = pd.concat([df, this_df])

            if not os.path.isfile(outfilename) or (args.clobber and index == 0):
                this_df.to_csv(outfilename, sep='\t', index=None, header='column_names')
                print(f'Wrote to catalog file {outfilename}')
            else:
                this_df.to_csv(outfilename, sep='\t', index=None, mode='a', header=False)
                print(f'Appended to catalog file {outfilename}')

        print(f'Completed id {args.id} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(args.id_arr) - index - 1} to go!')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
