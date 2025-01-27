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
             run make_diagnostic_maps.py --field Par28 --id 58,646,1457,1585,1588,2195,2343 --plot_radial_profiles --only_seg --plot_mappings
             run make_diagnostic_maps.py --field Par28 --id 255,1018,1069,1144,1707,2077,2360,2640,2657 --plot_radial_profiles --only_seg --plot_mappings
             run make_diagnostic_maps.py --field Par28 --id 58,646,1457,1585,1588,2195,2343 --plot_radial_profiles --only_seg --plot_mappings --vorbin --voronoi_snr 3
             run make_diagnostic_maps.py --field Par28 --id 58,646,1457,1585,1588,2195,2343 --plot_starburst --vorbin --voronoi_snr 3 --plot_radial_profile --only_seg
             run make_diagnostic_maps.py --field Par28 --id 58,646,1457,1585,1588,2195,2343 --plot_slope_vs_mass --vorbin --voronoi_snr 3 --only_seg
             run make_diagnostic_maps.py --field Par28 --id 58,646,1457,1588,2195,2343 --plot_metallicity --vorbin --voronoi_snr 3 --plot_radial_profile --only_seg
             run make_diagnostic_maps.py --field Par28 --id 58,646,1457,1585,1588,2195,2343 --plot_direct_filters --plot_radial_profiles --only_seg --vorbin --voronoi_snr 3
             run make_diagnostic_maps.py --field Par28 --id 646,1457,1585,1588,2195,2343 --plot_BPT --only_seg --vorbin --voronoi_snr 3 --plot_separately
             run make_diagnostic_maps.py --field Par28 --id 2343 --test_cutout
             run make_diagnostic_maps.py --do_all_fields --do_all_obj --plot_radial_profiles --only_seg --snr_cut 3 --write_file --clobber

             run make_diagnostic_maps.py --field Par28 --id 58,155,690,1200,1228,1395,1457,1585,1588 --plot_radial_profiles --only_seg --vorbin --voronoi_snr 3
             run make_diagnostic_maps.py --field Par28 --id 58,155,690,1200,1228,1395,1457,1585,1588 --plot_BPT --only_seg --vorbin --voronoi_snr 3 --plot_separately
             run make_diagnostic_maps.py --field Par28 --id 58,155,638,689,916,1457,1588 --plot_metallicity --plot_radial_profile --only_seg --vorbin --voronoi_snr 3
             run make_diagnostic_maps.py --field Par28 --id 58,155,638,689,916,1457,1588 --plot_starburst --plot_radial_profile --only_seg --vorbin --voronoi_snr 3

             run make_diagnostic_maps.py --field Par28 --id 1457 --plot_AGN_frac --only_seg --vorbin --voronoi_snr 3

             run make_diagnostic_maps.py --field Par28 --id 1332,1500,1565,1692,1697,192,68,754 --plot_slope_vs_mass --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --drv 0.5
             run make_diagnostic_maps.py --field Par28 --id 1332,1500,1565,1692,1697,192,68,754 --plot_metallicity --plot_radial_profile --plot_ion --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --drv 0.5
             run make_diagnostic_maps.py --field Par28 --id 192 --plot_metallicity --ignore_combined --plot_radial_profile --plot_ion --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --drv 0.5

             run make_diagnostic_maps.py --field Par28 --id 1303,1934,2734,2867,300,2903,2906 --plot_AGN_frac --plot_radial_profile --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --drv 0.5 --do_not_correct_pixel

   Afterwards, to make the animation: run /Users/acharyya/Work/astro/ayan_codes/animate_png.py --inpath /Volumes/Elements/acharyya_backup/Work/astro/passage/passage_output/Par028/all_diag_plots_wradprof_snr3.0_onlyseg/ --rootname Par028_*_all_diag_plots_wradprof_snr3.0_onlyseg.png --delay 0.1
'''

from header import *
from util import *
from matplotlib import cm as mpl_cm
import imageio
from get_field_stats import get_crossmatch_with_cosmos

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
        linefit = np.array([ufloat(linefit[0], np.sqrt(linecov[0][0])), ufloat(linefit[1], np.sqrt(linecov[1][1]))])
    except Exception:
        print(f'Could not fit radial profile in this case..')
        linefit = np.array([ufloat(np.nan, np.nan), ufloat(np.nan, np.nan)])

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
def plot_radial_profile(image, ax, args, label=None, ymin=None, ymax=None, hide_xaxis=False, hide_yaxis=False, image_err=None):
    '''
    Plots the average radial profile for a given 2D map in the given axis
    Returns the axis handle
    '''
    print(f'Plotting radial profile of {label}..')
    color = 'darkorange'

    distance_map = get_distance_map(np.shape(image), args)
    try: distance_map = np.ma.masked_where(image.mask, distance_map)
    except AttributeError: distance_map = np.ma.masked_where(False, distance_map)

    # ----making the dataframe before radial profile plot--------------
    xcol, ycol = 'radius', 'data'
    df = pd.DataFrame({xcol: np.ma.compressed(distance_map), ycol: np.ma.compressed(image)})
    if image_err is not None: df[ycol + '_err'] = np.ma.compressed(image_err)
    if len(df[df[ycol + '_err'] > 0]) == 0: image_err = None
    df = df[df[xcol] <= args.radius_max]
    df = df.sort_values(by=xcol)

    # --------processing the dataframe in case voronoi binning has been performed and there are duplicate data values------
    df_vorbinned = pd.DataFrame()
    df_vorbinned[xcol] = df.groupby([ycol], as_index=False).agg([(np.mean)])[xcol]['mean']
    df_vorbinned[ycol] = df.groupby([ycol], as_index=False).agg([(np.mean)])[ycol]
    if image_err is not None: df_vorbinned[ycol+ '_err'] = df.groupby([ycol], as_index=False).agg([(np.mean)])[ycol+ '_err']
    df = df_vorbinned

    # -------proceeding with plotting--------
    ax.scatter(df[xcol], df[ycol], c='grey', s=20 if args.vorbin else 1, alpha=1 if args.vorbin else 0.2)
    if image_err is not None: ax.errorbar(df[xcol], df[ycol], yerr=df[ycol + '_err'], c='grey', fmt='none', lw=2 if args.vorbin else 0.5, alpha=0.2 if args.vorbin else 0.1)

    ax.set_xlim(0, args.radius_max) # kpc
    ax.set_ylim(ymin, ymax)
    ax.set_box_aspect(1)

    if args.vorbin:
        print(f'Not radially binning {label} profile since already voronoi binned')
        # ----------to fit and plot the binned profile--------------
        try:
            linefit, linecov = np.polyfit(df[xcol], df[ycol], 1, cov=True, w=1. / (df[ycol + '_err']) ** 2 if image_err is not None else None)
            y_fitted = np.poly1d(linefit)(df[xcol])
            ax.plot(df[xcol], y_fitted, color=color, lw=1, ls='dashed')
            linefit = np.array([ufloat(linefit[0], np.sqrt(linecov[0][0])), ufloat(linefit[1], np.sqrt(linecov[1][1]))])
        except Exception:
            print(f'Could not fit radial profile in this case..')
            linefit = np.array([ufloat(np.nan, np.nan), ufloat(np.nan, np.nan)])
    else:
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
        try:
            image_err = unp.std_devs(image)
            image = unp.nominal_values(image)
        except:
            image_err = np.ma.masked_where(image.mask, unp.std_devs(image.data))
            image = np.ma.masked_where(image.mask, unp.nominal_values(image.data))

    if not np.ma.isMaskedArray(image): image = np.ma.masked_where(False, image)
    if cmap is None: cmap = 'cividis'
    print(f'Plotting 2D map of {label}..')

    if takelog:
        new_mask = image <= 0
        image[new_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        image_log = unp.log10(unp.uarray(image, image_err))
        image_log = np.ma.masked_where(new_mask | image.mask, image_log)
        image = unp.nominal_values(image_log)
        image_err = unp.std_devs(image_log)

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
        radprof_ax, radprof_fit = plot_radial_profile(image, radprof_ax, args, label=label.split(r'$_{\rm int}')[0], ymin=vmin, ymax=vmax, image_err=image_err)
    else:
        radprof_fit = [np.nan, np.nan] # dummy values for when the fit was not performed

    ax.contour(args.segmentation_map != args.id, levels=0, colors='w' if args.fortalk else 'k', extent=args.extent, linewidths=0.5) # demarcating the segmentation map zone

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
    for id in np.unique(np.ma.compressed(bin_IDs)):
        candidates = map[bin_IDs == id]
        candidates_wo_nan = np.ma.compressed(candidates)
        # print(f'Deb445: data, {id}, {len(candidates)}, {np.nanmean(candidates) if len(candidates_wo_nan) > 0 else np.nan}, {len(candidates_wo_nan)}, {np.mean(candidates_wo_nan)}')  ##
        if len(candidates_wo_nan) > 0: binned_map[bin_IDs == id] = np.mean(candidates_wo_nan)
        else: binned_map[bin_IDs == id] = np.nan

    if map_err is not None:
        binned_map_err = np.zeros(np.shape(map_err))
        for id in np.unique(np.ma.compressed(bin_IDs)):
            candidates = map_err[bin_IDs == id]
            candidates_wo_nan = np.ma.compressed(candidates)
            # print(f'Deb457: err, {id}, {len(candidates)}, {len(candidates_wo_nan)}, {np.sqrt(np.sum(candidates_wo_nan ** 2)) / len(candidates_wo_nan)}')  ##
            binned_map_err[bin_IDs == id] = np.sqrt(np.sum(candidates_wo_nan ** 2)) / len(candidates_wo_nan) # this is the appropriate error propagation for mean() operation (which the flux is undergoing above)

    if map_err is None: return binned_map
    else: return binned_map, binned_map_err

# --------------------------------------------------------------------------------------------------------------------
def make_combined_cmap(cmap1, cmap2, vmin, vmax, vcut):
    '''
    Combine cmap1 and cmap2 into a new colormap such that above vcut cmap1 is used and below vcut cmap2 is used
    Returns new colormap
    '''
    cutoff_frac = (vcut - vmin) / (vmax - vmin)
    if cutoff_frac < 0: cutoff_frac = 0.
    col_above = plt.get_cmap(cmap1)(np.linspace(cutoff_frac, 1, int(256 * (1 - cutoff_frac))))
    col_below = plt.get_cmap(cmap2)(np.linspace(0, cutoff_frac, int(256 * cutoff_frac)))
    colors = np.vstack((col_below, col_above))
    combined_cmap = mplcolors.LinearSegmentedColormap.from_list('combined_cmap', colors)

    return combined_cmap

# --------------------------------------------------------------------------------------------------------------------
def get_voronoi_bin_IDs(map, snr_thresh, plot=False, quiet=True, args=None):
    '''
    Compute the Voronoi bin IDs a given 2D map and corresponding uncertainty and SNR threshold
    Returns the 2D map (of same shape as input map) with just the IDs
    '''
    snr_cut_for_vorbin = 0
    x_size, y_size = np.shape(map)
    map_err = unp.std_devs(map)
    map = unp.nominal_values(map)

    map = np.ma.masked_where(map / map_err < snr_cut_for_vorbin, map)
    map_err = np.ma.masked_where(map / map_err < snr_cut_for_vorbin, map_err)

    if args is not None and args.debug_vorbin:
        fig, axes = plt.subplots(1, 7, figsize=(14, 3))
        fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.1, wspace=0.5)
        cmap = 'viridis'
        fig.text(0.05, 0.98, f'{args.field}: ID {args.id}: vorbin step', fontsize=args.fontsize, c='k', ha='left', va='top')
        plot_2D_map(map, axes[0], args, takelog=False, label=f'input {args.voronoi_line} map', cmap=cmap, hide_yaxis=False)
        plot_2D_map(map_err, axes[1], args, takelog=False, label=f'input {args.voronoi_line} err', cmap=cmap, hide_yaxis=True)

        snr = map / map_err
        combined_cmap = make_combined_cmap(cmap, 'Grays', np.min(snr), np.max(snr), snr_cut_for_vorbin)
        plot_2D_map(snr, axes[2], args, takelog=False, label=f'input {args.voronoi_line} SNR', cmap=combined_cmap, hide_yaxis=True)

    x_coords_grid = np.reshape(np.repeat(np.arange(x_size), y_size), (x_size, y_size))
    y_coords_grid = np.reshape(np.tile(np.arange(y_size), x_size), (x_size, y_size))
    x_coords_grid_masked = np.ma.masked_where(map.mask, x_coords_grid)
    y_coords_grid_masked = np.ma.masked_where(map.mask, y_coords_grid)

    map_array = np.ma.compressed(map)
    map_err_array = np.ma.compressed(map_err)
    x_coords_array = np.ma.compressed(x_coords_grid_masked)
    y_coords_array = np.ma.compressed(y_coords_grid_masked)

    binIDs, _, _, _, _, _, _, _ = voronoi_2d_binning(x_coords_array, y_coords_array, map_array, map_err_array, snr_thresh, plot=False, quiet=quiet, cvt=False, wvt=True)

    interp = NearestNDInterpolator(list(zip(x_coords_array, y_coords_array)), binIDs)
    binID_map = interp(x_coords_grid, y_coords_grid)
    binID_map = np.ma.masked_where(args.segmentation_map != args.id, binID_map)

    if args is not None and args.debug_vorbin:
        plot_2D_map(binID_map, axes[3], args, takelog=False, label='resultant bin IDs', cmap=cmap, hide_yaxis=True)

        map, map_err = bin_2D(map, binID_map, map_err=map_err)
        plot_2D_map(map, axes[4], args, takelog=False, label=f'binned {args.voronoi_line} map', cmap=cmap, hide_yaxis=False)
        plot_2D_map(map_err, axes[5], args, takelog=False, label=f'binned {args.voronoi_line} err', cmap=cmap, hide_yaxis=True)
        plot_2D_map(map / map_err, axes[6], args, takelog=False, label='binned SNR', cmap=cmap, hide_yaxis=True)
        plt.show(block=False)

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
    factor = 1.0

    # ----------deblending flux--------------------
    if not args.do_not_correct_flux:
        if line == 'OIII': # special treatment for OIII 5007 line, in order to account for and remove the OIII 4959 component
            ratio_5007_to_4959 = 2.98 # from grizli source code
            factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
            print(f'Correcting OIII for 4959 component, by factor of {factor:.3f}')
        elif line == 'SII': # special treatment for SII 6717 line, in order to account for and remove the SII 6731 component
            ratio_6717_to_6731 = 1. # from grizli source code
            factor = ratio_6717_to_6731 / (1 + ratio_6717_to_6731)
            print(f'Correcting SII for 6731 component, by factor of {factor:.3f}')
        elif line == 'Ha': # special treatment for Ha line, in order to account for and remove the NII component
            factor = 0.823 # from James et al. 2023?
            print(f'Correcting Ha for NII component, by factor of {factor:.3f}')

    line_map = line_map * factor
    line_map_err = line_map_err * factor

    line_map = trim_image(line_map, args)
    line_map_err = trim_image(line_map_err, args)

    if args.only_seg:
        line_map = cut_by_segment(line_map, args)
        line_map_err = cut_by_segment(line_map_err, args)

    if args.snr_cut is not None:
        snr_map = line_map / line_map_err
        mask = (~np.isfinite(snr_map)) | (snr_map < args.snr_cut)
        line_map = np.ma.masked_where(mask, line_map)
        line_map_err = np.ma.masked_where(mask, line_map_err)

    if args.vorbin and not for_vorbin:
        if args.voronoi_line is None: # No reference emission line specified, so Voronoi IDs need to be computed now
            bin_IDs = get_voronoi_bin_IDs(line_map, line_map_err, args.voronoi_snr, plot=args.debug_vorbin, quiet=not args.debug_vorbin)
        else: # Reference emission line specified for Voronoi binning, so bin IDs have been pre-computed
            bin_IDs = args.voronoi_bin_IDs

        if args.debug_vorbin:
            fig, axes = plt.subplots(1, 7, figsize=(14, 3))
            fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.1, wspace=0.5)
            cmap = 'viridis'
            fig.text(0.05, 0.8, f'{args.field}: ID {args.id}: line {line}', fontsize=args.fontsize, c='k', ha='left', va='top')
            plot_2D_map(np.log10(line_map), axes[0], args, takelog=False, label='map', cmap=cmap, hide_yaxis=True)
            plot_2D_map(np.log10(line_map_err), axes[1], args, takelog=False, label='map err', cmap=cmap, hide_yaxis=True)
            snr = line_map / line_map_err
            combined_cmap = make_combined_cmap(cmap, 'Grays', np.min(snr), np.max(snr), args.snr_cut)
            plot_2D_map(snr, axes[2], args, takelog=False, label='snr', cmap=combined_cmap, hide_yaxis=True)
            plot_2D_map(bin_IDs, axes[3], args, takelog=False, label='bin IDs', cmap=cmap, hide_yaxis=False)

        line_map, line_map_err = bin_2D(line_map, bin_IDs, map_err=line_map_err)

        if args.debug_vorbin:
            plot_2D_map(np.log10(line_map), axes[4], args, takelog=False, label='binned map', cmap=cmap, hide_yaxis=True)
            plot_2D_map(np.log10(line_map_err), axes[5], args, takelog=False, label='binned map err', cmap=cmap, hide_yaxis=True)
            snr = line_map / line_map_err
            combined_cmap = make_combined_cmap(cmap, 'Grays', np.nanmin(snr), np.nanmax(snr), args.snr_cut)
            plot_2D_map(snr, axes[6], args, takelog=False, label='binned snr', cmap=combined_cmap, hide_yaxis=True)
            plt.show(block=False)

    line_map = unp.uarray(line_map, line_map_err)

    # -----------getting the integrated flux value by summing the 2D map-----------------
    line_sum = np.sum(line_map) # ergs/s/cm^2
    print(f'Summed up {line} flux for object {args.id} is {line_sum/1e-17: .3f} x 10^-17 ergs/s/cm^2.')

    # -----------getting the integrated flux value from grizli-----------------
    line_int = full_hdu[0].header[f'FLUX{line_index + 1:03d}'] # ergs/s/cm^2
    line_int_err = full_hdu[0].header[f'ERR{line_index + 1:03d}'] # ergs/s/cm^2
    line_int = line_int * factor
    line_int_err = line_int_err * factor
    line_int = ufloat(line_int, line_int_err)

    # -----------getting the integrated EW value-----------------
    line_index_in_cov = int([item for item in list(full_hdu[2].header.keys()) if full_hdu[0].header[f'FLUX{line_index + 1:03d}'] == full_hdu[2].header[item]][0][5:])
    line_ew = full_hdu[2].header[f'EW50_{line_index_in_cov:03d}'] # rest-frame EW
    line_ew_err = full_hdu[2].header[f'EWHW_{line_index_in_cov:03d}'] # rest-frame EW uncertainty
    line_ew = ufloat(line_ew, line_ew_err)

    # -----------getting the dereddened flux value-----------------
    if dered:
        line_map = get_dereddened_flux(line_map, line_wave, args.EB_V)
        line_int = get_dereddened_flux(line_int, line_wave, args.EB_V)

    if not np.ma.isMaskedArray(line_map): line_map = np.ma.masked_where(False, line_map)

    return line_map, line_wave, line_int, line_ew

# --------------------------------------------------------------------------------------------------------------------
def plot_emission_line_map(line, full_hdu, ax, args, cmap='cividis', EB_V=None, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=False):
    '''
    Plots the emission map for a given line in the given axis
    Returns the axes handle
    '''

    line_map, line_wave, line_int, line_ew = get_emission_line_map(line, full_hdu, args, dered=False)
    ax, _ = plot_2D_map(line_map, ax, args, label=r'%s$_{\rm int}$ = %.1e' % (line, line_int.n), cmap=cmap, vmin=vmin, vmax=vmax, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar)
    ax.text(ax.get_xlim()[0] * 0.88, ax.get_ylim()[0] * 0.88, f'EW = {line_ew:.1f}', c='k', fontsize=args.fontsize, ha='left', va='bottom', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))

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
    lim, label = [-3, -1], 'SFR'
    SFR_map, SFR_int = get_SFR(full_hdu, args)
    ax, SFR_radfit = plot_2D_map(SFR_map, ax, args, label=r'log %s$_{\rm int}$ = %.1f' % (label, SFR_int.n), cmap='viridis', radprof_ax=radprof_ax, vmin=lim[0], vmax=lim[1], hide_yaxis=True)

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
    ax, Te_radfit = plot_2D_map(Te_map, ax, args, label=r'log %s$_{\rm int}$ = %.1e' % (label, Te_int.n), cmap='OrRd_r', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

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
def compute_Z_O3S2(OIII5007_flux, Hbeta_flux, SII6717_flux, Halpha_flux):
    '''
    Calculates and returns the O3S2 metallicity given observed line fluxes
    Conversion factor is from Curti+2019
    '''
    # -----handling masks separately because uncertainty package cannot handle masks---------
    if np.ma.isMaskedArray(OIII5007_flux):
        net_mask = SII6717_flux.mask | OIII5007_flux.mask | Hbeta_flux.mask | Halpha_flux.mask
        OIII5007_flux = OIII5007_flux.data
        Hbeta_flux = Hbeta_flux.data
        SII6717_flux = SII6717_flux.data
        Halpha_flux = Halpha_flux.data
    else:
        net_mask = False

    k = [0.191, -4.292, -2.538, 0.053, 0.332]  # c0-4 parameters from Table 2 of Curti+19 last row

    if hasattr(Hbeta_flux, "__len__"): # if it is an array
        # --------computing the ratio and appropriate errors------------
        new_mask1 = Hbeta_flux == 0
        Hbeta_flux[new_mask1] = -1 # arbitrary fill value to bypass unumpy's inability to handle math domain errors

        new_mask2 = Halpha_flux == 0
        Halpha_flux[new_mask2] = -1 # arbitrary fill value to bypass unumpy's inability to handle math domain errors

        new_mask3 = SII6717_flux == 0
        SII6717_flux[new_mask3] = -1 # arbitrary fill value to bypass unumpy's inability to handle math domain errors

        ratio = (OIII5007_flux / Hbeta_flux) / (SII6717_flux / Halpha_flux)
        ratio = np.ma.masked_where(new_mask1 | new_mask2 | new_mask3, ratio)

        # --------computing the log of the ratio and appropriate errors------------
        new_mask = ratio <= 0
        ratio[new_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        O3S2 = unp.log10(ratio.data)
        O3S2 = np.ma.masked_where(new_mask | ratio.mask, O3S2)

        # --------computing the polynomial and appropriate errors------------
        log_OH = []
        for this_O3S2 in O3S2.data.flatten():
            solution = [item.real for item in np.roots(np.hstack([k[::-1], -1 * unp.nominal_values(this_O3S2)])) if item.imag == 0][0]
            this_log_OH = solution + 8.69 # see Table 1 caption in Curti+19
            log_OH.append(ufloat(this_log_OH, 0.))
        log_OH = np.ma.masked_where(O3S2.mask | net_mask, np.reshape(log_OH, np.shape(O3S2)))

    else: # if it is scalar
        try:
            ratio = (OIII5007_flux / Hbeta_flux) / (SII6717_flux / Halpha_flux)
            O3S2 = unp.log10(ratio)
            solution = [item.real for item in np.roots(np.hstack([k[::-1], -1 * unp.nominal_values(O3S2)])) if item.imag == 0][0]
            log_OH = ufloat(solution + 8.69, 0.)  # see Table 1 caption in Curti+19
        except:
            log_OH = ufloat(np.nan, np.nan)

    return log_OH

# --------------------------------------------------------------------------------------------------------------------
def get_Z_O3S2(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated O3S2 metallicity from a given HDU
    '''
    OIII5007_map, line_wave, OIII5007_int, _ = get_emission_line_map('OIII', full_hdu, args)
    Hbeta_map, line_wave, Hbeta_int, _ = get_emission_line_map('Hb', full_hdu, args)
    SII6717_map, line_wave, SII6717_int, _ = get_emission_line_map('SII', full_hdu, args)
    Halpha_map, line_wave, Halpha_int, _ = get_emission_line_map('Ha', full_hdu, args)

    if not args.do_not_correct_flux:
        # special treatment for OIII 5007 line, in order to account for and ADD the OIII 4959 component back
        ratio_5007_to_4959 = 2.98  # from grizli source code
        factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
        print(f'Re-correcting OIII to include the 4959 component, for computing O3S2 metallicity, by factor of {factor:.3f}')
        OIII5007_map = np.ma.masked_where(OIII5007_map.mask, OIII5007_map.data / factor)
        OIII5007_int = OIII5007_int / factor

        # special treatment for SII 6717 line, in order to account for and ADD the SII 6731 component back
        ratio_6717_to_6731 = 1.  # from grizli source code
        factor = ratio_6717_to_6731 / (1 + ratio_6717_to_6731)
        print(f'Re-correcting SII to include 6731 component, for computing O3S2 metallicity, by factor of {factor:.3f}')
        SII6717_map = np.ma.masked_where(SII6717_map.mask, SII6717_map.data / factor)
        SII6717_int = SII6717_int / factor

    logOH_map = compute_Z_O3S2(OIII5007_map, Hbeta_map, SII6717_map, Halpha_map)
    logOH_int = compute_Z_O3S2(OIII5007_int, Hbeta_int, SII6717_int, Halpha_int)

    return logOH_map, logOH_int

# ----------------------------------------------------------------------------------------------------
def get_nearest(value, array):
    '''
    Finds the nearest neighbour to an item in a given list and returns that element of the list
    '''
    diffs = np.abs(np.array(array) - value)
    min_diff_index = np.where(diffs == np.min(diffs))[0][0]
    nearest_item = array[min_diff_index]

    return nearest_item

# --------------------------------------------------------------------------------------------------------------------
def compute_Z_R23(OII3727_flux, OIII5007_flux, Hbeta_flux, branch='low'):
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

    OH_coeff_dict = {5e6:[-27.0004, 6.03910, -0.327006], 1e7:[-31.2133, 7.15810, -0.399343], 2e7:[-36.0239, 8.44804, -0.483762], 4e7:[-40.9994, 9.78396, -0.571551], 8e7:[-44.7026, 10.8052, -0.640113], 1.5e8:[-46.1589, 11.2557, -0.672731], 3e8:[-45.6075, 11.2074, -0.674460]} # k0-2 parameters for all q from Table 3 of KD02 last row
    q_coeff = [0.0843640, 0.739315, 7.57817] # coeffciients from 3rd column of Table 2 of KD02

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

        logOH_Z94 = np.poly1d([-0.333, -0.207, -0.202, -0.33, 9.625])(R23)  # Eq 8 of KD02

        new_mask2 = OII3727_flux == 0
        OII3727_flux[new_mask2] = -1 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        ratio2 = OIII5007_flux / OII3727_flux

        new_mask2 = ratio2 <= 0
        ratio2[new_mask2] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        y = unp.log10(ratio2)
        y = np.ma.masked_where(new_mask2, y)
        logOH_M91 = 12 + np.poly1d([0.602, 0.767, -4.944])(R23) + y * np.poly1d([-0.331, 0.332, 0.29])(R23)  # Eq 9 of KD02

        logOH_avg = np.mean([logOH_Z94, logOH_M91], axis=0)
        logOH_R23 = np.ones(np.shape(Hbeta_flux)) * ufloat(np.nan, np.nan)

        for i in range(np.shape(Hbeta_flux)[0]):
            for j in range(np.shape(Hbeta_flux)[1]):
                if args.ignore_combined_method or logOH_avg[i][j] < 8.5:
                    # --------computing the polynomial and appropriate errors------------
                    q = 10 ** np.poly1d(q_coeff)(y[i][j])  # coeffciients from 3rd column of Table 2 of KD02
                    if np.isnan(unp.nominal_values(q)):
                        logOH_R23[i][j] = ufloat(np.nan, np.nan)
                        continue
                    nearest_q = get_nearest(unp.nominal_values(q), list(OH_coeff_dict.keys()))
                    k = OH_coeff_dict[nearest_q]
                    p = k[1] ** 2 - 4 * k[2] * (k[0] - R23[i][j])
                    if p >= 0:
                        if branch == 'low': logOH_R23[i][j] = (-k[1] + unp.sqrt(p))/(2 * k[2])
                        elif branch == 'high': logOH_R23[i][j] = (-k[1] - unp.sqrt(p))/(2 * k[2])
                        #print(f'Deb1141: i={i}, j={j}, range=({np.shape(Hbeta_flux)}), logOH_Z94={logOH_Z94[i][j]}, logOH_M91={logOH_M91[i][j]}, and the avg={logOH_avg[i][j]} which is < 8.5, so computing logOH_R23={logOH_R23[i][j]} instead (using q={q:.1e}, which is closest to branch {nearest_q:.1e}).')
                    else:
                        logOH_R23[i][j] = ufloat(np.nan, np.nan)
                        continue
                else:
                    #print(f'Deb1145: i={i}, j={j}, range=({np.shape(Hbeta_flux)}), logOH_Z94={logOH_Z94[i][j]}, logOH_M91={logOH_M91[i][j]}, and the avg={logOH_avg[i][j]} which is > 8.5, so using this as logOH_R23.')
                    logOH_R23[i][j] = logOH_avg[i][j]

        logOH_R23 = np.ma.masked_where(new_mask | R23.mask | net_mask, logOH_R23)

    else: # if it is scalar
        ratio = (OII3727_flux + OIII5007_flux) / Hbeta_flux
        R23 = unp.log10(ratio)
        logOH_Z94 = np.poly1d([-0.333, -0.207, -0.202, -0.33, 9.625])(R23) # Eq 8 of KD02

        ratio = OIII5007_flux / OII3727_flux
        y = unp.log10(ratio)
        logOH_M91 = 12 + np.poly1d([0.602, 0.767, -4.944])(R23) + y * np.poly1d([-0.331, 0.332, 0.29])(R23) # Eq 9 of KD02

        logOH_avg = np.mean([logOH_Z94, logOH_M91])

        if args.ignore_combined_method or logOH_avg < 8.5:
            q = 10 ** np.poly1d(q_coeff)(y) # coeffciients from 3rd column of Table 2 of KD02
            nearest_q = get_nearest(unp.nominal_values(q), list(OH_coeff_dict.keys()))
            k = OH_coeff_dict[nearest_q]
            p = k[1] ** 2 - 4 * k[2] * (k[0] - R23)
            if p >= 0:
                if branch == 'low': logOH_R23 = (-k[1] + unp.sqrt(p)) / (2 * k[2])
                elif branch == 'high': logOH_R23 = (-k[1] - unp.sqrt(p)) / (2 * k[2])
                print(f'Based on the integrated fluxes, logOH_Z94={logOH_Z94}, logOH_M91={logOH_M91}, and the avg={logOH_avg} which is < 8.5, so computing logOH_R23={logOH_R23} instead.')
            else:
                logOH_R23 = ufloat(np.nan, np.nan)
        else:
            print(f'Based on the integrated fluxes, logOH_Z94={logOH_Z94}, logOH_M91={logOH_M91}, and the avg={logOH_avg} which is > 8.5, so using this as logOH_R23.')
            logOH_R23 = logOH_avg

    return logOH_R23

# --------------------------------------------------------------------------------------------------------------------
def get_Z_R23(full_hdu, args, branch='low'):
    '''
    Computes and returns the spatially resolved as well as intregrated R23 metallicity from a given HDU
    '''
    OII3727_map, line_wave, OII3727_int, _ = get_emission_line_map('OII', full_hdu, args)
    OIII5007_map, line_wave, OIII5007_int, _ = get_emission_line_map('OIII', full_hdu, args)
    Hbeta_map, line_wave, Hbeta_int, _ = get_emission_line_map('Hb', full_hdu, args)

    if not args.do_not_correct_flux:
        # special treatment for OIII 5007 line, in order to account for and ADD the OIII 4959 component back
        ratio_5007_to_4959 = 2.98  # from grizli source code
        factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
        print(f'Re-correcting OIII to include the 4959 component, for computing R23 metallicity, by factor of {factor:.3f}')
        OIII5007_map = np.ma.masked_where(OIII5007_map.mask, OIII5007_map.data / factor)
        OIII5007_int = OIII5007_int / factor

    logOH_map = compute_Z_R23(OII3727_map, OIII5007_map, Hbeta_map, branch=branch)
    logOH_int = compute_Z_R23(OII3727_int, OIII5007_int, Hbeta_int, branch=branch)

    return logOH_map, logOH_int

# --------------------------------------------------------------------------------------------------------------------
def plot_Z_R23_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the R23 metallicity map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D metallicity map just produced
    '''
    lim, label = [7, 9], 'Z (R23)'
    logOH_map, logOH_int = get_Z_R23(full_hdu, args)
    ax, logOH_radfit = plot_2D_map(logOH_map, ax, args, takelog=False, label=r'%s$_{\rm int}$ = %.1f' % (label, logOH_int.n), cmap='viridis', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

    return ax, logOH_map, logOH_radfit, logOH_int

# --------------------------------------------------------------------------------------------------------------------
def compute_q_O32(OII3727_flux, OIII5007_flux):
    '''
    Calculates and returns the O32 ionisation parameter given observed line fluxes
    Conversion factor is from Kewley+2002
    '''
    # -----handling masks separately because uncertainty package cannot handle masks---------
    if np.ma.isMaskedArray(OIII5007_flux):
        net_mask = OII3727_flux.mask | OIII5007_flux.mask
        OIII5007_flux = OIII5007_flux.data
        OII3727_flux = OII3727_flux.data
    else:
        net_mask = False

    k = [0.0300472, -0.914212, 10.0581, -36.7948]  # k3-0 parameters from 3rd column of Table 2 of KD02

    if hasattr(OII3727_flux, "__len__"): # if it is an array
        # --------computing the ratio and appropriate errors------------
        new_mask = OII3727_flux == 0
        OII3727_flux[new_mask] = -1 # arbitrary fill value to bypass unumpy's inability to handle math domain errors

        ratio = (OIII5007_flux / OII3727_flux)
        ratio = np.ma.masked_where(new_mask, ratio)

        # --------computing the log of the ratio and appropriate errors------------
        new_mask = ratio <= 0
        ratio[new_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        O32 = unp.log10(ratio.data)
        O32 = np.ma.masked_where(new_mask | ratio.mask, O32)

        # --------computing the polynomial and appropriate errors------------
        logq = np.ones(np.shape(OII3727_flux)) * ufloat(np.nan, np.nan)
        for i in range(np.shape(OII3727_flux)[0]):
            for j in range(np.shape(OII3727_flux)[1]):
                coeff = np.hstack((k[:3], [k[3] - unp.nominal_values(O32[i][j])]))
                try:
                    roots = np.roots(coeff)
                    logq[i][j] = ufloat([item.real for item in roots if item.imag == 0][0], 0)
                except:
                    logq[i][j] = ufloat(np.nan, np.nan)
        logq = np.ma.masked_where(O32.mask | net_mask, np.reshape(logq, np.shape(O32)))

    else: # if it is scalar
        try:
            ratio = (OIII5007_flux / OII3727_flux)
            O32 = unp.log10(ratio)
            coeff = np.hstack((k[:3], [k[3] - unp.nominal_values(O32)]))
            roots = np.roots(coeff)
            logq = ufloat([item.real for item in roots if item.imag == 0][0], 0)
        except:
            logq = ufloat(np.nan, np.nan)

    return logq

# --------------------------------------------------------------------------------------------------------------------
def get_q_O32(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated O32 ionisation parameter from a given HDU
    '''
    OII3727_map, line_wave, OII3727_int, _ = get_emission_line_map('OII', full_hdu, args)
    OIII5007_map, line_wave, OIII5007_int, _ = get_emission_line_map('OIII', full_hdu, args)

    logq_map = compute_q_O32(OII3727_map, OIII5007_map)
    logq_int = compute_q_O32(OII3727_int, OIII5007_int)

    return logq_map, logq_int

# --------------------------------------------------------------------------------------------------------------------
def plot_q_O32_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the O32 iojisation parameter map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D metallicity map just produced
    '''
    lim, label = [7, 8], 'log(q) (O32)'
    logq_map, logq_int = get_q_O32(full_hdu, args)
    ax, logq_radfit = plot_2D_map(logq_map, ax, args, takelog=False, label=r'%s$_{\rm int}$ = %.1f' % (label, logq_int.n), cmap='cividis', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

    return ax, logq_map, logq_radfit, logq_int

# --------------------------------------------------------------------------------------------------------------------
def plot_this(data, ax):
    '''
    Tiny plotting routine for testing purposes
    Returns axis handle
    '''
    if data.dtype == 'object': data = unp.nominal_values(data)
    ax.imshow(np.log10(data))
    ax.set_title(np.shape(data))
    return ax

# --------------------------------------------------------------------------------------------------------------------
def get_cutout(filename, pos, size, target_header, args, plot_test_axes=None):
    '''
    Return a cutout from a given the filename of a fits image, around a given position within a given extent
    Optionally trims the cutout as per the segmentation map, and/or Voronoi binned
    Returns the 2D cutout as a 2D array
    '''
    data = fits.open(filename)

    image = data[0].data
    source_header = data[0].header
    wcs = pywcs.WCS(source_header)

    cutout = Cutout2D(image, pos, size, wcs=wcs)
    cutout_data = cutout.data * source_header['PHOTFNU'] * 1e6
    cutout_header = cutout.wcs.to_header()
    source_header.update(cutout_header)  # important step to update the header of the full mosaic to that of just the cut out

    cutout_data_hdu = fits.ImageHDU(cutout_data, header=source_header)
    cutout_data_rebinned, _ = reproject_interp(cutout_data_hdu, target_header)  # 125 x 125 pixels
    cutout_data_rebinned_trimmed = trim_image(cutout_data_rebinned, args)  # 50 x 50 pixels

    # --------test plots-------------------
    if plot_test_axes is not None:
        plot_test_axes[2] = plot_this(image, plot_test_axes[2])
        plot_test_axes[3] = plot_this(cutout_data, plot_test_axes[3])
        plot_test_axes[4] = plot_this(cutout_data_rebinned, plot_test_axes[4])
        plot_test_axes[5] = plot_this(cutout_data_rebinned_trimmed, plot_test_axes[5])
        plt.show(block=False)

    if args.only_seg:
        cutout_data_rebinned_trimmed = cut_by_segment(cutout_data_rebinned_trimmed, args)

    if args.vorbin:
        cutout_data_rebinned_trimmed = bin_2D(cutout_data_rebinned_trimmed, args.voronoi_bin_IDs)

    return cutout_data_rebinned_trimmed

# --------------------------------------------------------------------------------------------------------------------
def get_direct_image_per_filter(full_hdu, filter, target_header, args, plot_test_axes=None):
    '''
    Retrieve the direct image for a given filter from the HDU
    Returns the 2D image along with uncertainty
    '''
    # ---------detremining extents for cut out-------------
    size = args.arcsec_limit * 3 * u.arcsec # 3 arcsec so that it leads to shape 125 x 125 pixels
    pos = np.array([full_hdu[0].header['RA'], full_hdu[0].header['DEC']])
    pos = SkyCoord(*(pos * u.deg), frame='fk5')

    # ---------making cut outs-------------
    drizzled_image_filename = glob.glob(str(args.input_dir / args.drv / args.field / f'Products/{args.field}*{filter.lower()}_drz_sci.fits'))[0]
    direct_map = get_cutout(drizzled_image_filename, pos, size, target_header, args, plot_test_axes=plot_test_axes)
    direct_map_wht = get_cutout(drizzled_image_filename.replace('sci', 'wht'), pos, size, target_header, args)

    # -------------pixel offset----------------
    if not args.do_not_correct_pixel:
        ndelta_xpix, ndelta_ypix = 2, 0
        print(f'Correcting emission lines for pixel offset by {ndelta_xpix} on x and {ndelta_ypix} on y')
        direct_map = np.roll(direct_map, ndelta_xpix, axis=1)
        direct_map_wht = np.roll(direct_map_wht, ndelta_xpix, axis=1)
        direct_map = np.roll(direct_map, ndelta_ypix, axis=0)
        direct_map_wht = np.roll(direct_map_wht, ndelta_ypix, axis=0)

    # ---------computing uncertainty-------------
    direct_map_err = 1 / np.sqrt(direct_map_wht)
    direct_map = unp.uarray(direct_map, direct_map_err)

    return direct_map

# --------------------------------------------------------------------------------------------------------------------
def plot_direct_images_all_filters(full_hdu, args):
    '''
    Plots direct images from full_hdu vs corresponding cut outs from full field drizzled images
    '''
    filters_arr = args.filters

    # -------getting the full.fits image------------
    ext = 5
    direct_map = full_hdu[ext].data
    direct_map_trimmed = trim_image(direct_map, args)
    target_header = full_hdu[ext].header
    direct_map_err = 1. / np.sqrt(full_hdu[ext + 1].data)
    direct_map_err_trimmed = trim_image(direct_map_err, args)

    if args.only_seg:
        direct_map_trimmed = cut_by_segment(direct_map_trimmed, args)
        direct_map_err_trimmed = cut_by_segment(direct_map_err_trimmed, args)

    if args.vorbin:
        direct_map_trimmed, direct_map_err_trimmed = bin_2D(direct_map_trimmed, args.voronoi_bin_IDs, map_err=direct_map_err_trimmed)

    direct_map_trimmed = unp.uarray(direct_map_trimmed, direct_map_err_trimmed)
    if not np.ma.isMaskedArray(direct_map_trimmed): direct_map_trimmed = np.ma.masked_where(False, direct_map_trimmed)

    # ------------test plots---------------
    if args.test_cutout:
        fig, plot_test_axes = plt.subplots(1, 6, figsize=(14,3))
        plot_test_axes[0] = plot_this(direct_map, plot_test_axes[0])
        plot_test_axes[1] = plot_this(direct_map_trimmed, plot_test_axes[1])

    # ----------main plots-------------------
    if args.plot_direct_filters:
        nrows = 2
        if args.plot_radial_profiles: nrows += 2

        fig_size_dict = {2: [14, 7, 0.02, 0.98, 0.07, 0.95, 0.05, 0.2], \
                         4: [9, 7, 0.02, 0.98, 0.07, 0.95, 0.05, 0.3]}  # figsize_w, figsize_h, l, r, b, t, ws, hs

        fig, axes = plt.subplots(nrows, 1 + len(filters_arr), figsize=(fig_size_dict[nrows][0], fig_size_dict[nrows][1]))

        if nrows == 2:
            axes_direct = axes[0, :]
            axes_ratio = axes[1, :]
            radprof_axes_direct = np.tile(None, np.shape(axes_direct))
            radprof_axes_ratio = np.tile(None, np.shape(axes_ratio))
        else:
            axes_direct = axes[0, :]
            radprof_axes_direct = axes[1, :]
            axes_ratio = axes[2, :]
            radprof_axes_ratio = axes[3, :]

        fig.subplots_adjust(left=fig_size_dict[nrows][2], right=fig_size_dict[nrows][3], bottom=fig_size_dict[nrows][4], top=fig_size_dict[nrows][5], wspace=fig_size_dict[nrows][6], hspace=fig_size_dict[nrows][7])

        #fig, axes = plt.subplots(2, 1 + len(filters_arr), figsize=(8 + 2 * len(filters_arr), 6))
        #fig.subplots_adjust(left=0.06, right=0.97, bottom=0.12, top=0.9, wspace=0.2, hspace=0.1)
        limits = [-5, -2]
        cmap = 'viridis'

        axes_direct[0], _ = plot_2D_map(direct_map_trimmed, axes_direct[0], args, label='From *full.fits', vmin=limits[0], vmax=limits[1], cmap=cmap, radprof_ax=radprof_axes_direct[0])

    # ---------getting the direct image from cut out-------------
    filter_map_dict = {}
    for index, filter in enumerate(filters_arr):
        filter_map = get_direct_image_per_filter(full_hdu, filter, target_header, args, plot_test_axes=plot_test_axes if args.test_cutout and index == 0 else None)
        filter_map_dict.update({filter: filter_map})
        if args.plot_direct_filters:
            axes_direct[index + 1], _ = plot_2D_map(filter_map, axes_direct[index + 1], args, label=f'{filter}', vmin=limits[0], vmax=limits[1], cmap=cmap, hide_yaxis=True, hide_xaxis=True, radprof_ax=radprof_axes_direct[index + 1])

    # ---------plotting ratios of direct_images----------------------
    if args.plot_direct_filters:
        ratio_list = list(itertools.combinations(filters_arr, 2))
        for index in range(len(filters_arr)):
            num = ratio_list[index][0]
            den = ratio_list[index][1]
            num_map = filter_map_dict[num]
            den_map = filter_map_dict[den]

            try:
                new_mask = unp.nominal_values(den_map.data) == 0
                den_map[new_mask] = -1  # arbitrary fill value to bypass unumpy's inability to handle math domain errors
                ratio_map = num_map.data / den_map.data
                ratio_map = np.ma.masked_where(num_map.mask | den_map.mask | new_mask, ratio_map)
            except:
                new_mask = unp.nominal_values(den_map) == 0
                den_map[new_mask] = -1  # arbitrary fill value to bypass unumpy's inability to handle math domain errors
                ratio_map = num_map / den_map
                ratio_map = np.ma.masked_where(new_mask, ratio_map)

            if args.only_seg: ratio_map = cut_by_segment(ratio_map, args)
            axes_ratio[index + 1], _ = plot_2D_map(ratio_map, axes_ratio[index + 1], args, label=f'{num}/{den}', cmap=cmap, hide_yaxis=True, vmin=-1, vmax=1, radprof_ax=radprof_axes_ratio[index + 1])
        axes_ratio[0].set_visible(False)
        if args.plot_radial_profiles: radprof_axes_ratio[0].set_visible(False)

    plt.show(block=False)

    return fig

# --------------------------------------------------------------------------------------------------------------------
def plot_starburst_map(full_hdu, axes, args, radprof_axes=None, vorbin_Axes=None, snr_axes=None):
    '''
    Plots the Ha map, direct F115W map and their ratio (starbursty-ness map) in a given axis handle
    Returns the axis handle and the ratio map just produced
    '''
    # ---------getting the Ha map-------------
    ha_map, _, _, _ = get_emission_line_map('Ha', full_hdu, args, dered=True)

    # ---------getting the direct image-------------
    ext = 5
    filter = 'F115W'
    target_header = full_hdu[ext].header
    direct_map = get_direct_image_per_filter(full_hdu, filter, target_header, args)

    # ---------getting the ratio and seg maps-------------
    new_mask = unp.nominal_values(direct_map.data) == 0
    direct_map[new_mask] = 1 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
    direct_map = np.ma.masked_where(new_mask, direct_map)

    ratio_map = ha_map.data / direct_map.data
    ratio_map = np.ma.masked_where(ha_map.mask | direct_map.mask | new_mask, ratio_map)

    # --------making arrays for subplots-------------
    maps_dict = {'direct':direct_map, 'ha':ha_map, 'ratio':ratio_map}
    labels_dict = {'direct':filter, 'ha':r'H$\alpha$', 'ratio':r'H$\alpha$/' + filter}
    lims_dict = {'direct':[-4.5, -2], 'ha':[-20, -18], 'ratio':[-17, -14]}
    if len(np.atleast_1d(axes)) > 1: sequence_to_plot = ['direct', 'ha', 'ratio']
    else: sequence_to_plot = ['ratio']

    # ---------plotting-------------
    starburst_radfit = []
    for index, ax in enumerate(np.atleast_1d(axes)):
        map_err = np.ma.masked_where(maps_dict[sequence_to_plot[index]].mask, unp.std_devs(maps_dict[sequence_to_plot[index]].data))
        ax, radprof_fit = plot_2D_map(maps_dict[sequence_to_plot[index]], ax, args, label=labels_dict[sequence_to_plot[index]], cmap='viridis', vmin=lims_dict[sequence_to_plot[index]][0], vmax=lims_dict[sequence_to_plot[index]][1], radprof_ax=np.atleast_1d(radprof_axes)[index] if args.plot_radial_profiles else None, vorbin_ax=np.atleast_1d(vorbin_axes)[index] if args.plot_vorbin else None, snr_ax=np.atleast_1d(snr_axes)[index] if args.plot_snr else None, image_err=map_err if args.plot_snr else None, hide_yaxis=True)
        starburst_radfit.append(radprof_fit)

    return axes, ratio_map, starburst_radfit

# --------------------------------------------------------------------------------------------------------------------
def plot_starburst_figure(full_hdu, args):
    '''
    Plots the Ha map, direct F115W map and their ratio (starbursty-ness map) in a new figure
    Returns the figure handle and the ratio map just produced
    '''
    nrows = 1
    if args.plot_vorbin: nrows += 1
    if args.plot_snr: nrows += 1
    if args.plot_radial_profiles: nrows += 1

    fig_size_dict = {1: [14, 4, 0.05, 0.97, 0.15, 0.9, 0.2, 0.], \
                     2: [14, 7, 0.02, 0.98, 0.07, 0.95, 0.05, 0.2], \
                     3: [9, 7, 0.02, 0.98, 0.07, 0.95, 0.05, 0.3], \
                     4: [9, 7, 0.02, 0.98, 0.07, 0.95, 0.05, 0.3]}  # figsize_w, figsize_h, l, r, b, t, ws, hs

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

    axes, starburst_map, starburst_radfit = plot_starburst_map(full_hdu, axes, args, radprof_axes=radprof_axes, vorbin_axes=vorbin_axes, snr_axes=snr_axes)

    return fig, starburst_map, starburst_radfit

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_map(full_hdu, args):
    '''
    Plots the metallicity map, and optionally radial profile, in a new figure
    Returns the figure handle and the metallicity map just produced
    '''
    ncols = 2
    if args.plot_radial_profiles: ncols += 1
    if args.plot_ionisation_parameter: ncols += 1

    fig_size_dict = {2: [10, 4, 0.05, 0.97, 0.15, 0.9, 0.2, 0.], 3: [14, 4, 0.05, 0.97, 0.15, 0.9, 0.2, 0.], 4: [14, 3.5, 0.06, 0.99, 0.15, 0.95, 0.25, 0.]} # figsize_w, figsize_h, l, b, r, t, ws, hs

    fig, axes = plt.subplots(1, ncols, figsize=(fig_size_dict[ncols][0], fig_size_dict[ncols][1]))
    if ncols > 3:
        ip_ax = axes[0]
        axes = axes[1:]
    if ncols > 2:
        radprof_ax = axes[-1]
        axes = axes[:-1]
    else:
        radprof_ax = None
    fig.subplots_adjust(left=fig_size_dict[ncols][2], right=fig_size_dict[ncols][3], bottom=fig_size_dict[ncols][4], top=fig_size_dict[ncols][5], wspace=fig_size_dict[ncols][6], hspace=fig_size_dict[ncols][7])

    if all([line in args.available_lines for line in ['OIII', 'OII', 'Hb']]):
        # --------deriving the metallicity map-------------
        logOH_map, logOH_int = get_Z_R23(full_hdu, args, branch=args.Zbranch)
        #logOH_map, logOH_int = get_Z_O3S2(full_hdu, args)
        if args.plot_ionisation_parameter: logq_map, logq_int = get_q_O32(full_hdu, args)

        # ---------plotting-------------
        lim, label = [7.5, 9.0], 'log(O/H) (R23)'
        axes[0], logOH_radfit = plot_2D_map(logOH_map, axes[0], args, takelog=False, label=r'%s$_{\rm int}$ = %.1f $\pm$ %.1f' % (label, logOH_int.n, logOH_int.s), cmap='viridis', radprof_ax=radprof_ax, hide_yaxis=True if args.plot_ionisation_parameter else False, vmin=lim[0], vmax=lim[1])

        logOH_map_err = np.ma.masked_where(logOH_map.mask, unp.std_devs(logOH_map.data))
        logOH_map_snr = np.ma.masked_where(logOH_map.mask, unp.nominal_values(10 ** logOH_map.data)) / np.ma.masked_where(logOH_map.mask, unp.std_devs(10 ** logOH_map.data))
        axes[1], _ = plot_2D_map(logOH_map_snr, axes[1], args, takelog=False, hide_yaxis=True, label=r'%s SNR' % (label), cmap='cividis', vmin=0, vmax=6)
        if args.plot_ionisation_parameter: ip_ax, _ = plot_2D_map(logq_map, ip_ax, args, takelog=False, hide_yaxis=False, label=r'log q$_{\rm int}$ = %.1f $\pm$ %.1f' % (logq_int.n, logq_int.s), cmap='viridis', vmin=6.5, vmax=8.5)

    else:
        print(f'Not all lines out of OIII, OII and Hb are available, so cannot compute R23 metallicity')
        logOH_map. logOH_radfit = None, None

    return fig, logOH_map, logOH_radfit

# --------------------------------------------------------------------------------------------------------------------
def plot_BPT(full_hdu, ax, args, cmap='viridis', ax_inset=None):
    '''
    Plots spatially resolved BPT diagram based on fluxes from grizli, on an existing axis
    Then overplots theoretical lines
    Returns axis handle and the handle of the spatially reslved scatter plot
    '''
    print(f'Plotting BPT diagram..')
    if args.plot_separately:
        fig_indiv, ax_indiv = plt.subplots(1, figsize=(8, 6))

    # -----------getting the fluxes------------------
    OIII_map, OIII_wave, OIII_int, _ = get_emission_line_map('OIII', full_hdu, args)
    Hbeta_map, Hbeta_wave, Hbeta_int, _ = get_emission_line_map('Hb', full_hdu, args)
    SII_map, SII_wave, SII_int, _ = get_emission_line_map('SII', full_hdu, args)
    Halpha_map, Halpha_wave, Halpha_int, _ = get_emission_line_map('Ha', full_hdu, args)

    # ----------deblending flux--------------------
    if not args.do_not_correct_flux:
        # special treatment for SII 6717 line, in order to account for and ADD the SII 6731 component back
        ratio_6717_to_6731 = 1. # from grizli source code
        factor = ratio_6717_to_6731 / (1 + ratio_6717_to_6731)
        SII_map = np.ma.masked_where(SII_map.mask, SII_map.data / factor)
        SII_int = SII_int / factor
        print(f'Re-correcting SII to include the 6731 component, for computing BPT, by factor of {factor:.3f}')

    def func(x):
        return 1.3 + 0.72 / (x - 0.32)  # Eq 6 of K01

    try:
        # -----------integrated-----------------------
        color = mpl_cm.get_cmap(cmap)(0.5)

        y_ratio = unp.log10(OIII_int / Hbeta_int)
        x_ratio = unp.log10(SII_int / Halpha_int)

        p = ax.scatter(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), c=color, marker='o', s=200, lw=2, edgecolor='w' if args.fortalk else 'k', zorder=10)
        ax.errorbar(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), xerr=unp.std_devs(x_ratio), yerr=unp.std_devs(y_ratio), c=color, fmt='none', lw=2)

        if args.plot_separately:
            p = ax_indiv.scatter(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), c=color, marker='o', s=200, lw=2, edgecolor='w' if args.fortalk else 'k', zorder=10)
            ax_indiv.errorbar(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), xerr=unp.std_devs(x_ratio), yerr=unp.std_devs(y_ratio), c=color, fmt='none', lw=2)

    except ValueError:
        print(f'Galaxy {args.id} in {args.field} has a negative integrated flux in one of the following lines, hence skipping this.')
        print(f'OIII = {OIII_int}\nHb = {Hbeta_int}\nSII = {SII_int}\nHa = {Halpha_int}\n')
        scatter_plot_handle = None
        pass

    try:
        # -----------spatially_resolved-----------------------
        distance_map = get_distance_map(np.shape(OIII_map), args)
        distance_map = np.ma.masked_where(False, distance_map)
        if args.vorbin: distance_map = bin_2D(distance_map, args.voronoi_bin_IDs)

        y_ratio_mask = (OIII_map.data < 0) | (Hbeta_map.data <= 0) | (~np.isfinite(unp.nominal_values(OIII_map.data))) | (~np.isfinite(unp.nominal_values(Hbeta_map.data)))
        OIII_map[y_ratio_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        Hbeta_map[y_ratio_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        y_ratio = unp.log10(OIII_map.data / Hbeta_map.data)
        y_ratio = np.ma.masked_where(y_ratio_mask | OIII_map.mask | Hbeta_map.mask, y_ratio)

        x_ratio_mask = (SII_map.data < 0) | (Halpha_map.data <= 0)| (~np.isfinite(unp.nominal_values(SII_map.data))) | (~np.isfinite(unp.nominal_values(Halpha_map.data)))
        SII_map[x_ratio_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        Halpha_map[x_ratio_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        x_ratio = unp.log10(SII_map.data / Halpha_map.data)
        x_ratio = np.ma.masked_where(x_ratio_mask | SII_map.mask | Halpha_map.mask, x_ratio)

        net_mask = y_ratio.mask | x_ratio.mask

        sign_map = (unp.nominal_values(y_ratio.data) > func(unp.nominal_values(x_ratio.data))).astype(int)
        sign_map[sign_map == 0] = -1
        distance_from_K01_map = sign_map * get_distance_from_Kewley2001(unp.nominal_values(x_ratio.data), unp.nominal_values(y_ratio.data), args, x_num='SII')
        distance_from_K01_map = np.ma.masked_where(net_mask, distance_from_K01_map)

        x_ratio = np.ma.compressed(np.ma.masked_where(net_mask, x_ratio))
        y_ratio = np.ma.compressed(np.ma.masked_where(net_mask, y_ratio))
        distance_map = np.ma.compressed(np.ma.masked_where(net_mask, distance_map))
        distance_from_K01_arr = np.ma.compressed(distance_from_K01_map)
        dist_lim = max(np.abs(np.max(distance_from_K01_arr)), np.abs(np.min(distance_from_K01_arr)))

        df = pd.DataFrame({'log_sii/ha': unp.nominal_values(x_ratio).flatten(), 'log_sii/ha_err': unp.std_devs(x_ratio).flatten(), 'log_oiii/hb': unp.nominal_values(y_ratio).flatten(), 'log_oiii/hb_err':  unp.std_devs(y_ratio).flatten(), 'distance': distance_map.flatten(), 'distance_from_K01': distance_from_K01_arr.flatten()})
        df = df.sort_values(by='distance')
        df = df.drop_duplicates().reset_index(drop=True)

        scatter_plot_handle = ax.scatter(df['log_sii/ha'], df['log_oiii/hb'], c=df[args.colorcol], marker='o', s=50, lw=0, cmap=cmap, alpha=0.8, vmin=0 if args.colorcol == 'distance' else -dist_lim if args.colorcol =='distance_from_K01' else None, vmax=6 if args.colorcol == 'distance' else dist_lim if args.colorcol =='distance_from_K01' else None)
        ax.errorbar(df['log_sii/ha'], df['log_oiii/hb'], xerr=df['log_sii/ha_err'], yerr=df['log_oiii/hb_err'], c='gray', fmt='none', lw=0.5, alpha=0.1)

        if args.plot_AGN_frac and not args.plot_separately:
            if ax_inset is None: ax_inset = ax.inset_axes([0.55, 0.75, 0.3, 0.3])
            plot_2D_map(distance_from_K01_map, ax_inset, args, takelog=False, label='dist from K01', cmap=args.diverging_cmap, vmin=-dist_lim, vmax=dist_lim, hide_yaxis=False)

        if args.plot_separately:
            scatter_plot_handle_indiv = ax_indiv.scatter(df['log_sii/ha'], df['log_oiii/hb'], c=df[args.colorcol], marker='o', s=50, lw=0, cmap=cmap, alpha=0.8, vmin=0 if args.colorcol == 'distance' else -dist_lim if args.colorcol =='distance_from_K01' else None, vmax=6 if args.colorcol == 'distance' else dist_lim if args.colorcol =='distance_from_K01' else None)
            ax_indiv.errorbar(df['log_sii/ha'], df['log_oiii/hb'], xerr=df['log_sii/ha_err'], yerr=df['log_oiii/hb_err'], c='gray', fmt='none', lw=0.5, alpha=0.1)

            if args.plot_AGN_frac:
                if ax_inset is None: ax_inset = ax_indiv.inset_axes([0.55, 0.75, 0.3, 0.3])
                plot_2D_map(distance_from_K01_map, ax_inset, args, takelog=False, label='dist from K01', cmap=args.diverging_cmap, vmin=-dist_lim, vmax=dist_lim)

    except ValueError:
        print(f'Galaxy {args.id} in {args.field} has some negative spatially resolved fluxes, hence skipping this object.')
        scatter_plot_handle = None
        pass

    if args.plot_separately:
        # ---------annotate axes-------
        cbar = plt.colorbar(scatter_plot_handle_indiv)
        cbar.set_label('Distance (kpc)' if args.colorcol == 'distance' else 'Distance from K01' if args.colorcol == 'distance_from_K01' else '')

        ax_indiv.set_xlim(-2, 0.3)
        ax_indiv.set_ylim(-1, 2)
        ax_indiv.set_xlabel(f'log (SII 6717+31/Halpha)')
        ax_indiv.set_ylabel(f'log (OIII 5007/Hbeta)')

        # ---------adding literature lines from Kewley+2001 (https://iopscience.iop.org/article/10.1086/321545)----------
        x = np.linspace(ax_indiv.get_xlim()[0], ax_indiv.get_xlim()[1], 100)
        y = func(x)  # Eq 6 of K01
        ax_indiv.plot(x, y, c='w' if args.fortalk else 'k', ls='dashed', lw=2, label='Kewley+2001')
        plt.legend()
        fig_indiv.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.95)
        fig_indiv.text(0.15, 0.9, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')

        # -----------save figure----------------
        figname = fig_dir / f'{args.field}_{args.id:05d}_BPT{snr_text}{only_seg_text}{vorbin_text}.png'
        fig_indiv.savefig(figname, transparent=args.fortalk)
        print(f'Saved figure at {figname}')
        plt.show(block=False)

    # ---------annotate axes-------
    cbar = plt.colorbar(scatter_plot_handle)
    cbar.set_label('Distance (kpc)' if args.colorcol == 'distance' else 'Distance from K01' if args.colorcol == 'distance_from_K01' else '')

    ax.set_xlim(-2, 0.3)
    ax.set_ylim(-1, 2)
    ax.set_xlabel(f'log (SII 6717/Halpha)')
    ax.set_ylabel(f'log (OIII 5007/Hbeta)')

    # ---------adding literature lines from Kewley+2001 (https://iopscience.iop.org/article/10.1086/321545)----------
    x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
    y = 1.3 + 0.72 / (x - 0.32)  # Eq 6 of K01
    ax.plot(x, y, c='w' if args.fortalk else 'k', ls='dashed', lw=2, label='Kewley+2001')

    return ax, distance_from_K01_map

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.plot_slope_vs_mass: args.plot_starburst, args.plot_radial_profiles = True, True
    if args.colorcol == 'ez_z_phot': args.colorcol = 'distance'
    args.diverging_cmap = get_custom_cmap(args.diverging_cmap)

    # ---------determining filename suffixes-------------------------------
    radial_plot_text = '_wradprof' if args.plot_radial_profiles else ''
    snr_text = f'_snr{args.snr_cut}' if args.snr_cut is not None else ''
    only_seg_text = '_onlyseg' if args.only_seg else ''
    pixscale_text = '' if args.pixscale == 0.04 else f'_{args.pixscale}arcsec_pix'
    vorbin_text = '' if not args.vorbin else f'_vorbin_at_{args.voronoi_line}_SNR_{args.voronoi_snr}'
    description_text = f'all_diag_plots{radial_plot_text}{snr_text}{only_seg_text}'

    # ---------determining list of fields----------------
    if args.do_all_fields:
        field_list = [os.path.split(item[:-1])[1] for item in glob.glob(str(args.input_dir / args.drv / 'Par*') + '/')]
        field_list.sort(key=natural_keys)
    else:
        field_list = args.field_arr

    # --------loop over all fields------------------
    for index2, field in enumerate(field_list):
        start_time2 = datetime.now()
        args.field = f'Par{int(field[3:]):03}' if len(field) < 6 else field
        print(f'\n\nCommencing field {args.field} which is {index2 + 1} of {len(field_list)}..')

        product_dir = args.input_dir / args.drv / args.field / 'Products'
        output_dir = args.output_dir / args.field
        if args.re_extract: output_dir = output_dir / 're_extracted'
        output_dir.mkdir(parents=True, exist_ok=True)
        outfilename = output_dir / f'{args.field}_all_diag_results.csv'

        if os.path.exists(outfilename) and not args.clobber and args.write_file:
            print(f'result file for {args.field} already exists as {outfilename}, so skipping this field.')
            continue

        # ---------prep the catalog file--------------------
        catalog_file = product_dir / f'{args.field}_photcat.fits'
        if not os.path.exists(catalog_file):
            print(f'photcat file for {args.field} does not exist, please download from gdrive or run data reduction. Until then, skipping this field.')
            continue

        catalog = GTable.read(catalog_file)

        # --------determine which objects to loop over----------
        if args.do_all_obj:
            if args.re_extract: args.id_arr = ids_to_re_extract_dict[args.field]
            else: args.id_arr = catalog['id']
        else:
            args.id_arr = args.id

        if args.start_id: args.id_arr = args.id_arr[args.start_id - 1:]
        if len(args.id_arr) == 1: args.plot_separately = False
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
        all_lines_to_plot = ['OII', 'Hb', 'OIII', 'Ha', 'SII'] # from blue to red
        all_lines_to_save = args.line_list
        measured_quantities_to_plot = ['BPT', 'logq', 'logOH_R23', 'SFR', 'F115W/Ha'] # ['EB_V', 'SFR', 'Te', 'logOH_Te', 'logOH_R23']

        if args.write_file:
            basic_cols = ['field', 'objid', 'ra', 'dec', 'redshift']
            flag_cols = ['radfit_extent_kpc', 'snr_cut', 'flag_only_seg', 'flag_vorbin', 'vor_snr', 'vor_line']
            cols_in_df = np.hstack([basic_cols, flag_cols, np.hstack([[item + '_int', item + '_int_u', item + '_EW', item + '_EW_u'] for item in all_lines_to_save]), np.hstack([[item + '_int', item + '_cen', item + '_slope'] for item in measured_quantities_to_plot])])
            df = pd.DataFrame(columns=cols_in_df)

            # -------checking if about to write the same columns-----------
            if os.path.exists(str(outfilename).replace('.csv', '.txt')): os.remove(str(outfilename).replace('.csv', '.txt'))
            if os.path.exists(outfilename):
                if args.clobber:
                    os.remove(outfilename)
                    print(f'Deleting existing {outfilename}')
                else:
                    existing_df = pd.read_table(outfilename, delim_whitespace=True)
                    existing_cols = existing_df.columns
                    if set(existing_cols) != set(cols_in_df):
                        new_cols = set(cols_in_df) - set(existing_cols)
                        print(
                        f'Existing dataframe at {outfilename} has a different set of columns (the difference being {new_cols}) than currently trying to '
                        f'write it in. Either delete/double check the existing file, OR use --clobber to overwrite the '
                        f'existing file. Aborting.')
                        sys.exit()

        # ---------plotting fit vs stellar mass-----------------------------
        if args.plot_slope_vs_mass:
            df_starburst_slope = pd.DataFrame(columns=['field', 'objid', 'ra', 'dec', 'redshift', 'F115W_slope', 'Ha_slope', 'Ha/F115W_slope'])

        # ---------plotting spatially resolved BPT-----------------------------
        if args.plot_BPT:
            fig, ax = plt.subplots(1, figsize=(8, 6))
            cmap_arr = ['Reds_r', 'Greens_r', 'Purples_r', 'Greys_r', 'Oranges_r', 'Blues_r', 'YlGnBu_r', 'BuPu_r', 'GnBu_r', 'spring']

        # ---------lotting metallicity profiles and gradients----------------------
        if args.plot_metallicity:
            df_logOH_radfit = pd.DataFrame(columns=['field', 'objid', 'logOH_slope', 'logOH_slope_u', 'logOH_cen', 'logOH_cen_u'])

        # ------------looping over the provided object IDs-----------------------
        for index, args.id in enumerate(args.id_arr):
            start_time3 = datetime.now()
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
                args.voronoi_bin_IDs = get_voronoi_bin_IDs(line_map, args.voronoi_snr, plot=args.debug_vorbin, quiet=not args.debug_vorbin, args=args)
                if args.debug_vorbin:
                    print(f'Running in --debug_vorbin mode, hence not proceeding further.')
                    #continue

            # ---------------radial profile stuff---------------
            if args.plot_radial_profiles:
                seg_map = full_hdu['SEG'].data
                distance_map = get_distance_map(np.shape(seg_map), args)
                distance_map = np.ma.compressed(np.ma.masked_where(seg_map != args.id, distance_map))
                args.radius_max = np.max(distance_map)
            else:
                args.radius_max = np.nan

            # ---------------dust value---------------
            if all([line in args.available_lines for line in ['Ha', 'Hb']]) and not args.test_cutout:
                try: _, args.EB_V = get_EB_V(full_hdu, args, verbose=True)
                except: args.EB_V = 0.

            # ---------initialising the starburst figure------------------------------
            if args.test_cutout or args.plot_direct_filters:
                fig = plot_direct_images_all_filters(full_hdu, args)

                # ---------decorating and saving the figure------------------------------
                fig.text(0.05, 0.98, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
                main_text = 'direct_filter_images' if args.plot_direct_filters else 'cutout_tests'
                figname = fig_dir / f'{args.field}_{args.id:05d}_{main_text}{only_seg_text}{vorbin_text}.png'

            # ---------plotting spatially resolved BPT-----------------------------
            elif args.plot_BPT:
                ax, distance_from_K01_map = plot_BPT(full_hdu, ax, args, cmap=args.diverging_cmap if args.colorcol == 'distance_from_K01' else cmap_arr[index])

            # ---------initialising the starburst figure------------------------------
            elif args.plot_starburst:
                fig, ratio_map, starburst_radfit = plot_starburst_map(full_hdu, args)

                # ---------decorating and saving the figure------------------------------
                fig.text(0.05, 0.98, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
                figname = fig_dir / f'{args.field}_{args.id:05d}_starburst_maps{radial_plot_text}{snr_text}{only_seg_text}{vorbin_text}.png'

                # ---------plotting fit vs stellar mass-----------------------------
                if args.plot_slope_vs_mass:
                    this_row = np.hstack(([args.field, args.id, full_hdu[0].header['RA'], full_hdu[0].header['DEC'], args.z], np.array(starburst_radfit)[:, 0]))
                    this_df = pd.DataFrame(dict(map(lambda i, j: (i, [j]), df_starburst_slope.columns, this_row)))
                    df_starburst_slope = pd.concat([df_starburst_slope, this_df])

            # ---------initialising the metallicity figure------------------------------
            elif args.plot_metallicity:
                fig, logOH_map, logOH_radfit = plot_metallicity_map(full_hdu, args)
                df_logOH_radfit.loc[len(df_logOH_radfit)] = [args.field, args.id, logOH_radfit[0].n, logOH_radfit[0].s, logOH_radfit[1].n, logOH_radfit[1].s]

                # ---------decorating and saving the figure------------------------------
                fig.text(0.05, 0.98, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
                figname = fig_dir / f'{args.field}_{args.id:05d}_metallicity_maps{radial_plot_text}{snr_text}{only_seg_text}{vorbin_text}.png'

            # ---------initialising the full figure------------------------------
            else:
                nrow, ncol = 4 if args.plot_radial_profiles else 3, len(all_lines_to_plot)
                fig = plt.figure(figsize=(13/1., 9/1.) if args.plot_radial_profiles else (13, 6), layout='constrained')

                axis_dirimg = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 0), colspan=1)
                axis_1dspec = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 1), colspan=ncol - 1)
                ax_em_lines = [plt.subplot2grid(shape=(nrow, ncol), loc=(1, item), colspan=1) for item in np.arange(ncol)]  # H alpha, H beta, OII, OIII-4363, OIII
                ax_derived_quant_arr = [plt.subplot2grid(shape=(nrow, ncol), loc=(2, item), colspan=1) for item in np.arange(ncol)]  # SFR, E(B-V), Te, Z (Te), Z (R23)
                ax1, ax2, ax3, ax4, ax5 = ax_derived_quant_arr[0], ax_derived_quant_arr[1], ax_derived_quant_arr[2], ax_derived_quant_arr[3], ax_derived_quant_arr[4]
                if args.plot_radial_profiles:
                    ax_radprof_arr = [plt.subplot2grid(shape=(nrow, ncol), loc=(3, item), colspan=1) for item in np.arange(ncol)]  # SFR, E(B-V), Te, Z (Te), Z (R23)
                    rax1, rax2, rax3, rax4, rax5 = ax_radprof_arr[0], ax_radprof_arr[1], ax_radprof_arr[2], ax_radprof_arr[3], ax_radprof_arr[4]
                else:
                    [rax1, rax2, rax3, rax4, rax5] = np.tile(None, ncol)
                # ---------direct imaging------------------------------
                axis_dirimg = plot_direct_image(full_hdu, axis_dirimg, args)

                # ---------1D spectra------------------------------
                axis_1dspec = plot_1d_spectra(od_hdu, axis_1dspec, args)

                # -----------------emission line maps---------------
                for ind, line in enumerate(all_lines_to_plot):
                    if line in args.available_lines: ax_em_lines[ind] = plot_emission_line_map(line, full_hdu, ax_em_lines[ind], args, cmap='BuPu', vmin=-20, vmax=-18, hide_xaxis=True, hide_yaxis=ind > 0, hide_cbar=False)
                    else: fig.delaxes(ax_em_lines[ind])

                # # ---------------dust map---------------
                # if all([line in args.available_lines for line in ['Ha', 'Hb']]):
                #     try: ax2, EB_V_map, EB_V_radfit, EB_V_int = plot_EB_V_map(full_hdu, ax2, args, radprof_ax=rax2)
                #     except:
                #         EB_V_map, EB_V_radfit, EB_V_int = np.nan, [np.nan, np.nan], np.nan
                #         pass
                # else:
                #     fig.delaxes(ax2)
                #     if args.plot_radial_profiles: fig.delaxes(rax2)
                #     EB_V_map, EB_V_radfit, EB_V_int = np.nan, [np.nan, np.nan], np.nan
                #
                # # ---------------electron temperature map---------------
                # if all([line in args.available_lines for line in ['OIII-4363', 'OIII']]):
                #     try: ax3, Te_map, Te_radfit, Te_int = plot_Te_map(full_hdu, ax3, args, radprof_ax=rax3)
                #     except: pass
                # else:
                #     fig.delaxes(ax3)
                #     if args.plot_radial_profiles: fig.delaxes(rax3)
                #
                # # ---------------metallicity maps---------------
                # if all([line in args.available_lines for line in ['OIII', 'OII', 'Hb']]):
                #     if 'OIII-4363' in args.available_lines:
                #         try: ax4, logOH_Te_map, logOH_Te_radfit, logOH_Te_int = plot_Z_Te_map(full_hdu, ax4, args, radprof_ax=rax4)
                #         except: pass
                #     try: ax5, logOH_R23_map, logOH_R23_radfit, logOH_R23_int = plot_Z_R23_map(full_hdu, ax5, args, radprof_ax=rax5)
                #     except: pass
                # else:
                #     fig.delaxes(ax4)
                #     fig.delaxes(ax5)
                #     if args.plot_radial_profiles:
                #         fig.delaxes(rax4)
                #         fig.delaxes(rax5)

                # ---------------BPT map------------------
                if all([line in args.available_lines for line in ['OIII', 'Hb', 'SII', 'Ha']]):
                    ax1, BPT_map = plot_BPT(full_hdu, rax1, args, cmap='viridis', ax_inset=ax1)
                else:
                    fig.delaxes(ax1)
                    if args.plot_radial_profiles:
                        fig.delaxes(rax1)

                # ---------------logq map------------------
                if all([line in args.available_lines for line in ['OIII', 'OII']]):
                    ax2, logq_map, logq_radfit, logq_int = plot_q_O32_map(full_hdu, ax2, args, radprof_ax=rax2)
                else:
                    fig.delaxes(ax2)
                    if args.plot_radial_profiles:
                        fig.delaxes(rax2)

                # ---------------metallicity map---------------
                if all([line in args.available_lines for line in ['OIII', 'OII', 'Hb']]):
                    ax3, logOH_R23_map, logOH_R23_radfit, logOH_R23_int = plot_Z_R23_map(full_hdu, ax3, args, radprof_ax=rax3)
                else:
                    fig.delaxes(ax4)
                    if args.plot_radial_profiles:
                        fig.delaxes(rax4)

                # ---------------SFR map------------------
                if 'Ha' in args.available_lines:
                    ax4, SFR_map, SFR_radfit, SFR_int = plot_SFR_map(full_hdu, ax4, args, radprof_ax=rax4)
                else:
                    fig.delaxes(ax4)
                    if args.plot_radial_profiles: fig.delaxes(rax4)

                # ---------------F115W/Ha map------------------
                if 'Ha' in args.available_lines:
                    ax5, SB_map, SB_radfit = plot_starburst_map(full_hdu, ax5, args, radprof_axes=rax5)
                else:
                    fig.delaxes(ax5)
                    if args.plot_radial_profiles:
                        fig.delaxes(rax5)

                # ---------decorating and saving the figure------------------------------
                fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.95)
                fig.text(0.05, 0.98, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
                figname = fig_dir / f'{args.field}_{args.id:05d}_{description_text}{vorbin_text}.png'

            # --------for talk plots--------------
            if not args.plot_BPT:
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
            if args.write_file and not (args.test_cutout or args.plot_starburst or args.plot_metallicity):

                # -------collating all the integrated line fluxes from the HDU header----------
                line_properties = []
                for line in all_lines_to_save:
                    try:
                        line_index = np.where(args.available_lines == line)[0][0]
                        flux = full_hdu[0].header[f'FLUX{line_index + 1:03d}']
                        flux_err = full_hdu[0].header[f'ERR{line_index + 1:03d}']
                        line_properties += [flux, flux_err]
                    except IndexError:
                        line_properties += [np.nan, np.nan]
                    try:
                        line_index = np.where(args.available_lines == line)[0][0]
                        line_index_in_cov = int([item for item in list(full_hdu[2].header.keys()) if full_hdu[0].header[f'FLUX{line_index + 1:03d}'] == full_hdu[2].header[item]][0][5:])
                        ew = full_hdu[2].header[f'EW50_{line_index_in_cov:03d}'] # rest-frame EW
                        ew_err = full_hdu[2].header[f'EWHW_{line_index_in_cov:03d}']
                        line_properties += [ew, ew_err]
                    except IndexError:
                        line_properties += [np.nan, np.nan]

                # -------collating all the measured quantities----------
                measured_quants = []
                for quantity in measured_quantities_to_plot:
                    if quantity + '_int' in locals():
                        measured_quants += [locals()[quantity + '_int']]
                    else:
                        measured_quants += [np.nan]
                    if quantity + '_radfit' in locals():
                        try: measured_quants += [locals()[quantity + '_radfit'][0].n, locals()[quantity + '_radfit'][1].n]
                        except AttributeError: measured_quants += [np.nan, np.nan]
                    else:
                        measured_quants += [np.nan, np.nan]

                basic_data = [args.field, f'{args.id:05d}{pixscale_text}', full_hdu[0].header['RA'], full_hdu[0].header['DEC'], args.z]
                flag_data = [args.radius_max, args.snr_cut if args.snr_cut is not None else np.nan, args.only_seg, args.vorbin, args.voronoi_snr if args.vorbin else np.nan, args.voronoi_line if args.vorbin else np.nan]
                this_row = np.hstack([basic_data, flag_data, line_properties, measured_quants])
                this_df = pd.DataFrame(dict(map(lambda i, j: (i, [j]), cols_in_df, this_row)))
                this_df = this_df.replace('N/A', np.nan)
                df = pd.concat([df, this_df])

                if not os.path.isfile(outfilename) or (args.clobber and index == 0):
                    this_df.to_csv(outfilename, index=None, header='column_names')
                    print(f'Wrote to catalog file {outfilename}')
                else:
                    this_df.to_csv(outfilename, index=None, mode='a', header=False)
                    print(f'Appended to catalog file {outfilename}')

            print(f'Completed id {args.id} in {timedelta(seconds=(datetime.now() - start_time3).seconds)}, {len(args.id_arr) - index - 1} to go!')

        # ------------------writing out Z gradient fits, for making MZGR plot later--------------------------
        if args.plot_metallicity:
            outfilename = args.output_dir / 'catalogs' / f'logOHgrad_df{snr_text}{only_seg_text}{vorbin_text}.txt'
            df_logOH_radfit.to_csv(outfilename, index=None, mode='a', header=not os.path.exists(outfilename))
            print(f'Appended metallicity gradient fits to catalog file {outfilename}')

        # ---------plotting spatially resolved BPT-----------------------------
        if args.plot_BPT:
            plt.legend()
            fig.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.95)
            fig.text(0.15, 0.9, f'{args.field}: IDs {",".join(np.array(args.id_arr).astype(str))}', fontsize=args.fontsize, c='k', ha='left', va='top')

            # -----------save figure----------------
            figname = fig_dir / f'{args.field}_{",".join(np.array(args.id_arr).astype(str))}_BPT{snr_text}{only_seg_text}{vorbin_text}.png'
            fig.savefig(figname, transparent=args.fortalk)
            print(f'Saved figure at {figname}')
            plt.show(block=False)

        # -----------plotting F115W/Ha slope vs stellar mass--------------
        if args.plot_slope_vs_mass:
            print(f'Now plotting fit slopess vs stellar mass..')
            df_starburst_slope = get_crossmatch_with_cosmos(df_starburst_slope, args)
            df_starburst_slope['redshift'] = df_starburst_slope['redshift'].astype(np.float32)

            # -----------plotting stuff with the resultant intersecting dataframe--------
            fig, ax = plt.subplots(1, figsize=(8, 6))
            fig.subplots_adjust(left=0.1, right=0.85, bottom=0.1, top=0.95)
            figname = fig_dir / f'{args.field}_lp_mass_vs_ha-f115w_slope{snr_text}{only_seg_text}{vorbin_text}.png'

            # ---------SFMS from df-------
            p = ax.scatter(df_starburst_slope['lp_mass'].values, unp.nominal_values(df_starburst_slope['Ha/F115W_slope'].values), c=df_starburst_slope['redshift'].values, marker='s', s=100, lw=1, edgecolor='k')
            ax.errorbar(df_starburst_slope['lp_mass'].values, unp.nominal_values(df_starburst_slope['Ha/F115W_slope'].values), yerr=unp.std_devs(df_starburst_slope['Ha/F115W_slope'].values), c='gray', fmt='none', lw=1)
            cbar = plt.colorbar(p)
            cbar.set_label('Redshift')

            # ---------annotate axes and save figure-------
            ax.set_xlabel(r'log M$_*$/M$_{\odot}$')
            ax.set_ylabel(r'H$_\alpha$/F115W slope')

            ax.set_ylim(-1.0, 1.0)

            # --------for talk plots--------------
            if args.fortalk:
                mplcyberpunk.add_glow_effects()
                try: mplcyberpunk.make_lines_glow()
                except: pass
                try: mplcyberpunk.make_scatter_glow()
                except: pass

            fig.savefig(figname, transparent=args.fortalk)
            print(f'\nSaved figure as {figname}')
            plt.show(block=False)

        print(f'Completed field {field} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(field_list) - index2 - 1} to go!')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
