'''
    Filename: make_diagnostic_maps.py
    Notes: Plots integrated 1D spectra and 2D emission line maps (from existing .full.fits file), for a given object/s in a given field
    Author : Ayan
    Created: 17-07-24
    Example: run make_diagnostic_maps.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par50 --id 3667
             run make_diagnostic_maps.py --field Par50 --id 823 --plot_radial_profiles
             run make_diagnostic_maps.py --field Par50 --id 823 --pixscale 0.2 --vorbin --voronoi_line Ha --voronoi_snr 3 --plot_radial_profiles
             run make_diagnostic_maps.py --field Par50 --id 823 --plot_radial_profiles --only_seg --snr_cut 3 --plot_mappings
             run make_diagnostic_maps.py --field Par51 --do_all_obj --plot_radial_profiles --only_seg --snr_cut 3 --write_file
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
        image = cut_by_segment(image, full_hdu, args)
        ax.contour(image.mask, levels=0, colors='k', extent=args.extent, linewidths=0.5)

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

    col_arr = ['orange', 'orangered', 'firebrick'] # colors in order: for measured flux, fitted continuum, fitted continuum + line flux
    factor = 1e-19

    # -------plot 1D spectra for each filter-----------
    for index, filter in enumerate(filters):
        print(f'Plotting 1D spectra for filter {filter} which is {index+1} of {nfilters}..')
        table = Table(od_hdu[filter].data)
        table['rest_wave'] = table['wave'] / (1 + args.z)
        ax.plot(table['rest_wave'], table['flux'] / table['flat'] / factor, lw=0.5, c=col_arr[index], alpha=0.8) # need to divide all columns with 'flat' to get the right units (ergs/s/cm^2/A)
        ax.plot(table['rest_wave'], table['cont'] / table['flat'] / factor, lw=0.5, c='grey')
        ax.plot(table['rest_wave'], table['line'] / table['flat'] / factor, lw=1, c='indigo')
        ax.text(float(filters[0][1:-1]) * 1e2 * 0.85 / (1 + args.z), args.flam_max * 0.95 - index * 0.1, filter, c=col_arr[index], fontsize=args.fontsize, ha='left', va='top')

    ax.set_xlabel(r'Rest-frame wavelength ($\AA$)', fontsize=args.fontsize)
    ax.set_ylabel(r'f$_{\lambda}$ ' + '(%.0e ' % factor + r'ergs/s/cm$^2$/A)', fontsize=args.fontsize/1.2)

    ax.set_ylim(0, args.flam_max) # x factor

    # ---observed wavelength axis-------
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(ax.get_xticks())
    ax2.set_xticklabels(['%.1F' % (item * (1 + args.z) / 1e4) for item in ax2.get_xticks()], fontsize=args.fontsize)
    ax2.set_xlabel(r'Observed wavelength ($\mu$)', fontsize=args.fontsize)

    # ---vertical lines for emission line wavelengths------
    if args.plot_mappings: ax = plot_MAPPINGS_lines(ax)
    else: ax = plot_linelist(ax)

    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    return ax

# ---------------------------------------------------------------------------------
def plot_binned_profile(xdata, ydata, ax, color='darkorange', yerr=None):
    '''
    Function to overplot binned data on existing plot in a given axis
    Returns axis handle
    '''
    df = pd.DataFrame({'xcol': xdata.flatten(), 'ycol':ydata.flatten()})

    bin_edges = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 10)
    df['binned_cat'] = pd.cut(df['xcol'], bin_edges)

    if yerr is not None:
        df['weightcol'] = 1 / yerr.flatten() ** 2
        agg_func = lambda x: np.sum(x * df.loc[x.index, 'weightcol']) / np.sum(df.loc[x.index, 'weightcol']) # function to get weighted mean
        agg_u_func = lambda x: np.sqrt(((np.sum(df.loc[x.index, 'weightcol'] * x**2) / np.sum(df.loc[x.index, 'weightcol'])) - (np.sum(x * df.loc[x.index, 'weightcol']) / np.sum(df.loc[x.index, 'weightcol']))**2) * (np.sum(df.loc[x.index, 'weightcol']**2)) / (np.sum(df.loc[x.index, 'weightcol'])**2 - np.sum(df.loc[x.index, 'weightcol']**2))) # eq 6 of http://seismo.berkeley.edu/~kirchner/Toolkits/Toolkit_12.pdf
    else:
        agg_func, agg_u_func = np.mean, np.std

    y_binned = df.groupby('binned_cat', as_index=False).agg([('ycol', agg_func)])['ycol'].values.flatten()
    y_u_binned = df.groupby('binned_cat', as_index=False).agg([('ycol', agg_u_func)])['ycol'].values.flatten()
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
    pixscale_kpc = ((args.pixscale * u.arcsec).to(u.radian) * args.distance.to('kpc')).value # kpc
    center_pix = image_shape[0] / 2.
    distance_map = np.array([[np.sqrt((i - center_pix)**2 + (j - center_pix)**2) for j in range(image_shape[1])] for i in range(image_shape[0])]) * pixscale_kpc # kpc

    return distance_map
# --------------------------------------------------------------------------------------------------------------------
def plot_radial_profile(image, ax, args, label=None, ymin=None, ymax=None, hide_xaxis=False, hide_yaxis=False):
    '''
    Plots the average radial profile for a given 2D map in the given axis
    Returns the axis handle
    '''
    print(f'Plotting radial profile of {label}..')

    distance_map = get_distance_map(np.shape(image), args)
    ax.scatter(distance_map, image, c='grey', s=1, alpha=0.2)

    ax.set_xlim(0, 2 * ((args.arcsec_limit * u.arcsec).to(u.radian) * args.distance.to('kpc')).value) # kpc
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0, args.radius_max)
    ax.set_box_aspect(1)

    ax, linefit = plot_binned_profile(distance_map, image, ax)

    if hide_xaxis:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel('Distance (kpc)', fontsize=args.fontsize)
        #ax.set_xticklabels(['%d' % ((item * u.arcsec).to(u.radian) * args.distance.to('kpc')).value for item in ax.get_xticks()], fontsize=args.fontsize)
        ax.tick_params(axis='x', which='major', labelsize=args.fontsize)

    if hide_yaxis:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel(label, fontsize=args.fontsize)
        ax.tick_params(axis='y', which='major', labelsize=args.fontsize)

    return ax, linefit

# --------------------------------------------------------------------------------------------------------------------
def plot_2D_map(image, ax, args, label=None, cmap=None, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=False, radprof_ax=None):
    '''
    Plots the emission map for a given line in the given axis
    Returns the axis handle
    '''
    if cmap is None: cmap = 'cividis'
    print(f'Plotting 2D map of {label}..')

    p = ax.imshow(image, cmap=cmap, origin='lower', extent=args.extent, vmin=vmin, vmax=vmax)

    ax.set_xlim(-args.arcsec_limit, args.arcsec_limit) # arcsec
    ax.set_ylim(-args.arcsec_limit, args.arcsec_limit) # arcsec

    ax.text(ax.get_xlim()[0] * 0.88, ax.get_ylim()[1] * 0.88, label, c='k', fontsize=args.fontsize, ha='left', va='top', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
    ax.scatter(0, 0, marker='x', s=10, c='grey')

    if hide_xaxis:
        ax.set_xticklabels([])
    else:
        if args.plot_target_frame:
            ax.set_xlabel('Offset (kpc)', fontsize=args.fontsize)
            ax.set_xticklabels(['%d' % ((item * u.arcsec).to(u.radian) * args.distance.to('kpc')).value for item in ax.get_xticks()], fontsize=args.fontsize)
        else:
            ax.set_xlabel('RA (")', fontsize=args.fontsize)
        ax.tick_params(axis='x', which='major', labelsize=args.fontsize)

    if hide_yaxis:
        ax.set_yticklabels([])
    else:
        if args.plot_target_frame:
            ax.set_ylabel('Offset (kpc)', fontsize=args.fontsize)
            ax.set_yticklabels(['%d' % ((item * u.arcsec).to(u.radian) * args.distance.to('kpc')).value for item in ax.get_xticks()], fontsize=args.fontsize)
        else:
            ax.set_ylabel('Dec (")', fontsize=args.fontsize)
        ax.tick_params(axis='y', which='major', labelsize=args.fontsize)

    if not hide_cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad='2%')
        cbar = plt.colorbar(p, cax=cax, orientation='vertical')
        cbar.ax.tick_params(labelsize=args.fontsize)

    if args.plot_radial_profiles and radprof_ax is not None:
        radius_pix = ((args.radius_max / args.distance.to('kpc').value) * u.radian).to(u.arcsec)
        circle = plt.Circle((0, 0), radius_pix.value, color='k', fill=False, lw=0.5)
        ax.add_patch(circle)
        radprof_ax, radprof_fit = plot_radial_profile(image, radprof_ax, args, label=label.split(r'$_{\rm int}')[0], ymin=vmin, ymax=vmax)
    else:
        radprof_fit = [np.nan, np.nan] # dummy values for when the fit was not performed

    return ax, radprof_fit

# --------------------------------------------------------------------------------------------------------------------
def bin_2D(map, bin_IDs):
    '''
    Bin a given 2D map by given bin_IDs
    Returns the binned 2D map (of same shape as input map)
    '''
    binned_map = np.zeros(np.shape(map))
    for id in np.unique(bin_IDs):
        binned_map[bin_IDs == id] = map[bin_IDs == id].mean()

    return binned_map

# --------------------------------------------------------------------------------------------------------------------
def get_voronoi_bin_IDs(map, map_err, snr_thresh, plot=False, quiet=True):
    '''
    Compute the Voronoi bin IDs a given 2D map and corresponding uncertainty and SNR threshold
    Returns the 2D map (of same shape as input map) with just the IDs
    '''
    x_size, y_size = np.shape(map)
    x_coords = np.repeat(np.arange(x_size), y_size)
    y_coords = np.tile(np.arange(y_size), x_size)

    map = np.ma.masked_where(~np.isfinite(map_err), map)
    map_err = np.ma.masked_where(~np.isfinite(map_err), map_err)

    map = np.ma.masked_where(map < 0, map)
    map_err = np.ma.masked_where(map < 0, map_err)

    binIDs, _, _, _, _, _, _, _ = voronoi_2d_binning(x_coords, y_coords, map.flatten(), map_err.flatten(), snr_thresh, plot=plot, quiet=quiet, cvt=False, pixelsize=1)
    binID_map = binIDs.reshape(np.shape(map))

    return binID_map

# --------------------------------------------------------------------------------------------------------------------
def cut_by_segment(map, full_hdu, args):
    '''
    Mask a given 2D map according to the segmentation map from the HDU
    Returns the masked map
    '''
    segmentation_map = full_hdu['SEG'].data
    cut_map = np.ma.masked_where(segmentation_map != args.id, map)

    return cut_map

# --------------------------------------------------------------------------------------------------------------------
def get_emission_line_map(line, full_hdu, args, dered=True):
    '''
    Retrieve the emission map for a given line from the HDU
    Returns the 2D line image
    '''

    line_index = np.where(args.available_lines == line)[0][0]
    ext = 5 + 2 * args.ndfilt + 4 * line_index
    line_map = full_hdu[ext].data * 1e-17 # this gives 200 x 200 array; in units of ergs/s/cm^2
    line_wave = full_hdu[ext].header['RESTWAVE'] # in Angstrom
    line_map_err = 1e-17 / full_hdu[ext + 3].data  # this gives 200 x 200 array; 1/LINEWHT = flux uncertainty; in units of ergs/s/cm^2

    if args.only_seg:
        line_map = cut_by_segment(line_map, full_hdu, args)

    if args.snr_cut is not None:
        snr_map = line_map / line_map_err
        line_map = np.ma.masked_where(~np.isfinite(snr_map), line_map)
        line_map = np.ma.masked_where(snr_map < args.snr_cut, line_map)

    if args.vorbin:
        if args.voronoi_line is None: # No reference emission line specified, so Voronoi IDs need to be computed now
            bin_IDs = get_voronoi_bin_IDs(line_map, line_map_err, args.voronoi_snr)
        else: # Reference emission line specified for Voronoi binning, so bin IDs have been pre-computed
            bin_IDs = args.voronoi_bin_IDs
        line_map = bin_2D(line_map, bin_IDs)

    # -----------getting the integrated flux value-----------------
    line_int = full_hdu[0].header[f'FLUX{line_index + 1:03d}'] # ergs/s/cm^2

    # -----------getting the dereddened flux value-----------------
    if dered:
        line_map = get_dereddened_flux(line_map, line_wave, args.EB_V)
        line_int = get_dereddened_flux(line_int, line_wave, args.EB_V)

    return line_map, line_wave, line_int

# --------------------------------------------------------------------------------------------------------------------
def plot_emission_line_map(line, full_hdu, ax, args, cmap='cividis', EB_V=None, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=False):
    '''
    Plots the emission map for a given line in the given axis
    Returns the axes handle
    '''

    line_map, line_wave, line_int = get_emission_line_map(line, full_hdu, args, dered=False)
    ax, _ = plot_2D_map(np.log10(line_map), ax, args, label=r'%s$_{\rm int}$ = %.1e' % (line, line_int), cmap=cmap, vmin=vmin, vmax=vmax, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar)

    # ------to annotate with EW--------------------
    line_index = np.where(args.available_lines == line)[0][0]
    line_index_in_cov = int([item for item in list(full_hdu[2].header.keys()) if full_hdu[0].header[f'FLUX{line_index + 1:03d}'] == full_hdu[2].header[item]][0][5:])
    line_ew = full_hdu[2].header[f'EW50_{line_index_in_cov:03d}']

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
def compute_EB_V(Ha_flux, Hb_flux, args, verbose=False):
    '''
    Calculates and returns the color excess given observed H alpha and H beta fluxes
    Based on Eqn 4 of Dominguez+2013 (https://iopscience.iop.org/article/10.1088/0004-637X/763/2/145/pdf)
    '''
    theoretical_ratio = 2.86

    obs_ratio = Ha_flux / Hb_flux
    EB_V = 1.97 * np.log10(obs_ratio / theoretical_ratio)

    if type(EB_V) == np.float64:
        if obs_ratio < theoretical_ratio:
            EB_V = 0.
            if verbose: print(f'Based on integrated fluxes Ha = {Ha_flux:.2e}, Hb = {Hb_flux:.2e}, the observed ratio is LOWER than theoretical ratio, so E(B-V) would be unphysical, so just assuming {EB_V}')
        else:
            if verbose: print(f'Based on integrated fluxes Ha = {Ha_flux:.2e}, Hb = {Hb_flux:.2e}, determine E(B-V) = {EB_V:.2f}')
    else:
        EB_V[~np.isfinite(EB_V)] = 0.

    return EB_V

# --------------------------------------------------------------------------------------------------------------------
def get_EB_V(full_hdu, args, verbose=False):
    '''
    Computes and returns the spatially resolved as well as integrated dust extinction map from a given HDU
    Based on Eqn 4 of Dominguez+2013 (https://iopscience.iop.org/article/10.1088/0004-637X/763/2/145/pdf)
    '''

    Ha_map, Ha_wave, Ha_int = get_emission_line_map('Ha', full_hdu, args, dered=False) # do not need to deredden the lines when we are fetching the flux in order to compute reddening
    Hb_map, Hb_wave, Hb_int = get_emission_line_map('Hb', full_hdu, args, dered=False)

    EB_V_map = compute_EB_V(Ha_map, Hb_map, args)
    EB_V_int = compute_EB_V(Ha_int, Hb_int, args, verbose=verbose)

    return EB_V_map, EB_V_int

# --------------------------------------------------------------------------------------------------------------------
def plot_EB_V_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the dust extinction map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D E(B-V) map just produced
    '''
    lim, label = [0, 1], 'E(B-V)'

    EB_V_map, EB_V_int = get_EB_V(full_hdu, args)
    ax, EB_V_radfit = plot_2D_map(EB_V_map, ax, args, label=r'%s$_{\rm int}$ = %.1f' % (label, EB_V_int), cmap='YlOrBr', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

    return ax, EB_V_map, EB_V_radfit, EB_V_int

# --------------------------------------------------------------------------------------------------------------------
def compute_SFR(Ha_flux, distance):
    '''
    Calculates and returns the SFR given observed H alpha fluxes
    Conversion factor is from Kennicutt 1998 (Eq 2 of https://ned.ipac.caltech.edu/level5/Sept01/Rosa/Rosa3.html)
    '''
    Ha_flux = Ha_flux * 4 * np.pi * (distance.to('cm').value) ** 2 # converting to ergs/s
    sfr = Ha_flux * 7.9e-42 # line_map in args/s; SFR in Msun/yr

    return sfr

# --------------------------------------------------------------------------------------------------------------------
def get_SFR(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated SFR from a given HDU
    '''
    Ha_map, Ha_wave, Ha_int = get_emission_line_map('Ha', full_hdu, args)

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
    ax, SFR_radfit = plot_2D_map(np.log10(SFR_map), ax, args, label=r'%s$_{\rm int}$ = %.1f' % (label, SFR_int), cmap='Blues', radprof_ax=radprof_ax, vmin=lim[0], vmax=lim[1])

    return ax, SFR_map, SFR_radfit, SFR_int

# --------------------------------------------------------------------------------------------------------------------
def compute_Te(OIII4363_flux, OIII5007_flux, args):
    '''
    Calculates and returns the Te given observed line fluxes
    Conversion factor is from Nicholls+2017
    '''
    ratio = OIII4363_flux / OIII5007_flux
    logTe = np.poly1d([0., 9.18962, 3.30355])(np.log10(ratio)) / np.poly1d([-0.00935, -0.15034, 2.09136, 1.000])(np.log10(ratio))
    Te = 10 ** logTe

    return Te

# --------------------------------------------------------------------------------------------------------------------
def get_Te(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated Te from a given HDU
    '''
    OIII4363_map, OIII4363_wave, OIII4363_int = get_emission_line_map('OIII-4363', full_hdu, args)
    OIII5007_map, OIII5007_wave, OIII5007_int = get_emission_line_map('OIII', full_hdu, args)

    Te_map = compute_Te(OIII4363_map, OIII5007_map, args)
    Te_int = compute_Te(OIII4363_int, OIII5007_int, args)

    return Te_map, Te_int

# --------------------------------------------------------------------------------------------------------------------
def plot_Te_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the T_e map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D T_e map just produced
    '''
    lim, label = [1, 7], r'T$_e$'
    Te_map, Te_int = get_Te(full_hdu, args)
    ax, Te_radfit = plot_2D_map(np.log10(Te_map), ax, args, label=r'%s$_{\rm int}$ = %.1e' % (label, Te_int), cmap='OrRd_r', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

    return ax, Te_map, Te_radfit, Te_int

# --------------------------------------------------------------------------------------------------------------------
def compute_Z_Te(OII3727_flux, OIII5007_flux, Hbeta_flux, Te, args, ne=1e3):
    '''
    Calculates and returns the Te metallicity given observed line fluxes
    Conversion factor is from Nicholls+2017
    '''
    def poly(R, t, x, a, b, c, d, e):
        return np.log10(R) + a + b / t - c * np.log10(t) - d * t + np.log10(1 + e * x)  # eqn 3 I06 pattern

    t = Te * 1e-4
    x = 1e-4 * ne * np.sqrt(t)

    ratio1 = OII3727_flux / Hbeta_flux
    ratio2 = OIII5007_flux / Hbeta_flux
    log_O2H2 = poly(ratio1, t, x, 5.961, 1.676, 0.4, 0.034, 1.35) - 12  # coefficients from eqn 3 I06
    log_O3H2 = poly(ratio2, t, x, 6.200, 1.251, -0.55, -0.014, 0.0) - 12  # coefficients from eqn 5 I06
    log_OH = np.log10(10 ** log_O2H2 + 10 ** log_O3H2) + 12

    return log_OH

# --------------------------------------------------------------------------------------------------------------------
def get_Z_Te(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated Te metallicity from a given HDU
    '''
    OII3727_map, line_wave, OII3727_int = get_emission_line_map('OII', full_hdu, args)
    OIII5007_map, line_wave, OIII5007_int = get_emission_line_map('OIII', full_hdu, args)
    Hbeta_map, line_wave, Hbeta_int = get_emission_line_map('Hb', full_hdu, args)

    Te_map, Te_int = get_Te(full_hdu, args)

    logOH_map = compute_Z_Te(OII3727_map, OIII5007_map, Hbeta_map, Te_map, args)
    logOH_int = compute_Z_Te(OII3727_int, OIII5007_int, Hbeta_int, Te_int, args)

    return logOH_map, logOH_int

# --------------------------------------------------------------------------------------------------------------------
def plot_Z_Te_map(full_hdu, ax, args, Te_map, radprof_ax=None):
    '''
    Plots the T_e-based metallicity map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D metallicity map just produced
    '''
    lim, label = [6, 9], 'Z (Te)'
    logOH_map, logOH_int = get_Z_Te(full_hdu, args)
    ax, logOH_radfit = plot_2D_map(logOH_map, ax, args, label=r'%s$_{\rm int}$ = %.1f' % (label, logOH_int), cmap='Greens', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

    return ax, logOH_map, logOH_radfit, logOH_int

# --------------------------------------------------------------------------------------------------------------------
def compute_Z_R23(OII3727_flux, OIII5007_flux, Hbeta_flux, args):
    '''
    Calculates and returns the R23 metallicity given observed line fluxes
    Conversion factor is from Kewley+2002
    '''
    def poly(R23, k):
        return  (-k[1] + np.sqrt(k[1]**2 - 4*k[2]*(k[0] - R23)))/(2*k[2])

    ratio = (OII3727_flux + OIII5007_flux) / Hbeta_flux
    R23 = np.log10(ratio)
    log_OH = poly(R23, [-44.7026, 10.8052, -0.640113]) #k0-2 parameters for q=8e7 from Table 3 of KD02 last row for q=8e7

    return log_OH

# --------------------------------------------------------------------------------------------------------------------
def get_Z_R23(full_hdu, args, verbose=False):
    '''
    Computes and returns the spatially resolved as well as intregrated R23 metallicity from a given HDU
    '''
    OII3727_map, line_wave, OII3727_int = get_emission_line_map('OII', full_hdu, args)
    OIII5007_map, line_wave, OIII5007_int = get_emission_line_map('OIII', full_hdu, args)
    Hbeta_map, line_wave, Hbeta_int = get_emission_line_map('Hb', full_hdu, args)

    logOH_map = compute_Z_R23(OII3727_map, OIII5007_map, Hbeta_map, args)
    logOH_int = compute_Z_R23(OII3727_int, OIII5007_int, Hbeta_int, args)

    return logOH_map, logOH_int

# --------------------------------------------------------------------------------------------------------------------
def plot_Z_R23_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the R23 metallicity map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D metallicity map just produced
    '''
    lim, label = [6, 9], 'Z (R23)'
    logOH_map, logOH_int = get_Z_R23(full_hdu, args)
    ax, logOH_radfit = plot_2D_map(logOH_map, ax, args, label=r'%s$_{\rm int}$ = %.1f' % (label, logOH_int), cmap='Greens', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

    return ax, logOH_map, logOH_radfit, logOH_int

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    # ---------determining filename suffixes-------------------------------
    extract_dir = args.input_dir / args.field / 'Extractions'

    radial_plot_text = '_wradprof' if args.plot_radial_profiles else ''
    snr_text = f'_snr{args.snr_cut}' if args.snr_cut is not None else ''
    only_seg_text = '_onlyseg' if args.only_seg else ''
    pixscale_text = '' if args.pixscale == 0.04 else f'_{args.pixscale}arcsec_pix'

    outfilename = args.output_dir / args.field / f'{args.field}_all_diag_results.txt'

    catalog_file = extract_dir / f'{args.field}-ir.cat.fits'
    if os.path.exists(catalog_file): catalog = GTable.read(catalog_file)


    if args.do_all_obj:
        args.id_arr = catalog['NUMBER']
        if args.start_id: args.id_arr = args.id_arr[args.start_id - 1:]
    else:
        args.id_arr = args.id

    # ---------for diplay and amimations----------------
    if len(args.id_arr) > 10: args.hide = True # if too many plots, do not display them, just save them
    if len(args.id_arr) > 20: args.make_anim = True
    else: args.make_anim = False

    if args.make_anim:
        outputfile = args.output_dir / args.field / f'{args.field}_all_diag_plots{radial_plot_text}{snr_text}{only_seg_text}.mp4'
        duration_per_frame = 0.1 #sec
        writer = imageio.get_writer(outputfile, mode='I', fps=int(1. / duration_per_frame))

    # ----------------------initiliasing dataframe-------------------------------
    lines_to_plot = ['Ha', 'Hb', 'OII', 'OIII-4363', 'OIII']
    measured_quantities_to_plot = ['EB_V', 'SFR', 'Te', 'logOH_Te', 'logOH_R23']

    if args.write_file:
        basic_cols = ['field', 'objid', 'ra', 'dec', 'redshift']
        flag_cols = ['radfit_extent_kpc', 'snr_cut', 'flag_only_seg', 'flag_vorbin', 'vor_snr', 'vor_line']
        cols_in_df = np.hstack([basic_cols, flag_cols, [item + '_int' for item in lines_to_plot], np.hstack([[item + '_int', item + '_cen', item + '_slope'] for item in measured_quantities_to_plot])])
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
        output_subdir = args.output_dir / args.field / f'{args.id:05d}{pixscale_text}'
        if not args.do_all_obj: output_subdir.mkdir(parents=True, exist_ok=True)
        full_fits_file = f'{args.field}_{args.id:05d}.full.fits'

        if os.path.exists(output_subdir / full_fits_file): # if the fits files are in sub-directories for individual objects
            args.work_dir = output_subdir
        elif os.path.exists(extract_dir / full_fits_file): # if the fits files are in Extractions/
            args.work_dir = extract_dir
        else:
            print(f'Could not find {full_fits_file} for ID {args.id}, so skipping it.')
            continue

        # ------------read in fits files--------------------------------
        od_hdu = fits.open(args.work_dir / f'{args.field}_{args.id:05d}.1D.fits')
        full_hdu = fits.open(args.work_dir / f'{args.field}_{args.id:05d}.full.fits')

        # ----------determining global parameters------------
        args.available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
        args.z = full_hdu[0].header['REDSHIFT']
        args.distance = cosmo.comoving_distance(args.z)
        args.ndfilt = full_hdu[0].header['NDFILT']
        args.nlines = full_hdu[0].header['NUMLINES']
        args.pix_arcsec = full_hdu[7].header['PIXASEC']
        args.pa_arr = np.unique([full_hdu[0].header[item] for item in list(full_hdu[0].header.keys()) if 'PA00' in item])
        try: args.mag = catalog[catalog['NUMBER'] == args.id]['MAG_AUTO'].data.data[0]
        except: args.mag = np.nan

        line_wcs = pywcs.WCS(full_hdu['DSCI'].header)
        pix_size = utils.get_wcs_pscale(line_wcs)
        imsize_arcsec = full_hdu['DSCI'].data.shape[0] * pix_size
        dp = 0 # -0.5 * pix_size  # FITS reference is center of a pixel, array is edge
        args.extent = (-imsize_arcsec / 2. - dp, imsize_arcsec / 2. - dp, -imsize_arcsec / 2. - dp, imsize_arcsec / 2. - dp)
        args.EB_V = 0. # until gets over-written, if both H alpha and H beta lines are present

        if args.vorbin and args.voronoi_line is not None:
            line_index = np.where(args.available_lines == args.voronoi_line)[0][0]
            ext = 5 + 2 * args.ndfilt + 4 * line_index
            line_map = full_hdu[ext].data * 1e-17  # this gives 200 x 200 array; in units of ergs/s/cm^2
            line_map_err = 1e-17 / full_hdu[ext + 3].data   # this gives 200 x 200 array; in units of ergs/s/cm^2
            args.voronoi_bin_IDs = get_voronoi_bin_IDs(line_map, line_map_err, args.voronoi_snr)#, plot=True, quiet=False)

        if args.plot_radial_profiles:
            seg_map = full_hdu['SEG'].data
            distance_map = get_distance_map(np.shape(seg_map), args)
            distance_map = np.ma.compressed(np.ma.masked_where(seg_map != args.id, distance_map))
            args.radius_max = np.max(distance_map)
        else:
            args.radius_max = np.nan

        # ---------initialising the figure------------------------------
        if not args.keep: plt.close('all')
        nrow, ncol = 4 if args.plot_radial_profiles else 3, len(lines_to_plot)
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
        for ind, line in enumerate(lines_to_plot):
            if line in args.available_lines: ax_em_lines[ind] = plot_emission_line_map(line, full_hdu, ax_em_lines[ind], args, cmap='BuPu', vmin=-20, vmax=-18, hide_xaxis=True, hide_yaxis=ind > 0, hide_cbar=False) #ind != len(lines_to_plot) - 1) # line_map in ergs/s/cm^2
            else: fig.delaxes(ax_em_lines[ind])

        # ---------------dust map---------------
        if all([line in args.available_lines for line in ['Ha', 'Hb']]):
            _, args.EB_V = get_EB_V(full_hdu, args, verbose=True)
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
            ax_Z_Te, logOH_Te_map, logOH_Te_radfit, logOH_Te_int = plot_Z_Te_map(full_hdu, ax_Z_Te, args, Te_map, radprof_ax=rax_Z_Te)
            ax_Z_R23, logOH_R23_map, logOH_R23_radfit, logOH_R23_int = plot_Z_R23_map(full_hdu, ax_Z_R23, args, radprof_ax=rax_Z_R23)
        else:
            fig.delaxes(ax_Z_Te)
            fig.delaxes(ax_Z_R23)
            if args.plot_radial_profiles:
                fig.delaxes(rax_Z_Te)
                fig.delaxes(rax_Z_R23)

        # ---------decorating and saving the figure------------------------------
        fig.text(0.05, 0.98, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
        figdir = args.output_dir / args.field if args.do_all_obj else output_subdir
        figname = figdir / f'{args.field}_{args.id:05d}_all_diag_plots{radial_plot_text}{snr_text}{only_seg_text}.png'
        fig.savefig(figname)
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
            line_fluxes = []
            for line in lines_to_plot:
                try:
                    line_index = np.where(args.available_lines == line)[0][0]
                    flux = full_hdu[0].header[f'FLUX{line_index + 1:03d}']
                    line_fluxes.append(flux)
                except IndexError:
                    line_fluxes.append(np.nan)

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
            this_row = np.hstack([basic_data, flag_data, line_fluxes, measured_quants])
            this_df = pd.DataFrame(dict(map(lambda i, j: (i, [j]), cols_in_df, this_row)))
            df = pd.concat([df, this_df])

            if not os.path.isfile(outfilename) or (args.clobber and index == 0):
                this_df.to_csv(outfilename, sep='\t', index=None, header='column_names')
                print(f'Wrote to catalog file {outfilename}')
            else:
                this_df.to_csv(outfilename, sep='\t', index=None, mode='a', header=False)
                print(f'Appended to catalog file {outfilename}')

        print(f'Completed id {args.id} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(args.id_arr) - index - 1} to go!')

    os.chdir(args.code_dir)
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

run ~/Work/astro/ayan_codes/animate_png.py --inpath /Volumes/Elements/acharyya_backup/Work/as
     ...: tro/passage/passage_output/Par051/ --rootname Par051_*_all_diag_plots_wradprof_snr3.0_onlyseg
     ...: .png --delay 0.1