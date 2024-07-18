'''
    Filename: make_diagnostic_maps.py
    Notes: Plots integrated 1D spectra and 2D emission line maps (from existing .full.fits file), for a given object/s in a given field
    Author : Ayan
    Created: 17-07-24
    Example: run make_diagnostic_maps.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par50 --id 3667
             run make_diagnostic_maps.py --field Par50 --id 823 --vorbin --voronoi_line Ha --voronoi_snr 3 --plot_radial_profiles
'''

from header import *
from util import *
from matplotlib import cm as mpl_cm
from matplotlib.gridspec import GridSpec

start_time = datetime.now()

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
    ax.scatter(0, 0, marker='x', s=10, c='grey')
    #cbar = plt.colorbar(p)

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
        ax.text(lines_df.iloc[index]['wave'] + np.diff(ax.get_xlim())[0] * 0.01, ax.get_ylim()[1] * 0.9, lines_df.iloc[index]['species'], rotation=90, va='top', ha='center', fontsize=args.fontsize)

    return ax


# --------------------------------------------------------------------------------------------------------------------
def get_linelist(wave_lim=None):
    '''
    Reads in an emission line list
    Returns list of lines that are within the given wavelength limits, as a pandas dataframe
    '''
    line_list_file = HOME / 'Desktop/mage_plot/labframe.shortlinelist'

    lines_df = pd.read_table(line_list_file, comment='#', delim_whitespace=True)
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
        ax.text(lines_df.iloc[index]['restwave'] + np.diff(ax.get_xlim())[0] * 0.01, ax.get_ylim()[1] * 0.9, lines_df.iloc[index]['LineID'], rotation=90, va='top', ha='center', fontsize=args.fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_1d_spectra(od_hdu, ax, args):
    '''
    Plots the 1D spectra in the given axis
    Returns the axis handle
    '''
    nfilters = sum(['GRISM' in item for item in list(od_hdu[0].header.keys())])
    filters = [od_hdu[0].header[f'GRISM{item + 1:03d}'] for item in range(nfilters)]

    col_arr = ['gold', 'sandybrown', 'firebrick'] # colors in order: for measured flux, fitted continuum, fitted continuum + line flux
    factor = 1e-19

    # -------plot 1D spectra for each filter-----------
    for index, filter in enumerate(filters):
        print(f'Plotting 1D spectra for filter {filter} which is {index+1} of {nfilters}..')
        table = Table(od_hdu[filter].data)
        table['rest_wave'] = table['wave'] / (1 + args.z)
        ax.plot(table['rest_wave'], table['flux'] / table['flat'] / factor, lw=0.5, c=col_arr[index], alpha=0.5) # need to divide all columns with 'flat' to get the right units (ergs/s/cm^2/A)
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
    ax = plot_MAPPINGS_lines(ax)
    #ax = plot_linelist(ax)

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

    # ----------to plot mean binned y vs x profile--------------
    ax.errorbar(x_bin_centers, y_binned, c=color, yerr=y_u_binned, lw=1, ls='none', zorder=1)
    ax.scatter(x_bin_centers, y_binned, c=color, s=20, lw=0.2, ec='black', zorder=10)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_radial_profile(image, ax, args, label=None, cmap=None, ymin=None, ymax=None, hide_xaxis=False, hide_yaxis=False):
    '''
    Plots the average radial profile for a given 2D map in the given axis
    Returns the axis handle
    '''
    print(f'Plotting radial profile of {label}..')

    pixscale_kpc = ((args.pixscale * u.arcsec).to(u.radian) * args.distance.to('kpc')).value # kpc
    center_pix = np.shape(image)[0] / 2.
    distance_map = np.array([[np.sqrt((i - center_pix)**2 + (j - center_pix)**2) for j in range(np.shape(image)[0])] for i in range(np.shape(image)[0])]) * pixscale_kpc # kpc

    ax.scatter(distance_map, image, c='grey', s=1, alpha=0.2)

    ax.set_xlim(0, 2 * ((args.arcsec_limit * u.arcsec).to(u.radian) * args.distance.to('kpc')).value) # kpc
    ax.set_ylim(ymin, ymax)
    #ax.set_box_aspect(1)

    ax = plot_binned_profile(distance_map, image, ax)

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

    return ax

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
        fig = ax.figure
        pad = -0.1
        [ax_l, ax_b, ax_w, ax_h] = ax.get_position().bounds # left, bottom, width, height
        cax = fig.add_axes([ax_l + ax_w + pad, ax_b, ax_w * 0.1, ax_h]) # left, bottom, width, height
        cbar = plt.colorbar(p, cax=cax, orientation='vertical')
        cbar.ax.tick_params(labelsize=args.fontsize)

    if args.plot_radial_profiles and radprof_ax is not None: radprof_ax = plot_radial_profile(image, radprof_ax, args, label=label, ymin=vmin, ymax=vmax)

    return ax

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
def get_emission_line_map(line, full_hdu, args):
    '''
    Retrieve the emission map for a given line from the HDU
    Returns the 2D line image
    '''

    line_index = np.where(args.available_lines == line)[0][0]
    ext = 5 + 2 * args.ndfilt + 4 * line_index
    line_map = full_hdu[ext].data * 1e-17 # this gives 200 x 200 array; in units of ergs/s/cm^2
    line_wave = full_hdu[ext].header['RESTWAVE'] # in Angstrom

    if args.vorbin:
        if args.voronoi_line is None: # No reference emission line specified, so Voronoi IDs need to be computed now
            line_map_err = 1e-17 / full_hdu[ext + 3].data   # this gives 200 x 200 array; in units of ergs/s/cm^2
            bin_IDs = get_voronoi_bin_IDs(line_map, line_map_err, args.voronoi_snr)
        else: # Reference emission line specified for Voronoi binning, so bin IDs have been pre-computed
            bin_IDs = args.voronoi_bin_IDs
        line_map = bin_2D(line_map, bin_IDs)

    return line_map, line_wave

# --------------------------------------------------------------------------------------------------------------------
def plot_emission_line_map(line, full_hdu, ax, args, cmap='cividis', EB_V=None, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=False):
    '''
    Plots the emission map for a given line in the given axis
    Returns the axes handle
    '''

    line_map, line_wave = get_emission_line_map(line, full_hdu, args)
    dered_line_map = get_dereddened_flux(line_map, line_wave, EB_V)
    ax = plot_2D_map(np.log10(line_map), ax, args, label=line, cmap=cmap, vmin=vmin, vmax=vmax, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar)

    return dered_line_map, line_wave, ax

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
def calculate_EB_V(full_hdu, args):
    '''
    Calculates and returns the color excess given observed integrated H alpha and H beta fluxes
    Based on Eqn 4 of Dominguez+2013 (https://iopscience.iop.org/article/10.1088/0004-637X/763/2/145/pdf)
    '''
    theoretical_ratio = 2.86

    if all([line in args.available_lines for line in ['Ha', 'Hb']]):
        Ha_index = np.where(args.available_lines == 'Ha')[0][0]
        Ha_int_flux = full_hdu[0].header[f'FLUX{Ha_index + 1:03d}']

        Hb_index = np.where(args.available_lines == 'Hb')[0][0]
        Hb_int_flux = full_hdu[0].header[f'FLUX{Hb_index + 1:03d}']

        obs_ratio = Ha_int_flux / Hb_int_flux
        if obs_ratio < theoretical_ratio:
            EB_V = 0.
            print(f'Based on integrated fluxes Ha = {Ha_int_flux:.2e}, Hb = {Hb_int_flux:.2e}, the observed ratio is LOWER than theoretical ratio, so E(B-V) would be unphysical, so just assuming {EB_V}')
        else:
            EB_V = 1.97 * np.log10(obs_ratio / theoretical_ratio)
            print(f'Based on integrated fluxes Ha = {Ha_int_flux:.2e}, Hb = {Hb_int_flux:.2e}, determine E(B-V) = {EB_V:.2f}')
    else:
        EB_V = 0.
        print(f'Assumed E(B-V) = {EB_V}, due to absence of one or more of H alpha, H beta lines')

    return EB_V

# --------------------------------------------------------------------------------------------------------------------
def get_EB_V_map(full_hdu, args):
    '''
    Computes and returns the spatially resolved dust extinction map from a given HDU
    Based on Eqn 4 of Dominguez+2013 (https://iopscience.iop.org/article/10.1088/0004-637X/763/2/145/pdf)
    Returns the E(B-V) map
    '''

    Ha_map, Ha_wave = get_emission_line_map('Ha', full_hdu, args)
    Hb_map, Hb_wave = get_emission_line_map('Hb', full_hdu, args)

    obs_ratio_map = Ha_map / Hb_map
    EB_V_map = 1.97 * np.log10(obs_ratio_map / 2.86)

    return EB_V_map

# --------------------------------------------------------------------------------------------------------------------
def plot_dust_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the dust extinction map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D E(B-V) map just produced
    '''
    lim, label = [0, 1], 'E(B-V)'

    EB_V_map = get_EB_V_map(full_hdu, args)
    ax = plot_2D_map(EB_V_map, ax, args, label=label, cmap='YlOrBr', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

    return ax, EB_V_map

# --------------------------------------------------------------------------------------------------------------------
def plot_sfr_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the SFR map (and the emission line maps that go into it) in the given axes
    Conversion factor is from Kennicutt 1998 (Eq 2 of https://ned.ipac.caltech.edu/level5/Sept01/Rosa/Rosa3.html)
    Returns the axes handles and the 2D SFR density map just produced
    '''
    lim, label = [-4, -2], 'SFR'

    Ha_map, line_wave = get_emission_line_map('Ha', full_hdu, args)
    dered_Ha_map = get_dereddened_flux(Ha_map, line_wave, args.EB_V) # line_map in ergs/s/cm^2

    dered_Ha_map = dered_Ha_map * 4 * np.pi * (args.distance.to('cm').value) ** 2 # converting to ergs/s

    sfr_map = dered_Ha_map * 7.9e-42 # line_map in args/s; SFR in Msun/yr
    ax = plot_2D_map(np.log10(sfr_map), ax, args, label=label, cmap='Blues', radprof_ax=radprof_ax, vmin=lim[0], vmax=lim[1])

    return ax, sfr_map

# --------------------------------------------------------------------------------------------------------------------
def plot_Te_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the T_e map (and the emission line maps that go into it) in the given axes
    Conversion factor is from Nicholls+2017
    Returns the axes handles and the 2D T_e map just produced
    '''
    lim, label = [1, 7], r'T$_e$'

    OIII4363_map, line_wave = get_emission_line_map('OIII-4363', full_hdu, args)
    dered_OIII4363_map = get_dereddened_flux(OIII4363_map, line_wave, args.EB_V)

    OIII5007_map, line_wave = get_emission_line_map('OIII', full_hdu, args)
    dered_OIII5007_map = get_dereddened_flux(OIII5007_map, line_wave, args.EB_V)

    ratio_map = dered_OIII4363_map / dered_OIII5007_map
    logT_map = np.poly1d([0., 9.18962, 3.30355])(np.log10(ratio_map)) / np.poly1d([-0.00935, -0.15034, 2.09136, 1.000])(np.log10(ratio_map))
    Te_map = 10 ** logT_map
    ax = plot_2D_map(np.log10(Te_map), ax, args, label=label, cmap='OrRd_r', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

    return ax, Te_map

# --------------------------------------------------------------------------------------------------------------------
def plot_Z_Te_map(full_hdu, ax, args, Te_map, ne=1e3, radprof_ax=None):
    '''
    Plots the T_e-based metallicity map (and the emission line maps that go into it) in the given axes
    Conversion factor is from Nicholls+2017
    Returns the axes handles and the 2D metallicity map just produced
    '''
    lim, label = [6, 9], 'Z (Te)'
    t = Te_map * 1e-4
    x = 1e-4 * ne * np.sqrt(t)

    def poly(R, t, x, a, b, c, d, e):
        return np.log10(R) + a + b / t - c * np.log10(t) - d * t + np.log10(1 + e * x)  # eqn 3 I06 pattern

    OIII5007_map, line_wave = get_emission_line_map('OIII', full_hdu, args)
    dered_OIII5007_map = get_dereddened_flux(OIII5007_map, line_wave, args.EB_V)

    Hbeta_map, line_wave = get_emission_line_map('Hb', full_hdu, args)
    dered_Hbeta_map = get_dereddened_flux(Hbeta_map, line_wave, args.EB_V)

    OII3727_map, line_wave = get_emission_line_map('OII', full_hdu, args)
    dered_OII3727_map = get_dereddened_flux(OII3727_map, line_wave, args.EB_V)

    ratio1_map = dered_OII3727_map / dered_Hbeta_map
    ratio2_map = dered_OIII5007_map / dered_Hbeta_map
    log_O2H2_map = poly(ratio1_map, t, x, 5.961, 1.676, 0.4, 0.034, 1.35) - 12  # coefficients from eqn 3 I06
    log_O3H2_map = poly(ratio2_map, t, x, 6.200, 1.251, -0.55, -0.014, 0.0) - 12  # coefficients from eqn 5 I06
    log_OH_map = np.log10(10 ** log_O2H2_map + 10 ** log_O3H2_map) + 12

    ax = plot_2D_map(log_OH_map, ax, args, label=label, cmap='Greens', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

    return ax, log_OH_map

# --------------------------------------------------------------------------------------------------------------------
def plot_Z_R23_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the R23 metallicity map (and the emission line maps that go into it) in the given axes
    Conversion factor is from Kewley+2002
    Returns the axes handles and the 2D metallicity map just produced
    '''
    lim, label = [6, 9], 'Z (R23)'

    def poly(R23, k):
        return  (-k[1] + np.sqrt(k[1]**2 - 4*k[2]*(k[0] - R23)))/(2*k[2])

    OIII5007_map, line_wave = get_emission_line_map('OIII', full_hdu, args)
    dered_OIII5007_map = get_dereddened_flux(OIII5007_map, line_wave, args.EB_V)

    Hbeta_map, line_wave = get_emission_line_map('Hb', full_hdu, args)
    dered_Hbeta_map = get_dereddened_flux(Hbeta_map, line_wave, args.EB_V)

    OII3727_map, line_wave = get_emission_line_map('OII', full_hdu, args)
    dered_OII3727_map = get_dereddened_flux(OII3727_map, line_wave, args.EB_V)

    ratio_map = (dered_OII3727_map + dered_OIII5007_map) / dered_Hbeta_map
    R23 = np.log10(ratio_map)
    log_OH_map = poly(R23, [-44.7026, 10.8052, -0.640113]) #k0-2 parameters for q=8e7 from Table 3 of KD02 last row for q=8e7

    ax = plot_2D_map(log_OH_map, ax, args, label=label, cmap='Greens', radprof_ax=radprof_ax, hide_yaxis=True, vmin=lim[0], vmax=lim[1])

    return ax, log_OH_map

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    # ------------looping over the provided object IDs-----------------------
    for index, this_id in enumerate(args.id):
        start_time2 = datetime.now()
        print(f'\nCommencing ID {this_id} which is {index+1} of {len(args.id)}..')

        # ------determining directories---------
        extract_dir = args.input_dir / args.field / 'Extractions'
        pixscale_text = '' if args.pixscale == 0.04 else f'_{args.pixscale}arcsec_pix'
        output_subdir = args.output_dir / args.field / f'{this_id:05d}{pixscale_text}'
        output_subdir.mkdir(parents=True, exist_ok=True)
        full_fits_file = f'{args.field}_{this_id:05d}.full.fits'

        if os.path.exists(extract_dir / full_fits_file): # if the fits files are in Extractions/
            args.work_dir = extract_dir
        elif os.path.exists(output_subdir / full_fits_file): # if the fits files are in sub-directories for individual objects
            args.work_dir = output_subdir

        # ------------read in fits files--------------------------------
        od_hdu = fits.open(args.work_dir / f'{args.field}_{this_id:05d}.1D.fits')
        full_hdu = fits.open(args.work_dir / f'{args.field}_{this_id:05d}.full.fits')

        # ----------determining global parameters------------
        args.available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
        args.z = full_hdu[0].header['REDSHIFT']
        args.distance = cosmo.comoving_distance(args.z)
        args.ndfilt = full_hdu[0].header['NDFILT']
        args.nlines = full_hdu[0].header['NUMLINES']
        args.pix_arcsec = full_hdu[7].header['PIXASEC']

        line_wcs = pywcs.WCS(full_hdu['DSCI'].header)
        pix_size = utils.get_wcs_pscale(line_wcs)
        imsize_arcsec = full_hdu['DSCI'].data.shape[0] * pix_size
        dp = 0 # -0.5 * pix_size  # FITS reference is center of a pixel, array is edge
        args.extent = (-imsize_arcsec / 2. - dp, imsize_arcsec / 2. - dp, -imsize_arcsec / 2. - dp, imsize_arcsec / 2. - dp)
        args.EB_V = 0. # until gets over-written, if both H alpha and H beta lines are present

        if args.vorbin and args.voronoi_line is not None:
            line_index = np.where(np.array(args.available_lines.split(' ')) == args.voronoi_line)[0][0]
            ext = 5 + 2 * args.ndfilt + 4 * line_index
            line_map = full_hdu[ext].data * 1e-17  # this gives 200 x 200 array; in units of ergs/s/cm^2
            line_map_err = 1e-17 / full_hdu[ext + 3].data   # this gives 200 x 200 array; in units of ergs/s/cm^2
            args.voronoi_bin_IDs = get_voronoi_bin_IDs(line_map, line_map_err, args.voronoi_snr)#, plot=True, quiet=False)

        # ---------initialising the figure------------------------------
        if not args.keep: plt.close('all')
        lines_to_plot = ['Ha', 'Hb', 'OII', 'OIII-4363', 'OIII']
        nrow, ncol = 4 if args.plot_radial_profiles else 3, len(lines_to_plot)
        fig = plt.figure(layout='constrained', figsize=(12, 9) if args.plot_radial_profiles else (12, 6))
        '''
        groups = fig.add_gridspec(nrow, 1, hspace=0.1) # spacing between the groups
        group1 = groups[0].subgridspec(1, ncol, hspace=0, wspace=0.1) # spacing within the group
        group2 = groups[1].subgridspec(2, ncol, hspace=0, wspace=0.1)

        axis_dirimg = plt.subplot(group1.new_subplotspec((0, 0), colspan=1))
        axis_1dspec = plt.subplot(group1.new_subplotspec((0, 1), colspan=ncol - 1))
        [ax_em_lines, [ax_SFR, ax_EB_V, ax_Te, ax_Z_Te, ax_Z_R23]] = group2.subplots()

        if args.plot_radial_profiles:
            group3 = groups[2].subgridspec(1, ncol, hspace=0, wspace=0.1)
            [rax_SFR, rax_EB_V, rax_Te, rax_Z_Te, rax_Z_R23] = group3.subplots()
        '''
        axis_dirimg = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 0), colspan=1)
        axis_1dspec = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 1), colspan=ncol - 1)
        ax_em_lines = [plt.subplot2grid(shape=(nrow, ncol), loc=(1, item), colspan=1) for item in np.arange(ncol)]  # H alpha, H beta, OII, OIII-4363, OIII
        [ax_SFR, ax_EB_V, ax_Te, ax_Z_Te, ax_Z_R23] = [plt.subplot2grid(shape=(nrow, ncol), loc=(2, item), colspan=1) for item in np.arange(ncol)]  # SFR, E(B-V), Te, Z (Te), Z (R23)
        if args.plot_radial_profiles: [rax_SFR, rax_EB_V, rax_Te, rax_Z_Te, rax_Z_R23] = [plt.subplot2grid(shape=(nrow, ncol), loc=(3, item), colspan=1) for item in np.arange(ncol)]  # SFR, E(B-V), Te, Z (Te), Z (R23)

        # ---------direct imaging------------------------------
        axis_dirimg = plot_direct_image(full_hdu, axis_dirimg, args)

        # ---------1D spectra------------------------------
        axis_1dspec = plot_1d_spectra(od_hdu, axis_1dspec, args)

        # -----------------emission line maps---------------
        for ind, line in enumerate(lines_to_plot):
            if line in args.available_lines: _, _, ax_em_lines[ind] = plot_emission_line_map(line, full_hdu, ax_em_lines[ind], args, cmap='BuPu', vmin=-20, vmax=-18, hide_xaxis=True, hide_yaxis=ind > 0, hide_cbar=False) #ind != len(lines_to_plot) - 1) # line_map in ergs/s/cm^2
            else: fig.delaxes(ax_em_lines[ind])

        # ---------------dust map---------------
        if all([line in args.available_lines for line in ['Ha', 'Hb']]):
            args.EB_V = calculate_EB_V(full_hdu, args)
            ax_EB_V, dust_map = plot_dust_map(full_hdu, ax_EB_V, args, radprof_ax=rax_EB_V)
        else:
            fig.delaxes(ax_EB_V)
            if args.plot_radial_profiles: fig.delaxes(rax_EB_V)

        # ---------------SFR map------------------
        if 'Ha' in args.available_lines:
            ax_SFR, sfr_map = plot_sfr_map(full_hdu, ax_SFR, args, radprof_ax=rax_SFR)
        else:
            fig.delaxes(ax_SFR)
            if args.plot_radial_profiles: fig.delaxes(rax_SFR)

        # ---------------electron temperature map---------------
        if all([line in args.available_lines for line in ['OIII-4363', 'OIII']]):
            ax_Te, Te_map = plot_Te_map(full_hdu, ax_Te, args, radprof_ax=rax_Te)
        else:
            fig.delaxes(ax_Te)
            if args.plot_radial_profiles: fig.delaxes(rax_Te)

        # ---------------metallicity maps---------------
        if all([line in args.available_lines for line in ['OIII', 'OII', 'Hb']]):
            ax_Z_Te, log_OH_Te_map = plot_Z_Te_map(full_hdu, ax_Z_Te, args, Te_map, radprof_ax=rax_Z_Te)
            ax_Z_R23, log_OH_R23_map = plot_Z_R23_map(full_hdu, ax_Z_R23, args, radprof_ax=rax_Z_R23)
        else:
            fig.delaxes(ax_Z_Te)
            fig.delaxes(ax_Z_R23)
            if args.plot_radial_profiles:
                fig.delaxes(rax_Z_Te)
                fig.delaxes(rax_Z_R23)

        # ---------decorating and saving the figure------------------------------
        fig.text(0.05, 0.98, f'{args.field}: ID {this_id}', fontsize=args.fontsize, c='k', ha='left', va='top')

        figname = output_subdir / f'{args.field}_{this_id:05d}_all_diag_plots.png'
        fig.savefig(figname)
        print(f'Saved figure at {figname}')
        plt.show(block=False)

        print(f'Completed id {this_id} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(args.id) - index - 1} to go!')

    os.chdir(args.code_dir)
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
