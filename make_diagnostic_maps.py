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

             run make_diagnostic_maps.py --field Par28 --id 1303,1934,2734,2867,300,2903,2906 --plot_AGN_frac --plot_radial_profile --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --drv 0.5 --Zdiag O3S2 --keep
             run make_diagnostic_maps.py --field Par28 --id 1303 --plot_AGN_frac --mask_agn --plot_radial_profile --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --drv 0.5 --Zdiag O3S2 --keep
             run make_diagnostic_maps.py --field Par28 --id 1303 --plot_AGN_frac --plot_radial_profile --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --drv 0.5 --Zdiag O3S2 --plot_circle_at_arcsec 0.5
             run make_diagnostic_maps.py --field Par28 --id 1303 --plot_snr --plot_ratio_maps --plot_AGN_frac --plot_radial_profile --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --drv 0.5 --Zdiag O3S2 --plot_circle_at_arcsec 0.5
             run make_diagnostic_maps.py --field Par28 --id 2867 --plot_BPT --plot_AGN_frac --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --drv 0.5 --plot_circle_at_arcsec 0.25 --colorcol distance_from_AGN_line --AGN_diag H21
             run make_diagnostic_maps.py --field Par28 --id 1303 --plot_BPT --plot_AGN_frac --plot_radial_profiles --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --drv 0.5 --plot_circle_at_arcsec 0.5 --AGN_diag O2O3
             run make_diagnostic_maps.py --field Par28 --id 1303 --plot_BPT --plot_AGN_frac --plot_radial_profiles --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --drv 0.5 --plot_circle_at_arcsec 0.5 --use_variable_N2Ha
             run make_diagnostic_maps.py --field glass-a2744 --id 2928,5184 --plot_radial_profile --plot_AGN_frac --only_seg --vorbin --voronoi_line OIII --voronoi_snr 10 --drv 0.5 --Zdiag O3O2
             run make_diagnostic_maps.py --field Par28 --id 1303 --plot_BPT --plot_AGN_frac --plot_radial_profiles --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --drv 0.5 --plot_circle_at_arcsec 0.5 --AGN_diag H21 --plot_models --slice_at_quantity3 4,9

             run make_diagnostic_maps.py --field glass-a2744 --id 321,998,1694,1983,2355,2744,2938 --plot_snr --plot_ratio_maps --plot_radial_profiles --plot_AGN_frac --plot_radial_profiles --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --drv 0.5 --plot_circle_at_arcsec 0.5 --fontsize 5 --arcsec_limit 0.5
             run make_diagnostic_maps.py --field glass-a2744 --plot_snr --plot_ratio_maps --plot_radial_profiles --plot_AGN_frac --only_seg --vorbin --voronoi_line NeIII-3867 --voronoi_snr 2 --drv 0.5 --plot_circle_at_arcsec 0.25 --AGN_diag Ne3O2 --Zdiag R3 --mask_agn --fontsize 5 --id 5184,1407,629,1504,1622,2273,2926,3422,2891,562,1300,2928,600 --output_subdir diagnostics_OIII,OII,Hb,NeIII-3867,SNR>2.0
             run make_diagnostic_maps.py --field Par28 --id 1303 --plot_DIG --plot_radial_profile --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --drv 0.5 --Zdiag P25 --plot_circle_at_arcsec 0.5 --plot_snr

             run make_diagnostic_maps.py --field Par28 --id 300,1303,1849,2171,2727,2867 --plot_ratio_maps --plot_snr --plot_AGN_frac --plot_radial_profile --only_seg --vorbin --voronoi_line NeIII-3867 --voronoi_snr 2 --drv 0.5 --Zdiag R3 --AGN_diag Ne3O2 --mask_agn --fontsize 5
             run make_diagnostic_maps.py --field Par28 --id 300,1303,1849,2171,2727,2867 --plot_metallicity --plot_radial_profile --only_seg --vorbin --voronoi_line NeIII-3867 --voronoi_snr 4 --drv 0.5 --Zdiag R3,R2,R23,O3O2,P25,NB --AGN_diag Ne3O2 --mask_agn --exclude_line SII
   Afterwards, to make the animation: run /Users/acharyya/Work/astro/ayan_codes/animate_png.py --inpath /Volumes/Elements/acharyya_backup/Work/astro/passage/passage_output/Par028/all_diag_plots_wradprof_snr3.0_onlyseg/ --rootname Par028_*_all_diag_plots_wradprof_snr3.0_onlyseg.png --delay 0.1
'''

from header import *
from util import *
from matplotlib import cm as mpl_cm
import imageio
from get_field_stats import get_crossmatch_with_cosmos
from plot_mappings_grid import plot_ratio_grid

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def annotate_PAs(pa_arr, ax, fontsize=10, color='k'):
    '''
    Annotates a given plot with the PAs of all available filters in the given axis
    Returns the axis handle
    '''
    x_cen = ax.get_xlim()[1] - 0.15 * np.diff(ax.get_xlim())[0]
    y_cen = ax.get_ylim()[1] - 0.15 * np.diff(ax.get_ylim())[0]
    len = np.diff(ax.get_xlim())[0] * 0.1

    for pa in pa_arr:
        x_comp = len * np.sin(pa * np.pi / 180)
        y_comp = len * np.cos(pa * np.pi / 180)
        ax.plot([x_cen, x_cen - x_comp], [y_cen, y_cen + y_comp], lw=1, c=color)
        ax.text(x_cen - 1.2 * x_comp, y_cen + 1.2 * y_comp, r'%d$^\circ$' % pa, color=color, fontsize=fontsize, ha='center', va='center', rotation=pa)

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

        filter_hdu = full_hdu['DSCI', f'{filt.upper()}']
        image = filter_hdu.data
        image = trim_image(image, args, skip_re_trim=True)

        extent = (-args.arcsec_limit, args.arcsec_limit, -args.arcsec_limit, args.arcsec_limit)
        p = ax.imshow(image, cmap=cmap_arr[index], origin='lower', extent=extent, alpha=1)#, vmin=0, vmax=0.03)

        ax.set_xlim(-args.arcsec_limit, args.arcsec_limit)  # arcsec
        ax.set_ylim(-args.arcsec_limit, args.arcsec_limit)  # arcsec

        textcolor = mpl_cm.get_cmap(cmap_arr[index])(0.9)
        ax.text(ax.get_xlim()[0] * 0.9, ax.get_ylim()[1] * 0.7 - index * 0.1, filt, c=textcolor, fontsize=args.fontsize, ha='left', va='top')

   # ----------annotate axis---------------
    if args.re_limit is not None:
        limit = args.re_limit * args.re_arcsec
        ax.add_patch(plt.Rectangle((-limit, -limit), 2 * limit, 2 * limit, lw=0.5, color='r', fill=False))
        #ax.add_patch(plt.Circle((0, 0), args.re_arcsec, color='k', fill=False, lw=0.5))

    ax.text(ax.get_xlim()[0] * 0.9, ax.get_ylim()[1] * 0.95, f'z={args.z:.2f}', c='k', fontsize=args.fontsize, ha='left', va='top')
    ax.text(ax.get_xlim()[1] * 0.95, ax.get_ylim()[0] * 0.95, f'Mag={args.mag:.1f}', c='k', fontsize=args.fontsize, ha='right', va='bottom')
    ax.scatter(0, 0, marker='x', s=10, c='grey')
    ax = annotate_PAs(args.pa_arr, ax, fontsize=args.fontsize)
    #cbar = plt.colorbar(p)

    if args.only_seg:
        ax.contour(args.segmentation_map != args.id, levels=0, colors='k', extent=extent, linewidths=0.5)

    if args.plot_circle_at_arcsec is not None:
        ax.add_patch(plt.Circle((0, 0), args.plot_circle_at_arcsec, color='r', fill=False, lw=0.5)) # additional circle for debugging purpose

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
    line_list_file = HOME / 'Work/astro/Mappings/P_spherical/sp_P70_a05modelfiles/Q700/spec0003.csv'

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
def get_linelist(wave_lim=None, line_list_file=None):
    '''
    Reads in an emission line list
    Returns list of lines that are within the given wavelength limits, as a pandas dataframe
    '''
    if line_list_file is None: line_list_file = HOME / 'Work/astro/Mappings/labframe.shortlinelist'

    lines_df = pd.read_table(line_list_file, comment='#', delim_whitespace=True)
    lines_df = lines_df[~lines_df['LineID'].str.contains('Fe')] # not interested in the numerous Fe lines
    if wave_lim is not None: lines_df = lines_df[lines_df['restwave'].between(wave_lim[0], wave_lim[1])]

    print(f'\nFound {len(lines_df)} lines in this wavelength regime from {line_list_file}')

    return lines_df

# --------------------------------------------------------------------------------------------------------------------
def plot_linelist(ax, fontsize=10, line_list_file=None):
    '''
    Plots a list of emission line wavelengths on the given axis
    Returns axis handle
    '''
    lines_df = get_linelist(wave_lim=ax.get_xlim(), line_list_file=line_list_file)
    print(f'over-plotting lines now..')

    for index in range(len(lines_df)):
        ax.axvline(lines_df.iloc[index]['restwave'], c='cornflowerblue', lw=1)
        xpos = lines_df.iloc[index]['restwave'] + np.diff(ax.get_xlim())[0] * 0.01
        ypos = ax.get_ylim()[1] * 0.98 if index % 2 else 0.02 + ax.get_ylim()[0] * 1.02
        ax.text(xpos, ypos, lines_df.iloc[index]['LineID'].strip(), rotation=90, va='top' if index % 2 else 'bottom', ha='left', fontsize=fontsize)

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
    ax.set_ylabel(r'f$_{\lambda}$ ' + '(%.0e ' % norm_factor + r'ergs/s/cm$^2$/A)', fontsize=args.fontsize)
    ax.set_ylim(0, args.flam_max) # flam_max should be in units of 1e-19 ergs/s/cm^2/A

    # ---observed wavelength axis-------
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(ax.get_xticks())
    ax2.set_xticklabels(['%.2F' % (item * (1 + args.z) / 1e4) for item in ax2.get_xticks()], fontsize=args.fontsize)
    ax2.set_xlabel(r'Observed wavelength ($\mu$)', fontsize=args.fontsize)

    # ---vertical lines for emission line wavelengths------
    if args.plot_mappings: ax = plot_MAPPINGS_lines(ax)
    else: ax = plot_linelist(ax, fontsize=args.fontsize)

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
def get_distance_map(image_shape, args, for_distmap=False):
    '''
    Get map of distances from the center, in target rest-frame, on a given 2D grid
    Returns 2D distance map
    '''
    if not args.vorbin or for_distmap:
        pixscale_kpc = args.pix_size_arcsec/ cosmo.arcsec_per_kpc_proper(args.z).value # kpc
        center_pix = image_shape[0] / 2.

        if args.use_elliptical_bins:
            y, x = np.indices(image_shape)
            x_shift = x - center_pix 
            y_shift = y - center_pix

            # Rotate coordinates to align with major axis
            x_rot = x_shift * np.cos(args.pa) + y_shift * np.sin(args.pa) # args.pa is in radians
            y_rot = -x_shift * np.sin(args.pa) + y_shift * np.cos(args.pa)
            
            distance_map = np.sqrt(x_rot**2 + (y_rot / args.q)**2) * pixscale_kpc # kpc
        else:
            distance_map = np.array([[np.sqrt((i - center_pix)**2 + (j - center_pix)**2) for j in range(image_shape[1])] for i in range(image_shape[0])]) * pixscale_kpc # kpc
    else:
        distance_map = args.voronoi_bin_distances

    return distance_map

# --------------------------------------------------------------------------------------------------------------------
def plot_radial_profile(image, ax, args, label=None, ymin=None, ymax=None, hide_xaxis=False, hide_yaxis=False, image_err=None, metallicity_multi_color=False, xcol='radius', ycol='data', color='darkorange', fontsize=None):
    '''
    Plots the average radial profile for a given 2D map in the given axis
    Returns the axis handle
    '''
    if fontsize is None: fontsize = args.fontsize
    print(f'Plotting radial profile of {label}..')
    if type(image) == pd.core.frame.DataFrame:
        df  = image
    else:
        distance_map = get_distance_map(np.shape(image), args)
        try: distance_map = np.ma.masked_where(image.mask, distance_map)
        except AttributeError: distance_map = np.ma.masked_where(False, distance_map)

        # ----making the dataframe before radial profile plot--------------
        df = pd.DataFrame({xcol: np.ma.compressed(distance_map), ycol: np.ma.compressed(image)})
        if args.re_limit is not None: df[xcol] = df[xcol] / args.re_kpc # converting from kpc to Re units
        if image_err is not None: df[ycol + '_err'] = np.ma.compressed(image_err)
        if metallicity_multi_color: df['agn_dist'] = np.ma.compressed(np.ma.masked_where(image.mask, args.distance_from_AGN_line_map.data))
        if len(df[df[ycol + '_err'] > 0]) == 0: image_err = None

        # --------processing the dataframe in case voronoi binning has been performed and there are duplicate data values------
        if args.vorbin:
            df['bin_ID'] = np.ma.compressed(np.ma.masked_where(image.mask, args.voronoi_bin_IDs.data))
            df = df.groupby('bin_ID', as_index=False).agg(np.mean)
        
        if args.radius_max is not None: df = df[df[xcol] <= args.radius_max]
        df = df.sort_values(by=xcol)

    # -------proceeding with plotting--------
    ax.scatter(df[xcol], df[ycol], c=df['agn_dist'] if metallicity_multi_color else 'grey', cmap=args.diverging_cmap if metallicity_multi_color else None, s=20 if args.vorbin else 1, alpha=1 if args.vorbin else 0.2)
    if ycol + '_err' in df: ax.errorbar(df[xcol], df[ycol], yerr=df[ycol + '_err'], c='grey', fmt='none', lw=2 if args.vorbin else 0.5, alpha=0.2 if args.vorbin else 0.1)
    for index, row in df.iterrows(): ax.text(row[xcol], row[ycol], row['bin_ID'].astype(int), c='k', fontsize=args.fontsize/1.5)

    if args.radius_max is not None: ax.set_xlim(0, args.radius_max) # kpc
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
        ax.set_xlabel('Distance (kpc)' if args.re_limit is None else r'Galactocentric distance (R/R$_e$)', fontsize=fontsize)
        #ax.set_xticklabels(['%d' % item / cosmo.arcsec_per_kpc_proper(args.z).value for item in ax.get_xticks()], fontsize=fontsize)
        ax.tick_params(axis='x', which='major', labelsize=fontsize)

    if hide_yaxis:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel(label, fontsize=args.fontsize)
        ax.tick_params(axis='y', which='major', labelsize=fontsize)

    return ax, linefit

# --------------------------------------------------------------------------------------------------------------------
def plot_2D_map(image, ax, args, takelog=True, label=None, cmap=None, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=False, radprof_ax=None, vorbin_ax=None, snr_ax=None, image_err=None, metallicity_multi_color=False, fontsize=None, plot_circles_at=[]):
    '''
    Plots the emission map for a given line in the given axis
    Returns the axis handle
    '''
    if fontsize is None: fontsize = args.fontsize
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

    orig_mask = image.mask
    if takelog:
        new_mask = image <= 0
        image[new_mask] = 1e-9 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
        image_log = unp.log10(unp.uarray(image, image_err))
        image_log = np.ma.masked_where(new_mask | image.mask, image_log)
        image = unp.nominal_values(image_log)
        image_err = unp.std_devs(image_log)

    if metallicity_multi_color:
        image_sfr = np.ma.masked_where(args.distance_from_AGN_line_map > 0, image)
        image_agn = np.ma.masked_where(args.distance_from_AGN_line_map < 0, image)
        p = ax.imshow(image_agn, cmap='pink', origin='lower', extent=args.extent, vmin=vmin, vmax=vmax)
        p = ax.imshow(image_sfr, cmap='cool', origin='lower', extent=args.extent, vmin=vmin, vmax=vmax)
    else:
        p = ax.imshow(image, cmap=cmap, origin='lower', extent=args.extent, vmin=vmin, vmax=vmax)

    if not args.no_text_on_plot: ax.text(ax.get_xlim()[0] * 0.88, ax.get_ylim()[1] * 0.88, label, c='k', fontsize=fontsize if args.arcsec_limit >= 1 else fontsize/1.5, ha='left', va='top', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
    ax.scatter(0, 0, marker='x', s=10, c='grey')

    offset_units = 'kpc' if args.plot_target_frame else r'R$_e$' if args.re_limit is not None else '"'
    if hide_xaxis:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel(f'Offset ({offset_units})', fontsize=fontsize)
        ax.tick_params(axis='x', which='major', labelsize=fontsize)

    if hide_yaxis:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel(f'Offset ({offset_units})', fontsize=fontsize)
        ax.tick_params(axis='y', which='major', labelsize=fontsize)

    if not hide_cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad='2%')
        cbar = plt.colorbar(p, cax=cax, orientation='vertical')
        cbar.ax.tick_params(labelsize=fontsize)

    if args.plot_radial_profiles and radprof_ax is not None:
        radprof_ax, radprof_fit = plot_radial_profile(image, radprof_ax, args, label=label.split(r'$_{\rm int}')[0], ymin=vmin, ymax=vmax, image_err=image_err, metallicity_multi_color=metallicity_multi_color, fontsize=fontsize)
        radius_max_arcsec = radprof_ax.get_xlim()[0] * cosmo.arcsec_per_kpc_proper(args.z).value # converting args.radius_max (in kpc) to arcsec
        if radius_max_arcsec <= args.arcsec_limit:
            circle = plt.Circle((0, 0), radius_max_arcsec, color='k', fill=False, lw=0.5)
            ax.add_patch(circle)
    else:
        radprof_fit = [np.nan, np.nan] # dummy values for when the fit was not performed

    if args.plot_circle_at_arcsec is not None: ax.add_patch(plt.Circle((0, 0), args.plot_circle_at_arcsec, color='r', fill=False, lw=0.5)) # additional circle for debugging purpose
    if 're_arcsec' in args: ax.add_patch(plt.Circle((0, 0), args.re_arcsec if args.re_limit is None else 1., color='w', fill=False, lw=0.5))
    if 'segmentation_map' in args: ax.contour(args.segmentation_map != args.id, levels=0, colors='w' if args.fortalk else 'k', extent=args.extent, linewidths=0.5) # demarcating the segmentation map zone
    #for this_circle in plot_circles_at: ax.add_patch(plt.Circle((0, 0), this_circle, color='sienna', fill=False, lw=1))
    for this_circle in plot_circles_at: ax.add_patch(matplotlib.patches.Ellipse((0, 0), width = 2 * this_circle, height = 2 * this_circle * args.q, angle=np.degrees(args.pa), color='sienna', fill=False, lw=1))
    
    if args.vorbin and args.plot_vorbin and vorbin_ax is not None:
        vorbin_IDs = args.voronoi_bin_IDs
        vorbin_IDs = np.ma.masked_where(image.mask, vorbin_IDs)
        _, _ = plot_2D_map(vorbin_IDs, vorbin_ax, args, takelog=False, label=label + ' vorbin', cmap='rainbow', fontsize=fontsize)

    if args.plot_snr and image_err is not None and snr_ax is not None:
        if takelog:
            quant = 10 ** unp.uarray(image.data, image_err.data) # undoing the log step that was done initially
            image = unp.nominal_values(quant)
            image_err = unp.std_devs(quant)
        snr_map = image / image_err
        snr_map = np.ma.masked_where(orig_mask, snr_map)
        _, _ = plot_2D_map(snr_map, snr_ax, args, takelog=False, label=label.split(r'$_{\rm int}')[0] + ' SNR', cmap='cividis', vmin=0, vmax=2 if '/' in label else 8, hide_xaxis=hide_xaxis, hide_yaxis=True, hide_cbar=hide_cbar, fontsize=fontsize)

    return ax, radprof_fit

# --------------------------------------------------------------------------------------------------------------------
def bin_fluxes(fluxes, errors=None):
    '''
    Bin the fluxes (and errors if provided) in the pixels that belong to a single voronoi bin (either taking mean or weighted mean)
    Returns binned flux
    '''
    fluxes = np.atleast_1d(fluxes)
    if errors is None: errors = np.zeros(len(fluxes))
    else: errors = np.atleast_1d(errors)
    
    #print(f'Deb475: fluxes={fluxes}, errors={errors}') ##
    binned_flux = np.average(unp.uarray(fluxes, errors))#, weights = 1 / (errors ** 2))
    #print(f'Deb477: result={binned_flux}') ##
 
    return binned_flux

# --------------------------------------------------------------------------------------------------------------------
def bin_2D(map, bin_IDs, map_err=None, debug_vorbin=False, ax=None):
    '''
    Bin a given 2D map by given bin_IDs
    Returns the binned 2D map (of same shape as input map)
    '''
    binned_map = np.zeros(np.shape(map))
    if map_err is not None: binned_map_err = np.zeros(np.shape(map_err))

    unique_IDs = np.unique(np.ma.compressed(bin_IDs))
    col_arr = plt.cm.viridis(np.linspace(0, 1, len(unique_IDs)))

    for index, id in enumerate(unique_IDs):
        all_pixel_fluxes = map[bin_IDs == id] # this is in ergs/s/cm^2
        good_pixel_fluxes = np.ma.compressed(all_pixel_fluxes)
        all_pixel_fluxes = all_pixel_fluxes.data

        if map_err is not None:
            all_pixel_errors = map_err[bin_IDs == id] # this is in ergs/s/cm^2
            good_pixel_errors = np.ma.compressed(all_pixel_errors)
            all_pixel_errors = all_pixel_errors.data
        else:
            good_pixel_errors, all_pixel_errors = None, None

        all_pixels_binned = bin_fluxes(all_pixel_fluxes, errors=all_pixel_errors)
        good_pixels_binned = bin_fluxes(good_pixel_fluxes, errors=good_pixel_errors)
        chosen_data = all_pixels_binned # choose from [all_pixels_binned, good_pixels_binned] #
        if debug_vorbin:
            print(f'Deb445: id {int(id)}/{len(unique_IDs)}: all pixels: {len(all_pixel_fluxes)}, {all_pixels_binned:.1e}, {all_pixels_binned.n / all_pixels_binned.s:.1f}; good pixels: {len(good_pixel_fluxes)}, {good_pixels_binned:.1e}, {good_pixels_binned.n / good_pixels_binned.s:.1f}; snr = {chosen_data.n / chosen_data.s:.1f}') ##
            if ax is not None:
                all_pixels_snr = all_pixel_fluxes / all_pixel_errors
                hist, bin_edges = np.histogram(all_pixels_snr, bins=50, density=True, range=(-5, 5))
                bin_centers = bin_edges[:-1] + np.diff(bin_edges)
                ax.step(bin_centers, hist + index * 1., lw=0.5, where='mid', color=col_arr[index])

        binned_map[bin_IDs == id] = chosen_data.n
        if map_err is not None: binned_map_err[bin_IDs == id] = chosen_data.s

    binned_map = np.ma.masked_where(bin_IDs.mask, binned_map)
    if map_err is not None: binned_map_err = np.ma.masked_where(bin_IDs.mask, binned_map_err)

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
    elif cutoff_frac > 1: cutoff_frac = 1.
    #col_above = plt.get_cmap(cmap1)(np.linspace(cutoff_frac, 1, int(256 * (1 - cutoff_frac))))
    #col_below = plt.get_cmap(cmap2)(np.linspace(0, cutoff_frac, int(256 * cutoff_frac)))
    col_above = plt.get_cmap(cmap1)(np.linspace(0, 1, int(256 * (1 - cutoff_frac))))
    col_below = plt.get_cmap(cmap2)(np.linspace(0, 1, int(256 * cutoff_frac)))
    colors = np.vstack((col_below, col_above))
    combined_cmap = mplcolors.LinearSegmentedColormap.from_list('combined_cmap', colors)

    return combined_cmap

# --------------------------------------------------------------------------------------------------------------------
def get_voronoi_bin_distances(full_hdu, line, args):
    '''
    Compute the galactocentric distances of Voronoi bins as the luminosity-weighted average of the pixels distances, based on the given line map
    Returns the 2D map (of same shape as input args.voronoi_bin_IDs) with just the distances
    '''
    print(f'Computing Voronoi bin galactocentric distances using {line}-map as weight..')
    line_map, _, _, _, _ = get_emission_line_map(line, full_hdu, args, for_vorbin=True)

    pixel_distance_map = get_distance_map(np.shape(line_map), args, for_distmap=True)
    try: pixel_distance_map = np.ma.masked_where(line_map.mask, pixel_distance_map)
    except AttributeError: pixel_distance_map = np.ma.masked_where(False, pixel_distance_map)

    # ----making the dataframe and computing luminosity-weighted bin distances--------------
    df = pd.DataFrame({'pixel_distance': np.ma.compressed(pixel_distance_map), 'flux': unp.nominal_values(np.ma.compressed(line_map)), 'flux_u': unp.std_devs(np.ma.compressed(line_map)), 'bin_ID': np.ma.compressed(np.ma.masked_where(line_map.mask, args.voronoi_bin_IDs.data))})
    df['snr'] = df['flux'] / df['flux']
    
    weight_col = 'snr' # 'snr' or 'flux'
    wm = lambda x: np.average(x[df[weight_col] >= 0], weights=df.loc[x[df[weight_col] >= 0].index, weight_col]) # weighting the pixel distance by the SNR of the pixel
    df2 = df.groupby(['bin_ID']).agg(bin_distance=('pixel_distance', wm)).reset_index()
    binID_distance_dict = dict(zip(df2['bin_ID'], df2['bin_distance']))

    bin_distance_map = np.zeros(np.shape(line_map))
    for xind in range(np.shape(line_map)[0]):
        for yind in range(np.shape(line_map)[1]):
            try: bin_distance_map[xind][yind] = binID_distance_dict[args.voronoi_bin_IDs.data[xind][yind]]
            except KeyError: bin_distance_map[xind][yind] = np.nan

    bin_distance_map = np.ma.masked_where(args.voronoi_bin_IDs.mask, bin_distance_map)

    return bin_distance_map

# --------------------------------------------------------------------------------------------------------------------
def get_voronoi_bin_IDs(full_hdu, snr_thresh, plot=False, quiet=True, args=None):
    '''
    Compute the Voronoi bin IDs a given 2D map and corresponding uncertainty and SNR threshold
    Returns the 2D map (of same shape as input map) with just the IDs
    '''
    # ---------getting the spatially resolved line flux map----------------
    try:
        line_hdu = full_hdu['LINE', args.voronoi_line]
        initial_map_err = 1e-17 / (full_hdu['LINEWHT', args.voronoi_line].data ** 0.5)  # ERR = 1/sqrt(LINEWHT) = flux uncertainty; in units of ergs/s/cm^2
    except KeyError:
        if args.voronoi_line == 'OIII':
            line_hdu = full_hdu['LINE', 'OIII-5007']
            initial_map_err = 1e-17 / (full_hdu['LINEWHT', 'OIII-5007'].data ** 0.5)  # ERR = 1/sqrt(LINEWHT) = flux uncertainty; in units of ergs/s/cm^2

    initial_map = line_hdu.data * 1e-17 # in units of ergs/s/cm^2

    # -------------pixel offset to true center----------------
    print(f'Correcting emission lines for pixel offset by {args.ndelta_xpix} on x and {args.ndelta_ypix} on y')
    initial_map = np.roll(initial_map, args.ndelta_xpix, axis=0)
    initial_map_err = np.roll(initial_map_err, args.ndelta_xpix, axis=0)
    initial_map = np.roll(initial_map, args.ndelta_ypix, axis=1)
    initial_map_err = np.roll(initial_map_err, args.ndelta_ypix, axis=1)

    # ----------getting a smaller cutout around the object center-----------
    initial_map = trim_image(initial_map, args)
    initial_map_err = trim_image(initial_map_err, args)
    initial_snr = initial_map / initial_map_err
    if args.only_seg: seg_mask = args.segmentation_map != args.id
    else: seg_mask = False
    initial_map = np.ma.masked_where(seg_mask, initial_map)
    initial_map_err = np.ma.masked_where(seg_mask, initial_map_err)
    initial_snr = np.ma.masked_where(seg_mask, initial_snr)

    # -------debugging checks: initial maps-----------
    if args is not None and args.debug_vorbin:
        fig, axes_all = plt.subplots(3, 4, figsize=(11, 7))
        fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.07, wspace=0.3, hspace=0.1)
        cmap = 'viridis'
        axes_top = axes_all[0]
        fig.text(0.05, 0.98, f'{args.field}: ID {args.id}: Voronoi binning diagnostics', fontsize=args.fontsize, c='k', ha='left', va='top')

        snr_min, snr_max = np.nanmin(initial_snr.data), np.nanmax(initial_snr.data)
        combined_cmap = get_combined_cmap([snr_min, 0, snr_thresh, snr_max], ['Greys', 'Reds', cmap], new_name='my_colormap')
        plot_2D_map(initial_map, axes_top[0], args, takelog=False, label=f'initial {args.voronoi_line} map', cmap=cmap, hide_yaxis=False, hide_xaxis=True)
        plot_2D_map(initial_map_err, axes_top[1], args, takelog=False, label=f'initial {args.voronoi_line} err', cmap=cmap, hide_yaxis=True, hide_xaxis=True)
        plot_2D_map(initial_snr, axes_top[2], args, takelog=False, label='initial SNR', cmap=combined_cmap, hide_yaxis=True, vmin=snr_min, vmax=snr_max, hide_xaxis=True)

    # ------smoothing the map before voronoi binning--------
    smoothing_kernel = Box2DKernel(args.kernel_size, mode=args.kernel_mode)
    smoothed_map = np.ma.masked_where(initial_map.mask, convolve(initial_map.data, smoothing_kernel))
    smoothed_map_err = np.ma.masked_where(initial_map_err.mask, convolve(initial_map_err.data, smoothing_kernel))

    # -------assigning input maps for voronoi binning------------
    #input_map, input_map_err = initial_map, initial_map_err
    input_map, input_map_err = smoothed_map, smoothed_map_err
    input_snr = input_map / input_map_err
    
    # -------applying snr mask before sending arrays to vorbin-----------
    snr_cut_for_vorbin = 0
    bad_mask = (input_snr < snr_cut_for_vorbin) | seg_mask
    #bad_mask = seg_mask
    input_map = np.ma.masked_where(bad_mask, input_map)
    input_map_err = np.ma.masked_where(bad_mask, input_map_err)

    # ------------debugging checks: input maps---------------------- 
    if args is not None and args.debug_vorbin:
        axes_mid = axes_all[1]

        snr_min, snr_max = np.nanmin(input_snr.data), np.nanmax(input_snr.data)
        combined_cmap = get_combined_cmap([snr_min, 0, snr_thresh, snr_max], ['Greys', 'Reds', cmap], new_name='my_colormap')
        plot_2D_map(input_map, axes_mid[0], args, takelog=False, label=f'input {args.voronoi_line} map', cmap=cmap, hide_yaxis=False, hide_xaxis=True)
        plot_2D_map(input_map_err, axes_mid[1], args, takelog=False, label=f'input {args.voronoi_line} err', cmap=cmap, hide_yaxis=True, hide_xaxis=True)
        plot_2D_map(input_snr, axes_mid[2], args, takelog=False, label='input SNR', cmap=combined_cmap, hide_yaxis=True, vmin=snr_min, vmax=snr_max, hide_xaxis=True)

    # -------making appropriate x and y coords for sending into to vorbin-----------
    x_size, y_size = np.shape(input_map)
    x_coords_grid = np.reshape(np.repeat(np.arange(x_size), y_size), (x_size, y_size))
    y_coords_grid = np.reshape(np.tile(np.arange(y_size), x_size), (x_size, y_size))
    x_coords_grid_masked = np.ma.masked_where(input_map.mask, x_coords_grid)
    y_coords_grid_masked = np.ma.masked_where(input_map.mask, y_coords_grid)

    input_map_array = np.ma.compressed(input_map)
    input_map_err_array = np.ma.compressed(input_map_err)
    x_coords_array = np.ma.compressed(x_coords_grid_masked)
    y_coords_array = np.ma.compressed(y_coords_grid_masked)

    if len(x_coords_array) == 0:
        print(f'No unmasked pixels to Voronoi bin on..so skipping this galaxy.')
        return None

    # -------actually calling vorbin with non-negative fluxes and errors-----------
    binIDs, _, _, _, _, _, _, _ = voronoi_2d_binning(x_coords_array, y_coords_array, input_map_array, input_map_err_array, snr_thresh, plot=False, quiet=quiet, cvt=True, wvt=False)

    # -------interpolating binIDs derived from vorbin----------
    interp = NearestNDInterpolator(list(zip(x_coords_array, y_coords_array)), binIDs)
    binID_map = interp(x_coords_grid, y_coords_grid)
    binID_map[seg_mask] = -1
    binID_map = np.ma.masked_where(seg_mask, binID_map)

    # --------extracting only those bins with SNR above a certain threshold for ALL relevant lines--------------
    relevant_lines = [args.voronoi_line] #  ['OII', 'Hb', 'OIII'] # 
    final_snr_cut = snr_thresh # snr_thresh OR 0 OR -100 (some very low number)
    bad_mask = seg_mask
    for index, line in enumerate(relevant_lines):
        line_map2, _, _, _, _ = get_emission_line_map(line, full_hdu, args, for_vorbin=True)
        map2 = np.ma.masked_where(line_map2.mask, unp.nominal_values(line_map2.data))
        map2_err = np.ma.masked_where(line_map2.mask, unp.std_devs(line_map2.data))
        map2, map2_err = bin_2D(map2, binID_map, map_err=map2_err, debug_vorbin=False)    
        snr = map2 / map2_err
        bad_mask = bad_mask | (snr < final_snr_cut)
    binID_map = np.ma.masked_where(bad_mask, binID_map.data)
    if args is not None and args.debug_vorbin: print(f'Eventually {len(np.unique(np.ma.compressed(binID_map)))} out of {len(np.unique(binIDs))} vorbins are above SNR={snr_thresh} for all {len(relevant_lines)} lines of interest')

    # -------debugging checks: final binned maps-----------
    if args is not None and args.debug_vorbin:
        axes_bottom = axes_all[2]

        ax = axes_mid[3]
        p = ax.scatter(y_coords_array, x_coords_array, c=binIDs, cmap=random_cmap, marker='s', s=10, vmin=-1, vmax=np.max(binIDs))
        ax.set_xlim(0, np.shape(input_map)[0])
        ax.set_ylim(0, np.shape(input_map)[1])
        ax.set_aspect('equal')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad='2%')
        cbar = plt.colorbar(p, cax=cax, orientation='vertical')
        ax.text(ax.get_xlim()[0] + np.diff(ax.get_xlim())[0] * 0.1, ax.get_ylim()[1] - np.diff(ax.get_ylim())[0] * 0.05, 'resultant bin IDs', c='k', ha='left', va='top', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))

        plot_2D_map(binID_map, axes_bottom[3], args, takelog=False, label='interpolated bin IDs', cmap=random_cmap, hide_yaxis=True, vmin=-1, vmax=np.max(binIDs))
        
        map, map_err = bin_2D(input_map, binID_map, map_err=input_map_err, debug_vorbin=args.debug_vorbin, ax=axes_top[3])
        axes_top[3].text(axes_top[3].get_xlim()[0] + np.diff(axes_top[3].get_xlim())[0] * 0.05, axes_top[3].get_ylim()[1] - np.diff(axes_top[3].get_ylim())[0] * 0.05, 'SNR of pixels in each bin', c='k', ha='left', va='top', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
        axes_top[3].axvline(0, ls='dotted', c='grey', lw=0.5)
        snr = map / map_err
        snr_min, snr_max = np.nanmin(snr.data), np.nanmax(snr.data)
        combined_cmap = get_combined_cmap([snr_min, 0, snr_thresh, snr_max], ['Greys', 'Reds', cmap], new_name='my_colormap')
        plot_2D_map(map, axes_bottom[0], args, takelog=False, label=f'binned {args.voronoi_line} map', cmap=cmap, hide_yaxis=False)
        plot_2D_map(map_err, axes_bottom[1], args, takelog=False, label=f'binned {args.voronoi_line} err', cmap=cmap, hide_yaxis=True)
        plot_2D_map(snr, axes_bottom[2], args, takelog=False, label='binned SNR', cmap=combined_cmap, hide_yaxis=True, vmin=snr_min, vmax=snr_max)
        print(f'\nDeb677: {len(np.unique(np.ma.compressed(np.ma.masked_where(snr.mask | (snr < snr_thresh), snr))))} out of {len(np.unique(binIDs))} vorbins are above SNR={snr_thresh}')

        figname = fig_dir / f'{args.field}_{args.id:05d}_vorbin_debug{snr_text}{vorbin_text}.png'
        fig.savefig(figname)
        print(f'Saved figure to {figname}')
        plt.show(block=False)
        sys.exit(f'Exiting here because of --debug_vorbin mode; if you want to run the full code as usual then remove the --debug_vorbin option and re-run')

    return binID_map

# --------------------------------------------------------------------------------------------------------------------
def bin_2D_radial(map, bin_IDs, map_err=None, debug_vorbin=False, ax=None):
    '''
    Binning a given 2D array according to a given bin IDs map, to get the mean flux in each bin
    '''
    valid_mask = bin_IDs.data.ravel() >= 0
    map_flat = map.ravel()[valid_mask]
    bin_IDs_flat = bin_IDs.data.ravel()[valid_mask]
    n_bins = len(np.unique(bin_IDs_flat))

    # mean fluxes per bin
    map_sum = np.bincount(bin_IDs_flat, weights=map_flat, minlength=n_bins)
    bin_counts = np.bincount(bin_IDs_flat, minlength=n_bins)
    map_mean = map_sum / bin_counts

    binned_map = np.full_like(map, np.nan, dtype=float)
    valid_pixels = bin_IDs.data >= 0
    binned_map[valid_pixels] = map_mean[bin_IDs.data[valid_pixels]]
    binned_map = np.ma.masked_where(bin_IDs.mask, binned_map)

    # Sum errors in quadrature per bin: sqrt(sum(err^2))
    if map_err is not None:
        map_err_flat = map_err.ravel()[valid_mask]
        map_err_squared_sum = np.bincount(bin_IDs_flat, weights=map_err_flat**2, minlength=n_bins)
        map_err_rms = np.sqrt(map_err_squared_sum) / bin_counts

        binned_map_err = np.full_like(map, np.nan, dtype=float)
        binned_map_err[valid_pixels] = map_err_rms[bin_IDs.data[valid_pixels]]
        binned_map_err = np.ma.masked_where(bin_IDs.mask, binned_map_err)

        # --------computing how many bins with SNR above a certain threshold for the given line--------------
        snr = binned_map / binned_map_err
        print(f'Upon radial binning, out of {len(np.unique(np.ma.compressed(bin_IDs)))} radbins {len(np.unique(np.ma.compressed(np.ma.masked_where(bin_IDs.mask | (snr <= 0), bin_IDs.data))))} have SNR>0 and {len(np.unique(np.ma.compressed(np.ma.masked_where(bin_IDs.mask | (snr < 3), bin_IDs.data))))} have SNR>=3\n')

    if debug_vorbin:
        col_arr = plt.cm.viridis(np.linspace(0, 1, n_bins))
        for index, this_bin_ID in enumerate(np.unique(bin_IDs_flat)):
            this_bin_mask = bin_IDs_flat == this_bin_ID
            all_pixel_fluxes = map_flat[this_bin_mask]
            all_pixel_errors = map_err_flat[this_bin_mask]
            all_pixels_snr = all_pixel_fluxes / all_pixel_errors

            good_fluxes_mask = all_pixels_snr > 0
            good_pixel_fluxes = all_pixel_fluxes[good_fluxes_mask]

            print(f'Deb763: id {int(this_bin_ID)}/{n_bins}: {len(all_pixel_fluxes)} pixels; {len(good_pixel_fluxes)} good pixels') ##
            if ax is not None:
                hist, bin_edges = np.histogram(all_pixels_snr, bins=50, density=True, range=(-5, 5))
                bin_centers = bin_edges[:-1] + np.diff(bin_edges)
                ax.step(bin_centers, hist + index * 1., lw=0.5, where='mid', color=col_arr[index])

    if map_err is None: return binned_map
    else: return binned_map, binned_map_err

# --------------------------------------------------------------------------------------------------------------------
def get_radial_bin_IDs(full_hdu, snr_thresh=3, plot=False, quiet=True, args=None):
    '''
    Compute the radial bin IDs a given 2D map and corresponding uncertainty
    Returns the 2D map (of same shape as input map) with just the IDs
    '''

    # -------debugging checks: initial maps-----------
    if args is not None and args.debug_vorbin:
        line_map, _, _, _, _ = get_emission_line_map('OII', full_hdu, args, dered=False, for_vorbin=True)
        input_map = np.ma.masked_where(line_map.mask, unp.nominal_values(line_map))
        input_map_err = np.ma.masked_where(line_map.mask, unp.std_devs(line_map))
        input_snr = input_map / input_map_err

        fig, axes_all = plt.subplots(2, 4, figsize=(11, 5))
        fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.07, wspace=0.3, hspace=0.1)
        cmap = 'viridis'
        axes_top = axes_all[0]
        fig.text(0.05, 0.98, f'{args.field}: ID {args.id}: Voronoi binning diagnostics', fontsize=args.fontsize, c='k', ha='left', va='top')

        snr_min, snr_max = np.nanmin(input_snr.data), np.nanmax(input_snr.data)
        combined_cmap = get_combined_cmap([snr_min, 0, snr_thresh, snr_max], ['Greys', 'Reds', cmap], new_name='my_colormap')
        plot_2D_map(input_map, axes_top[0], args, takelog=False, label=f'input {args.voronoi_line} map', cmap=cmap, hide_yaxis=False, hide_xaxis=True)
        plot_2D_map(input_map_err, axes_top[1], args, takelog=False, label=f'input {args.voronoi_line} err', cmap=cmap, hide_yaxis=True, hide_xaxis=True)
        plot_2D_map(input_snr, axes_top[2], args, takelog=False, label='input SNR', cmap=combined_cmap, hide_yaxis=True, vmin=snr_min, vmax=snr_max, hide_xaxis=True)

    # -------getting the radial distance map-----------
    distance_map = get_distance_map(np.shape(args.segmentation_map), args, for_distmap=True) # kpc

    # -------radially binning-----------
    if args.re_limit is None:
        radial_bin_edges = custom_spaced_grid(0,  5, args.nbins + 1, power=1.5)
    else:
        radial_bin_edges = custom_spaced_grid(0,  args.re_limit, args.nbins + 1, power=1.5) # user power = 1 for linspace
        distance_map /= args.re_kpc # converting from kpc to Re units

    print(f'Using {args.nbins} bins with bin edges {radial_bin_edges} {"kpc" if args.re_limit is None else "Re"} for radial binning..')
    binID_map = np.digitize(distance_map, radial_bin_edges)
    if args.only_seg: seg_mask = args.segmentation_map != args.id
    else: seg_mask = False
    binID_map = np.ma.masked_where(seg_mask | binID_map == args.nbins + 1, binID_map)
    if args.re_limit is None: radial_bin_edges *= cosmo.arcsec_per_kpc_proper(args.z).value # converting from kpc to arcsec units

    # -------debugging checks: final binned maps-----------
    if args is not None and args.debug_vorbin:
        axes_bottom = axes_all[1]

        plot_2D_map(binID_map, axes_bottom[3], args, takelog=False, label='radial bin IDs', cmap=random_cmap, hide_yaxis=True, vmin=-1, vmax=np.max(binID_map), plot_circles_at=radial_bin_edges[1:])
        
        map, map_err = bin_2D_radial(input_map, binID_map, map_err=input_map_err, debug_vorbin=args.debug_vorbin, ax=axes_top[3])
        axes_top[3].text(axes_top[3].get_xlim()[0] + np.diff(axes_top[3].get_xlim())[0] * 0.05, axes_top[3].get_ylim()[1] - np.diff(axes_top[3].get_ylim())[0] * 0.05, 'SNR of pixels in each bin', c='k', ha='left', va='top', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))
        axes_top[3].axvline(0, ls='dotted', c='grey', lw=0.5)
        snr = map / map_err
        snr_min, snr_max = np.nanmin(snr.data), np.nanmax(snr.data)
        combined_cmap = get_combined_cmap([snr_min, 0, snr_thresh, snr_max], ['Greys', 'Reds', cmap], new_name='my_colormap')
        plot_2D_map(map, axes_bottom[0], args, takelog=False, label=f'binned {args.voronoi_line} map', cmap=cmap, hide_yaxis=False, plot_circles_at=radial_bin_edges[1:])
        plot_2D_map(map_err, axes_bottom[1], args, takelog=False, label=f'binned {args.voronoi_line} err', cmap=cmap, hide_yaxis=True, plot_circles_at=radial_bin_edges[1:])
        plot_2D_map(snr, axes_bottom[2], args, takelog=False, label='binned SNR', cmap=combined_cmap, hide_yaxis=True, vmin=snr_min, vmax=snr_max, plot_circles_at=radial_bin_edges[1:])

        figname = fig_dir / f'{args.field}_{args.id:05d}_vorbin_debug{snr_text}{vorbin_text}.png'
        fig.savefig(figname)
        print(f'Saved figure to {figname}')
        plt.show(block=False)
        sys.exit(f'Exiting here because of --debug_vorbin mode; if you want to run the full code as usual then remove the --debug_vorbin option and re-run')

    return binID_map

# --------------------------------------------------------------------------------------------------------------------
def custom_spaced_grid(start, stop, num, power=1.5):
    linear = np.linspace(0, 1, num) # Generate linearly spaced numbers between 0 and 1
    nonlinear = linear ** power # Apply a non-linear transformation
    
    return start + (stop - start) * nonlinear

# --------------------------------------------------------------------------------------------------------------------
def cut_by_segment(map, args):
    '''
    Mask a given 2D map according to the segmentation map from the HDU
    Returns the masked map
    '''
    cut_map = np.ma.masked_where(args.segmentation_map != args.id, map)

    return cut_map

# --------------------------------------------------------------------------------------------------------------------
def get_emission_line_int(line, full_hdu, args, dered=True, silent=False):
    '''
    Retrieve the integrated flux for a given emission line from the HDU
    Returns the 2D line image
    '''
    try:
        line_hdu = full_hdu['LINE', line]
    except KeyError:
        if line == 'OIII': line_hdu = full_hdu['LINE', 'OIII-5007']
    line_wave = line_hdu.header['RESTWAVE'] # in Angstrom

    line_index = np.where(args.available_lines == line)[0][0]
    try:
        line_int = full_hdu[0].header[f'FLUX{line_index + 1:03d}'] # ergs/s/cm^2
        line_int_err = full_hdu[0].header[f'ERR{line_index + 1:03d}'] # ergs/s/cm^2
    except KeyError:
        line_int, line_int_err = np.nan, np.nan

    # ----------deblending flux--------------------
    factor = 1.
    if not args.do_not_correct_flux:
        if line == 'OIII': # special treatment for OIII 5007 line, in order to account for and remove the OIII 4959 component
            ratio_5007_to_4959 = 2.98 # from grizli source code
            factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
            if not silent: print(f'Correcting OIII for 4959 component, by factor of {factor:.3f}')
        elif line == 'Ha': # special treatment for Ha line, in order to account for and remove the NII component
            factor = 0.823 # from James et al. 2023?
            if not silent: print(f'Correcting Ha for NII component, by factor of {factor:.3f}')

    line_int = line_int * factor
    line_int_err = line_int_err * factor
    line_int = ufloat(line_int, line_int_err)
    if dered: line_int = get_dereddened_flux(line_int, line_wave, args.EB_V)
    if not silent: print(f'Integrated {line} flux for object {args.id} is {line_int} ergs/s/cm^2.')

    return line_int

# --------------------------------------------------------------------------------------------------------------------
def get_emission_line_map(line, full_hdu, args, dered=True, for_vorbin=False, silent=False):
    '''
    Retrieve the emission map for a given line from the HDU
    Returns the 2D line image
    '''
    # -----------getting the integrated flux value from grizli-----------------
    line_int = get_emission_line_int(line, full_hdu, args, dered=dered and not for_vorbin, silent=silent)

    # -----------getting the integrated EW value-----------------
    line_index = np.where(args.available_lines == line)[0][0]
    try:
        line_index_in_cov = int([item for item in list(full_hdu[2].header.keys()) if full_hdu[0].header[f'FLUX{line_index + 1:03d}'] == full_hdu[2].header[item]][0][5:])
        line_ew = full_hdu[2].header[f'EW50_{line_index_in_cov:03d}'] # rest-frame EW
        line_ew_err = full_hdu[2].header[f'EWHW_{line_index_in_cov:03d}'] # rest-frame EW uncertainty
        line_ew = ufloat(line_ew, line_ew_err)
    except KeyError:
        line_ew = ufloat(np.nan, np.nan)

    # ---------getting the spatially resolved line flux map----------------
    try:
        line_hdu = full_hdu['LINE', line]
        line_map_err = 1e-17 / (full_hdu['LINEWHT', line].data ** 0.5)  # ERR = 1/sqrt(LINEWHT) = flux uncertainty; in units of ergs/s/cm^2
    except KeyError:
        if line == 'OIII':
            line_hdu = full_hdu['LINE', 'OIII-5007']
            line_map_err = 1e-17 / (full_hdu['LINEWHT', 'OIII-5007'].data ** 0.5)  # ERR = 1/sqrt(LINEWHT) = flux uncertainty; in units of ergs/s/cm^2

    line_map = line_hdu.data * 1e-17 # in units of ergs/s/cm^2
    line_wave = line_hdu.header['RESTWAVE'] # in Angstrom
    factor = 1.0

    # ----------deblending flux--------------------
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

    # -------------pixel offset to true center----------------
    line_map = np.roll(line_map, args.ndelta_xpix, axis=0)
    line_map_err = np.roll(line_map_err, args.ndelta_xpix, axis=0)
    line_map = np.roll(line_map, args.ndelta_ypix, axis=1)
    line_map_err = np.roll(line_map_err, args.ndelta_ypix, axis=1)

    # ----------getting a smaller cutout around the object center-----------
    line_map = trim_image(line_map, args)
    line_map_err = trim_image(line_map_err, args)

    # -------getting segmentation map cutout---------
    if args.only_seg:
        line_map = cut_by_segment(line_map, args)
        line_map_err = cut_by_segment(line_map_err, args)

    seg_mask = line_map.mask if np.ma.isMaskedArray(line_map) else False

    # -----------getting the dereddened flux value-----------------
    if dered and not for_vorbin:
        line_map_quant = get_dereddened_flux(unp.uarray(line_map, line_map_err), line_wave, args.EB_V)
        line_map = unp.nominal_values(line_map_quant)
        line_map_err = unp.std_devs(line_map_quant)

    # -----------getting the integrated flux value by summing the 2D map-----------------
    line_sum = np.sum(np.ma.masked_where(seg_mask | ~np.isfinite(line_map.data), unp.uarray(line_map.data, line_map_err.data))) # ergs/s/cm^2
    if not silent: print(f'Summed up {line} flux for object {args.id} is {line_sum:.1e} ergs/s/cm^2.')

    if args.only_integrated:
        line_map = None
    else:
        # -----------voronoi binning the flux and err maps-----------------
        if args.radbin and not for_vorbin:
            line_map, line_map_err = bin_2D_radial(line_map, args.voronoi_bin_IDs, map_err=line_map_err)
            
        elif args.vorbin and not for_vorbin:

            # -----------discarding low-snr pixels BEFORE vorbin, if any-----------------
            if args.snr_cut is not None:
                snr_map = line_map / line_map_err
                snr_mask = (~np.isfinite(snr_map)) | (snr_map < args.snr_cut)
                line_map = np.ma.masked_where(snr_mask, line_map.data)
                line_map_err = np.ma.masked_where(snr_mask, line_map_err.data)

            # ------getting vorornoi bin IDs------------------
            if args.voronoi_line is None: # No reference emission line specified, so Voronoi IDs need to be computed now
                bin_IDs = get_voronoi_bin_IDs(unp.uarray(line_map, line_map_err), args.voronoi_snr, plot=args.debug_vorbin, quiet=not args.debug_vorbin, args=args)
            else: # Reference emission line specified for Voronoi binning, so bin IDs have been pre-computed
                bin_IDs = args.voronoi_bin_IDs

            # ------smoothing the map before voronoi binning--------
            smoothing_kernel = Box2DKernel(args.kernel_size, mode=args.kernel_mode)
            line_map = np.ma.masked_where(line_map.mask, convolve(line_map.data, smoothing_kernel))
            line_map_err = np.ma.masked_where(line_map_err.mask, convolve(line_map_err.data, smoothing_kernel))

            # ------debugging plots-------------------
            if args.debug_vorbin:
                fig, axes = plt.subplots(1, 7, figsize=(14, 3))
                fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.1, wspace=0.5)
                cmap = 'viridis'
                fig.text(0.05, 0.8, f'{args.field}: ID {args.id}: line {line}', fontsize=args.fontsize, c='k', ha='left', va='top')
                plot_2D_map(line_map, axes[0], args, takelog=False, label='map', cmap=cmap, hide_yaxis=True)
                plot_2D_map(line_map_err, axes[1], args, takelog=False, label='map err', cmap=cmap, hide_yaxis=True)
                snr = line_map / line_map_err
                combined_cmap = make_combined_cmap(cmap, 'Grays', np.min(snr), np.max(snr), args.snr_cut)
                plot_2D_map(snr, axes[2], args, takelog=False, label='snr', cmap=combined_cmap, hide_yaxis=True)
                plot_2D_map(bin_IDs, axes[3], args, takelog=False, label='bin IDs', cmap=random_cmap, hide_yaxis=False, hide_cbar=True)

            line_map, line_map_err = bin_2D(line_map, bin_IDs, map_err=line_map_err)
        
        # ---------converting from flux to surface brightness units------------
        pixscale_kpc = args.pix_size_arcsec/ cosmo.arcsec_per_kpc_proper(args.z).value # kpc
        line_map /= (pixscale_kpc ** 2) # this is now in ergs/s/cm^2/kpc^2 (i.e. Surface Brightness)
        line_map_err /= (pixscale_kpc ** 2) # this is now in ergs/s/cm^2/kpc^2 (i.e. Surface Brightness)

        # -----------discarding low-snr pixels AFTER vorbin, if any-----------------
        snr_map = line_map / line_map_err
        if args.snr_cut is not None:
            snr_mask = (~np.isfinite(snr_map)) | (snr_map < args.snr_cut)
        else:
            snr_mask = ~np.isfinite(snr_map)

        line_map = np.ma.masked_where(seg_mask | snr_mask, unp.uarray(line_map.data, line_map_err.data))

        if args.vorbin and not for_vorbin and args.debug_vorbin:
            map = np.ma.masked_where(line_map.mask, unp.nominal_values(line_map.data))
            map_err = np.ma.masked_where(line_map.mask, unp.std_devs(line_map.data))

            plot_2D_map(map, axes[4], args, takelog=False, label='binned map', cmap=cmap, hide_yaxis=True)
            plot_2D_map(map_err, axes[5], args, takelog=False, label='binned map err', cmap=cmap, hide_yaxis=True)
            snr = map / map_err
            combined_cmap = make_combined_cmap(cmap, 'Grays', np.min(snr), np.max(snr), args.snr_cut)
            plot_2D_map(snr, axes[6], args, takelog=False, label='binned snr', cmap=combined_cmap, hide_yaxis=True)
            plt.show(block=False)

    return line_map, line_wave, line_int, line_sum, line_ew

# --------------------------------------------------------------------------------------------------------------------
def plot_emission_line_map(line, full_hdu, ax, args, cmap='cividis', EB_V=None, vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=False, snr_ax=None, radprof_ax=None):
    '''
    Plots the emission map for a given line in the given axis
    Returns the axes handle
    '''

    line_map, line_wave, line_int, _, line_ew = get_emission_line_map(line, full_hdu, args, dered=True)
    ax, _ = plot_2D_map(line_map, ax, args, label=r'%s$_{\rm int}$ = %.1e' % (line, line_int.n), cmap=cmap, vmin=vmin, vmax=vmax, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar, snr_ax=snr_ax, radprof_ax=radprof_ax)
    if args.arcsec_limit >= 1 and not args.no_text_on_plot: ax.text(ax.get_xlim()[0] * 0.88, ax.get_ylim()[0] * 0.88, f'EW = {line_ew:.1f}', c='k', fontsize=args.fontsize, ha='left', va='bottom', bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_line_ratio_map(line_num, line_den, full_hdu, ax, args, cmap='cividis', vmin=None, vmax=None, hide_xaxis=False, hide_yaxis=False, hide_cbar=False, snr_ax=None, radprof_ax=None):
    '''
    Plots the emission line ratio map for a given pair of lines in the given axis
    Returns the axes handle
    '''

    num_map, _, num_int, _, _ = get_emission_line_map('Ha' if line_num == 'NII' else line_num, full_hdu, args, dered=True)
    den_map, _, den_int, _, _ = get_emission_line_map('Ha' if line_den == 'NII' else line_den, full_hdu, args, dered=True)

    # ----------deblending flux--------------------
    if not args.do_not_correct_flux:
        # special treatment for H-alpha line, in order to account for NII 6584 component
        # factor = factor by which Ha+NII complex was originally corrected to get Ha; need to divide Ha map by this factor to get back the Ha+NII complex
        factor = 0.823  # from grizli source code
    else:
        factor = 1.

    if line_num == 'NII': # converting original Ha map to NII map by simple scaling
        num_map = np.ma.masked_where(num_map.mask, num_map.data * (1 - 0.823) / factor)
        num_int = num_int * (1 - 0.823) / factor

    if line_den == 'NII': # converting original Ha map to NII map by simple scaling
        den_map = np.ma.masked_where(den_map.mask, den_map.data * (1 - 0.823) / factor)
        den_int = den_int * (1 - 0.823) / factor

    ratio_map = take_safe_log_ratio(num_map, den_map)

    try: ratio_int = unp.log10(num_int / den_int)
    except: ratio_int = ufloat(np.nan, np.nan)

    ax, _ = plot_2D_map(ratio_map, ax, args, takelog=False, label=r'%s/%s$_{\rm int}$ = %.1f$\pm$%.1f' % (line_num, line_den, unp.nominal_values(ratio_int), unp.std_devs(ratio_int)), cmap=cmap, vmin=vmin, vmax=vmax, hide_xaxis=hide_xaxis, hide_yaxis=hide_yaxis, hide_cbar=hide_cbar, snr_ax=snr_ax, radprof_ax=radprof_ax)

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
def get_EB_V(full_hdu, args, verbose=False, silent=False, for_vorbin=False):
    '''
    Computes and returns the spatially resolved as well as integrated dust extinction map from a given HDU
    Based on Eqn 4 of Dominguez+2013 (https://iopscience.iop.org/article/10.1088/0004-637X/763/2/145/pdf)
    '''

    Ha_map, Ha_wave, Ha_int, Ha_sum, _ = get_emission_line_map('Ha', full_hdu, args, dered=False, silent=silent, for_vorbin=for_vorbin) # do not need to deredden the lines when we are fetching the flux in order to compute reddening
    Hb_map, Hb_wave, Hb_int, Hb_sum, _ = get_emission_line_map('Hb', full_hdu, args, dered=False, silent=silent, for_vorbin=for_vorbin)

    EB_V_map = compute_EB_V(Ha_map, Hb_map)
    EB_V_int = compute_EB_V(Ha_int, Hb_int, verbose=verbose)
    EB_V_sum = compute_EB_V(Ha_sum, Hb_sum, verbose=verbose)

    return EB_V_map, EB_V_int, EB_V_sum

# --------------------------------------------------------------------------------------------------------------------
def plot_EB_V_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the dust extinction map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D E(B-V) map just produced
    '''
    lim, label = [0, 1], 'E(B-V)'

    EB_V_map, EB_V_int, EB_V_sum = get_EB_V(full_hdu, args)
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

    Ha_lum = Ha_flux * 4 * np.pi * (distance.to('cm').value) ** 2 # converting to ergs/s (luminosity)
    #sfr = Ha_lum * 7.9e-42 # luminosity in ergs/s; SFR in Msun/yr; from Kennicutt+1998
    #sfr = Ha_lum * 10**(-41.67) # luminosity in ergs/s; SFR in Msun/yr; from Reddy+2022
    sfr = Ha_lum * ufloat(7.5, 1.3) * 10 ** (-42) # luminosity in ergs/s; SFR in Msun/yr; from Shivaei+2015

    if hasattr(Ha_flux, "__len__"): # if it is an array
        sfr = np.ma.masked_where(net_mask, sfr)

    return sfr

# --------------------------------------------------------------------------------------------------------------------
def get_SFR(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated SFR from a given HDU
    '''
    Ha_map, Ha_wave, Ha_int, Ha_sum, _ = get_emission_line_map('Ha', full_hdu, args)

    SFR_map = compute_SFR(Ha_map, args.distance)
    SFR_int = compute_SFR(Ha_int, args.distance)
    SFR_sum = compute_SFR(Ha_sum, args.distance)

    return SFR_map, SFR_int, SFR_sum

# --------------------------------------------------------------------------------------------------------------------
def plot_SFR_map(full_hdu, ax, args, radprof_ax=None, snr_ax=None):
    '''
    Plots the SFR map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D SFR density map just produced
    '''
    lim, label = [-3, -1], 'SFR'
    SFR_map, SFR_int, SFR_sum = get_SFR(full_hdu, args)
    ax, SFR_radfit = plot_2D_map(SFR_map, ax, args, label=r'log %s$_{\rm int}$ = %.1f' % (label, SFR_int.n), cmap='viridis', radprof_ax=radprof_ax, snr_ax=snr_ax, vmin=lim[0], vmax=lim[1], hide_yaxis=True, hide_xaxis=True)

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
        log_O3 = take_safe_log_ratio(OIII4363_flux, OIII5007_flux)
        Te = []
        for this_log_O3 in log_O3.data.flatten():
            try:
                this_logTe = np.poly1d([0., 9.18962, 3.30355])(this_log_O3) / np.poly1d([-0.00935, -0.15034, 2.09136, 1.000])(this_log_O3)
                this_Te = 10 ** this_logTe
                Te.append(this_Te)
            except:
                Te.append(ufloat(np.nan, np.nan))
        Te = np.ma.masked_where(log_O3.mask, np.reshape(Te, np.shape(log_O3)))
        Te = np.ma.masked_where(~np.isfinite(unp.nominal_values(Te.data)), Te)

    else: # if it is scalar
        try:
            log_O3 = unp.log10(OIII4363_flux / OIII5007_flux)
            logTe = np.poly1d([0., 9.18962, 3.30355])(log_O3) / np.poly1d([-0.00935, -0.15034, 2.09136, 1.000])(log_O3)
            Te = 10 ** logTe
        except:
            Te = ufloat(np.nan, np.nan)

    return Te

# --------------------------------------------------------------------------------------------------------------------
def get_Te(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated Te from a given HDU
    '''
    OIII4363_map, OIII4363_wave, OIII4363_int, OIII4363_sum, _ = get_emission_line_map('OIII-4363', full_hdu, args)
    OIII5007_map, OIII5007_wave, OIII5007_int, OIII5007_sum, _ = get_emission_line_map('OIII', full_hdu, args)

    Te_map = compute_Te(OIII4363_map, OIII5007_map)
    Te_int = compute_Te(OIII4363_int, OIII5007_int)
    Te_sum = compute_Te(OIII4363_sum, OIII5007_sum)

    return Te_map, Te_int, Te_sum

# --------------------------------------------------------------------------------------------------------------------
def plot_Te_map(full_hdu, ax, args, radprof_ax=None):
    '''
    Plots the T_e map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D T_e map just produced
    '''
    lim, label = [1, 7], r'T$_e$'
    Te_map, Te_int, Te_sum = get_Te(full_hdu, args)
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

    if hasattr(Hbeta_flux, "__len__"): # if it is an array
        O2H2 = take_safe_log_ratio(OII3727_flux, Hbeta_flux, skip_log=True)
        O3H2 = take_safe_log_ratio(OIII5007_flux, Hbeta_flux, skip_log=True)

        log_OH = []
        for (this_O2H2, this_O3H2, this_Te) in zip(O2H2.data.flatten(), O3H2.data.flatten(), Te.flatten()):
            try:
                t =this_Te * 1e-4
                x = 1e-4 * ne * unp.sqrt(t)

                this_log_O2H2 = poly(this_O2H2, t, x, 5.961, 1.676, 0.4, 0.034, 1.35) - 12  # coefficients from eqn 3 I06 # in case 'ratio1' happens to be < 0
                this_log_O3H2 = poly(this_O3H2, t, x, 6.200, 1.251, -0.55, -0.014,0.0) - 12  # coefficients from eqn 5 I06 # in case 'ratio1' happens to be < 0

                this_log_OH = unp.log10(10 ** this_log_O2H2 + 10 ** this_log_O3H2) + 12
                log_OH.append(this_log_OH)
            except:
                log_OH.append(ufloat(np.nan, np.nan))
        log_OH = np.ma.masked_where(O2H2.mask | O3H2.mask, np.reshape(log_OH, np.shape(O2H2)))
        log_OH = np.ma.masked_where(~np.isfinite(unp.nominal_values(log_OH.data)), log_OH)

    else: # if it is scalar
        try:
            t = Te * 1e-4
            x = 1e-4 * ne * unp.sqrt(t)

            O2H2 = OII3727_flux / Hbeta_flux # in case 'Hbeta_flux' happens to be 0
            O3H2 = OIII5007_flux / Hbeta_flux

            log_O2H2 = poly(O2H2, t, x, 5.961, 1.676, 0.4, 0.034, 1.35) - 12  # coefficients from eqn 3 I06 # in case 'ratio1' happens to be < 0
            log_O3H2 = poly(O3H2, t, x, 6.200, 1.251, -0.55, -0.014, 0.0) - 12  # coefficients from eqn 5 I06 # in case 'ratio1' happens to be < 0

            log_OH = unp.log10(10 ** log_O2H2 + 10 ** log_O3H2) + 12
        except:
            log_OH = ufloat(np.nan, np.nan)

    return log_OH

# --------------------------------------------------------------------------------------------------------------------
def get_Z_Te(full_hdu, args):
    '''
    Computes and returns the spatially resolved as well as intregrated Te metallicity from a given HDU
    '''
    OII3727_map, line_wave, OII3727_int, OII3727_sum, _ = get_emission_line_map('OII', full_hdu, args)
    OIII5007_map, line_wave, OIII5007_int, OIII5007_sum, _ = get_emission_line_map('OIII', full_hdu, args)
    Hbeta_map, line_wave, Hbeta_int, Hbeta_sum, _ = get_emission_line_map('Hb', full_hdu, args)

    Te_map, Te_int, Te_sum = get_Te(full_hdu, args)

    logOH_map = compute_Z_Te(OII3727_map, OIII5007_map, Hbeta_map, Te_map)
    logOH_int = compute_Z_Te(OII3727_int, OIII5007_int, Hbeta_int, Te_int)
    logOH_sum = compute_Z_Te(OII3727_sum, OIII5007_sum, Hbeta_sum, Te_sum)

    return logOH_map, logOH_int, logOH_sum

# --------------------------------------------------------------------------------------------------------------------
def compute_Z_P25_SFR(O3S2, N2S2):
    '''
    Calculates and returns the SFR metallicity given observed O3S2 and N2S2 line ratios
    Conversion factor is from Eq 1 of Peluso+2025
    '''
    return 8.78 + 0.97 * N2S2 - 0.11 * O3S2 - 0.39 * N2S2 ** 2 + 0.09 * N2S2 * O3S2 + 0.02 * O3S2 ** 2

# --------------------------------------------------------------------------------------------------------------------
def compute_Z_P25_AGN(O3S2, N2S2):
    '''
    Calculates and returns the AGN metallicity given observed O3S2 and N2S2 line ratios
    Conversion factor is from Eq 2 of Peluso+2025
    '''
    return 8.85 + 1.10 * N2S2 - 0.04 * O3S2

# --------------------------------------------------------------------------------------------------------------------
def compute_Z_P25_Comp(O3S2, N2S2):
    '''
    Calculates and returns the Composite region metallicity given observed O3S2 and N2S2 line ratios
    Conversion factor is from Eq 3 of Peluso+2025
    '''
    return 8.83 + 1.07 * N2S2 + 0.10 * O3S2

# --------------------------------------------------------------------------------------------------------------------
def compute_Z_P25(OIII5007_flux, NII6584_flux, SII6717_flux, AGN_map):
    '''
    Calculates and returns the SF/AGN metallicity given observed line fluxes
    Conversion factor is from Peluso+2025
    '''
    # -----handling masks separately because uncertainty package cannot handle masks---------
    if np.ma.isMaskedArray(OIII5007_flux):
        net_mask = SII6717_flux.mask | OIII5007_flux.mask | NII6584_flux.mask
        OIII5007_flux = OIII5007_flux.data
        SII6717_flux = SII6717_flux.data
        NII6584_flux = NII6584_flux.data
    else:
        net_mask = False

    if hasattr(SII6717_flux, "__len__"): # if it is an array
        # --------computing the ratio and appropriate errors------------
        O3S2 = take_safe_log_ratio(OIII5007_flux, SII6717_flux)
        N2S2 = take_safe_log_ratio(NII6584_flux, SII6717_flux)

        # --------computing the polynomial and appropriate errors------------
        log_OH = []
        if args.debug_Zdiag:
            fig, ax = plt.subplots(2, 1, figsize=(6, 8), sharex=True)
            ax[0].set_ylabel('O3S2')
            ax[1].set_ylabel('N2S2')
            ax[1].set_xlabel('log(O/H)+12')
            ax[0].set_xlim(7.5, 9.5)

        O3S2_flat = O3S2.data.flatten()
        N2S2_flat = N2S2.data.flatten()
        AGN_map_flat = AGN_map.flatten()
        for index in range(len(O3S2_flat)):
            try:
                if AGN_map_flat[index] > 0: this_log_OH = compute_Z_P25_AGN(O3S2_flat[index], N2S2_flat[index])
                else: this_log_OH = compute_Z_P25_SFR(O3S2_flat[index], N2S2_flat[index])
                log_OH.append(this_log_OH)
                if args.debug_Zdiag:
                    ax[0].scatter(unp.nominal_values(this_log_OH), unp.nominal_values(O3S2_flat[index]), lw=0, s=50, c='r' if AGN_map_flat[index] > 0 else 'b')
                    ax[0].errorbar(unp.nominal_values(this_log_OH), unp.nominal_values(O3S2_flat[index]), xerr=unp.std_devs(this_log_OH), yerr=unp.std_devs(O3S2_flat[index]), fmt='none', lw=0.5, c='r' if AGN_map_flat[index] > 0 else 'b')
                    ax[1].scatter(unp.nominal_values(this_log_OH), unp.nominal_values(N2S2_flat[index]), lw=0, s=50, c='r' if AGN_map_flat[index] > 0 else 'b')
                    ax[1].errorbar(unp.nominal_values(this_log_OH), unp.nominal_values(N2S2_flat[index]), xerr=unp.std_devs(this_log_OH), yerr=unp.std_devs(N2S2_flat[index]), fmt='none', lw=0.5, c='r' if AGN_map_flat[index] > 0 else 'b')
            except:
                log_OH.append(ufloat(np.nan, np.nan))
        log_OH = np.ma.masked_where(O3S2.mask | N2S2.mask, np.reshape(log_OH, np.shape(O3S2)))

    else: # if it is scalar
        try:
            O3S2 = unp.log10(OIII5007_flux / SII6717_flux)
            N2S2 = unp.log10(NII6584_flux / SII6717_flux)
            if AGN_map > 0: log_OH = compute_Z_P25_AGN(O3S2, N2S2)
            else: log_OH = compute_Z_P25_SFR(O3S2, N2S2)
        except:
            log_OH = ufloat(np.nan, np.nan)

    return log_OH

# --------------------------------------------------------------------------------------------------------------------
def get_Z_P25(full_hdu, args, branch='low'):
    '''
    Computes and returns the spatially resolved as well as intregrated Peluso+2025 metallicity from a given HDU
    '''
    SII6717_map, line_wave, SII6717_int, SII6717_sum, _ = get_emission_line_map('SII', full_hdu, args)
    OIII5007_map, line_wave, OIII5007_int, OIII5007_sum, _ = get_emission_line_map('OIII', full_hdu, args)
    Halpha_map, line_wave, Halpha_int, Halpha_sum, _ = get_emission_line_map('Ha', full_hdu, args)

    if not args.do_not_correct_flux:
        # special treatment for H-alpha line, in order to account for NII 6584 component
        factor = 0.823  # from grizli source code
    else:
        factor = 1.

    NII6584_map = np.ma.masked_where(Halpha_map.mask, Halpha_map.data * (1 - 0.823) / factor)
    NII6584_int = Halpha_int * (1 - 0.823) / factor
    NII6584_sum = Halpha_sum * (1 - 0.823) / factor

    logOH_map = compute_Z_P25(OIII5007_map, NII6584_map, SII6717_map, args.distance_from_AGN_line_map)
    logOH_int = compute_Z_P25(OIII5007_int, NII6584_int, SII6717_int, args.distance_from_AGN_line_int)
    logOH_sum = compute_Z_P25(OIII5007_sum, NII6584_sum, SII6717_sum, args.distance_from_AGN_line_int)

    return logOH_map, logOH_int, logOH_sum

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
def compute_Z_KD02_R23(OII3727_flux, OIII5007_flux, Hbeta_flux, branch='low'):
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
        R23 = take_safe_log_ratio(OIII5007_flux + OII3727_flux, Hbeta_flux)
        y = take_safe_log_ratio(OIII5007_flux, OII3727_flux)

        logOH_Z94 = np.poly1d([-0.333, -0.207, -0.202, -0.33, 9.625])(R23)  # Eq 8 of KD02
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
                    else:
                        logOH_R23[i][j] = ufloat(np.nan, np.nan)
                        continue
                else:
                    logOH_R23[i][j] = logOH_avg[i][j]

        logOH_R23 = np.ma.masked_where(R23.mask | net_mask, logOH_R23)

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
def get_Z_KD02_R23(full_hdu, args, branch='low'):
    '''
    Computes and returns the spatially resolved as well as intregrated R23 metallicity from a given HDU
    '''
    OII3727_map, line_wave, OII3727_int, OII3727_sum, _ = get_emission_line_map('OII', full_hdu, args)
    OIII5007_map, line_wave, OIII5007_int, OIII5007_sum, _ = get_emission_line_map('OIII', full_hdu, args)
    Hbeta_map, line_wave, Hbeta_int, Hbeta_sum, _ = get_emission_line_map('Hb', full_hdu, args)

    if not args.do_not_correct_flux:
        # special treatment for OIII 5007 line, in order to account for and ADD the OIII 4959 component back
        ratio_5007_to_4959 = 2.98  # from grizli source code
        factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
        print(f'Un-correcting OIII to include the 4959 component, for computing R23 metallicity, by factor of {factor:.3f}')
        OIII5007_map = np.ma.masked_where(OIII5007_map.mask, OIII5007_map.data / factor)
        OIII5007_int = OIII5007_int / factor
        OIII5007_sum = OIII5007_sum  / factor

    logOH_map = compute_Z_KD02_R23(OII3727_map, OIII5007_map, Hbeta_map, branch=branch)
    logOH_int = compute_Z_KD02_R23(OII3727_int, OIII5007_int, Hbeta_int, branch=branch)
    logOH_sum = compute_Z_KD02_R23(OII3727_sum, OIII5007_sum, Hbeta_sum, branch=branch)

    return logOH_map, logOH_int, logOH_sum

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
def compute_Z_C19(ratio, coeff, ax=None, branch='high'):
    '''
    Calculates and returns the metallicity given observed line fluxes ratio and coefficient, according to Curti+2019
    '''
    # -----handling turnover situations, where measured ratio is beyond model peak ratio---------
    metallicity_offset = 0 if args.use_C25 else 8.69
    reasonable_Z_limit = [7.0, 8.4] if args.use_C25 else [7.6, 8.9] # Z limits within which each calibration is valid
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
    print(f'\nDeb1691: ratio at turnover = {ratio_turnover}') ##
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
    logOH_map = compute_Z_C19(ratio_map, coeff, ax=ax, branch=args.Zbranch) if not args.only_integrated else None
    logOH_int = compute_Z_C19(ratio_int, coeff, branch=args.Zbranch)
    logOH_int = ufloat(unp.nominal_values(np.atleast_1d(logOH_int))[0], unp.std_devs(np.atleast_1d(logOH_int))[0])
    logOH_sum = compute_Z_C19(ratio_sum, coeff, branch=args.Zbranch)
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
    output_fitsname = args.output_dir / 'catalogs' / f'{args.field}_{args.id:05d}_logOH_map_upto_{extent_text}{snr_text}{only_seg_text}{vorbin_text}_Zdiag_{args.Zdiag}{Zbranch_text}{AGN_diag_text}{NB_text}{exclude_text}{C25_text}{dered_text}.fits'
    
    # ------checking if the outputfile already exists--------------
    if os.path.exists(output_fitsname) and not args.clobber:
        print(f'\nReading metallicity map from existing fits file {output_fitsname}')
        hdul = fits.open(output_fitsname)
        hdu = hdul['log_OH']
        logOH_map = np.ma.masked_where(np.isnan(hdu.data), unp.uarray(hdu.data, hdul['log_OH_u'].data))
        if hdu.header['log_oh_int'] is None: logOH_int = ufloat(np.nan, np.nan)
        else: logOH_int = ufloat(hdu.header['log_oh_int'], hdu.header['log_oh_int_err'])
        if hdu.header['log_oh_sum'] is None: logOH_sum = ufloat(np.nan, np.nan)
        else: logOH_sum = ufloat(hdu.header['log_oh_sum'], hdu.header['log_oh_sum_err'])
        if args.Zdiag == 'NB': line_label_array = hdu.header['LINES'].split(',')

    else:
        print(f'\nNot found{output_fitsname}: making new one..')
        # --------deriving the metallicity map-------------
        if args.Zdiag == 'NB':
            logOH_map, logOH_int, logOH_sum, line_label_array = get_Z_NB(full_hdu, args)
        elif args.Zdiag == 'KD02_R23' and all([line in args.available_lines for line in ['OIII', 'OII', 'Hb']]):
            logOH_map, logOH_int, logOH_sum = get_Z_KD02_R23(full_hdu, args, branch=args.Zbranch)
        elif args.Zdiag == 'Te' and all([line in args.available_lines for line in ['OIII', 'OIII-4363', 'OII', 'Hb']]):
            logOH_map, logOH_int, logOH_sum = get_Z_Te(full_hdu, args)
        elif args.Zdiag == 'P25' and all([line in args.available_lines for line in ['OIII', 'Ha', 'SII']]):
            logOH_map, logOH_int, logOH_sum = get_Z_P25(full_hdu, args)
        else:
            logOH_map, logOH_int, logOH_sum = get_Z_C19(full_hdu, args)

        # ---------saving the metallicity maps as fits files-------------
        if logOH_map is not None:
            logOH_map_val = np.where(logOH_map.mask, np.nan, unp.nominal_values(logOH_map.data))
            logOH_map_err = np.where(logOH_map.mask, np.nan, unp.std_devs(logOH_map.data))

            distance_map = get_distance_map(np.shape(logOH_map_val), args)
            if args.vorbin:  # getting how many unique IDs present, so that NB does not have to run unnecessary repeats
                bin_IDs_map = np.where(args.voronoi_bin_IDs.mask, np.nan, args.voronoi_bin_IDs.data)
                distance_map = np.where(args.voronoi_bin_IDs.mask, np.nan, distance_map.data)
            else:
                bin_IDs_map = np.arange(np.shape(logOH_map_val)[0] * np.shape(logOH_map_val)[1]).reshape(np.shape(logOH_map_val))
            print(f'\nDeb1935: id {args.id}, Zdiag={args.Zdiag}, logOH_map shape={np.shape(logOH_map)}, logOH_map_val shape={np.shape(logOH_map_val)} distance_map shape={np.shape(distance_map)}, args.re_limit={args.re_limit}, args.arcsec_limit={args.arcsec_limit}, re={args.re_arcsec} arcsec\n') ##

        # --------processing the dataframe in case voronoi binning has been performed and there are duplicate data values------
            df = pd.DataFrame({'distance': distance_map.flatten(), 'log_OH': logOH_map_val.flatten(),'log_OH_u': logOH_map_err.flatten(), 'bin_ID': bin_IDs_map.flatten()})
            
            if 'distance_from_AGN_line_map' in args and args.distance_from_AGN_line_map is not None: df['agn_dist'] = args.distance_from_AGN_line_map.data.flatten()
            df = df.dropna().reset_index(drop=True)
            df['bin_ID'] = df['bin_ID'].astype(int)
            df = df.groupby(['bin_ID'], as_index=False).agg(np.mean)

            hdr1, hdr2 = fits.Header(), fits.Header()
            hdr1['field'] = args.field
            hdr1['object'] = args.id
            hdr1['redshift'] = args.z
            hdr1['re_arcsec'] = args.re_arcsec
            primary_hdu = fits.PrimaryHDU(header=hdr1)

            hdr2['z_diag'] = args.Zdiag
            hdr2['zbranch'] = args.Zbranch
            hdr2['agn_diag'] = args.AGN_diag
            hdr2['vorbin_line'] = args.voronoi_line if args.vorbin else None
            hdr2['vorbin_snr'] = args.voronoi_snr if args.vorbin else None
            hdr2['log_oh_int'] = None if ~np.isfinite(logOH_int.n) else logOH_int.n
            hdr2['log_oh_int_err'] = None if ~np.isfinite(logOH_int.s) else logOH_int.s
            hdr2['log_oh_sum'] = None if ~np.isfinite(logOH_sum.n) else logOH_sum.n
            hdr2['log_oh_sum_err'] = None if ~np.isfinite(logOH_sum.s) else logOH_sum.s
            hdr2['nb_old_grid'] = True if args.use_original_NB_grid and args.Zdiag == 'NB' else False
            if 'NB' in args.Zdiag:
                hdr2['lines'] = ','.join(line_label_array)
                hdr2['nlines'] = len(line_label_array)

            logOH_val_hdu = fits.ImageHDU(data=logOH_map_val, name='log_OH', header=hdr2)
            logOH_err_hdu = fits.ImageHDU(data=logOH_map_err, name='log_OH_u')
            binID_hdu = fits.ImageHDU(data=bin_IDs_map, name='bin_ID')
            df_hdu = fits.BinTableHDU(Table.from_pandas(df), name='tab')
            hdul = fits.HDUList([primary_hdu, logOH_val_hdu, logOH_err_hdu, binID_hdu, df_hdu])
            hdul.writeto(output_fitsname, overwrite=True)

            print(f'Saved metallicity maps in {output_fitsname}')

    if logOH_map is not None and args.mask_agn: logOH_map = np.ma.masked_where((args.distance_from_AGN_line_map > 0) | logOH_map.mask, logOH_map)

    if 'line_label_array' in locals(): return logOH_map, logOH_int, logOH_sum, line_label_array
    else: return logOH_map, logOH_int, logOH_sum

# --------------------------------------------------------------------------------------------------------------------
def plot_Z_map(full_hdu, ax, args, radprof_ax=None, snr_ax=None):
    '''
    Plots the metallicity map in the given axes
    Returns the axes handles and the 2D metallicity map just produced
    '''
    if args.Zdiag == 'NB': logOH_map, logOH_int, logOH_sum, line_label_array = get_Z(full_hdu, args)
    else: logOH_map, logOH_int, logOH_sum = get_Z(full_hdu, args)

    if logOH_map is not None:
        lim = [7, 9]
        ax, logOH_radfit = plot_2D_map(logOH_map, ax, args, takelog=False, label=r'Z (%s)$_{\rm int}$ = %.1f' % (args.Zdiag, logOH_int.n), cmap='viridis', radprof_ax=radprof_ax, snr_ax=snr_ax, hide_yaxis=True, hide_xaxis=True, vmin=lim[0], vmax=lim[1], metallicity_multi_color=args.Zdiag == 'P25')
    else:
        fig = ax.figure
        fig.delaxes(ax)
        if args.plot_radial_profiles: fig.delaxes(radprof_ax)
        if args.plot_snr: fig.delaxes(snr_ax)
        logOH_radfit = np.nan

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

        O32 = take_safe_log_ratio(OIII5007_flux, OII3727_flux)

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
    OII3727_map, line_wave, OII3727_int, _, _ = get_emission_line_map('OII', full_hdu, args)
    OIII5007_map, line_wave, OIII5007_int, _, _ = get_emission_line_map('OIII', full_hdu, args)

    logq_map = compute_q_O32(OII3727_map, OIII5007_map)
    logq_int = compute_q_O32(OII3727_int, OIII5007_int)

    return logq_map, logq_int

# --------------------------------------------------------------------------------------------------------------------
def plot_q_O32_map(full_hdu, ax, args, radprof_ax=None, snr_ax=None):
    '''
    Plots the O32 iojisation parameter map (and the emission line maps that go into it) in the given axes
    Returns the axes handles and the 2D metallicity map just produced
    '''
    lim, label = [7, 8], 'log(q) (O32)'
    logq_map, logq_int = get_q_O32(full_hdu, args)
    ax, logq_radfit = plot_2D_map(logq_map, ax, args, takelog=False, label=r'%s$_{\rm int}$ = %.1f' % (label, logq_int.n), cmap='viridis', radprof_ax=radprof_ax, snr_ax=snr_ax, hide_yaxis=True, hide_xaxis=True, vmin=lim[0], vmax=lim[1])

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
def get_cutout(filename, pos, size, target_header, args, plot_test_axes=None, skip_re_trim=False):
    '''
    Return a cutout from a given filename of a fits image, around a given position within a given extent,
    and then reproject it on to a given target header parameters
    Optionally trims the cutout as per the segmentation map, and/or Voronoi bins it
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
    cutout_data_rebinned_trimmed = trim_image(cutout_data_rebinned, args, skip_re_trim=skip_re_trim)  # 50 x 50 pixels

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
    direct_image_path = args.input_dir / args.field / 'Products'
    drizzled_image_filename = glob.glob(str(direct_image_path / f'{args.field}*{filter.lower()}_drz_sci.fits'))
    if len(drizzled_image_filename) == 0:
        print(f'\nNo direct image found in {direct_image_path}, therefore cannot plot starburst map')
        return None
    else:
        drizzled_image_filename = drizzled_image_filename[0]
        direct_map = get_cutout(drizzled_image_filename, pos, size, target_header, args, plot_test_axes=plot_test_axes)
        direct_map_wht = get_cutout(drizzled_image_filename.replace('sci', 'wht'), pos, size, target_header, args)

        # -------------pixel offset----------------
        print(f'Correcting direct image maps for pixel offset by {args.ndelta_xpix} on x and {args.ndelta_ypix} on y')
        direct_map = np.roll(direct_map, args.ndelta_xpix, axis=0)
        direct_map_wht = np.roll(direct_map_wht, args.ndelta_xpix, axis=0)
        direct_map = np.roll(direct_map, args.ndelta_ypix, axis=1)
        direct_map_wht = np.roll(direct_map_wht, args.ndelta_ypix, axis=1)

        # ---------computing uncertainty-------------
        direct_map_err = 0.1 * np.abs(direct_map) # 1 / np.sqrt(direct_map_wht)
        direct_map = unp.uarray(direct_map, direct_map_err)

        return direct_map

# --------------------------------------------------------------------------------------------------------------------
def plot_direct_images_all_filters(full_hdu, args):
    '''
    Plots direct images from full_hdu vs corresponding cut outs from full field drizzled images
    '''
    filters_arr = args.filters

    # -------getting the full.fits image------------
    filter = 'f115w'
    filter_hdu = full_hdu['DSCI', f'{filter.upper()}-CLEAR']
    direct_map = filter_hdu.data
    direct_map_trimmed = trim_image(direct_map, args)
    target_header = filter_hdu.header
    direct_map_err = 1. / np.sqrt(full_hdu['DWHT', f'{filter.upper()}-CLEAR'].data)
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
def correct_ha_C20(N2_plus_Ha_map, logOH_map):
    '''
    Extract the H-alpha flux from the N2+Halpha compund, following C20 metallicity calibration
    '''
    #log_N2Ha_logOH_poly_coeff = [1, -10] # from approx fit to MAPPINGS models
    log_N2Ha_logOH_poly_coeff = [-0.489, 1.513, -2.554, -5.293, -2.867][::-1] # from Table 2 N2 row of Curti+2019

    if hasattr(N2_plus_Ha_map, "__len__"): # if it is an array
        logOH_arr = logOH_map.data.flatten()
        log_N2Ha_arr = []
        for logOH in logOH_arr:
            try:
                #log_N2Ha = np.poly1d(log_N2Ha_logOH_poly_coeff)(logOH)
                log_N2Ha = np.poly1d(log_N2Ha_logOH_poly_coeff)(logOH - 8.69) # because using C19 N2 calibration
                log_N2Ha_arr.append(log_N2Ha)
            except:
                log_N2Ha_arr.append(np.nan)
        log_N2Ha_map = np.ma.masked_where(logOH_map.mask, np.reshape(log_N2Ha_arr, np.shape(logOH_map)))
        N2Ha_map = np.ma.masked_where(log_N2Ha_map.mask, 10 ** log_N2Ha_map.data)

        Ha_data = N2_plus_Ha_map.data / (1 + N2Ha_map.data)
        Ha_map = np.ma.masked_where(N2_plus_Ha_map.mask | N2Ha_map.mask, Ha_data)
    else: # for single values
        log_N2Ha_int = np.poly1d(log_N2Ha_logOH_poly_coeff)(logOH_map - 8.69)
        N2Ha_int = 10 ** log_N2Ha_int
        Ha_map = N2_plus_Ha_map / (1 + N2Ha_int)

    return Ha_map

# --------------------------------------------------------------------------------------------------------------------
def get_direct_image(full_hdu, filter, args, skip_vorbin=False):
    '''
    Loads the direct image for a given filter for a given object
    Returns the image
    '''
    try:
        hdu = full_hdu['DSCI', filter.upper()]
        image = hdu.data
        image_wht = None
        exptime = 1
    except:
        try:
            hdu = full_hdu['DSCI', f'{filter}-CLEAR']
            image = hdu.data
            image_wht = None
            exptime = full_hdu[0].header[f'T_{filter.upper()}']
        except:
            try:
                hdu = full_hdu['DSCI', f'{filter.upper()}-{filter.upper()}-CLEAR']
                image = hdu.data
                image_wht = None
                exptime = full_hdu[0].header[f'T_{filter.upper()}']
            except:
                full_field_filename = args.input_dir / f'{args.field}' / 'Products' / f'{args.field}_{filter.lower()}-clear_drz_sci.fits'
                print(f'{filter.upper()} not found in full_hdu extension. Therefore trying to get cutout from full field image {full_field_filename}')
            
                exptime = fits.open(full_field_filename)[0].header['EXPTIME']
                pos = SkyCoord(full_hdu[0].header['RA'], full_hdu[0].header['DEC'], unit = 'deg')
                size = 2 * args.arcsec_limit * u.arcsec
                target_header = full_hdu['DSCI', 'F140W'].header
                
                temp1, temp2 = args.only_seg, args.vorbin
                args.only_seg, args.vorbin = False, False
                image = get_cutout(full_field_filename, pos, size, target_header, args)
                image_wht = get_cutout(str(full_field_filename).replace('sci', 'wht'), pos, size, target_header, args)
                args.only_seg, args.vorbin = temp1, temp2

    image = np.roll(image, args.ndelta_xpix, axis=0)
    image = np.roll(image, args.ndelta_ypix, axis=1)

    image = trim_image(image, args, skip_re_trim=False)

    if image_wht is not None:
        image_wht = np.roll(image_wht, args.ndelta_xpix, axis=0)
        image_wht = np.roll(image_wht, args.ndelta_ypix, axis=1)

        image_wht = trim_image(image_wht, args, skip_re_trim=False)
        image_err = 0.1 * np.abs(image_wht) # 1 / np.sqrt(image_wht)
    else:
        image_err = None

    if args.only_seg: seg_mask = args.segmentation_map != args.id
    else: seg_mask = False
    image = np.ma.masked_where(seg_mask, image)
    if image_err is not None: image_err = np.ma.masked_where(seg_mask, image_err)        

    if args.vorbin and not skip_vorbin:
        # -----------discarding low-snr pixels BEFORE vorbin, if any-----------------
        if args.snr_cut is not None:
            snr_map = image / image_err
            snr_mask = (~np.isfinite(snr_map)) | (snr_map < args.snr_cut)
            image = np.ma.masked_where(snr_mask, image.data)
            if image_err is not None: image_err = np.ma.masked_where(snr_mask, image_err.data)

        # ------smoothing the map before voronoi binning--------
        smoothing_kernel = Box2DKernel(args.kernel_size, mode=args.kernel_mode)
        image = np.ma.masked_where(image.mask, convolve(image.data, smoothing_kernel))
        if image_err is not None: 
            image_err = np.ma.masked_where(image_err.mask, convolve(image_err.data, smoothing_kernel))
            image, image_err = bin_2D(image, args.voronoi_bin_IDs, map_err=image_err)
        else:
            image = bin_2D(image, args.voronoi_bin_IDs)

    if True: ##image_err is None:
        image_err = np.zeros(np.shape(image))
    image = unp.uarray(image, image_err)

    return image, exptime

# --------------------------------------------------------------------------------------------------------------------
def plot_starburst_map(full_hdu, axes, args, radprof_axes=None, vorbin_axes=None, snr_axes=None):
    '''
    Plots the Ha map, direct F115W map and their ratio (starbursty-ness map) in a given axis handle
    Returns the axis handle and the ratio map just produced
    '''
    # ---------getting the Ha map-------------
    ha_map, _, _, _, _ = get_emission_line_map('Ha', full_hdu, args, dered=True)
    logOH_map, _, _ = get_Z(full_hdu, args)
    ha_map = correct_ha_C20(ha_map, logOH_map)

    # ---------getting the direct image-------------
    filter = 'F115W'
    direct_map, _ = get_direct_image(full_hdu, filter, args)
    if direct_map is None:
        fig = plt.gcf()
        for index, ax in enumerate(np.atleast_1d(axes)):
            fig.delaxes(ax)
            if args.plot_radial_profiles: fig.delaxes(np.atleast_1d(radprof_axes)[index])
            if args.plot_vorbin: fig.delaxes(np.atleast_1d(vorbin_axes)[index])
            if args.plot_snr: fig.delaxes(np.atleast_1d(snr_axes)[index])
        return axes, None, None

    if args.vorbin:
        direct_map
    # ---------getting the ratio and seg maps-------------
    new_mask = unp.nominal_values(direct_map.data) == 0
    direct_map[new_mask] = 1 # arbitrary fill value to bypass unumpy's inability to handle math domain errors
    direct_map = np.ma.masked_where(new_mask, direct_map)

    ratio_map = ha_map.data / direct_map.data
    ratio_map = np.ma.masked_where(ha_map.mask | direct_map.mask | new_mask, ratio_map)

    # --------making arrays for subplots-------------
    maps_dict = {'direct':direct_map, 'ha':ha_map, 'ratio':ratio_map}
    labels_dict = {'direct':filter, 'ha':r'H$\alpha$', 'ratio':r'H$\alpha$/' + filter}
    #lims_dict = {'direct':[-2.7, -1.2], 'ha':[-19, -17], 'ratio':[-17, -14]}
    lims_dict = {'direct':[None, None], 'ha':[-19, -17], 'ratio':[None, None]}
    if len(np.atleast_1d(axes)) > 1: sequence_to_plot = ['direct', 'ha', 'ratio']
    else: sequence_to_plot = ['ratio']

    # ---------plotting-------------
    starburst_radfit = []
    for index, ax in enumerate(np.atleast_1d(axes)):
        map_err = np.ma.masked_where(maps_dict[sequence_to_plot[index]].mask, unp.std_devs(maps_dict[sequence_to_plot[index]].data))
        ax, radprof_fit = plot_2D_map(maps_dict[sequence_to_plot[index]], ax, args, takelog=True, label=labels_dict[sequence_to_plot[index]], cmap='viridis', vmin=lims_dict[sequence_to_plot[index]][0], vmax=lims_dict[sequence_to_plot[index]][1], radprof_ax=np.atleast_1d(radprof_axes)[index] if args.plot_radial_profiles else None, vorbin_ax=np.atleast_1d(vorbin_axes)[index] if args.plot_vorbin else None, snr_ax=np.atleast_1d(snr_axes)[index] if args.plot_snr else None, image_err=map_err if args.plot_snr else None, hide_yaxis=True, hide_xaxis=False)
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

    axes, starburst_map, starburst_radfit = plot_starburst_map(full_hdu, axes, args, radprof_axes=radprof_axes if 'radprof_axes' in locals() else None, vorbin_axes=vorbin_axes if 'vorbin_axes' in locals() else None, snr_axes=snr_axes if 'snr_axes' in locals() else None)

    return fig, starburst_map, starburst_radfit

# --------------------------------------------------------------------------------------------------------------------
def plot_metallicity_fig(full_hdu, args):
    '''
    Plots the metallicity map, and optionally radial profile, in a new figure
    Saves the 2D metallicity map (along with associated uncertainty and voronoi bins if applicable) as fits file
    Returns the figure handle and the metallicity map just produced
    '''
    # ---------lotting metallicity profiles and gradients----------------------
    df_logOH_radfit = pd.DataFrame(columns=['field', 'objid', 'logOH_int', 'logOH_int_u', 'logOH_sum', 'logOH_sum_u', 'logOH_slope', 'logOH_slope_u', 'logOH_cen', 'logOH_cen_u', 'logOH_diagnostic', 'logOH_branch', 'AGN_diag', 'extent_value', 'extent_unit'])

    # -------looping over all Zdiagnostics-----------
    for index, args.Zdiag in enumerate(args.Zdiag_arr):
        print(f'\n\nCommencing Zdiag {args.Zdiag} which is {index + 1} of {len(args.Zdiag_arr)}..')

        # ---------setting up figure--------------
        ncols = 1
        if args.plot_snr: ncols += 1
        if args.plot_radial_profiles: ncols += 1
        if args.plot_ionisation_parameter: ncols += 1

        fig_size_dict = {2: [10, 4, 0.05, 0.97, 0.15, 0.9, 0.2, 0.], 3: [14, 4, 0.05, 0.97, 0.15, 0.9, 0.2, 0.], 4: [14, 3.5, 0.06, 0.99, 0.15, 0.95, 0.25, 0.]} # figsize_w, figsize_h, l, b, r, t, ws, hs

        fig, axes = plt.subplots(1, ncols, figsize=(fig_size_dict[ncols][0], fig_size_dict[ncols][1]))
        if ncols > 3:
            ip_ax = axes[0]
            axes = axes[1:]

        ax = axes[0]
        axes = axes[1:]

        if args.plot_snr:
            snr_ax = axes[0]
            axes = axes[1:]
        else:
            snr_ax = None
        if args.plot_radial_profiles:
            radprof_ax = axes[0]
        else:
            radprof_ax = None

        fig.subplots_adjust(left=fig_size_dict[ncols][2], right=fig_size_dict[ncols][3], bottom=fig_size_dict[ncols][4], top=fig_size_dict[ncols][5], wspace=fig_size_dict[ncols][6], hspace=fig_size_dict[ncols][7])

        # --------deriving the metallicity map-------------
        if args.Zdiag == 'NB': logOH_map, logOH_int, logOH_sum, line_label_array = get_Z(full_hdu, args)
        else: logOH_map, logOH_int, logOH_sum = get_Z(full_hdu, args)

        # --------deriving the ionisation parameter map-------------
        if args.plot_ionisation_parameter: logq_map, logq_int = get_q_O32(full_hdu, args)

        # ---------plotting-------------
        if logOH_map is not None:
            lim = [7.0, 9.2]
            ax, logOH_radfit = plot_2D_map(logOH_map, ax, args, takelog=False, label=r'Z (%s)$_{\rm int}$ = %.1f $\pm$ %.1f' % (args.Zdiag, logOH_int.n, logOH_int.s), cmap='viridis', radprof_ax=radprof_ax, hide_yaxis=True if args.plot_ionisation_parameter else False, vmin=lim[0], vmax=lim[1], metallicity_multi_color=args.Zdiag == 'P25')
            if args.plot_snr:
                OH_map = 10 ** logOH_map.data
                logOH_map_snr = np.ma.masked_where(logOH_map.mask, unp.nominal_values(OH_map) / unp.std_devs(OH_map))
                snr_ax, _ = plot_2D_map(logOH_map_snr, snr_ax, args, takelog=False, hide_yaxis=True, label=r'Z (%s) SNR' % (args.Zdiag), cmap='cividis', vmin=0, vmax=6)
            if args.plot_ionisation_parameter:
                ip_ax, _ = plot_2D_map(logq_map, ip_ax, args, takelog=False, hide_yaxis=False, label=r'log q$_{\rm int}$ = %.1f $\pm$ %.1f' % (logq_int.n, logq_int.s), cmap='viridis', vmin=6.5, vmax=8.5)
        else:
            logOH_radfit = [ufloat(np.nan, np.nan), ufloat(np.nan, np.nan)]
            plt.close(fig)

        # -----------saving the figure-------------
        if not args.only_integrated:
            # ---------decorating and saving the figure------------------------------
            fig.text(0.05, 0.98, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
            Zbranch_text = '' if args.Zdiag in ['NB', 'P25', 'Te'] else f'-{args.Zbranch}'
            exclude_text = f'_without_{",".join(args.exclude_lines)}' if len(args.exclude_lines) > 0 and args.Zdiag == 'NB' else ''
            figname = fig_dir / f'{args.field}_{args.id:05d}_metallicity_maps{radial_plot_text}{snr_text}{only_seg_text}{vorbin_text}_Zdiag_{args.Zdiag}{Zbranch_text}{exclude_text}.png'
            fig.savefig(figname, transparent=args.fortalk)
            print(f'Saved figure at {figname}')
            plt.show(block=False)

            # ------------------appending to Zgrad dataframe--------------------------
            df_logOH_radfit.loc[len(df_logOH_radfit)] = [args.field, args.id, logOH_int.n, logOH_int.s, logOH_sum.n, logOH_sum.s, logOH_radfit[0].n, logOH_radfit[0].s, logOH_radfit[1].n, logOH_radfit[1].s, args.Zdiag, args.Zbranch, args.AGN_diag, args.arcsec_limit if args.re_limit is None else args.re_limit, 'arcsec' if args.re_limit is None else 're']

    # ------------------writing out Z gradient fits, for making MZGR plot later--------------------------
    if not args.only_integrated:
        outfilename = args.output_dir / 'catalogs' / f'logOHgrad_df{snr_text}{only_seg_text}{vorbin_text}.csv'
        df_logOH_radfit.to_csv(outfilename, index=None, mode='a', header=not os.path.exists(outfilename))
        print(f'Appended metallicity gradient fits to catalog file {outfilename}')

    return fig, logOH_map, logOH_int, logOH_sum, logOH_radfit

# --------------------------------------------------------------------------------------------------------------------
def get_AGN_func_methods(args):
    '''
    Determine which AGN demarcation method/s and line ratio combinations to be used in the BPT, given the BPT diagnostic requested
    Returns array of strings
    '''
    if args.AGN_diag in ['VO87', 'K01', 'S24']:
        args.ynum_line, args.yden_line, args.xnum_line, args.xden_line, theoretical_lines = 'OIII', 'Hb', 'SII', 'Ha', ['K01', 'S24']
    elif args.AGN_diag in ['H21', 'B22']:
        args.ynum_line, args.yden_line, args.xnum_line, args.xden_line, theoretical_lines = 'OIII', 'Hb', 'SII', 'Ha', ['H21', 'B22_S2Ha']
    elif args.AGN_diag == 'O2O3':
        args.ynum_line, args.yden_line, args.xnum_line, args.xden_line, theoretical_lines = 'OIII', 'Hb', 'OII', 'OIII', []
    elif args.AGN_diag == 'O2Hb':
        args.ynum_line, args.yden_line, args.xnum_line, args.xden_line, theoretical_lines = 'OII', 'Hb', 'OIII', 'Hb', []
    elif args.AGN_diag == 'Ne3O2':
        args.ynum_line, args.yden_line, args.xnum_line, args.xden_line, theoretical_lines = 'OIII', 'Hb', 'NeIII-3867', 'OII', ['NB', 'F24', 'B22_Ne3O2']
    else:
        sys.exit('Choose AGN_diag to be one among VO87,H21,O2O3,O2Hb,Ne3O2')

    label_dict = {'K01': 'Kewley+2001', 'S24': 'Schultz+2024', 'H21':'Henry+2021', 'B22_S2Ha':'Backhaus+2022', 'B22_Ne3O2':'Backhaus+2022', 'MAPPINGS':'MAPPINGS', 'F24':'Feuillet+2024', 'NB':'This work'}
    line_labels = [label_dict[item] for item in theoretical_lines]

    return theoretical_lines, line_labels

# --------------------------------------------------------------------------------------------------------------------
def AGN_func(x, theoretical_line):
    '''
    Equation for AGN demarcation line on R3-S2 BPT, from different literature sources
    '''
    if theoretical_line == 'K01': # Eq 6 of Kewley+2001 (https://iopscience.iop.org/article/10.1086/321545)
        y = 1.3 + 0.72 / (x - 0.32)
    elif theoretical_line == 'S24': # Eq 2 of Schultz+2024 (https://arxiv.org/abs/2311.18731), parameters from Table 3 for S2
        y = np.piecewise(x, [x >= -0.92, x < -0.92], [lambda x: (0.78 / (x - 0.34)) + 1.36, lambda x: -0.91 - 1.79 * x])
    elif theoretical_line == 'H21': # Eq 2 of Henry+2021 (https://iopscience.iop.org/article/10.3847/1538-4357/ac1105/pdf)
        y = np.piecewise(x, [x < -0.14, x >= -0.14], [lambda x: 1.27 + (0.28 / (x + 0.14)), lambda x: -np.inf])
    elif theoretical_line == 'B22_S2Ha':  # Eq 2 of Backhaus+2022 (https://iopscience.iop.org/article/10.3847/1538-4357/ac3919/pdf)
        y = np.piecewise(x, [x < -0.11, x >= -0.11], [lambda x: 1.3 + (0.48 / (1.09 * x + 0.12)), lambda x: -np.inf])
    elif theoretical_line == 'B22_Ne3O2':  # Eq 3 of Backhaus+2022 (https://iopscience.iop.org/article/10.3847/1538-4357/ac3919/pdf)
        y = np.piecewise(x, [x < 0.286, x >= 0.286], [lambda x: 0.64 + (0.35 / (2.8 * x - 0.8)), lambda x: -np.inf])
    elif theoretical_line == 'F24':  # Eq 6 of Feuillet+2024 (https://arxiv.org/abs/2312.17381)
        y = np.sqrt((x + 3.4) / 1.8) - 0.55
    elif theoretical_line == 'MAPPINGS':  # new eqn based on MAPPINGS models
        y = np.poly1d([0.05, -0.25, 0.15, 0.89])(x)
    elif theoretical_line == 'NB':  # new eqn based on NebulaBayes models
        y = np.poly1d([0.09, -0.03, 0.01, 0.70])(x) #without accounting for HeI, so x-axis is NeIII/OII
    else:
        sys.exit(f'Requested theoreitcal line {theoretical_line} should be one of K01,S24,H21,B@@_S2Ha,B22_Ne3O2')
    return y

# --------------------------------------------------------------------------------------------------------------------
def overplot_AGN_line_on_BPT(ax, theoretical_line, label, color='k', fontsize=10, lw=2, ls='dashed'):
    '''
    Overplots a given AGN demarcation line on R3 vs S2 ratio BPT, on an existing axis
    Returns axis handle
    '''
    x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
    y = AGN_func(x, theoretical_line)
    ax.plot(x, y, c=color, ls=ls, lw=lw, label=label)
    if label is not None: ax.legend(loc='lower left', fontsize=fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def take_safe_log_ratio(num_map, den_map, skip_log=False):
    '''
    Takes the log of ratio of two 2D masked arrays by properly accounting for bad values so as to avoid math errors
    Returns 2D masked array
    '''
    if np.ma.isMaskedArray(num_map):
        net_mask = num_map.mask | den_map.mask
        num_map = num_map.data
        den_map = den_map.data
    else:
        net_mask = False

    bad_mask = (unp.nominal_values(num_map) <= 0) | (unp.nominal_values(den_map) <= 0) | (~np.isfinite(unp.nominal_values(num_map))) | (~np.isfinite(unp.nominal_values(den_map))) | (~np.isfinite(unp.std_devs(num_map))) | (~np.isfinite(unp.std_devs(den_map)))
    num_map[bad_mask] = 1e-9  # arbitrary fill value to bypass unumpy's inability to handle math domain errors
    den_map[bad_mask] = 1e-9  # arbitrary fill value to bypass unumpy's inability to handle math domain errors
    ratio_map = num_map / den_map
    if not skip_log: ratio_map = unp.log10(ratio_map)
    ratio_map[bad_mask | net_mask] = -99.
    ratio_map = np.ma.masked_where(bad_mask | net_mask, ratio_map)

    return ratio_map

# --------------------------------------------------------------------------------------------------------------------
def take_safe_log_sum(map1, map2, skip_log=False):
    '''
    Takes the log of the sum of two 2D masked arrays by properly accounting for bad values so as to avoid math errors
    Returns 2D masked array
    '''
    if np.ma.isMaskedArray(map1):
        net_mask = map1.mask | map2.mask
        map1 = map1.data
        map2 = map2.data
    else:
        net_mask = False

    bad_mask = (unp.nominal_values(map1) <= 0) | (unp.nominal_values(map2) <= 0) | (~np.isfinite(unp.nominal_values(map1))) | (~np.isfinite(unp.nominal_values(map2))) | (~np.isfinite(unp.std_devs(map1))) | (~np.isfinite(unp.std_devs(map2)))
    map1[bad_mask] = 1e-9  # arbitrary fill value to bypass unumpy's inability to handle math domain errors
    map2[bad_mask] = 1e-9  # arbitrary fill value to bypass unumpy's inability to handle math domain errors
    sum_map = map1 + map2
    if not skip_log: sum_map = unp.log10(sum_map)
    sum_map[bad_mask | net_mask] = -99.
    sum_map = np.ma.masked_where(bad_mask | net_mask, sum_map)

    return sum_map

# --------------------------------------------------------------------------------------------------------------------
def annotate_BPT_axes(scatter_plot_handle, ax, args, theoretical_lines=[], line_labels=[], color_label=''):
    '''
    Annotates the axis labels, limits etc for a BPT diagram in a given axis handle
    Returns the axis handle
    '''

    line_labels_dict = {'OII':'O II', 'SII':'S II 6717+31', 'OIII':'O III', 'Hb':'H beta', 'NeIII-3867':'Ne III', 'Ha':'N II + H alpha' if args.AGN_diag in ['H21', 'B22'] else 'H alpha'}
    ratios_limits_dict = {'OIII/Hb':[-1, 2], 'SII/Ha':[-2, 0.3], 'OII/OIII':[-2, 1], 'OII/Hb':[-1, 1],'NeIII-3867/OII':[-1.6, 0.5]}

    # ---------overplot MAPPINGS models-------
    if args.plot_models:
        args2 = copy.deepcopy(args)
        if args2.AGN_diag in ['H21', 'B22']: args2.xden_line = 'Ha,NII'
        print(f'Overplotting MAPPINGS models for {args2.ynum_line}/{args2.yden_line} vs {args2.xnum_line}/{args2.xden_line}..')
        geom_path_dict = {'s':'sp', 'p':'pp'}  # to determine directory structures based on geometry and iso parameters
        grid_filename = args2.mappings_dir / 'grids' / f'mappings_grid_{geom_path_dict[args2.geometry]}_iso_{args2.iso}.txt'
        df_ratios = pd.read_table(grid_filename, delim_whitespace=True)
        ax, _, _ = plot_ratio_grid(df_ratios, ax, args2)

    # ---------annotate axes-------
    cbar = plt.colorbar(scatter_plot_handle)
    cbar.set_label(color_label, fontsize=args.fontsize)
    cbar.ax.tick_params(labelsize=args.fontsize)

    ax.set_xlim(ratios_limits_dict[f'{args.xnum_line}/{args.xden_line}'])
    ax.set_ylim(ratios_limits_dict[f'{args.ynum_line}/{args.yden_line}'])
    ax.set_xlabel(f'log ({line_labels_dict[args.xnum_line]}/{line_labels_dict[args.xden_line]})', fontsize=args.fontsize)
    ax.set_ylabel(f'log ({line_labels_dict[args.ynum_line]}/{line_labels_dict[args.yden_line]})', fontsize=args.fontsize)
    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    # ---------adding literature AGN demarcation lines----------
    color_arr = ['brown', 'darkgreen', 'dodgerblue', 'cyan', 'sienna']
    for index, (theoretical_line, line_label) in enumerate(zip(theoretical_lines, line_labels)): overplot_AGN_line_on_BPT(ax, theoretical_line=theoretical_line, label=line_label, color=color_arr[index], fontsize=args.fontsize, lw=0.5 if index else 1, ls='dashed' if index else 'solid')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_BPT(full_hdu, ax, args, cmap='viridis', ax_inset=None, hide_plot=False, index=0):
    '''
    Plots spatially resolved BPT diagram based on fluxes from grizli, on an existing axis
    Then overplots theoretical lines
    Returns axis handle and the handle of the spatially resolved scatter plot
    '''
    print(f'Plotting BPT diagram..')
    if args.plot_separately: fig_indiv, ax_indiv = plt.subplots(1, figsize=(8, 6))

    theoretical_lines, line_labels = get_AGN_func_methods(args)
    dist_method = theoretical_lines[0] if len(theoretical_lines) > 0 else None

    # -----------getting the fluxes------------------
    try:
        ynum_map, _, ynum_int, ynum_sum, _ = get_emission_line_map(args.ynum_line, full_hdu, args)
        yden_map, _, yden_int, yden_sum, _ = get_emission_line_map(args.yden_line, full_hdu, args)

        xnum_map, _, xnum_int, xnum_sum, _ = get_emission_line_map(args.xnum_line, full_hdu, args)
        xden_map, _, xden_int, xden_sum, _ = get_emission_line_map(args.xden_line, full_hdu, args)
    except:
        print(f'Required emission lines not available for {args.id} with {args.AGN_diag} AGN diagnostic. So skipping this object')
        return ax, None, None

    if not args.do_not_correct_flux and args.AGN_diag in ['H21', 'B22'] and args.xden_line == 'Ha': # special treatment for H-alpha line, in order to add the NII 6584 component back
        factor = 0.823  # from grizli source code
        print(f'Adding the NII component back to Ha, i.e. dividing by factor {factor} because using Henry+21 for AGN-SF separation')
        xden_map = np.ma.masked_where(xden_map.mask, xden_map.data / factor)
        xden_int = xden_int / factor
        xden_sum = xden_sum / factor

    # -----------integrated-----------------------
    try:
        color = mpl_cm.get_cmap(cmap)(0.5)

        y_ratio = unp.log10(ynum_int / yden_int)
        x_ratio = unp.log10(xnum_int / xden_int)

        # ------distance of integrated value from AGN line--------------
        if dist_method is not None:
            sign = (unp.nominal_values(y_ratio) > AGN_func(unp.nominal_values(x_ratio), dist_method)).astype(int)
            if sign == 0: sign = -1
            distance_from_AGN_line_int = sign * get_distance_from_line(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), AGN_func, dist_method)
        else:
            print(f'\nFor the given combination of BPT line ratios {args.ynum_line}/{args.yden_line} vs {args.xnum_line}/{args.xden_line}, no AGN demarcation line method was found, so distance from AGN line will not be computed or plotted.\n')
            distance_from_AGN_line_int = None

        if not hide_plot:
            p = ax.scatter(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), c=color, marker='o', s=200 / args.fig_scale_factor, lw=2, edgecolor='w' if args.fortalk else 'k', zorder=10)
            ax.errorbar(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), xerr=unp.std_devs(x_ratio), yerr=unp.std_devs(y_ratio), c=color, fmt='none', lw=2)

            if args.plot_separately:
                p = ax_indiv.scatter(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), c=color, marker='o', s=200 / args.fig_scale_factor, lw=2, edgecolor='w' if args.fortalk else 'k', zorder=10)
                ax_indiv.errorbar(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), xerr=unp.std_devs(x_ratio), yerr=unp.std_devs(y_ratio), c=color, fmt='none', lw=2)

    except ValueError:
        print(f'Galaxy {args.id} in {args.field} has a negative integrated flux in one of the following lines, hence skipping this.')
        print(f'{args.ynum_line} = {ynum_int}\n{args.yden_line} = {yden_int}\n{args.xnum_line} = {xnum_int}\n{args.xden_line} = {xden_int}\n')
        distance_from_AGN_line_int = None
        pass

    # -----------summed-----------------------
    try:
        color = mpl_cm.get_cmap(cmap)(0.5)

        y_ratio_sum = unp.log10(ynum_sum / yden_sum)
        x_ratio_sum = unp.log10(xnum_sum / xden_sum)

        # ------distance of integrated value from AGN line--------------
        if dist_method is not None:
            sign = (unp.nominal_values(y_ratio_sum) > AGN_func(unp.nominal_values(x_ratio_sum), dist_method)).astype(int)
            if sign == 0: sign = -1
            distance_from_AGN_line_sum = sign * get_distance_from_line(unp.nominal_values(x_ratio_sum), unp.nominal_values(y_ratio_sum), AGN_func, dist_method)
        else:
            print(f'\nFor the given combination of BPT line ratios {args.ynum_line}/{args.yden_line} vs {args.xnum_line}/{args.xden_line}, no AGN demarcation line method was found, so distance from AGN line will not be computed or plotted.\n')
            distance_from_AGN_line_sum = None

        if not hide_plot:
            p = ax.scatter(unp.nominal_values(x_ratio_sum), unp.nominal_values(y_ratio_sum), c=color, marker='s', s=200 / args.fig_scale_factor, lw=2, edgecolor='w' if args.fortalk else 'k', zorder=10)
            ax.errorbar(unp.nominal_values(x_ratio_sum), unp.nominal_values(y_ratio_sum), xerr=unp.std_devs(x_ratio_sum), yerr=unp.std_devs(y_ratio_sum), c=color, fmt='none', lw=2)

            if args.plot_separately:
                p = ax_indiv.scatter(unp.nominal_values(x_ratio_sum), unp.nominal_values(y_ratio_sum), c=color, marker='s', s=200 / args.fig_scale_factor, lw=2, edgecolor='w' if args.fortalk else 'k', zorder=10)
                ax_indiv.errorbar(unp.nominal_values(x_ratio_sum), unp.nominal_values(y_ratio_sum), xerr=unp.std_devs(x_ratio_sum), yerr=unp.std_devs(y_ratio_sum), c=color, fmt='none', lw=2)

    except ValueError:
        print(f'Galaxy {args.id} in {args.field} has a negative summed flux in one of the following lines, hence skipping this.')
        print(f'{args.ynum_line} = {ynum_sum}\n{args.yden_line} = {yden_sum}\n{args.xnum_line} = {xnum_sum}\n{args.yden_line} = {xden_sum}\n')
        distance_from_AGN_line_sum = None
        pass

    # -----------spatially_resolved-----------------------
    distance_map = get_distance_map(np.shape(ynum_map), args)

    if args.use_variable_N2Ha and args.xden_line == 'Ha':
        if not args.do_not_correct_flux:  # special treatment for H-alpha line, in order to add the NII 6584 component back
            factor = 0.823  # from grizli source code
            print(
                f'Adding the NII component back to Ha, i.e. dividing by factor {factor} because using Henry+21 for AGN-SF separation')
            xden_map = np.ma.masked_where(xden_map.mask, xden_map.data / factor)

        print(f'Correcting Ha for NII component, by factor varying with distance')
        factor = 0.5 + 0.1 * distance_map  # such that in the center (high Z), Ha/(NII + Ha) = 0.5 and at 4 kpc (low Z),  Ha/(NII + Ha) = 0.9
        xden_map = np.ma.masked_where(xden_map.mask, xden_map.data * factor)

    y_ratio = take_safe_log_ratio(ynum_map, yden_map)
    x_ratio = take_safe_log_ratio(xnum_map, xden_map)
    # if args.xnum_line == 'SII' and args.xden_line == 'Ha' and not args.AGN_diag in ['H21', 'B22']: x_ratio = np.ma.masked_where(x_ratio_mask | xnum_map.mask | xden_map.mask | (x_ratio.data > np.log10(0.29)), x_ratio) # using a SII/Ha based DIG cut-off from Petrodjojo+19

    net_mask = y_ratio.mask | x_ratio.mask

    # ------distance of pixels from AGN line--------------
    if dist_method is not None:
        sign_map = (unp.nominal_values(y_ratio.data) > AGN_func(unp.nominal_values(x_ratio.data), dist_method)).astype(int)
        sign_map[sign_map == 0] = -1
        distance_from_AGN_line_map = sign_map * get_distance_from_line(unp.nominal_values(x_ratio.data), unp.nominal_values(y_ratio.data), AGN_func, dist_method)
        distance_from_AGN_line_map = np.ma.masked_where(net_mask, distance_from_AGN_line_map)
    else:
        distance_from_AGN_line_map = None

    # -----------plotting-----------------------
    if not hide_plot:
        x_ratio = np.ma.compressed(np.ma.masked_where(net_mask, x_ratio))
        y_ratio = np.ma.compressed(np.ma.masked_where(net_mask, y_ratio))
        distance_map = np.ma.compressed(np.ma.masked_where(net_mask, distance_map))
        if distance_from_AGN_line_map is not None:
            distance_from_AGN_line_arr = np.ma.compressed(distance_from_AGN_line_map)
            if len(distance_from_AGN_line_arr) > 0:
                dist_lim = max(np.abs(np.max(distance_from_AGN_line_arr)), np.abs(np.min(distance_from_AGN_line_arr)))
            else:
                #distance_from_AGN_line_arr = np.zeros(len(x_ratio.flatten()))
                dist_lim = None
        else:
            distance_from_AGN_line_arr = np.zeros(len(x_ratio.flatten()))

        df = pd.DataFrame({'log_xratio': unp.nominal_values(x_ratio).flatten(), 'log_xratio_err': unp.std_devs(x_ratio).flatten(), 'log_yratio': unp.nominal_values(y_ratio).flatten(), 'log_yratio_err': unp.std_devs(y_ratio).flatten(), 'distance': distance_map.flatten(), 'distance_from_AGN_line': distance_from_AGN_line_arr.flatten()})
        df = df.sort_values(by='distance')
        df = df.drop_duplicates().reset_index(drop=True)

        scatter_plot_handle = ax.scatter(df['log_xratio'], df['log_yratio'], c=df[args.colorcol], marker='o', s=50 / args.fig_scale_factor, lw=0, cmap=args.diverging_cmap if args.colorcol == 'distance_from_AGN_line' else cmap, alpha=0.8, vmin=0 if args.colorcol == 'distance' else -dist_lim if args.colorcol == 'distance_from_AGN_line' and dist_lim is not None else None, vmax=6 if args.colorcol == 'distance' else dist_lim if args.colorcol == 'distance_from_AGN_line' else None)
        ax.errorbar(df['log_xratio'], df['log_yratio'], xerr=df['log_xratio_err'], yerr=df['log_yratio_err'], c='gray', fmt='none', lw=0.5, alpha=0.5 if args.fortalk else 0.5, zorder=-10)

        if args.plot_AGN_frac and distance_from_AGN_line_map is not None and not args.plot_separately and not (
                len(args.id_arr) > 1 and args.plot_BPT):
            if ax_inset is None: ax_inset = ax.inset_axes([0.55, 0.75, 0.3, 0.3])
            # plot_2D_map(factor, ax_inset, args, takelog=False, label='Ha/(NII+Ha)', cmap=args.diverging_cmap, hide_yaxis=not args.plot_BPT, hide_xaxis=not args.plot_BPT)
            plot_2D_map(distance_from_AGN_line_map, ax_inset, args, takelog=False, label=f'Dist from {dist_method}', cmap=args.diverging_cmap, vmin=-dist_lim if dist_lim is not None else None, vmax=dist_lim, hide_yaxis=not args.plot_BPT, hide_xaxis=not args.plot_BPT)

        if args.plot_separately:
            scatter_plot_handle_indiv = ax_indiv.scatter(df['log_xratio'], df['log_yratio'], c=df[args.colorcol], marker='o', s=50 / args.fig_scale_factor, lw=0, cmap=cmap, alpha=0.8, vmin=0 if args.colorcol == 'distance' else -dist_lim if args.colorcol == 'distance_from_AGN_line' else None, vmax=6 if args.colorcol == 'distance' else dist_lim if args.colorcol == 'distance_from_AGN_line' else None)
            ax_indiv.errorbar(df['log_xratio'], df['log_yratio'], xerr=df['log_xratio_err'], yerr=df['log_yratio_err'], c='gray', fmt='none', lw=0.5, alpha=0.1)

            if args.plot_AGN_frac and distance_from_AGN_line_map is not None:
                if ax_inset is None: ax_inset = ax_indiv.inset_axes([0.55, 0.75, 0.3, 0.3])
                plot_2D_map(distance_from_AGN_line_map, ax_inset, args, takelog=False, label=f'Dist from {dist_method}', cmap=args.diverging_cmap, vmin=-dist_lim, vmax=dist_lim)

    # -----------annotating axes-----------------------
    if not hide_plot:
        color_label = 'Distance (kpc)' if args.colorcol == 'distance' else 'Distance from ' + theoretical_lines[0] if args.colorcol == 'distance_from_AGN_line' else ''
        if args.plot_separately:
            ax_indiv = annotate_BPT_axes(scatter_plot_handle_indiv, ax_indiv, args, color_label=color_label, theoretical_lines=theoretical_lines, line_labels=line_labels)
            fig_indiv.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.95)
            fig_indiv.text(0.15, 0.9, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')

            # -----------save figure----------------
            figname = fig_dir / f'{args.field}_{args.id:05d}_BPT{snr_text}{only_seg_text}{vorbin_text}_AGN_diag_{args.AGN_diag}.png'
            fig_indiv.savefig(figname, transparent=args.fortalk)
            print(f'Saved figure at {figname}')
            plt.show(block=False)

        if index == 0:
            ax = annotate_BPT_axes(scatter_plot_handle, ax, args,  color_label=color_label, theoretical_lines=theoretical_lines, line_labels=line_labels)

    return ax, distance_from_AGN_line_map, distance_from_AGN_line_int

# --------------------------------------------------------------------------------------------------------------------
def plot_DIG_maps(full_hdu, axes, args, radprof_axes=None, snr_axes=None):
    '''
    Plots spatially resolved DIG contribution based on Tomicic+21, on an existing axis
    Returns axis handle and the handle of the spatially resolved scatter plot
    '''
    print(f'Plotting DIG diagnostics..')
    # ---------getting the Ha and SII maps-------------
    ha_map, _, _, _, _ = get_emission_line_map('Ha', full_hdu, args, dered=True)
    sii_map, _, _, _, _ = get_emission_line_map('SII', full_hdu, args, dered=True)

    # ---------getting the S2Ha ratio map-------------
    S2Ha_map = take_safe_log_ratio(sii_map, ha_map, skip_log=True)

    # ---------getting the metallicity map-------------
    logOH_map, logOH_int = get_Z(full_hdu, args)

    # -----correcting S2Ha by metallicity, following Tomicic+21----------------
    Z_MW = 8.69 # T21 (and Franchetto+20) assumed MW abundance = solar neighbourhood abundance
    Z_ratio = 10 ** (logOH_map.data - Z_MW)
    S2Ha_map_corr = np.ma.masked_where(S2Ha_map.mask | logOH_map.mask, (S2Ha_map.data / Z_ratio))

    # ------converting units of Ha map from ergs/s/cm^2 to ergs/s/kpc^2---------
    quant = (unp.nominal_values(ha_map.data) * u.erg / u.second / u.cm ** 2).to(u.erg / u.second / (u.kpc ** 2)).value
    quant_u = (unp.std_devs(ha_map.data) * u.erg / u.second / u.cm ** 2).to(u.erg / u.second / (u.kpc ** 2)).value
    ha_map = np.ma.masked_where(ha_map.mask, unp.uarray(quant, quant_u))

    # --------making arrays for subplots-------------
    maps_dict = {'ha':ha_map, 's2ha':S2Ha_map, 'metal':logOH_map, 's2ha_corr':S2Ha_map_corr}
    labels_dict = {'ha':r'H$\alpha$', 's2ha':r'SII/H$\alpha$', 'metal':f'log(O/H)+12 ({args.Zdiag})', 's2ha_corr':r'SII/H$\alpha$ corr'}
    lims_dict = {'ha':[23.5, 25.5], 's2ha':[0, 2], 'metal':[7.5, 9.2], 's2ha_corr':[0, 2]}
    takelog_dict = {'ha':True, 's2ha':False, 'metal':False, 's2ha_corr':False}
    cmap_dict = {'ha':'BuPu_r', 's2ha':'RdPu_r', 'metal':'viridis', 's2ha_corr':'RdPu_r'}
    if len(np.atleast_1d(axes)) > 1: sequence_to_plot = ['ha', 's2ha', 'metal', 's2ha_corr']
    else: sequence_to_plot = ['s2ha_corr']

    # ---------plotting the 2D maps-------------
    axes = np.atleast_1d(axes)
    for index, quant in enumerate(sequence_to_plot):
        map_err = np.ma.masked_where(maps_dict[quant].mask, unp.std_devs(maps_dict[quant].data))
        axes[index], _ = plot_2D_map(maps_dict[quant], axes[index], args, takelog=takelog_dict[quant], label=labels_dict[quant], cmap=cmap_dict[quant], vmin=lims_dict[quant][0], vmax=lims_dict[quant][1], radprof_ax=np.atleast_1d(radprof_axes)[index] if args.plot_radial_profiles else None, snr_ax=np.atleast_1d(snr_axes)[index] if args.plot_snr else None, image_err=map_err if args.plot_snr else None, hide_yaxis=True, hide_xaxis=False, hide_cbar=False, metallicity_multi_color=args.Zdiag == 'P25' and quant == 'metal')

    # ---------plotting S2Ha_corr vs Ha SB-------------
    bad_mask = S2Ha_map_corr.mask | ha_map.mask
    S2Ha_corr_array = np.ma.compressed(np.ma.masked_where(bad_mask, S2Ha_map_corr))
    ha_array = np.ma.compressed(np.ma.masked_where(bad_mask, ha_map))
    log_ha_array = unp.log10(ha_array)
    logOH_array = np.ma.compressed(np.ma.masked_where(bad_mask, logOH_map))

    ax = axes[-2]
    ax.errorbar(unp.nominal_values(log_ha_array), unp.nominal_values(S2Ha_corr_array), xerr=unp.std_devs(log_ha_array), yerr=unp.std_devs(S2Ha_corr_array), c='grey', fmt='none', lw=1 if args.vorbin else 0.5, alpha=0.2 if args.vorbin else 0.1)
    p = ax.scatter(unp.nominal_values(log_ha_array), unp.nominal_values(S2Ha_corr_array), s=30, lw=0.5, edgecolors='k', c=unp.nominal_values(logOH_array), cmap=args.diverging_cmap if args.Zdiag == 'P25' else 'viridis', vmin=7.5, vmax=9.2)

    ax.set_xlim(23.5, 25.5)
    ax.set_ylim(0, 2)

    ax.set_xlabel(r'log H$\alpha$ (ergs/s/kpc^2)', fontsize=args.fontsize)
    ax.set_ylabel(r'SII/H$\alpha$ corr', fontsize=args.fontsize)

    ax.tick_params(axis='both', which='major', labelsize=args.fontsize)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad='2%')
    cbar = plt.colorbar(p, cax=cax, orientation='vertical')
    cbar.ax.tick_params(labelsize=args.fontsize)
    cbar.set_label(f'log(O/H)+12 ({args.Zdiag})', fontsize=args.fontsize)

    fig = ax.figure
    if args.plot_radial_profiles: fig.delaxes(np.atleast_1d(radprof_axes)[-2])
    if args.plot_snr: fig.delaxes(np.atleast_1d(snr_axes)[-2])

    # ---------computing spatially resolved C_DIG, following Tomicic+21 S3.1-------------
    df = pd.DataFrame({'S2Ha_corr': unp.nominal_values(S2Ha_corr_array), 'S2Ha_corr_u': unp.std_devs(S2Ha_corr_array), \
                       'log_Ha_SB':unp.nominal_values(log_ha_array), 'log_Ha_SB_u':unp.std_devs(log_ha_array), \
                       'logOH':unp.nominal_values(logOH_array), 'logOH_u':unp.std_devs(logOH_array)})
    df['pct_Ha'] = df['log_Ha_SB'].rank(pct=True)
    S2Ha_corr_DIG = np.median(df[df['pct_Ha'] <= 0.05]['S2Ha_corr'])
    S2Ha_corr_dense = np.median(df[df['pct_Ha'] >= 0.95]['S2Ha_corr'])
    c_dig_array = (S2Ha_corr_dense - S2Ha_corr_array) / (S2Ha_corr_dense - S2Ha_corr_DIG)

    if np.isnan(S2Ha_corr_dense) or np.isnan(S2Ha_corr_DIG):
        print(f'Not enough pixels to compute 5 and 95 percentile of Ha surface brightness')
        fig.delaxes(axes[-1])
        if args.plot_radial_profiles: fig.delaxes(np.atleast_1d(radprof_axes)[-1])
        return axes, None

    # ---------fitting C_DIG vs Ha SB------------------
    def cdig_func(x, *popt): return np.piecewise(x, [x <= popt[0], x > popt[0]], [1, lambda x: (popt[0]/x) ** popt[1]])
    popt, pcov = curve_fit(cdig_func, unp.nominal_values(ha_array), unp.nominal_values(c_dig_array), p0=[1e24, 0.7], sigma=unp.std_devs(c_dig_array), absolute_sigma=True)

    # ---------deriving the fitted C_DIG map and uncertainty------------------
    c_dig_map = cdig_func(unp.nominal_values(ha_map), *popt)
    c_dig_map_lowlim = cdig_func(unp.nominal_values(ha_map) + unp.std_devs(ha_map), *popt)
    c_dig_map_uplim = cdig_func(unp.nominal_values(ha_map) - unp.std_devs(ha_map), *popt)
    c_dig_map_u = np.mean(np.array([c_dig_map_lowlim, c_dig_map_uplim]), axis=0)
    c_dig_map = np.ma.masked_where(S2Ha_map_corr.mask, unp.uarray(c_dig_map, c_dig_map_u))

    # -------plotting the fit-------------
    fit_color = 'g'
    x_array = 10 ** np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 10)
    c_dig_array = cdig_func(x_array, *popt)
    S2Ha_corr_array_fit = (1 - c_dig_array) * S2Ha_corr_dense + c_dig_array * S2Ha_corr_DIG
    ax.plot(np.log10(x_array), S2Ha_corr_array_fit, lw=2, c=fit_color)
    ax.text(0.99, 0.99, f'log(f0)={np.log10(popt[0]):.2f}, ' + r'$\beta$' + f'={popt[1]:.2f}', c=fit_color, ha='right', va='top', fontsize=args.fontsize, transform=ax.transAxes)

    # ---------plotting the C_DIG map-------------
    dig_cmap = get_combined_cmap([0., 0.3, 0.7, 1.], ['binary', 'spring', 'hot']) # making C_DIG colormap, broken at 0.3 and 0.7
    map_err = np.ma.masked_where(c_dig_map.mask, unp.std_devs(c_dig_map.data))
    axes[-1], _ = plot_2D_map(c_dig_map, axes[-1], args, takelog=False, label=r'C$_{\mathrm{DIG}}$', cmap=dig_cmap, vmin=0, vmax=1, radprof_ax=np.atleast_1d(radprof_axes)[-1] if args.plot_radial_profiles else None, snr_ax=np.atleast_1d(snr_axes)[-1] if args.plot_snr else None, image_err=map_err if args.plot_snr else None, hide_yaxis=True, hide_xaxis=False, hide_cbar=False)
    if args.plot_radial_profiles:
        ax = np.atleast_1d(radprof_axes)[-1]
        for val in [0.3, 0.7]: ax.axhline(val, lw=0.5, ls='--', c='k')

    return axes, c_dig_map

# --------------------------------------------------------------------------------------------------------------------
def plot_DIG_figure(full_hdu, args):
    '''
    Plots the DIG diagnostics (Ha SB map, S2Ha map, metallicity map) in a new figure
    Returns the figure handle and the ratio map just produced
    '''
    nrows, ncols = 1, 6
    if args.plot_snr: nrows += 1
    if args.plot_radial_profiles: nrows += 1

    fig_size_dict = {1: [14, 4, 0.05, 0.97, 0.15, 0.9, 0.2, 0.], \
                     2: [14, 5, 0.05, 0.98, 0.05, 0.95, 0.2, 0.2], \
                     3: [14, 7, 0.05, 0.98, 0.05, 0.95, 0.2, 0.2], \
                     4: [9, 7, 0.02, 0.98, 0.07, 0.95, 0.05, 0.3]}  # figsize_w, figsize_h, l, r, b, t, ws, hs

    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_size_dict[nrows][0], fig_size_dict[nrows][1]))
    if nrows > 1:
        extra_axes = axes[1:,:]
        axes = axes[0,:]
    if args.plot_snr:
        snr_axes = extra_axes[0,:].flatten()
        extra_axes = extra_axes[1:,:]
    if args.plot_radial_profiles:
        radprof_axes = extra_axes[0:, :].flatten()
        extra_axes = extra_axes[1:, :]

    fig.subplots_adjust(left=fig_size_dict[nrows][2], right=fig_size_dict[nrows][3], bottom=fig_size_dict[nrows][4], top=fig_size_dict[nrows][5], wspace=fig_size_dict[nrows][6], hspace=fig_size_dict[nrows][7])

    axes, cdig_map = plot_DIG_maps(full_hdu, axes, args, radprof_axes=radprof_axes if 'radprof_axes' in locals() else None, snr_axes=snr_axes if 'snr_axes' in locals() else None)

    return fig, cdig_map

# --------------------------------------------------------------------------------------------------------------------
def myimshow(data, ax, contour=None, re_pix=None, label='', cmap='viridis'):
    '''
    Utility function to plot a 2D data array on to ax, and plot center and segmentation map
    Returns ax
    '''
    ax.imshow(data, origin='lower', cmap=cmap)
    cen_x, cen_y = int(np.shape(data)[0] / 2), int(np.shape(data)[1] / 2)
    ax.scatter(cen_y, cen_x, marker='x', c='k')
    if contour is not None: ax.contour(contour, levels=0, colors='w', linewidths=0.5)
    ax.text(0.9, 0.9, label, c='w', fontsize=args.fontsize, ha='right', va='top', transform=ax.transAxes)
    if re_pix is not None: ax.add_patch(plt.Circle((cen_y, cen_x), re_pix, color='r', fill=False, lw=0.5))

    return ax

# --------------------------------------------------------------------------------------------------------------------
def myradprof(radius, flux, ax, args, re_kpc, label=''):
    '''
    Utility function to plot the flux radial profile on to ax, and plot the half light radius
    Returns ax
    '''
    ax.scatter(radius, np.cumsum(flux))
    ax.axvline(re_kpc, c='k', ls='dotted')
    ax.axhline(0.5 * np.sum(flux), c='k', ls='dashed')
    ax.text(0.9, 0.9, label, c='k', fontsize=args.fontsize, ha='right', va='top', transform=ax.transAxes)
    ax.text(0.9, 0.1, f're = {re_kpc:.2f} kpc', c='k', fontsize=args.fontsize, ha='right', va='top', transform=ax.transAxes)
    ax.set_xlabel('Radius (kpc)')
    ax.set_ylabel(f'Cumulative flux')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def get_re(full_hdu, args, filter='F200W'):
    '''
    Computes the half-light radius in the direct image with the given filter
    Returns the radius
    '''
    # ----------getting the PSF--------------
    niriss = webbpsf.NIRISS()
    niriss.filter = filter
    niriss.pixelscale = args.pix_size_arcsec
    psf = niriss.calc_psf(fov_arcsec=1 if 'glass' in args.field else 0.5, oversample=1)
    psf_array = psf[0].data
    psf_array /= psf_array.sum()

    # -------reading in direct image----------
    dir_img, filter = read_direct_image(full_hdu, filter=filter)

    # -------compute re----------
    if dir_img is not None:
        dir_img = trim_image(dir_img, args=args, skip_re_trim=True)
        segmentation_map = full_hdu['SEG'].data
        segmentation_map = trim_image(segmentation_map, args, skip_re_trim=True)
        
        dir_img_shifted = np.roll(dir_img, args.ndelta_xpix, axis=0)
        dir_img_shifted = np.roll(dir_img_shifted, args.ndelta_ypix, axis=1)
        segmentation_map_shifted = np.roll(segmentation_map, args.ndelta_xpix, axis=0)
        segmentation_map_shifted = np.roll(segmentation_map_shifted, args.ndelta_ypix, axis=1)

        dir_img_deconvolved = richardson_lucy(dir_img_shifted, psf_array, num_iter=20)

        distance_map = get_distance_map(np.shape(dir_img_shifted), args, for_distmap=True)
        radius = np.ma.compressed(np.ma.masked_where(segmentation_map_shifted != args.id, distance_map)).flatten()
        
        flux_shifted = np.ma.compressed(np.ma.masked_where(segmentation_map_shifted != args.id, dir_img_shifted)).flatten()
        flux_deconvolved = np.ma.compressed(np.ma.masked_where(segmentation_map_shifted != args.id, dir_img_deconvolved)).flatten()

        flux_shifted = [x for _, x in sorted(zip(radius, flux_shifted), key=lambda pair: pair[0])]
        flux_deconvolved = [x for _, x in sorted(zip(radius, flux_deconvolved), key=lambda pair: pair[0])]
        radius = sorted(radius)

        re_kpc_shifted = radius[np.where(np.cumsum(flux_shifted) >= 0.5 * np.sum(flux_shifted))[0][0]]
        re_kpc_deconvolved = radius[np.where(np.cumsum(flux_deconvolved) >= 0.5 * np.sum(flux_deconvolved))[0][0]]
        re_kpc = re_kpc_deconvolved
        
        re_arcsec = re_kpc * cosmo.arcsec_per_kpc_proper(args.z).value # arcsec
        re_pix = re_arcsec / args.pix_size_arcsec
        print(f'For {args.field}:{args.id}: Determined half-light radius from {filter} direct image = {re_kpc:.2f} kpc (={re_arcsec:.2f} arcsec = {re_pix:.1f} pixel)')

        if args.debug_re:
            fig, axes = plt.subplots(2, 3, figsize=(10, 6))
            fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.1, wspace=0.3, hspace=0.1)
            cmap = 'viridis'
            fig.text(0.05, 0.98, f'{args.field}: ID {args.id}: Effective radius diagnostics', fontsize=args.fontsize, c='k', ha='left', va='top')

            axes[0][0] = myimshow(dir_img, axes[0][0], contour=segmentation_map != args.id, re_pix=re_pix, label='Original', cmap=cmap)
            axes[0][1] = myimshow(dir_img_shifted, axes[0][1], contour=segmentation_map_shifted != args.id, re_pix=re_pix, label='Shifted', cmap=cmap)
            axes[1][0] = myimshow(psf_array, axes[1][0], re_pix=re_pix, label='PSF', cmap=cmap)
            axes[1][1] = myimshow(dir_img_deconvolved, axes[1][1], contour=segmentation_map_shifted != args.id, re_pix=re_pix, label='Deconvolved', cmap=cmap)

            axes[0][2] = myradprof(radius, flux_shifted, axes[0][2], args, re_kpc_shifted, label='Light profile (shifted)')
            axes[1][2] = myradprof(radius, flux_deconvolved, axes[1][2], args, re_kpc_deconvolved, label='Light profile (deconvolved)')

            plt.show(block=False)
            sys.exit(f'Exiting here because of --debug_re mode; if you want to run the full code as usual then remove the --debug_re option and re-run')

    else:
        re_kpc, re_arcsec = np.nan, np.nan
        print(f'Could not find image for filter {filter} in {args.field}:{args.id}, so forcing reto be {re}')

    return re_kpc, re_arcsec

# --------------------------------------------------------------------------------------------------------------------
def read_direct_image(full_hdu, filter='F150W'):
    '''
    Reads in the direct image of an object in the given filter, by trying out a few combinations of how the filter name might be in the header
    Returns direct image and the final filter name it found in the header
    '''
    dummy_filter = 'F140W'
    try:
        dir_img = full_hdu['DSCI', f'{filter}-{filter}-CLEAR'].data
    except:
        try:
            dir_img = full_hdu['DSCI', f'{filter}-CLEAR'].data
        except:
            try: 
                dir_img = full_hdu['DSCI', filter].data
            except:
                try:
                    dir_img = full_hdu['DSCI', dummy_filter].data
                    filter = dummy_filter
                except:
                    dir_img = None

    return dir_img, filter

# --------------------------------------------------------------------------------------------------------------------
def get_offsets_from_center(full_hdu, args, filter='F200W'):
    '''
    Computes the offset from the original center of image to the brightest pixel in the direct image with the given filter
    Returns two integers (offset in x and y axes)
    '''
    dir_img, filter = read_direct_image(full_hdu, filter=filter)

    if dir_img is not None:
        dir_img = trim_image(dir_img, args=args, skip_re_trim=True)
        segmentation_map = full_hdu['SEG'].data
        segmentation_map = trim_image(segmentation_map, args, skip_re_trim=True)

        smoothing_kernel = Box2DKernel(5, mode=args.kernel_mode)
        dir_img_smoothed = convolve(dir_img, smoothing_kernel)
        brightest_coords = np.where(dir_img_smoothed == dir_img_smoothed.max())
        brightest_x, brightest_y = brightest_coords[0][0], brightest_coords[1][0]
        cen_x, cen_y = int(np.shape(dir_img)[0] / 2), int(np.shape(dir_img)[1] / 2)
        ndelta_xpix = cen_x - brightest_x
        ndelta_ypix = cen_x - brightest_y
        #ndelta_xpix, ndelta_ypix = 0, 0 ##
        print(f'For {args.field}:{args.id}: Determined x and y offsets from {filter} direct image = {ndelta_xpix}, {ndelta_ypix}')

        if args.debug_offset:
            print(f'Deb2934: original shape = {np.shape(dir_img)}') ##
            print(f'Deb2935: original center = {cen_x}, {cen_y}') ##
            print(f'Deb2934: brightest pixel = {brightest_x}, {brightest_y}') ##

            fig, axes = plt.subplots(1, 3, figsize=(14, 4))
            fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.07, wspace=0.3, hspace=0.1)
            cmap = 'viridis'
            fig.text(0.05, 0.98, f'{args.field}: ID {args.id}: Centering offset diagnostics', fontsize=args.fontsize, c='k', ha='left', va='top')

            axes[0].imshow(dir_img, origin='lower', cmap=cmap)
            axes[0].scatter(cen_y, cen_x, marker='x', c='k')
            axes[0].scatter(brightest_y, brightest_x, marker='x', c='r')
            axes[0].contour(segmentation_map != args.id, levels=0, colors='w', linewidths=0.5)
            axes[0].text(0.9, 0.9, 'Original', c='w', fontsize=args.fontsize, ha='right', va='top', transform=axes[0].transAxes)

            axes[1].imshow(dir_img_smoothed, origin='lower', cmap=cmap)
            axes[1].scatter(cen_y, cen_x, marker='x', c='k')
            axes[1].scatter(brightest_y, brightest_x, marker='x', c='r')
            axes[1].contour(segmentation_map != args.id, levels=0, colors='w', linewidths=0.5)
            axes[1].text(0.9, 0.9, 'Smoothed', c='w', fontsize=args.fontsize, ha='right', va='top', transform=axes[1].transAxes)

            dir_img_shifted = np.roll(dir_img_smoothed, ndelta_xpix, axis=0)
            dir_img_shifted = np.roll(dir_img_shifted, ndelta_ypix, axis=1)
            
            segmentation_map = np.roll(segmentation_map, ndelta_xpix, axis=0)
            segmentation_map = np.roll(segmentation_map, ndelta_ypix, axis=1)

            new_cen_x, new_cen_y = int(np.shape(dir_img_shifted)[0] / 2), int(np.shape(dir_img_shifted)[1] / 2)
            
            axes[2].imshow(dir_img_shifted, origin='lower', cmap=cmap)
            axes[2].scatter(new_cen_y, new_cen_x, marker='x', c='k')
            axes[2].contour(segmentation_map != args.id, levels=0, colors='w', linewidths=0.5)
            axes[2].text(0.9, 0.9, 'Shifted', c='w', fontsize=args.fontsize, ha='right', va='top', transform=axes[2].transAxes)

            plt.show(block=False)
            sys.exit(f'Exiting here because of --debug_offset mode; if you want to run the full code as usual then remove the --debug_offset option and re-run')
    
    else:
        ndelta_xpix = 0
        ndelta_ypix = 0
        print(f'Could not find image for filter {filter} in {args.field}:{args.id}, so forcing center offsets to be {ndelta_xpix}, {ndelta_ypix}')

    return ndelta_xpix, ndelta_ypix

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    if args.fontsize == 10 and args.plot_BPT: args.fontsize = 15
    args.fig_scale_factor = 1.0 if args.plot_direct_filters or args.plot_BPT or args.plot_starburst or args.plot_metallicity else 1.6
    args.fontsize /= args.fig_scale_factor
    if args.plot_slope_vs_mass: args.plot_starburst, args.plot_radial_profiles = True, True
    if args.colorcol == 'ez_z_phot': args.colorcol = 'distance'

    # ---------determining filename suffixes-------------------------------
    radial_plot_text = '_wradprof' if args.plot_radial_profiles else ''
    snr_text = f'_snr{args.snr_cut}' if args.snr_cut is not None else ''
    only_seg_text = '_onlyseg' if args.only_seg else ''
    pixscale_text = '' if args.pixscale == 0.04 else f'_{args.pixscale}arcsec_pix'
    vorbin_text = f'_vorbin_at_{args.voronoi_line}_SNR_{args.voronoi_snr}' if args.vorbin else f'_ellbin_{args.nbins}bins' if args.radbin and args.use_elliptical_bins else f'_radbin_{args.nbins}bins' if args.radbin else ''
    description_text = f'all_diag_plots{radial_plot_text}{snr_text}{only_seg_text}'

    # ---------determining list of fields----------------
    if args.do_all_fields:
        field_list = [os.path.split(item[:-1])[1] for item in glob.glob(str(args.input_dir / 'Par*') + '/')]
        field_list.sort(key=natural_keys)
    else:
        field_list = args.field_arr

    id_arr = args.id

    # --------loop over all fields------------------
    for index2, field in enumerate(field_list):
        start_time2 = datetime.now()
        args.field = f'Par{int(field[3:]):03}' if len(field) < 6 else field
        print(f'\n\nCommencing field {args.field} which is {index2 + 1} of {len(field_list)}..')

        product_dir = args.input_dir / args.field / 'Products'
        output_dir = args.output_dir / args.field
        if args.re_extract: output_dir = output_dir / 're_extracted'
        if args.output_subdir is not None: output_dir = output_dir / args.output_subdir
        output_dir.mkdir(parents=True, exist_ok=True)
        outfilename = output_dir / f'{args.field}_all_diag_results.csv'

        if os.path.exists(outfilename) and not args.clobber and args.write_file:
            print(f'result file for {args.field} already exists as {outfilename}, so skipping this field.')
            continue

        # ---------prep the catalog file--------------------
        catalog_file = product_dir / f'{args.field}_photcat.fits'
        if os.path.exists(catalog_file):
            catalog = GTable.read(catalog_file)
        elif args.do_all_obj:
            print(f'photcat file for {args.field} does not exist, please download from gdrive or run data reduction. Until then, skipping this field.')
            continue

        # --------determine which objects to loop over----------
        if args.do_all_obj:
            if args.re_extract: args.id_arr = ids_to_re_extract_dict[args.field]
            else: args.id_arr = catalog['id']
        else:
            args.id_arr = id_arr

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
        all_lines_to_plot = ['OII', 'Hb', 'OIII', 'Ha', 'NeIII-3867' if args.AGN_diag == 'Ne3O2' else 'SII'] # from blue to red
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
            cmap_arr = np.tile(['Reds_r', 'Greens_r', 'Purples_r', 'Greys_r', 'Oranges_r', 'Blues_r', 'YlGnBu_r', 'BuPu_r', 'GnBu_r', 'spring'], 10)

        # ------------looping over the provided object IDs-----------------------
        for index, args.id in enumerate(args.id_arr):
            start_time3 = datetime.now()
            print(f'\nCommencing ID {args.id} which is {index+1} of {len(args.id_arr)}..')

            # ------determining directories---------
            output_subdir = output_dir / f'{args.id:05d}{pixscale_text}'
            full_fits_file1 = output_subdir / f'{args.field}_{args.id:05d}.full.fits'
            full_fits_file2 = product_dir / 'full' / f'{args.field}_{args.id:05d}.full.fits'
            maps_fits_file = product_dir / 'maps' / f'{args.field}_{args.id:05d}.maps.fits'

            if os.path.exists(full_fits_file1): # if the fits files are in sub-directories for individual objects
                full_filename = full_fits_file1
                od_filename = output_subdir/ f'{args.field}_{args.id:05d}.1D.fits'

            elif os.path.exists(maps_fits_file): # if the fits files are in Products/
                full_filename = maps_fits_file
                od_filename = product_dir / 'spec1D' / f'{args.field}_{args.id:05d}.1D.fits'

            elif os.path.exists(full_fits_file2): # if the fits files are in Products/
                full_filename = full_fits_file2
                od_filename = product_dir / 'spec1D' / f'{args.field}_{args.id:05d}.1D.fits'

            else:
                print(f'Could not find {full_fits_file1} or {full_fits_file2} or {maps_fits_file} for ID {args.id}, so skipping it.')
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
            args.available_lines = np.array(['OIII' if item == 'OIII-5007' else item for item in args.available_lines]) # replace 'OIII-5007' with 'OIII'
            args.z = full_hdu[0].header['REDSHIFT']
            args.ndfilt = full_hdu[0].header['NDFILT']
            args.nlines = full_hdu[0].header['NUMLINES']
            args.distance = cosmo.luminosity_distance(args.z)
            args.pix_arcsec = full_hdu[5].header['PIXASEC']
            args.pa_arr = np.unique([full_hdu[0].header[item] for item in list(full_hdu[0].header.keys()) if 'PA00' in item])
            try:
                obj = catalog[catalog['id']==args.id][0]
                args.mag = obj['mag_auto']
            except:
                print(f'ID {args.id} not found in phot catalog')
                args.mag = np.nan

            if args.use_elliptical_bins:
                args.q = obj['b_image'] / obj['a_image']
                args.pa = obj['theta_image'] # radians

            line_wcs = pywcs.WCS(full_hdu['DSCI'].header)
            args.pix_size_arcsec = utils.get_wcs_pscale(line_wcs)
            imsize_arcsec = full_hdu['DSCI'].data.shape[0] * args.pix_size_arcsec
            offset = args.pix_size_arcsec / 2 # half a pixel offset to make sure cells in 2D plot are aligned with centers and not edges
            args.extent = (-args.arcsec_limit - offset, args.arcsec_limit - offset, -args.arcsec_limit - offset, args.arcsec_limit - offset)
            args.EB_V = 0. # until it gets over-written, if both H alpha and H beta lines are present

            # --------determining true center of object---------------------
            args.ndelta_xpix, args.ndelta_ypix = get_offsets_from_center(full_hdu, args, filter='F150W')

            # --------determining effective radius of object---------------------
            args.re_kpc, args.re_arcsec = get_re(full_hdu, args, filter='F150W')
            if args.re_limit is not None:
                args.pix_size_re = args.pix_size_arcsec / args.re_arcsec
                offset = args.pix_size_re / 2 # half a pixel offset to make sure cells in 2D plot are aligned with centers and not edges
                args.extent = (-args.re_limit - offset, args.re_limit - offset, -args.re_limit - offset, args.re_limit - offset)
            
            # ---------------segmentation map---------------
            segmentation_map = full_hdu['SEG'].data
            segmentation_map = np.roll(segmentation_map, args.ndelta_xpix, axis=0)
            segmentation_map = np.roll(segmentation_map, args.ndelta_ypix, axis=1)
            args.segmentation_map = trim_image(segmentation_map, args)

            # ---------------voronoi binning stuff---------------
            if args.vorbin and args.voronoi_line is not None and not args.only_integrated:
                if args.voronoi_line in args.available_lines:
                    args.voronoi_bin_IDs = get_voronoi_bin_IDs(full_hdu, args.voronoi_snr, plot=args.debug_vorbin, quiet=not args.debug_vorbin, args=args)
                    args.voronoi_bin_distances = get_voronoi_bin_distances(full_hdu, 'OIII', args)
                else:
                    print(f'Requested line for voronoi binning {args.voronoi_line} not available, therefore skipping this object..')
                    continue

                if args.voronoi_bin_IDs is None:
                    print(f'Voronoi binning was not possible, therefore skipping this object..')
                    continue

            # ---------------radially binning stuff---------------
            if args.radbin and not args.only_integrated:
                args.voronoi_bin_IDs = get_radial_bin_IDs(full_hdu, snr_thresh=args.voronoi_snr, plot=args.debug_vorbin, quiet=not args.debug_vorbin, args=args)
                args.voronoi_bin_distances = get_voronoi_bin_distances(full_hdu, 'OIII', args)
                args.vorbin = True

                if args.voronoi_bin_IDs is None:
                    print(f'Voronoi binning was not possible, therefore skipping this object..')
                    continue

            # ---------------dust value---------------
            if all([line in args.available_lines for line in ['Ha', 'Hb']]) and not args.test_cutout:
                try: 
                    EB_V_map, EB_V_int, EB_V_sum = get_EB_V(full_hdu, args, verbose=True, silent=args.only_integrated)
                    args.EB_V = EB_V_sum
                except:
                    args.EB_V = 0.
                    print(f'Could not properly compute EB-V so assigning E(B-V)={args.EB_V}')

            # ---------initialising the starburst figure------------------------------
            if args.test_cutout or args.plot_direct_filters:
                fig = plot_direct_images_all_filters(full_hdu, args)

                # ---------decorating and saving the figure------------------------------
                fig.text(0.05, 0.98, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
                main_text = 'direct_filter_images' if args.plot_direct_filters else 'cutout_tests'
                figname = fig_dir / f'{args.field}_{args.id:05d}_{main_text}{only_seg_text}{vorbin_text}.png'

            # ---------plotting spatially resolved BPT-----------------------------
            elif args.plot_BPT:
                ax, args.distance_from_AGN_line_map, args.distance_from_AGN_line_int = plot_BPT(full_hdu, ax, args, cmap=args.diverging_cmap if args.colorcol == 'distance_from_AGN_line' else cmap_arr[index], index=index)

            # ---------initialising the starburst figure------------------------------
            elif args.plot_starburst:
                fig, ratio_map, starburst_radfit = plot_starburst_figure(full_hdu, args)

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
                if not args.only_integrated and (args.mask_agn or args.Zdiag == 'P25'):
                    _, args.distance_from_AGN_line_map, args.distance_from_AGN_line_int = plot_BPT(full_hdu, None, args, cmap=None, hide_plot=True)  # just to get the distance_from_AGN_line map, without actually plotting the BPT diagram
                else:
                    args.distance_from_AGN_line_map, args.distance_from_AGN_line_int = None, None
                fig, logOH_map, logOH_int, logOH_sum, logOH_radfit = plot_metallicity_fig(full_hdu, args)

            # ---------initialising the metallicity figure------------------------------
            elif args.plot_DIG:
                if all([line in args.available_lines for line in ['OIII', 'Hb', 'SII', 'Ha']]):
                    if args.mask_agn or args.Zdiag == 'P25': _, args.distance_from_AGN_line_map, args.distance_from_AGN_line_int = plot_BPT(full_hdu, None, args, cmap=None, hide_plot=True) # just to get the distance_from_AGN_line map, without actually plotting the BPT diagram
                    fig, dig_map = plot_DIG_figure(full_hdu, args)

                    # ---------decorating and saving the figure------------------------------
                    fig.text(0.05, 0.98, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
                    figname = fig_dir / f'{args.field}_{args.id:05d}_DIG_maps{radial_plot_text}{only_seg_text}{vorbin_text}_Zdiag_{args.Zdiag}.png'
                else:
                    print(f'Necessary lines not present to be able to do DIG diagnostics. So skipping this object..')
                    continue

            # ---------initialising the full figure------------------------------
            else:
                nrow = len(all_lines_to_plot) + 1
                nquant_col = 2
                if args.plot_ratio_maps: nquant_col += 1
                ncol = nquant_col
                if args.plot_snr: ncol += nquant_col
                if args.plot_radial_profiles: ncol += nquant_col

                figsize = (18, 12) if ncol == 9 else (13, 12) if ncol == 6 else (10, 2) if ncol == 3 else (6, 2)
                figsize = (figsize[0] / args.fig_scale_factor, figsize[1] / args.fig_scale_factor)
                fig = plt.figure(figsize=figsize, layout='constrained')

                axis_dirimg = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 0), colspan=1)
                axis_1dspec = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 1), colspan=ncol - 1)
                col_loc = 0

                # ---------setting up emission line map axes----------------
                ax_em_lines = [plt.subplot2grid(shape=(nrow, ncol), loc=(item, col_loc), colspan=1) for item in np.arange(1, nrow)] # OII, H beta, OIII, H alpha, SII maps
                if args.plot_snr:
                    col_loc += 1
                    ax_em_lines_snr = [plt.subplot2grid(shape=(nrow, ncol), loc=(item, col_loc), colspan=1) for item in np.arange(1, nrow)] # OII, H beta, OIII, H alpha, SII SNRs
                if args.plot_radial_profiles:
                    col_loc += 1
                    ax_em_lines_radprof = [plt.subplot2grid(shape=(nrow, ncol), loc=(item, col_loc), colspan=1) for item in np.arange(1, nrow)] # OII, H beta, OIII, H alpha, SII radial profiles

                # ---------setting up line ratio map axes----------------
                if args.plot_ratio_maps:
                    col_loc += 1
                    ax_ratio_maps = [plt.subplot2grid(shape=(nrow, ncol), loc=(item, col_loc), colspan=1) for item in np.arange(1, nrow)] # O3Hb, S2Ha, O3S2, N2S2, O3O2 ratio maps
                    ax_o3hb, ax_s2ha, ax_o3s2, ax_ne3o2, ax_o3o2 = ax_ratio_maps

                    if args.plot_snr:
                        col_loc += 1
                        ax_ratio_maps_snr = [plt.subplot2grid(shape=(nrow, ncol), loc=(item, col_loc), colspan=1) for item in np.arange(1, nrow)] # O3Hb, S2Ha, O3S2, N2S2, O3O2 ratio SNRs
                        ax_o3hb_snr, ax_s2ha_snr, ax_o3s2_snr, ax_ne3o2_snr, ax_o3o2_snr = ax_ratio_maps_snr
                    if args.plot_radial_profiles:
                        col_loc += 1
                        ax_ratio_maps_radprof = [plt.subplot2grid(shape=(nrow, ncol), loc=(item, col_loc), colspan=1) for item in np.arange(1, nrow)] # O3Hb, S2Ha, O3S2, N2S2, O3O2 ratio radial profiles
                        ax_o3hb_radprof, ax_s2ha_radprof, ax_o3s2_radprof, ax_ne3o2_radprof, ax_o3o2_radprof = ax_ratio_maps_radprof

                # ---------setting up measured quantity axes----------------
                col_loc += 1
                ax_quant = [plt.subplot2grid(shape=(nrow, ncol), loc=(item, col_loc), colspan=1) for item in np.arange(1, nrow)] # BPT, log(q), Z, SFR, Ha/F115W maps
                ax1, ax2, ax3, ax4, ax5 = ax_quant
                if args.plot_snr:
                    col_loc += 1
                    ax_quant_snr = [plt.subplot2grid(shape=(nrow, ncol), loc=(item, col_loc), colspan=1) for item in np.arange(1, nrow)] # BPT, log(q), Z, SFR, Ha/F115W SNRs
                    ax1_snr, ax2_snr, ax3_snr, ax4_snr, ax5_snr = ax_quant_snr
                if args.plot_radial_profiles:
                    col_loc += 1
                    ax_quant_radprof = [plt.subplot2grid(shape=(nrow, ncol), loc=(item, col_loc), colspan=1) for item in np.arange(1, nrow)]  # BPT, log(q), Z, SFR, Ha/F115W radial profiles
                    ax1_radprof, ax2_radprof, ax3_radprof, ax4_radprof, ax5_radprof = ax_quant_radprof

                # ---------direct imaging------------------------------
                axis_dirimg = plot_direct_image(full_hdu, axis_dirimg, args)

                # ---------1D spectra------------------------------
                axis_1dspec = plot_1d_spectra(od_hdu, axis_1dspec, args)

                # -----------------emission line maps---------------
                for ind, line in enumerate(all_lines_to_plot):
                    if line in args.available_lines:
                        ax_em_lines[ind] = plot_emission_line_map(line, full_hdu, ax_em_lines[ind], args, cmap='BuPu_r', vmin=-20, vmax=-18, hide_yaxis=False, hide_xaxis=ind < nrow - 2, hide_cbar=False, snr_ax=ax_em_lines_snr[ind] if args.plot_snr else None, radprof_ax=ax_em_lines_radprof[ind] if args.plot_radial_profiles else None)
                    else:
                        fig.delaxes(ax_em_lines[ind])
                        if args.plot_snr: fig.delaxes(ax_em_lines_snr[ind])
                        if args.plot_radial_profiles: fig.delaxes(ax_em_lines_radprof[ind])

                # -----------------emission line ratio maps---------------
                if args.plot_ratio_maps:
                    cmap_ratio = 'RdPu_r'
                    # -----------------emission line ratio maps: O3Hb---------------
                    if all([line in args.available_lines for line in ['OIII', 'Hb']]):
                        ax_o3hb = plot_line_ratio_map('OIII', 'Hb', full_hdu, ax_o3hb, args, cmap=cmap_ratio, vmin=-1, vmax=1, hide_xaxis=True, hide_yaxis=True, hide_cbar=False, snr_ax=ax_o3hb_snr if args.plot_snr else None, radprof_ax=ax_o3hb_radprof if args.plot_radial_profiles else None)
                    else:
                        fig.delaxes(ax_o3hb)
                        if args.plot_snr: fig.delaxes(ax_o3hb_snr)
                        if args.plot_radial_profiles: fig.delaxes(ax_o3hb_radprof)

                    # -----------------emission line ratio maps: S2Ha---------------
                    if all([line in args.available_lines for line in ['SII', 'Ha']]):
                        ax_s2ha = plot_line_ratio_map('SII', 'Ha', full_hdu, ax_s2ha, args, cmap=cmap_ratio, vmin=-1.5, vmax=0.5, hide_xaxis=True, hide_yaxis=True, hide_cbar=False, snr_ax=ax_s2ha_snr if args.plot_snr else None, radprof_ax=ax_s2ha_radprof if args.plot_radial_profiles else None)
                    else:
                        fig.delaxes(ax_s2ha)
                        if args.plot_snr: fig.delaxes(ax_s2ha_snr)
                        if args.plot_radial_profiles: fig.delaxes(ax_s2ha_radprof)

                    # -----------------emission line ratio maps: O3S2---------------
                    if all([line in args.available_lines for line in ['OIII', 'SII']]):
                        ax_o3s2 = plot_line_ratio_map('OIII', 'SII', full_hdu, ax_o3s2, args, cmap=cmap_ratio, vmin=-0.5, vmax=1.5, hide_xaxis=True, hide_yaxis=True, hide_cbar=False, snr_ax=ax_o3s2_snr if args.plot_snr else None, radprof_ax=ax_o3s2_radprof if args.plot_radial_profiles else None)
                    else:
                        fig.delaxes(ax_o3s2)
                        if args.plot_snr: fig.delaxes(ax_o3s2_snr)
                        if args.plot_radial_profiles: fig.delaxes(ax_o3s2_radprof)

                    # -----------------emission line ratio maps: N2S2---------------
                    if all([line in args.available_lines for line in ['NeIII-3867', 'OII']]):
                        ax_ne3o2 = plot_line_ratio_map('NeIII-3867', 'OII', full_hdu, ax_ne3o2, args, cmap=cmap_ratio, vmin=-1, vmax=1, hide_xaxis=True, hide_yaxis=True, hide_cbar=False, snr_ax=ax_ne3o2_snr if args.plot_snr else None, radprof_ax=ax_ne3o2_radprof if args.plot_radial_profiles else None)
                    else:
                        fig.delaxes(ax_ne3o2)
                        if args.plot_snr: fig.delaxes(ax_ne3o2_snr)
                        if args.plot_radial_profiles: fig.delaxes(ax_ne3o2_radprof)

                    # -----------------emission line ratio maps: O3O2---------------
                    if all([line in args.available_lines for line in ['OIII', 'OII']]):
                        ax_o3o2 = plot_line_ratio_map('OIII', 'OII', full_hdu, ax_o3o2, args, cmap=cmap_ratio, vmin=-0.5, vmax=0.5, hide_xaxis=False, hide_yaxis=True, hide_cbar=False, snr_ax=ax_o3o2_snr if args.plot_snr else None, radprof_ax=ax_o3o2_radprof if args.plot_radial_profiles else None)
                    else:
                        fig.delaxes(ax_o3o2)
                        if args.plot_snr: fig.delaxes(ax_o3o2_snr)
                        if args.plot_radial_profiles: fig.delaxes(ax_o3o2_radprof)

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
                #     try: ax5, logOH_map, logOH_radfit, logOH_int = plot_Z_R23_map(full_hdu, ax5, args, radprof_ax=rax5)
                #     except: pass
                # else:
                #     fig.delaxes(ax4)
                #     fig.delaxes(ax5)
                #     if args.plot_radial_profiles:
                #         fig.delaxes(rax4)
                #         fig.delaxes(rax5)

                # ---------------BPT map------------------
                if (args.AGN_diag in ['VO87', 'H21'] and all([line in args.available_lines for line in ['OIII', 'Hb', 'SII', 'Ha']])) or \
                    (args.AGN_diag in ['O2O3', 'O2Hb'] and all([line in args.available_lines for line in ['OIII', 'Hb', 'OII']])) or \
                    (args.AGN_diag in ['Ne3O2'] and all([line in args.available_lines for line in ['OIII', 'Hb', 'OII', 'NeIII-3867']])):
                    ax1, args.distance_from_AGN_line_map, args.distance_from_AGN_line_int = plot_BPT(full_hdu, ax1_radprof, args, cmap='viridis', ax_inset=ax1)
                else:
                    args.distance_from_AGN_line_map = None
                    fig.delaxes(ax1)
                    if args.plot_radial_profiles: fig.delaxes(ax1_radprof)
                if args.plot_snr: fig.delaxes(ax1_snr)

                # ---------------logq map------------------
                if all([line in args.available_lines for line in ['OIII', 'OII']]):
                    ax2, logq_map, logq_radfit, logq_int = plot_q_O32_map(full_hdu, ax2, args, radprof_ax=ax2_radprof if args.plot_radial_profiles else None, snr_ax=ax2_snr if args.plot_snr else None)
                else:
                    fig.delaxes(ax2)
                    if args.plot_radial_profiles: fig.delaxes(ax2_radprof)
                    if args.plot_snr: fig.delaxes(ax2_snr)

                # ---------------metallicity map---------------
                ax3, logOH_map, logOH_radfit, logOH_int = plot_Z_map(full_hdu, ax3, args, radprof_ax=ax3_radprof if args.plot_radial_profiles else None, snr_ax=ax3_snr if args.plot_snr else None)

                # ---------------SFR map------------------
                if 'Ha' in args.available_lines:
                    ax4, SFR_map, SFR_radfit, SFR_int = plot_SFR_map(full_hdu, ax4, args, radprof_ax=ax4_radprof if args.plot_radial_profiles else None, snr_ax=ax4_snr if args.plot_snr else None)
                else:
                    fig.delaxes(ax4)
                    if args.plot_radial_profiles: fig.delaxes(ax4_radprof)
                    if args.plot_snr: fig.delaxes(ax4_snr)

                # ---------------F115W/Ha map------------------
                if 'Ha' in args.available_lines:
                    ax5, SB_map, SB_radfit = plot_starburst_map(full_hdu, ax5, args, radprof_axes=ax5_radprof if args.plot_radial_profiles else None, snr_axes=ax5_snr if args.plot_snr else None)
                else:
                    fig.delaxes(ax5)
                    if args.plot_radial_profiles: fig.delaxes(ax5_radprof)
                    if args.plot_snr: fig.delaxes(ax5_snr)

                # ---------decorating and saving the figure------------------------------
                fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.95)
                fig.text(0.05, 0.98, f'{args.field}: ID {args.id}', fontsize=args.fontsize, c='k', ha='left', va='top')
                figname = fig_dir / f'{args.field}_{args.id:05d}_{description_text}{vorbin_text}.png'

            # --------for talk plots--------------
            if not args.plot_BPT and not args.plot_metallicity:
                if args.fortalk:
                    mplcyberpunk.add_glow_effects()
                    try: mplcyberpunk.make_lines_glow()
                    except: pass
                    try: mplcyberpunk.make_scatter_glow()
                    except: pass

                fig.savefig(figname, transparent=args.fortalk, dpi=1000 if args.fontsize <= 5 else 200)
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

        # ---------plotting spatially resolved BPT-----------------------------
        if args.plot_BPT:
            fig.subplots_adjust(left=0.15, right=0.99, bottom=0.1, top=0.95)
            fig.text(0.15, 0.9, f'{args.field}: IDs {",".join(np.array(args.id_arr).astype(str))}', fontsize=args.fontsize, c='k', ha='left', va='top')

            # -----------save figure----------------
            figname = fig_dir / f'{args.field}_{",".join(np.array(args.id_arr).astype(str))}_BPT{snr_text}{only_seg_text}{vorbin_text}_AGN_diag_{args.AGN_diag}.png'
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
