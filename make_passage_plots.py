'''
    Filename: make_passage_plots.py
    Notes: Plots various quantities from a given passage dataframe (produced by get_field_stats.py by intersecting Venn diagrams of various properties)
    Author : Ayan
    Created: 10-09-24
    Example: run make_passage_plots.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/
             run make_passage_plots.py --line_list OIII,Ha --plot_conditions EW,mass,PA --xcol lp_mass --ycol lp_SFR --colorcol OIII_EW
             run make_passage_plots.py --line_list OIII,Ha --plot_conditions EW,mass,PA --xcol lp_mass --ycol logOH_slope --colorcol redshift
             run make_passage_plots.py --line_list OIII,Ha --plot_conditions EW,mass,PA --xcol redshift --ycol logOH_slope --foggie_comp
             run make_passage_plots.py --line_list OIII,Ha --plot_conditions EW,mass,PA,a_image --plot_BPT
             run make_passage_plots.py --plot_flux_vs_mag
             run make_passage_plots.py --line_list OIII,Ha --plot_conditions EW,mass,PA --xcol log_SFR_int --ycol lp_SFR --colorcol redshift

             run make_passage_plots.py --line_list OIII,Ha --plot_conditions EW,mass,PA --xcol log_mass_bgp --ycol lp_mass --colorcol redshift --run narrow_z
             run make_passage_plots.py --line_list OIII,Ha --plot_conditions EW,mass,PA --xcol log_mass_bgp --ycol ez_mass --colorcol redshift --run narrow_z
             run make_passage_plots.py --line_list OIII,Ha --plot_conditions EW,mass,PA --xcol z_bgp --ycol lp_zBEST --colorcol redshift --run narrow_z
             run make_passage_plots.py --line_list OIII,Ha --plot_conditions EW,mass,PA --xcol log_sfr_bgp --ycol log_SFR_int --colorcol redshift --run narrow_z
             run make_passage_plots.py --line_list OIII,Ha --plot_conditions EW,mass,PA --xcol log_mass_bgp --ycol log_SFR_int --colorcol OIII_EW --run narrow_z

             run make_passage_plots.py --line_list Ha --plot_conditions SNR,mass,F115W,F150W,F200W --xcol log_mass_bgp --ycol log_SFR_int --colorcol OIII_EW --run narrow_z
             run make_passage_plots.py --line_list Ha --plot_conditions SNR,mass,F115W,F150W,F200W --xcol log_mass_bgp --ycol log_SFR_int --colorcol Ha_EW --run narrow_z_narrow_mass
             run make_passage_plots.py --line_list Ha --plot_conditions SNR,mass,F115W,F150W,F200W --xcol log_mass_bgp --ycol log_SFR_int --colorcol Ha_EW --run only_st_bands
             run make_passage_plots.py --line_list Ha --plot_conditions SNR,mass,F115W,F150W,F200W --xcol log_mass_bgp --ycol log_SFR_int --colorcol redshift --run including_nircam

             run make_passage_plots.py --line_list Ha --plot_conditions SNR,mass,F115W,F150W,F200W --plot_BPT --colorcol distance_from_K01 --run narrow_z_narrow_mass
             run make_passage_plots.py --line_list Ha --plot_conditions SNR,mass,F115W,F150W,F200W --xcol log_mass_bgp --ycol log_SFR_int --colorcol distance_from_K01 --run narrow_z_narrow_mass
'''

from header import *
from util import *
from make_diagnostic_maps import compute_EB_V, get_dereddened_flux

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_Whitaker14(ax, redshift, color='olivegreen'):
    '''
    Overplots fitted SFMS based on Whitaker+14 (https://iopscience.iop.org/article/10.1088/0004-637X/795/2/104/pdf) Table 1, eq 2
    Then overplots this on a given existing axis handle
    Returns axis handle
    '''
    log_mass_low_lim = 9.2 # lower limit of mass they fitted up to
    if redshift >= 0.5 and redshift < 1.0:
        a, b, c, z1, z2 = -27.40, 5.02, -0.22, 0.5, 1.0 # polynomial coefficients from Table 1, for this z range
    elif redshift >= 1.0 and redshift < 1.5:
        a, b, c, z1, z2 = -26.03, 4.62, -0.19, 1.0, 1.5 # polynomial coefficients from Table 1, for this z range
    elif redshift >= 1.5 and redshift < 2.0:
        a, b, c, z1, z2 = -24.04, 4.17, -0.16, 1.5, 2.0 # polynomial coefficients from Table 1, for this z range
    elif redshift >= 2.0 and redshift < 2.5:
        a, b, c, z1, z2 = -19.99, 3.44, -0.13, 2.0, 2.5 # polynomial coefficients from Table 1, for this z range
    else:
        print(f'Provided redshift {redshift} should be between 0.5 and 2.5, otherwise cannot plot')
        return ax

    log_mass1 = np.linspace(ax.get_xlim()[0], log_mass_low_lim, 20)
    log_SFR1 = np.poly1d([c, b, a])(log_mass1)

    log_mass2 = np.linspace(log_mass_low_lim, ax.get_xlim()[1], 20)
    log_SFR2 = np.poly1d([c, b, a])(log_mass2)

    ax.plot(log_mass1, log_SFR1, ls='dashed', c=color, lw=2)
    ax.plot(log_mass2, log_SFR2, ls='solid', c=color, lw=2, label=f'Whitaker+14: {z1}<z<{z2}')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_Shivaei15(ax, color='salmon'):
    '''
    Overplots fitted SFMS based on Shivaei+15 (https://iopscience.iop.org/article/10.1088/0004-637X/815/2/98/pdf) Table 1
    Then overplots this on a given existing axis handle
    Returns axis handle
    '''
    a1, b1, c1 = 0.65, -5.40, 0.40 # slope, intercept, scatter from Table 1, values for SFR(H alpha), z range 1.37-2.61
    a2, b2, c2 = 0.58, -4.65, 0.36 # slope, intercept, scatter from Table 1, values for SFR(H alpha), z range 2.09-2.61
    log_mass_low_lim = 9.5 # lower limit of mass they fitted up to

    log_mass1 = np.linspace(ax.get_xlim()[0], log_mass_low_lim, 20)
    log_SFR1 = a1 * log_mass1 + b1

    log_mass2 = np.linspace(log_mass_low_lim, ax.get_xlim()[1], 20)
    log_SFR2 = a1 * log_mass2 + b1

    ax.plot(log_mass1, log_SFR1, ls='dashed', c=color, lw=2)
    ax.fill_between(log_mass1, log_SFR1 - c1/2, log_SFR1 + c1/2, alpha=0.3, facecolor=color)
    ax.plot(log_mass2, log_SFR2, ls='solid', c=color, lw=2, label=f'Shivaei+15: z~2')
    ax.fill_between(log_mass2, log_SFR2 - c1/2, log_SFR2 + c1/2, alpha=0.3, facecolor=color)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_Popesso22(ax, redshift, color='cornflowerblue'):
    '''
    Computes an empirical SFMS based on Popesso+22 (https://arxiv.org/abs/2203.10487) Eq 10, for given redshift
    Then overplots this on a given existing axis handle
    Returns axis handle
    '''
    a0, a1, b0, b1, b2 = 0.2, -0.034, -26.134, 4.722, -0.1925  # Table 2, Eq 10
    log_mass_low_lim = 8.7 # lower limit of mass they fitted up to

    age_at_z = cosmo.age(redshift).value # Gyr

    log_mass1 = np.linspace(ax.get_xlim()[0], log_mass_low_lim, 20)
    log_SFR1 = (a1 * age_at_z + b1) * log_mass1 + b2 * (log_mass1) ** 2 + b0 + a0 * age_at_z

    log_mass2 = np.linspace(log_mass_low_lim, ax.get_xlim()[1], 20)
    log_SFR2 = (a1 * age_at_z + b1) * log_mass2 + b2 * (log_mass2) ** 2 + b0 + a0 * age_at_z

    ax.plot(log_mass1, log_SFR1, ls='dashed', c=color, lw=2)
    ax.plot(log_mass2, log_SFR2, ls='solid', c=color, lw=2, label=f'Popesso+22: z~{redshift}')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_MZR_literature(ax):
    '''
    Computes empirical MZRs based on several studies
    Then overplots them on a given existing axis handle
    Returns axis handle
    '''
    log_mass = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 20)

    # ---------Andrews+13-----------
    log_mass_TO, logOH_asym, gamma = 8.901, 8.798, 0.640 # First row of Table 3 in https://iopscience.iop.org/article/10.1088/0004-637X/765/2/140
    log_OH = logOH_asym - np.log10(1 + (log_mass_TO / log_mass)**gamma) # Eq 5 of A13
    ax.plot(log_mass, log_OH, c='cornflowerblue', lw=2, label=f'Andrews+13')

    # ---------Zahid+14-----------
    log_mass_TO, logOH_asym, gamma = 9.08, 10.06, 0.61 # COSMOS row of Table 2 in https://iopscience.iop.org/article/10.1088/0004-637X/791/2/130#apj498704t2
    log_OH = logOH_asym + np.log10(1 - np.exp(-(log_mass/log_mass_TO)**gamma)) # Eq 5 of Z14
    ax.plot(log_mass, log_OH, c='blue', lw=2, label=f'Zahid+14:COSMOS')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def get_integrated_line_flux(line, full_hdu, args, dered=True):
    '''
    To obtain the integrated flux for a given line, given a full_hdu object
    Returns the flux along with uncertainty
    '''
    line_index = np.where(args.available_lines == line)[0][0]
    ext = 5 + 2 * args.ndfilt + 4 * line_index
    line_wave = full_hdu[ext].header['RESTWAVE'] # in Angstrom
    line_int = full_hdu[0].header[f'FLUX{line_index + 1:03d}'] # ergs/s/cm^2
    line_int_err = full_hdu[0].header[f'ERR{line_index + 1:03d}'] # ergs/s/cm^2
    line_int = ufloat(line_int, line_int_err)
    if dered: line_int = get_dereddened_flux(line_int, line_wave, args.EB_V)

    return line_int

# --------------------------------------------------------------------------------------------------------------------
def distance(x, y, x0, y0):
    """
    Return distance between point
    P[x0,y0] and a curve (x,y)
    """
    d_x = x - x0
    d_y = y - y0
    dis = np.sqrt(d_x ** 2 + d_y ** 2)
    return dis

# --------------------------------------------------------------------------------------------------------------------
def min_distance(x, y, P, precision=5):
    """
    Compute minimum/a distance/s between
    a point P[x0,y0] and a curve (x,y)
    rounded at `precision`.

    ARGS:
        x, y      (array)
        P         (tuple)
        precision (int)

    Returns min indexes and distances array.
    """
    # compute distance
    d = distance(x, y, P[0], P[1])
    d = np.round(d, precision)
    # find the minima
    glob_min_idxs = np.argwhere(d == np.min(d)).ravel()
    return glob_min_idxs, d

# --------------------------------------------------------------------------------------------------------------------
def get_distance_from_Kewley2001(xdata, ydata, args, x_num='SII'):
    '''
    Computes distance of each object in the given xdata and ydata (line ratios) arrays, from the Kewley+2011 AGN-SF line
    Returns the distance as an array
    '''
    print(f'Computing distance form Kewley+2001 line on the BPT diagram..')
    x = np.linspace(-2, 0, 100)

    if x_num == 'NII':
        y = 1.19 + 0.61 / (x - 0.47) # Eq 5 of K01
    elif x_num == 'SII':
        y = 1.3 + 0.72 / (x - 0.32) # Eq 6 of K01

    min_dist_arr = []
    for P in zip(xdata, ydata):
        min_idxs, distances = min_distance(x, y, P)
        if len(min_idxs) > 0: min_dist = distances[min_idxs[0]]
        else: min_dist = np.nan
        min_dist_arr.append(min_dist)

    return np.array(min_dist_arr)

# --------------------------------------------------------------------------------------------------------------------
def plot_BPT(df, ax, args):
    '''
    Plots BPT diagram based on integrated fluxes from grizli
    Then overplots theoretical lines
    Returns axis handle
    '''
    print(f'Plotting integrated BPT diagram..')
    objid_arr = (df['field'].astype(str) + '-' + df['objid'].astype(str)).values  # making a unique combination of field and object id
    df['distance_from_K01'] = np.ones(len(df)) * np.nan

    x_lim, y_lim = [-2, 0], [-0.5, 1.5]
    y_num, y_den = 'OIII', 'Hb'
    x_num, x_den = 'SII', 'Ha'
    skip_count = 0

    if ('distance' in args.colorcol and 'K01' in args.colorcol) or args.colorcol == 'distance_from_K01':
        cmap = 'BrBG'
        cmin, cmax = -0.8, 0.8
    else:
        cmap = 'rainbow'
        cmin, cmax = 0, 6
    scalarMap = mpl_cm.ScalarMappable(norm=mplcolors.Normalize(vmin=cmin, vmax=cmax), cmap=plt.get_cmap(cmap))

    # ---------adding literature lines from Kewley+2001 (https://iopscience.iop.org/article/10.1086/321545)----------
    x = np.linspace(x_lim[0], x_lim[1], 100)
    if y_num == 'OIII' and y_den == 'Hb' and x_num == 'NII' and x_den == 'Ha':
        def func(x): return 1.19 + 0.61 / (x - 0.47) # Eq 5 of K01
    elif y_num == 'OIII' and y_den == 'Hb' and x_num == 'SII' and x_den == 'Ha':
        def func(x): return 1.3 + 0.72 / (x - 0.32) # Eq 6 of K01
    y = func(x)

    ax.plot(x, y, c='khaki' if args.fortalk else 'brown', ls='dashed', lw=2, label='Kewley+2001')
    plt.legend()

    # --------looping through the objects-----------------
    for index, objid in enumerate(objid_arr):
        print(f'Doing galaxy {index + 1} of {len(objid_arr)}..')
        args.field = objid.split('-')[0]
        args.id = int(objid.split('-')[1])

        full_fits_file = args.output_dir / args.field / f'{args.id:05d}' / f'{args.field}_{args.id:05d}.full.fits'
        if not os.path.exists(full_fits_file): # if the fits files are actually maps.fits
            full_fits_file = args.input_dir / f'{args.field}/Products/maps/{args.field}_{args.id:05d}.maps.fits'
        full_hdu = fits.open(full_fits_file)

        args.available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
        args.ndfilt = full_hdu[0].header['NDFILT']

        try:
            Halpha = get_integrated_line_flux('Ha', full_hdu, args, dered=False)
            Hbeta = get_integrated_line_flux('Hb', full_hdu, args, dered=False)
            args.EB_V = compute_EB_V(Halpha, Hbeta)
        except:
            args.EB_V = ufloat(np.nan, np.nan)

        try:
            y_num_flux = get_integrated_line_flux(y_num, full_hdu, args)
            y_den_flux = get_integrated_line_flux(y_den, full_hdu, args)
            x_num_flux = get_integrated_line_flux(x_num, full_hdu, args)
            x_den_flux = get_integrated_line_flux(x_den, full_hdu, args)
            z = full_hdu[0].header['REDSHIFT']
        except IndexError as e:
            skip_count += 1
            print(f'Skipping {index + 1}: ID {objid} due to {e}')
            continue

        try:
            y_ratio = unp.log10(y_num_flux / y_den_flux)
            x_ratio = unp.log10(x_num_flux / x_den_flux)

            # -------to plot minimum distances to the K01 line-----------
            min_idxs, distances = min_distance(x, y, (unp.nominal_values(x_ratio), unp.nominal_values(y_ratio)))
            if len(min_idxs) > 0:
                sign = 1 if unp.nominal_values(y_ratio) > func(unp.nominal_values(x_ratio)) else -1
                min_dist = sign * distances[min_idxs[0]]
                df.at[index, 'distance_from_K01'] = min_dist

            if ('distance' in args.colorcol and 'K01' in args.colorcol) or args.colorcol == 'distance_from_K01':
                if len(min_idxs) > 0:
                    color = scalarMap.to_rgba(min_dist)
                    #ax.plot([unp.nominal_values(x_ratio), x[min_idxs[0]]], [unp.nominal_values(y_ratio), y[min_idxs[0]]], lw=1, linestyle='dotted', c=color)
                else:
                    color = np.nan
            else:
                color = scalarMap.to_rgba(z)

            p = ax.scatter(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), c=color, marker='o', s=100, lw=1, edgecolor='k')
            ax.errorbar(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), xerr=unp.std_devs(x_ratio), yerr=unp.std_devs(y_ratio), c='gray', fmt='none', lw=1)

        except ValueError:
            skip_count += 1
            print(f'Galaxy {args.id} in {args.field} has a negative flux in one of the following, hence skipping this.')
            print(f'{y_num} = {y_num_flux}\n{y_den} = {y_den_flux}\n{x_num} = {x_num_flux}\n{x_den} = {x_den_flux}\n')
            pass

    # ---------annotate axes and save figure-------
    if args.fontsize == 10: args.fontsize = 15
    cax, kw = matplotlib.colorbar.make_axes_gridspec(ax, orientation='vertical', pad=0.01, fraction=0.1, shrink=1, aspect=35)
    cbar = plt.colorbar(scalarMap, cax=cax,  orientation='vertical')
    cbar.set_label('Redshift' if args.colorcol == 'ez_z_phot' else label_dict[args.colorcol] if args.colorcol in label_dict else args.colorcol, fontsize=args.fontsize)

    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)

    ax.set_xlabel(f'log ({x_num}/{x_den})', fontsize=args.fontsize)
    ax.set_ylabel(f'log ({y_num}/{y_den})', fontsize=args.fontsize)

    print(f'Eventually managed to plot {len(objid_arr) - skip_count} out of {len(objid_arr)} galaxies.')

    return ax, df

# --------------------------------------------------------------------------------------------------------------------
def plot_flux_vs_mag(ax, args):
    '''
    Plots line flux vs magnitude based on integrated fluxes from grizli
    Returns axis handle
    '''
    df_filename = args.output_dir / 'catalogs' / f'all_fields_diag_results.txt'
    df = pd.read_table(df_filename)

    xcol, ycol1, ycol2, colorcol = 'mag', 'OIII_int', 'Ha_int', 'redshift'
    size = 5
    limit = -10.5

    df[ycol1[:-4] + '_limit'] = df[xcol] + 2.5 * np.log10(df[ycol1])
    df[ycol2[:-4] + '_limit'] = df[xcol] + 2.5 * np.log10(df[ycol2])
    df = df.dropna(subset=[xcol, ycol1, ycol2, colorcol], axis='rows')
    #df = df[df['field'] == 'Par061']

    df1 = df[df[ycol1[:-4] + '_limit'] < limit]
    #p = ax.scatter(df1[xcol], np.log10(df1[ycol1]), c=df1[colorcol], cmap='Greens', s=size, lw=0, label=f'{ycol1[:-4]}, total {len(df1)}')

    df2 = df[df[ycol2[:-4] + '_limit'] < limit]
    p = ax.scatter(df2[xcol], np.log10(df2[ycol2]), c=df2[colorcol], cmap='Reds', s=size, lw=0, label=f'{ycol2[:-4]}, total {len(df2)}')

    plt.legend()
    cbar = plt.colorbar(p)
    cbar.set_label(colorcol)

    try:
        plt.tricontour(df[xcol], np.log10(df[ycol1]), df[colorcol], 15, lw=2, colors='g')
        plt.tricontourf(df[xcol], np.log10(df[ycol2]), df[colorcol], 15, lw=2, colors='r')
    except:
        pass

    ax.set_xlim(32, 17)
    ax.set_ylim(-20, -12)

    ax.set_xlabel('Magnitude AUTO (from photcat)')
    ax.set_ylabel('log (Line flux ergs/s/cm^2)')

    return ax

# -------------------------------------global dictionaries-------------------------------------------------------------------------------
label_dict = {'lp_mass': r'log M$_*$/M$_{\odot}$ (LePhare)', 'ez_mass': r'log M$_*$/M$_{\odot}$ (EAZY)', 'log_mass_bgp': r'log M$_*$/M$_{\odot}$ (Bagpipes)', \
              'lp_SFR': r'log SFR (M$_{\odot}$/yr) (LePhare)', 'ez_sfr': r'log SFR (M$_{\odot}$/yr) (EAZY)', 'log_sfr_bgp': r'log SFR (M$_{\odot}$/yr) (Bagpipes)', 'log_SFR_int': r'log SFR (M$_{\odot}$/yr) (Grizli)', \
              'lp_zBEST': 'Redshift (LePhare)', 'ez_z_phot': 'Redshift (EAZY)', 'z_bgp': 'Redshift (Bagpipes)', 'redshift': 'Redshift (Grizli)', \
              'logOH_slope':r'log $\nabla$Z$_r$ (dex/kpc)'}
bounds_dict = {'lp_mass': (6, 9), 'ez_mass': (6, 9), 'log_mass_bgp': (6.5, 10.5), \
               'lp_SFR': (-3, 1), 'ez_sfr': (-3, 1), 'log_sfr_bgp': (-3, 2), 'log_SFR_int': (-3, 2.5), \
               'ez_z_phot': (0, 3), 'lp_zBEST': (0, 3), 'z_bgp': (0, 3), 'redshift': (0.5, 2.2), \
               'logOH_slope': (-0.4, 0.1)}
colormap_dict = defaultdict(lambda: 'viridis', ez_z_phot='plasma', distance_from_K01='BrBG')

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ----------to initialise axes limits and colormaps-----------------------
    # args.cmin, args.cmax = bounds_dict[args.colorcol]
    # args.colormap = colormap_dict[args.colorcol]

   # -----------plotting stuff with the resultant intersecting dataframe--------
    if args.xcol == 'redshift' and args.ycol == 'logOH_slope' and args.foggie_comp:  # special case, to match the FOGGIE plot size
        fig, ax = plt.subplots(1, figsize=(12, 6))
        fig.subplots_adjust(top=0.95, bottom=0.12, left=0.12, right=0.97)
    else:
        fig, ax = plt.subplots(1, figsize=(8, 6))
        fig.subplots_adjust(left=0.1, right=0.98, bottom=0.1, top=0.95)

    # ---------flux vs mag for full sample------
    if args.plot_flux_vs_mag:
        figname = args.output_dir / 'plots' / f'allpar_flux_vs_mag.png'
        ax = plot_flux_vs_mag(ax, args)

    else:
        # -------reading in dataframe produced by get_field_stats.py or by compute_stellar_masses.py----------------
        plot_conditions_text = ','.join(args.line_list) + ',' + ','.join(args.plot_conditions)
        plot_conditions_text = plot_conditions_text.replace('SNR', f'SNR>{args.SNR_thresh}').replace('EW', f'EW>{args.EW_thresh}').replace('a_image', f'a>{args.a_thresh}')
        df_infilename = args.output_dir / 'catalogs' / f'allpar_venn_{plot_conditions_text}_df_withSED_{args.run}.csv'
        df = pd.read_csv(df_infilename)
        print(f'Reading in main df from {df_infilename}')

        # -------combing with metallicity dataframe if it exists----------------
        logOHgrad_filename = args.output_dir / 'catalogs' / f'logOHgrad_df_onlyseg_vorbin_at_Ha_SNR_3.0.txt'
        if os.path.exists(logOHgrad_filename):
            print(f'Reading in and merging logOH gradient df: {logOHgrad_filename}')
            df_logOHgrad = pd.read_csv(logOHgrad_filename)
            df = pd.merge(df, df_logOHgrad, on=['field', 'objid'], how='outer')

        # -------making the dsired plots----------------
        if args.plot_BPT:
            colorby_text = f'_colorby_z' if args.colorcol == 'ez_z_phot' else f'_colorby_{args.colorcol}'
            figname = args.output_dir / 'plots' / f'allpar_venn_{plot_conditions_text}_run_{args.run}_BPT{colorby_text}.png'
            ax, df = plot_BPT(df, ax, args)

            # ------writing out distance from K01 AGN-SF line--------------------------
            if ('distance' in args.colorcol and 'K01' in args.colorcol) or args.colorcol == 'distance_from_K01':
                df.to_csv(df_infilename, index=None)
                print(f'\nAdded distance_from_K01 column to df and saved in {df_infilename}.')

        else:
            figname = args.output_dir / 'plots' / f'allpar_venn_{plot_conditions_text}_run_{args.run}_df_{args.xcol}_vs_{args.ycol}_colorby_{args.colorcol}.png'
            if df['SFR_int'].dtype == object: # accompanied by uncertainty in the same column
                quant_arr = []
                for item in df['SFR_int']:
                    if 'e' in item:
                        pow = float(item[item.find('e')+1:])
                        base = np.array(item[1:item.find('e')-1].split('+/-')).astype(np.float64)
                        quant = base * 10 ** pow
                    else:
                        quant = np.array(item.split('+/-')).astype(np.float64)
                    quant_arr.append(ufloat(quant[0], quant[1]))

                df['SFR_int'] = unp.nominal_values(quant_arr)
                df['SFR_int_u']= unp.std_devs(quant_arr)
                df['log_SFR_int_u'] = unp.std_devs(unp.log10(quant_arr))
            df['log_SFR_int'] = np.log10(df['SFR_int'])


            # ---------SFMS from df-------
            p = ax.scatter(df[args.xcol], df[args.ycol], c='gold' if args.foggie_comp else df[args.colorcol], plotnonfinite=True, marker='*' if args.foggie_comp else 'o', s=1000 if args.foggie_comp else 100, lw=1, edgecolor='w' if args.fortalk else 'k', vmin=bounds_dict[args.colorcol][0] if args.colorcol in bounds_dict else None, vmax=bounds_dict[args.colorcol][1] if args.colorcol in bounds_dict else None, cmap=colormap_dict[args.colorcol])
            if args.ycol + '_u' in df and not args.foggie_comp: # if uncertainty column exists
                ax.errorbar(df[args.xcol], df[args.ycol], yerr=df[args.ycol + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)
            if args.xcol + '_u' in df and not args.foggie_comp: # if uncertainty column exists
                ax.errorbar(df[args.xcol], df[args.ycol], xerr=df[args.xcol + '_u'], c='gray', fmt='none', lw=1, alpha=0.5)

            if not args.foggie_comp:
                cbar = plt.colorbar(p)
                cbar.set_label(label_dict[args.colorcol] if args.colorcol in label_dict else args.colorcol)

            # ---------SFMS from literature-------
            if 'mass' in args.xcol and 'sfr' in args.ycol.lower():
                ax = plot_SFMS_Popesso22(ax, 0.5, color='darkblue')
                ax = plot_SFMS_Popesso22(ax, 1.2, color='blue')
                ax = plot_SFMS_Popesso22(ax, 2.0, color='cornflowerblue')
                ax = plot_SFMS_Shivaei15(ax, color='salmon')
                ax = plot_SFMS_Whitaker14(ax, 0.6, color='forestgreen')
                ax = plot_SFMS_Whitaker14(ax, 1.2, color='limegreen')
                ax = plot_SFMS_Whitaker14(ax, 1.8, color='yellowgreen')

            # ---------MZR from literature-------
            if args.xcol == 'lp_mass' and 'logOH' in args.ycol and 'slope' not in args.ycol:
                ax = plot_MZR_literature(ax)

            # ---------annotate axes and save figure-------
            if len(ax.get_legend_handles_labels()[0]) > 0: plt.legend()
            ax.set_xlabel(label_dict[args.xcol] if args.xcol in label_dict else args.xcol)
            ax.set_ylabel(label_dict[args.ycol] if args.ycol in label_dict else args.ycol)

            if args.xcol in bounds_dict: ax.set_xlim(bounds_dict[args.xcol][0], bounds_dict[args.xcol][1])
            if args.ycol in bounds_dict: ax.set_ylim(bounds_dict[args.ycol][0], bounds_dict[args.ycol][1])

            # ---------comparing SFRs-------
            if 'sfr' in args.xcol.lower() and 'ssfr' not in args.xcol.lower() and 'sfr' in args.ycol.lower() and 'ssfr' not in args.ycol.lower() or ('mass' in args.xcol and 'mass' in args.ycol) or ('z' in args.xcol and 'z' in args.ycol):
                line = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 10)
                ax.plot(line, line, ls='dashed', c='k', lw=1)

            if args.xcol == 'redshift': # redshift axis should be reversed in values
                ax.set_xlim(ax.get_xlim()[1], ax.get_xlim()[0])

            if args.xcol == 'redshift' and args.ycol == 'logOH_slope' and args.foggie_comp: # special case, to match the FOGGIE plot limits
                ax.set_xlim(4, 0.5)
                ax.set_ylim(-0.5, 0.4)

    # --------for talk plots--------------
    if args.fortalk:
        #mplcyberpunk.add_glow_effects()
        #try: mplcyberpunk.make_lines_glow()
        #except: pass
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    fig.savefig(figname, transparent=args.fortalk)
    print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
