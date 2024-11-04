'''
    Filename: make_passage_plots.py
    Notes: Plots various quantities from a given passage dataframe (produced by get_field_stats.py by intersecting Venn diagrams of various properties)
    Author : Ayan
    Created: 10-09-24
    Example: run make_passage_plots.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/
             run make_passage_plots.py --plot_conditions EW,mass,PA --xcol lp_mass --ycol lp_SFR --colorcol OIII_EW
             run make_passage_plots.py --plot_conditions EW,mass,PA,a_image --plot_BPT
             run make_passage_plots.py --plot_flux_vs_mag
'''

from header import *
from util import *
from make_diagnostic_maps import compute_EB_V, get_dereddened_flux

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def convert_redshift_to_time(redshift):
    '''
    Computes time (age of the Universe) for a given redshift, and returns the time (in Gyr)
    Currently not complete
    '''
    if redshift == 1.8: time = 10 # Gyr
    elif redshift == 1.2: time = 11 # Gyr

    return time

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_Popesso22(ax, redshift, color='cornflowerblue'):
    '''
    Computes an empirical SFMS based on Popesso+22 (https://arxiv.org/abs/2203.10487) Eq 10, for given redshift
    Then overplots this on a given existing axis handle
    Returns axis handle
    '''
    a0, a1, b0, b1, b2 = 0.2, -0.034, -26.134, 4.722, -0.1925  # Table 2, Eq 10

    log_mass = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 20)
    time = convert_redshift_to_time(redshift) # Gyr
    log_SFR = (a1 * time + b1) * log_mass + b2 * (log_mass) ** 2 + b0 + a0 * time

    ax.plot(log_mass, log_SFR, c=color, lw=2, label=f'Popesso+22: z~{redshift}')

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
def plot_BPT(objid_arr, ax, args):
    '''
    Plots BPT diagram based on integrated fluxes from grizli
    Then overplots theoretical lines
    Returns axis handle
    '''
    print(f'Plotting integrated BPT diagram..')
    y_num, y_den = 'OIII', 'Hb'
    x_num, x_den = 'SII', 'Ha'

    # --------looping through the objects-----------------
    for index, objid in enumerate(objid_arr):
        print(f'Doing galaxy {index + 1} of {len(objid_arr)}..')
        args.field = objid.split('-')[0]
        args.id = int(objid.split('-')[1])

        full_fits_file = args.output_dir / args.field / f'{args.id:05d}' / f'{args.field}_{args.id:05d}.full.fits'
        full_hdu = fits.open(full_fits_file)

        args.available_lines = np.array(full_hdu[0].header['HASLINES'].split(' '))
        args.ndfilt = full_hdu[0].header['NDFILT']

        try:
            Halpha = get_integrated_line_flux('Ha', full_hdu, args, dered=False)
            Hbeta = get_integrated_line_flux('Hb', full_hdu, args, dered=False)
            args.EB_V = compute_EB_V(Halpha, Hbeta)
        except:
            args.EB_V = ufloat(np.nan, np.nan)

        y_num_flux = get_integrated_line_flux(y_num, full_hdu, args)
        y_den_flux = get_integrated_line_flux(y_den, full_hdu, args)
        x_num_flux = get_integrated_line_flux(x_num, full_hdu, args)
        x_den_flux = get_integrated_line_flux(x_den, full_hdu, args)
        z = full_hdu[0].header['REDSHIFT']

        try:
            y_ratio = unp.log10(y_num_flux / y_den_flux)
            x_ratio = unp.log10(x_num_flux / x_den_flux)

            p = ax.scatter(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), c=z, marker='o', s=100, lw=1, edgecolor='k', cmap='rainbow')
            ax.errorbar(unp.nominal_values(x_ratio), unp.nominal_values(y_ratio), xerr=unp.std_devs(x_ratio), yerr=unp.std_devs(y_ratio), c='gray', fmt='none', lw=1)
        except ValueError:
            print(f'Galaxy {args.id} in {args.field} has a negative flux in one of the following, hence skipping this.')
            print(f'{y_num} = {y_num_flux}\n{y_den} = {y_den_flux}\n{x_num} = {x_num_flux}\n{x_den} = {x_den_flux}\n')
            pass

    # ---------adding literature lines from Kewley+2001 (https://iopscience.iop.org/article/10.1086/321545)----------
    x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
    if y_num == 'OIII' and y_den == 'Hb' and x_num == 'NII' and x_den == 'Ha':
        y = 1.19 + 0.61 / (x - 0.47) # Eq 5 of K01
    elif y_num == 'OIII' and y_den == 'Hb' and x_num == 'SII' and x_den == 'Ha':
        y = 1.3 + 0.72 / (x - 0.32) # Eq 6 of K01

    ax.plot(x, y, c='khaki', ls='solid', lw=0.5, label='Kewley+2001')
    plt.legend()

    # ---------annotate axes and save figure-------
    cbar = plt.colorbar(p)
    cbar.set_label('Redshift')

    ax.set_xlabel(f'log ({x_num}/{x_den})')
    ax.set_ylabel(f'log ({y_num}/{y_den})')

    return ax
# --------------------------------------------------------------------------------------------------------------------
def plot_flux_vs_mag(ax, args):
    '''
    Plots line flux vs magnitude based on integrated fluxes from grizli
    Returns axis handle
    '''
    df_filename = args.output_dir / f'all_fields_diag_results.txt'
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


# --------------------------------------------------------------------------------------------------------------------
label_dict = {'lp_mass': r'log M$_*$/M$_{\odot}$', 'lp_SFR': r'log SFR (M$_{\odot}$/yr)', 'ez_z_phot': 'Redshift', 'redshift': 'Redshift'}
bounds_dict = {'lp_mass': (6, 12), 'lp_SFR': (-3, 3), 'ez_z_phot': (0, 3), 'redshift': (0, 3)}
colormap_dict = defaultdict(lambda: 'viridis', ez_z_phot='plasma')

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ----------to initialise axes limits and colormaps-----------------------
    # args.xmin, args.xmax = bounds_dict[args.xcol]
    # args.ymin, args.ymax = bounds_dict[args.ycol]
    # args.cmin, args.cmax = bounds_dict[args.colorcol]
    # args.colormap = colormap_dict[args.colorcol]

   # -----------plotting stuff with the resultant intersecting dataframe--------
    fig, ax = plt.subplots(1, figsize=(8, 6))
    fig.subplots_adjust(left=0.1, right=0.98 if args.plot_flux_vs_mag else 0.85, bottom=0.1, top=0.95)

    # ---------flux vs mag for full sample------
    if args.plot_flux_vs_mag:
        figname = args.output_dir / f'allpar_flux_vs_mag.png'
        ax = plot_flux_vs_mag(ax, args)

    else:
        # -------reading in dataframe produced by get_field_stats.py----------------
        df_infilename = args.output_dir / f'allpar_venn_{",".join(args.plot_conditions)}_df.txt'
        df = pd.read_csv(df_infilename)

        if args.plot_BPT:
            figname = args.output_dir / f'allpar_venn_{",".join(args.plot_conditions)}_BPT.png'
            df['par_obj'] = df['field'].astype(str) + '-' + df['objid'].astype(str)  # making a unique combination of field and object id
            objid_arr = df['par_obj'].values
            ax = plot_BPT(objid_arr, ax, args)

        else:
            figname = args.output_dir / f'allpar_venn_{",".join(args.plot_conditions)}_df_{args.xcol}_vs_{args.ycol}_colorby_{args.colorcol}.png'

            # ---------SFMS from df-------
            p = ax.scatter(df[args.xcol], df[args.ycol], c=df[args.colorcol], marker='s', s=100, lw=1, edgecolor='k')
            cbar = plt.colorbar(p)
            cbar.set_label(label_dict[args.colorcol] if args.colorcol in label_dict else args.colorcol)

            # ---------SFMS from literature-------
            if args.xcol == 'lp_mass' and args.ycol == 'lp_SFR':
                ax = plot_SFMS_Popesso22(ax, 1.8, color='cornflowerblue')
                ax = plot_SFMS_Popesso22(ax, 1.2, color='darkblue')

            # ---------MZR from literature-------
            if args.xcol == 'lp_mass' and 'logOH' in args.ycol:
                ax = plot_MZR_literature(ax)

            # ---------annotate axes and save figure-------
            plt.legend()
            ax.set_xlabel(label_dict[args.xcol] if args.xcol in label_dict else args.xcol)
            ax.set_ylabel(label_dict[args.ycol] if args.ycol in label_dict else args.ycol)

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
