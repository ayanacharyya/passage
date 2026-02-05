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

             run make_passage_plots.py  --do_field Par028 --plot_conditions all_match --xcol log_mass_bgp --ycol log_SFR_int --colorcol redshift --run Par028_including_nircam
             run make_passage_plots.py  --do_field Par052 --plot_conditions all_match --xcol log_mass_bgp --ycol log_SFR_int --colorcol redshift --run Par052_including_nircam

             run make_passage_plots.py --do_field Par028 --drv 0.5 --plot_conditions SNR --line_list OIII,Ha,OII,Hb,SII --SNR_thresh 2 --xcol log_mass_bgp --ycol log_SFR_int --colorcol redshift --run including_nircam --use_only_good
             run make_passage_plots.py --do_field Par028 --drv 0.5 --plot_conditions SNR --line_list OIII,Ha,OII,Hb,SII --SNR_thresh 2 --xcol redshift --ycol logOH_slope --foggie_comp

             run make_passage_plots.py --do_field Par028 --drv 0.5 --SNR_thresh 3 --plot_full_BPT --colorcol sn_Hb --log_colorcol

             run make_passage_plots.py --do_field Par028 --drv 0.5 --plot_conditions SNR,mass --line_list OIII,Ha,OII,Hb,SII --SNR_thresh 2 --plot_mass_excitation --colorcol redshift
             run make_passage_plots.py --do_field Par028 --drv 0.5 --plot_conditions SNR,mass --line_list OIII,Ha,OII,Hb,SII --SNR_thresh 2 --xcol lp_mass --ycol log_SFR_int --colorcol redshift
             run make_passage_plots.py --do_field Par028 --drv 0.5 --plot_conditions SNR,mass --line_list OIII,Ha,OII,Hb,SII --SNR_thresh 2 --only_seg --vorbin --voronoi_line Ha --voronoi_snr 5 --Zdiag R3 --xcol lp_mass --ycol logOH_int --colorcol redshift
'''

from header import *
from util import *
from make_diagnostic_maps import compute_EB_V, get_dereddened_flux

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_SFMS_Whitaker14(log_mass_min, log_mass_max, redshift, nbins=40):
    '''
    Overplots fitted SFMS based on Whitaker+14 (https://iopscience.iop.org/article/10.1088/0004-637X/795/2/104/pdf) Table 1, eq 2
    Then returns a two-part log SFR array based on an input minimum and maximum log mass; the two part are based on the lower limit of this empirical relation and any extrapolation below it
    Returns three tuples: (log_mass1, log_SFR1) and (log_mass2, log_SFR2) and (z1, z2)
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

    log_mass1 = np.linspace(log_mass_min, log_mass_low_lim, nbins//2)
    log_SFR1 = np.poly1d([c, b, a])(log_mass1)

    log_mass2 = np.linspace(log_mass_low_lim, log_mass_max, nbins//2)
    log_SFR2 = np.poly1d([c, b, a])(log_mass2)

    return (log_mass1, log_SFR1), (log_mass2, log_SFR2), (z1, z2)

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_Whitaker14(ax, redshift, color='olivegreen'):
    '''
    Overplots fitted SFMS based on Whitaker+14 (https://iopscience.iop.org/article/10.1088/0004-637X/795/2/104/pdf) Table 1, eq 2
    Then overplots this on a given existing axis handle
    Returns axis handle
    '''

    (log_mass1, log_SFR1), (log_mass2, log_SFR2), (z1, z2) = get_SFMS_Whitaker14(ax.get_xlim()[0], ax.get_xlim()[1], redshift)

    ax.plot(log_mass1, log_SFR1, ls='dashed', c=color, lw=2)
    ax.plot(log_mass2, log_SFR2, ls='solid', c=color, lw=2, label=f'Whitaker+14: {z1}<z<{z2}')

    return ax

# --------------------------------------------------------------------------------------------------------------------
def get_SFMS_Shivaei15(log_mass_min, log_mass_max, nbins=40):
    '''
    Overplots fitted SFMS based on Shivaei+15 (https://iopscience.iop.org/article/10.1088/0004-637X/815/2/98/pdf) Table 1
    Then returns a two-part log SFR array based on an input minimum and maximum log mass; the two part are based on the lower limit of this empirical relation and any extrapolation below it
    Returns two tuples and a float which is the scatter in the relation: (log_mass1, log_SFR1) and (log_mass2, log_SFR2), c
    '''
    a1, b1, c1 = 0.65, -5.40, 0.40 # slope, intercept, scatter from Table 1, values for SFR(H alpha), z range 1.37-2.61
    a2, b2, c2 = 0.58, -4.65, 0.36 # slope, intercept, scatter from Table 1, values for SFR(H alpha), z range 2.09-2.61
    log_mass_low_lim = 9.5 # lower limit of mass they fitted up to

    a, b, c = a1, b1, c1 # choosing the first set of coefficients
    #a, b, c = a2, b2, c2 # choosing the second set of coefficients

    log_mass1 = np.linspace(log_mass_min, log_mass_low_lim, nbins//2)
    log_SFR1 = a * log_mass1 + b

    log_mass2 = np.linspace(log_mass_low_lim, log_mass_max, nbins//2)
    log_SFR2 = a * log_mass2 + b

    return (log_mass1, log_SFR1), (log_mass2, log_SFR2), c

# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_Shivaei15(ax, color='salmon'):
    '''
    Overplots fitted SFMS based on Shivaei+15 (https://iopscience.iop.org/article/10.1088/0004-637X/815/2/98/pdf) Table 1
    Then overplots this on a given existing axis handle
    Returns axis handle
    '''

    (log_mass1, log_SFR1), (log_mass2, log_SFR2), scatter = get_SFMS_Shivaei15(ax.get_xlim()[0], ax.get_xlim()[1])

    ax.plot(log_mass1, log_SFR1, ls='dashed', c=color, lw=2)
    ax.fill_between(log_mass1, log_SFR1 - scatter/2, log_SFR1 + scatter/2, alpha=0.3, facecolor=color)
    
    ax.plot(log_mass2, log_SFR2, ls='solid', c=color, lw=2, label=f'Shivaei+15: z~2')
    ax.fill_between(log_mass2, log_SFR2 - scatter/2, log_SFR2 + scatter/2, alpha=0.3, facecolor=color)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def get_SFMS_Popesso23(log_mass_min, log_mass_max, redshift, nbins=40):
    '''
    Computes an empirical SFMS based on Popesso+23 (https://arxiv.org/abs/2203.10487) Eq 10, for given redshift
    Then returns a two-part log SFR array based on an input minimum and maximum log mass; the two part are based on the lower limit of this empirical relation and any extrapolation below it
    Returns two tuples: (log_mass1, log_SFR1) and (log_mass2, log_SFR2)
    '''
    a0, a1, b0, b1, b2 = ufloat(0.2,0.02), ufloat(-0.034, 0.002), ufloat(-26.134,0.015), ufloat(4.722, 0.012), ufloat(-0.1925, 0.0011)  # Table 2, Eq 10
    log_mass_low_lim = 8.7 # lower limit of mass they fitted up to

    age_at_z = cosmo.age(redshift).value # Gyr

    log_mass1 = np.linspace(log_mass_min, log_mass_low_lim, nbins//2)
    log_SFR1 = (a1 * age_at_z + b1) * log_mass1 + b2 * (log_mass1) ** 2 + b0 + a0 * age_at_z

    log_mass2 = np.linspace(log_mass_low_lim, log_mass_max, nbins//2)
    log_SFR2 = (a1 * age_at_z + b1) * log_mass2 + b2 * (log_mass2) ** 2 + b0 + a0 * age_at_z

    return (log_mass1, log_SFR1), (log_mass2, log_SFR2)
   
# --------------------------------------------------------------------------------------------------------------------
def plot_SFMS_Popesso23(ax, redshift, color='cornflowerblue'):
    '''
    Computes an empirical SFMS based on Popesso+23 (https://arxiv.org/abs/2203.10487) Eq 10, for given redshift
    Then overplots this on a given existing axis handle
    Returns axis handle
    '''
    (log_mass1, log_SFR1), (log_mass2, log_SFR2) = get_SFMS_Popesso23(ax.get_xlim()[0], ax.get_xlim()[1], redshift)

    ax.plot(log_mass1, unp.nominal_values(log_SFR1), ls='dashed', c=color, lw=2)
    ax.fill_between(log_mass1, unp.nominal_values(log_SFR1) - unp.std_devs(log_SFR1)/2, unp.nominal_values(log_SFR1) + unp.std_devs(log_SFR1)/2, alpha=0.3, facecolor=color)

    ax.plot(log_mass2, unp.nominal_values(log_SFR2), ls='solid', c=color, lw=2, label=f'Popesso+23: z = {redshift}')
    ax.fill_between(log_mass2, unp.nominal_values(log_SFR2) - unp.std_devs(log_SFR2)/2, unp.nominal_values(log_SFR2) + unp.std_devs(log_SFR2)/2, alpha=0.3, facecolor=color)

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
def plot_mass_excitation(df, ax, args, mass_col='lp_mass'):
    '''
    Plots mass-excitation diagram based on integrated fluxes from grizli read in from the df provided,
    Then overplots theoretical lines
    Returns axis handle
    '''
    print(f'Plotting integrated MEx diagram..')

    df = df[(df['Hb_SNR'] > args.SNR_thresh) & (df['OIII_SNR'] > args.SNR_thresh)]

    OIII_factor = 2.98 / (1 + 2.98)
    df['OIII_int'] *= OIII_factor
    df['OIII_int_u'] *= OIII_factor

    df['O3Hb'] = np.log10(df['OIII_int'] / df['Hb_int'])

    colorcol = np.log10(df[args.colorcol]) if args.log_colorcol else df[args.colorcol]
    vmin = np.percentile(colorcol[np.isfinite(colorcol)], 5.)
    vmax = np.percentile(colorcol[np.isfinite(colorcol)], 95.)

    p = ax.scatter(df[mass_col], df['O3Hb'], c=colorcol, s=50, lw=0, cmap='cividis', vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(p)

    # ---------annotate axes and save figure-------
    if args.fontsize == 10: args.fontsize = 15
    ax.set_xlim(6.5, 11.5)
    ax.set_ylim(-1, 1.5)

    ax.set_xlabel(f'log (M*/Msun)', fontsize=args.fontsize)
    ax.set_ylabel(f'log (OIII/Hb)', fontsize=args.fontsize)

    colorby_text = f'log({args.colorcol})' if args.log_colorcol else args.colorcol
    cbar.set_label(colorby_text, fontsize=args.fontsize)

    plt.title(f'{args.do_field}', fontsize=args.fontsize)
    plt.tight_layout()

    # ---------adding literature lines from Juneau+2014 (https://iopscience.iop.org/article/10.1088/0004-637X/788/1/88/pdf)----------
    x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
    y_up = np.piecewise(x, [x <= 10, x > 10], [lambda x: (0.375 / (x - 10.5)) + 1.14, lambda x: np.poly1d([410.24, -109.333, 9.71731, -0.288244][::-1])(x)]) # J14 eq 1
    y_lo = np.piecewise(x, [x <= 9.6, x > 9.6], [lambda x: (0.375 / (x - 10.5)) + 1.14, lambda x: np.poly1d([352.066, -93.8249, 8.32651, -0.246416][::-1])(x)]) # J14 eq 2
    ax.plot(x, y_up, c='k', ls='dashed', lw=2, label='Juneau+2014')
    ax.plot(x, y_lo, c='brown', ls='dashed', lw=2)
    plt.legend()

    return ax, df

# --------------------------------------------------------------------------------------------------------------------
def plot_BPT(df, ax, args):
    '''
    Plots BPT diagram based on integrated fluxes from grizli read in from the df provided,
    Then overplots theoretical lines
    Returns axis handle
    '''
    print(f'Plotting integrated BPT diagram..')

    df = df[(df['Ha_SNR'] > args.SNR_thresh) & (df['Hb_SNR'] > args.SNR_thresh) & (df['OIII_SNR'] > args.SNR_thresh) & (df['SII_SNR'] > args.SNR_thresh)]

    OIII_factor = 2.98 / (1 + 2.98)
    df['OIII_int'] *= OIII_factor
    df['OIII_int_u'] *= OIII_factor

    Ha_factor = 0.823  # choose between 0.823, 0.754, 0.695
    df['Ha_int'] *= Ha_factor
    df['Ha_int_u'] *= Ha_factor

    df['O3Hb'] = np.log10(df['OIII_int'] / df['Hb_int'])
    df['S2Ha'] = np.log10(df['SII_int'] / df['Ha_int'])

    colorcol = np.log10(df[args.colorcol]) if args.log_colorcol else df[args.colorcol]
    vmin = np.percentile(colorcol[np.isfinite(colorcol)], 5.)
    vmax = np.percentile(colorcol[np.isfinite(colorcol)], 95.)

    p = ax.scatter(df['S2Ha'], df['O3Hb'], c=colorcol, s=50, lw=0, cmap='cividis', vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(p)

    # ---------annotate axes and save figure-------
    if args.fontsize == 10: args.fontsize = 15
    ax.set_xlim(-1.5, 0.25)
    ax.set_ylim(-1, 1.5)

    ax.set_xlabel(f'log (SII/Ha)', fontsize=args.fontsize)
    ax.set_ylabel(f'log (OIII/Hb)', fontsize=args.fontsize)

    colorby_text = f'log({args.colorcol})' if args.log_colorcol else args.colorcol
    cbar.set_label(colorby_text, fontsize=args.fontsize)

    plt.title(f'{args.do_field}', fontsize=args.fontsize)
    plt.tight_layout()

    # ---------adding literature lines from Kewley+2001 (https://iopscience.iop.org/article/10.1086/321545)----------
    x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
    y_K01 = 1.3 + 0.72 / (x - 0.32)  # Eq 6 of K01
    y_S24 = np.piecewise(x, [x >= -0.92, x < -0.92], [lambda x: (0.78 / (x - 0.34)) + 1.36, lambda x: -0.91 - 1.79 * x])
    ax.plot(x, y_K01, c='k', ls='dashed', lw=2, label='Kewley+2001')
    ax.plot(x, y_S24, c='brown', ls='dashed', lw=2, label='Schultz+2024')
    plt.legend()

    return ax, df

# --------------------------------------------------------------------------------------------------------------------
def plot_BPT_from_speccat(speccat_file, args, ax, mass_col='lp_mass'):
    '''
    Plots BPT diagram based on integrated fluxes from grizli, read in from the spec cat file
    Then overplots theoretical lines
    Returns axis handle
    '''
    tab = Table.read(speccat_file)
    columns = ['id', 'flux_Ha', 'err_Ha', 'sn_Ha', 'flux_OIII', 'err_OIII', 'sn_OIII', 'flux_SII', 'err_SII', 'sn_SII', 'flux_Hb', 'err_Hb', 'sn_Hb']
    df = tab[columns].to_pandas()

    if args.colorcol in tab.columns:
        df[args.colorcol] = tab[[args.colorcol]].to_pandas()[args.colorcol]
    else:
        photcat_file = speccat_file.parent / f'{args.do_field}_photcat.fits'
        if photcat_file.is_file():
            tab2 = Table.read(photcat_file)
            if args.colorcol in tab2.columns:
                print(f'{args.colorcol} not available in speccat, so obtaining it from photcat')
                dfp = tab2[['id', args.colorcol]].to_pandas()
                df  = pd.merge(df, dfp, on='id', how='inner')
            else:
                print(f'{args.colorcol} not available in either speccat or photcat')
                print(f'\nspeccat columns =', tab.columns)
                print(f'\nphotcat columns =', tab2.columns)

    columns_new = ['id', 'Ha_int', 'Ha_int_u', 'Ha_SNR', 'OIII_int', 'OIII_int_u', 'OIII_SNR', 'SII_int', 'SII_int_u', 'SII_SNR', 'Hb_int', 'Hb_int_u', 'Hb_SNR']
    df = df.rename(columns={columns[index]:columns_new[index] for index in range(len(columns))})
    if args.colorcol in columns: args.colorcol = dict(zip(columns, columns_new))[args.colorcol]

    if args.plot_full_BPT: ax, df = plot_BPT(df, ax, args)
    elif args.plot_full_mass_excitation: ax, df = plot_mass_excitation(df, ax, args, mass_col=mass_col)

    return df, ax

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

# --------------------------------------------------------------------------------------------------------------------
def plot_fluxcomp(field, colx, coly_arr, ylabel, args, photcat_file=None):
    '''
    Plots fluxes in different apertures for comparison
    Returns figure handle
    '''
    def label_dict(col):
        dict = {'aper_3':'0.5"', 'aper_5':'1.0"', 'aper_7':'2.0"'}
        text = [item for item in list(dict.keys()) if item in col]
        if len(text) == 0: return col
        else: return dict[text[0]]

    if photcat_file is None: photcat_file = args.input_dir / f'v0.5/data/{field}/DATA/DIRECT_GRISM/{field}_photcat.fits'
    df = Table.read(photcat_file).to_pandas()

    fig, ax = plt.subplots(1)
    for coly in coly_arr: ax.scatter(np.log10(df[colx]), np.log10(df[coly]), label=label_dict(coly), s=5)
    ax.plot([-4 ,2], [-4, 2], ls='--', c='k')
    ax.legend()
    ax.set_ylim(-4, 2)
    ax.set_xlim(-4, 2)
    ax.set_xlabel(colx)
    ax.set_ylabel(ylabel)
    plt.title(field)

    fig.savefig(photcat_file.parent / f'{field}_fluxcomp.png')
    return fig

# --------------------------------------------------------------------------------------------------------------------
def break_column_into_uncertainty(df, col, make_log=False):
    '''
    If a given column in a dataframe is of type ufloat (with uncertainy) then separate it into different columsn for values and errors
    and also make corresponding log columns
    Returns the modified dataframe
    '''

    if df[col].dtype == object:  # accompanied by uncertainty in the same column
        counter = 0
        for index, item in enumerate(df[col]):
            if type(item) != str and np.isnan(item):
                quant = [np.nan, np.nan]
            elif 'e' in item:
                pow = float(item[item.find('e') + 1:])
                base = np.array(item[1:item.find('e') - 1].split('+/-')).astype(np.float64)
                quant = base * 10 ** pow
            else:
                quant = np.array(item.split('+/-')).astype(np.float64)
            df.at[index, col] = quant[0]
            df.at[index, col + '_u'] = quant[1]
            if make_log:
                try:
                    log_quant = unp.log10(ufloat(quant[0], quant[1]))
                    df.at[index, 'log_' + col] = unp.nominal_values(log_quant)
                    df.at[index, 'log_' + col + '_u'] = unp.std_devs(log_quant)
                except ValueError:
                    df.at[index, 'log_' + col] = np.nan
                    df.at[index, 'log_' + col + '_u'] = np.nan
                    counter += 1
                    pass
        if counter > 0: print(f'\n{counter} out of {len(df)} had -ve {col}!')
    elif make_log:
        if col + '_u' not in df: df[col + '_u'] = 0.
        quant = unp.uarray(df[col], df[col + '_u'])
        log_quant = unp.log10(quant)
        df['log_' + col] = unp.nominal_values(log_quant)
        df['log_' + col + '_u'] = unp.std_devs(log_quant)

    return df

# -------------------------------------global dictionaries-------------------------------------------------------------------------------
label_dict = {'lp_mass': r'log M$_*$/M$_{\odot}$ (LePhare)', 'ez_mass': r'log M$_*$/M$_{\odot}$ (EAZY)', 'log_mass_bgp': r'log M$_*$/M$_{\odot}$ (Bagpipes)', \
              'lp_SFR': r'log SFR (M$_{\odot}$/yr) (LePhare)', 'ez_sfr': r'log SFR (M$_{\odot}$/yr) (EAZY)', 'log_sfr_bgp': r'log SFR (M$_{\odot}$/yr) (Bagpipes)', 'log_SFR_int': r'log SFR (M$_{\odot}$/yr) (Grizli)', \
              'lp_zBEST': 'Redshift (LePhare)', 'ez_z_phot': 'Redshift (EAZY)', 'z_bgp': 'Redshift (Bagpipes)', 'redshift': 'Redshift (Grizli)', \
              'logOH_slope':r'log $\nabla$Z$_r$ (dex/kpc)'}
bounds_dict = {'lp_mass': (6.5, 11), 'ez_mass': (6, 9), 'log_mass_bgp': (6.5, 10.5), \
               'lp_SFR': (-3, 1), 'ez_sfr': (-3, 1), 'log_sfr_bgp': (-3, 2), 'log_SFR_int': (-3, 2.5), \
               'ez_z_phot': (0, 3), 'lp_zBEST': (0, 3), 'z_bgp': (0, 3), 'redshift': (1.7, 3.1), \
               'logOH_slope': (-0.5, 0.5)}
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

    # ---------BPT for for full sample------
    speccat_file = args.input_dir / f'{args.drv}/{args.do_field}/Products/{args.do_field}_speccat.fits'
    if args.plot_full_BPT:
        if args.colorcol == 'ez_z_phot': args.colorcol = 'redshift'
        df, ax = plot_BPT_from_speccat(speccat_file, args, ax)
        colorby_text = f'log_{args.colorcol}' if args.log_colorcol else args.colorcol
        figname = speccat_file.parent / f'{args.do_field}_BPT_colorby_{colorby_text}.png'

    # ---------mass-excitation diagram for for full sample------
    elif args.plot_full_mass_excitation:
        if args.colorcol == 'ez_z_phot': args.colorcol = 'redshift'
        df, ax = plot_BPT_from_speccat(speccat_file, args, ax)
        colorby_text = f'log_{args.colorcol}' if args.log_colorcol else args.colorcol
        figname = speccat_file.parent / f'{args.do_field}_mass_excitation_colorby_{colorby_text}.png'

    # ---------flux vs mag for full sample------
    elif args.plot_flux_vs_mag:
        figname = args.output_dir / 'plots' / f'allpar_flux_vs_mag.png'
        ax = plot_flux_vs_mag(ax, args)

    else:
        # -------reading in dataframe produced by get_field_stats.py or by compute_stellar_masses.py----------------
        if args.do_field is None:
            plot_conditions_text = ','.join(args.line_list) + ',' + ','.join(args.plot_conditions)
            plot_conditions_text = plot_conditions_text.replace('SNR', f'SNR>{args.SNR_thresh}').replace('EW', f'EW>{args.EW_thresh}').replace('a_image', f'a>{args.a_thresh}')
            args.field_set_plot_conditions_text = f'allpar_{args.drv}_venn_{plot_conditions_text}'
            if len(args.run) > 0 and args.run[0] != '_': args.run = '_' + args.run
            df_infilename = args.output_dir / 'catalogs' / f'{args.field_set_plot_conditions_text}_df_withSED{args.run}.csv'
        elif args.plot_conditions == 'all_match':
            args.field = f'Par{int(args.do_field.split("Par")[1]):03d}'
            args.field_set_plot_conditions_text = f'{args.field}_{args.drv}_allmatch'
            df_infilename = args.output_dir / args.field / f'{args.field}_all_diag_results_withSED_{args.run}.csv'
        else:
            args.field = f'Par{int(args.do_field.split("Par")[1]):03d}'
            plot_conditions_text = ','.join(args.line_list) + ',' + ','.join(args.plot_conditions)
            plot_conditions_text = plot_conditions_text.replace('SNR', f'SNR>{args.SNR_thresh}').replace('EW', f'EW>{args.EW_thresh}').replace('a_image', f'a>{args.a_thresh}')
            args.field_set_plot_conditions_text = f'{args.field}_{args.drv}_venn_{plot_conditions_text}'
            if len(args.run) > 0 and args.run[0] != '_': args.run = '_' + args.run
            df_infilename = args.output_dir / 'catalogs' / f'{args.field_set_plot_conditions_text}_df_withSED{args.run}.csv'

        if os.path.exists(df_infilename):
            print(f'Reading in main df from {df_infilename}')
        else:
            print(f'Could not find {df_infilename},')
            df_infilename = Path(str(df_infilename)[:str(df_infilename).find('withSED') - 1] + '.txt')
            if os.path.exists(df_infilename):
                print(f'Loading pre-SED df from {df_infilename} instead')
            else:
                print(f'Could not find {df_infilename},')
                df_infilename = Path(str(df_infilename).replace(f'_{args.drv}', ''))
                if os.path.exists(df_infilename):
                    print(f'Loading df from {df_infilename} instead')
                else:
                    sys.exit(f'{df_infilename} does not exist')

        df = pd.read_csv(df_infilename)
        if args.use_only_good and args.drv == 'v0.5' and 'SNR' in args.plot_conditions:
            if set(args.line_list) == set(['OIII', 'Ha', 'OII', 'Hb', 'SII']): df = df[df['objid'].isin([1303,1934,2734,2867,300,2903])].reset_index(drop=True) # only choosing the pre-determined good galaxies
            elif set(args.line_list) == set(['OIII', 'OII', 'Hb', 'NeIII-3867']): df = df[df['objid'].isin([300,1303,2171,2727,2867])].reset_index(drop=True) # only choosing the pre-determined good galaxies
            print(f'\nUsing only the pre-determined good galaxies, and there are {len(df)} of them..')

        # -------combing with metallicity dataframe if it exists----------------
        snr_text = f'_snr{args.snr_cut}' if args.snr_cut is not None else ''
        only_seg_text = '_onlyseg' if args.only_seg else ''
        vorbin_text = '' if not args.vorbin else f'_vorbin_at_{args.voronoi_line}_SNR_{args.voronoi_snr}'
        logOHgrad_filename = args.output_dir / 'catalogs' / f'logOHgrad_df{snr_text}{only_seg_text}{vorbin_text}.txt'

        if os.path.exists(logOHgrad_filename):
            print(f'Reading in and merging logOH gradient df: {logOHgrad_filename}')
            df_logOHgrad = pd.read_csv(logOHgrad_filename)
            df_logOHgrad = df_logOHgrad.drop_duplicates(subset=['field', 'objid', 'logOH_diagnostic'], keep='last')
            df_logOHgrad = df_logOHgrad[(df_logOHgrad['logOH_diagnostic'] == args.Zdiag) & (df_logOHgrad['logOH_branch'] == args.Zbranch)]
            df = pd.merge(df, df_logOHgrad, on=['field', 'objid'], how='outer')

        # -------making the mass excitation plot----------------
        if args.plot_mass_excitation:
            colorby_text = f'_colorby_z' if args.colorcol == 'ez_z_phot' else f'_colorby_{args.colorcol}'
            figname = args.output_dir / 'plots' / f'{args.field_set_plot_conditions_text}_run_{args.run}_mass_excitation{colorby_text}.png'
            ax, df = plot_mass_excitation(df, ax, args, mass_col='lp_mass')

        # -------making the BPT plot----------------
        elif args.plot_BPT:
            colorby_text = f'_colorby_z' if args.colorcol == 'ez_z_phot' else f'_colorby_{args.colorcol}'
            figname = args.output_dir / 'plots' / f'{args.field_set_plot_conditions_text}_run_{args.run}_BPT{colorby_text}.png'
            ax, df = plot_BPT(df, ax, args)

            # ------writing out distance from K01 AGN-SF line--------------------------
            if ('distance' in args.colorcol and 'K01' in args.colorcol) or args.colorcol == 'distance_from_K01':
                df.to_csv(df_infilename, index=None)
                print(f'\nAdded distance_from_K01 column to df and saved in {df_infilename}.')

        # -------making the desired plots----------------
        else:
            figname = args.output_dir / 'plots' / f'{args.field_set_plot_conditions_text}_run_{args.run}_df_{args.xcol}_vs_{args.ycol}_colorby_{args.colorcol}.png'
            df = break_column_into_uncertainty(df, 'SFR_int', make_log=True)

            for thiscol in [item for item in df.columns if 'logOH' in item and 'int' in item]:
                df = break_column_into_uncertainty(df, thiscol, make_log=False)

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
                ax = plot_SFMS_Popesso23(ax, 0.5, color='darkblue')
                ax = plot_SFMS_Popesso23(ax, 1.2, color='blue')
                ax = plot_SFMS_Popesso23(ax, 2.0, color='cornflowerblue')
                ax = plot_SFMS_Shivaei15(ax, color='salmon')
                ax = plot_SFMS_Whitaker14(ax, 0.6, color='forestgreen')
                ax = plot_SFMS_Whitaker14(ax, 1.2, color='limegreen')
                ax = plot_SFMS_Whitaker14(ax, 1.8, color='yellowgreen')

            # ---------MZR from literature-------
            if 'mass' in args.xcol and 'logOH' in args.ycol and 'slope' not in args.ycol:
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

            if args.ycol == 'logOH_slope':
                ax.axhline(0, ls='--', c='k', lw=0.5)
                if args.xcol == 'redshift' and args.foggie_comp: # special case, to match the FOGGIE plot limits
                    ax.set_xlim(4, 0.5)
                    ax.set_ylim(-0.5, 0.4)

            if 'mass' in args.xcol and 'logOH' in args.ycol and 'slope' not in args.ycol:
                ax.set_xlim(6.5, 12.0)
                ax.set_ylim(6.0, 10.0)

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
