'''
    Filename: make_diagnostic_maps.py
    Notes: Plots integrated 1D spectra and 2D emission line maps (from existing .full.fits file), for a given object/s in a given field
    Author : Ayan
    Created: 17-07-24
    Example: run make_diagnostic_maps.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par50 --id 3667
             run make_diagnostic_maps.py --Par 50 --id 823
'''

from header import *
from util import *
from matplotlib import cm as mpl_cm

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def plot_direct_image(full_hdu, ax, args):
    '''
    Plots the combined direct image of all available filters in the given axis
    Returns the axis handle
    '''
    cmap_arr = ['Reds', 'Blues', 'Greens']

    # -------plot direct image for each filter---------
    for index in range(args.ndfilt):
        filt = full_hdu[0].header[f'DFILT{index:02d}']
        print(f'Plotting direct image for filter {filt} which is {index+1} of {args.ndfilt}..')

        ext = 5 + index * 2
        image = full_hdu[ext].data

        ax.imshow(image, cmap=cmap_arr[index], origin='lower', extent=args.extent, alpha=0.7)

        textcolor = mpl_cm.get_cmap(cmap_arr[index])(0.9)
        ax.text(0.1, 0.9 - index * 0.1, filt, c=textcolor, fontsize=args.fontsize/2, ha='left', va='top')

    ax.set_xlabel('RA', fontsize=args.fontsize)
    ax.set_ylabel('Dec', fontsize=args.fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_1d_spectra(od_hdu, ax, args):
    '''
    Plots the 1D spectra in the given axis
    Returns the axis handle
    '''
    nfilters = sum(['GRISM' in item for item in list(od_hdu[0].header.keys())])
    filters = [od_hdu[0].header[f'GRISM{item + 1:03d}'] for item in range(nfilters)]

    col_arr = ['gray', 'black', 'salmon'] # colors in order: for measured flux, fitted continuum, fitted continuum + line flux

    # -------plot 1D spectra for each filter-----------
    for filter in filters:
        table = Table(od_hdu[filter].data)
        table['rest_wave'] = table['wave'] / (1 + args.z)
        ax.plot(table['rest_wave'], table['flux'] / table['flat'], lw=0.5, c=col_arr[0], alpha=0.5) # need to divide all columns with 'flat' to get the right units (ergs/s/cm^2/A)
        ax.plot(table['rest_wave'], table['cont'] / table['flat'], lw=0.5, c=col_arr[1])
        ax.plot(table['rest_wave'], table['line'] / table['flat'], lw=1, c=col_arr[2])

    ax.set_xlabel(r'Rest-frame wavelength ($\AA$)', fontsize=args.fontsize)
    ax.set_ylabel(r'f$_{\lambda}$ (ergs/s/cm$^2$/A)', fontsize=args.fontsize)

    # ---observed wavelength axis-------
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(ax.get_xticks())
    ax2.set_xticklabels(['%.1F' % item * (1 + args.z) / 1e4 for item in ax2.get_xticks()], fontsize=args.fontsize)
    ax2.set_xlabel(r'Observed wavelength ($\mu$)', fontsize=args.fontsize)

    # ---vertical lines for emission line wavelengths------
    lines_df = pd.read_csv(HOME / 'Desktop/Lisa_UV_diag/P_spherical/sp_P70_a05modelfiles/Q700/ spec0003.csv', skiprows=55, names=['wave', 'ev', 'flux', 'species', 'kind', 'acc'])
    lines_df = lines_df[lines_df['wave'].between(table['rest_wave'].min(), table['rest_wave'].max())]
    lines_df = lines_df[lines_df['flux'] > 0.003]
    lines_df = lines_df[lines_df['kind'].isin(['CM', 'RCAB'])]

    max_flux = lines_df['flux'].max() * 1.01
    min_flux = lines_df['flux'].min() * 0.09

    for index in range(len(lines_df)):
        ax.axvline(lines_df.loc[index]['wave'], c='k', alpha=(lines_df.lox[index]['flux'] - min_flux) / (max_flux - min_flux))
        ax.text(lines_df.loc[index]['wave'] * 1.01, ax.get_ylim()[1] * 0.9, lines_df.lox[index]['species'], rotation=90, fontsize=args.fontsize / 2)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_2D_map(image, ax, args, label=None, cmap=None, vmin=None, vmax=None):
    '''
    Plots the emission map for a given line in the given axes
    Returns the axes handle
    '''
    if cmap is None: cmap = 'cividis'

    ax.imshow(image, cmap=cmap, origin='lower', extent=args.extent, vmin=vmin, vmax=vmax)

    ax.text(0.1, 0.9, label, c='k', fontsize=args.fontsize/2, ha='left', va='top')
    ax.set_xlabel('RA', fontsize=args.fontsize)
    ax.set_ylabel('Dec', fontsize=args.fontsize)

    return ax

# --------------------------------------------------------------------------------------------------------------------
def plot_emission_line_map(line, full_hdu, ax, args):
    '''
    Plots the emission map for a given line in the given axes
    Returns the axes handle
    '''

    line_index = np.where(np.array(args.available_lines.split(' ')) == line)[0][0]
    ext = 5 + 2 * args.ndfilt + 4 * line_index
    line_map = full_hdu[ext].data

    ax = plot_2D_map(line_map, ax, args, label=line, cmap=args.emline_cmap, vmin=None, vmax=None)

    return line_map, ax

# --------------------------------------------------------------------------------------------------------------------
def plot_sfr_map(full_hdu, axes, args):
    '''
    Plots the SFR map (and the emission line maps that go into it) in the given axes
    Returns the axes handles
    '''

    line_map, axes[0] = plot_emission_line_map('Ha', full_hdu, axes[0], args)

    sfr_map = line_map * 
    axes[1] = plot_2D_map(sfr_map, axes[1], args, label='SFR', cmap='Blues', vmin=None, vmax=None)

    return axes

# --------------------------------------------------------------------------------------------------------------------
def plot_dust_map(full_hdu, axes, args):
    '''
    Plots the dust extinction map (and the emission line maps that go into it) in the given axes
    Returns the axes handles
    '''

    return axes

# --------------------------------------------------------------------------------------------------------------------
def plot_Z_R23_map(full_hdu, axes, args):
    '''
    Plots the R23 metallicity map (and the emission line maps that go into it) in the given axes
    Returns the axes handles
    '''

    return axes

# --------------------------------------------------------------------------------------------------------------------
def plot_Z_Te_map(full_hdu, axes, args):
    '''
    Plots the Te metallicity map (and the emission line maps that go into it) in the given axes
    Returns the axes handles
    '''

    return axes

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()

    # ------------looping over the provided object IDs-----------------------
    for index, this_id in enumerate(args.id):
        start_time2 = datetime.now()
        print(f'\nCommencing ID {this_id} which is {index+1} of {len(args.id)}..')

        # ------determining directories---------
        extract_dir = args.input_dir / args.field / 'Extractions'
        output_subdir = args.output_dir / args.field / f'{this_id:05d}'
        output_subdir.mkdir(parents=True, exist_ok=True)
        full_fits_file = f'{args.field}_{this_id:05d}.full.fits'

        if os.path.exists(extract_dir / full_fits_file): # if the fits files are in Extractions/
            args.work_dir = extract_dir
        elif os.path.exists(output_subdir / full_fits_file): # if the fits files are in sub-directories for individual objects
            args.work_dir = output_subdir
        os.chdir(args.work_dir)

        # ------------read in fits files--------------------------------
        od_hdu = fits.open(args.work_dir / f'{args.field}_{this_id:05d}.1D.fits')
        full_hdu = fits.open(args.work_dir / f'{args.field}_{this_id:05d}.full.fits')

        # ----------determining global parameters------------
        args.available_lines = full_hdu[0].header['HASLINES']
        args.z = full_hdu[0].header['REDSHIFT']
        args.ndfilt = full_hdu[0].header['NDFILT']
        args.nlines = full_hdu[0].header['NUMLINES']
        args.emline_cmap = 'cividis'

        line_wcs = pywcs.WCS(full_hdu['DSCI'].header)
        pix_size = utils.get_wcs_pscale(line_wcs)
        imsize_arcsec = full_hdu['DSCI'].data.shape[0] * pix_size
        dp = 0 # -0.5 * pix_size  # FITS reference is center of a pixel, array is edge
        args.extent = (-imsize_arcsec / 2. - dp, imsize_arcsec / 2. - dp, -imsize_arcsec / 2. - dp, imsize_arcsec / 2. - dp)

        # ---------initialising the figure------------------------------
        nrow, ncol = 5, 3
        fig = plt.figure(figsize=(10,8))
        axis_dirimg = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 0), colspan=1)
        axis_1dspec = plt.subplot2grid(shape=(nrow, ncol), loc=(0, 1), colspan=ncol - 1)
        axes_sfr = [plt.subplot2grid(shape=(nrow, ncol), loc=(1, item), colspan=2) for item in np.arange(2)]
        axes_dust = [plt.subplot2grid(shape=(nrow, ncol), loc=(2, item), colspan=2) for item in np.arange(2)]
        axes_Z_R23 = [plt.subplot2grid(shape=(nrow, ncol), loc=(3, item), colspan=3) for item in np.arange(3)]
        axes_Z_Te = [plt.subplot2grid(shape=(nrow, ncol), loc=(4, item), colspan=2) for item in np.arange(2)]
        fig.tight_layout()
        fig.subplots_adjust(top=0.98, bottom=0.07, left=0.1, right=0.87, wspace=2.0, hspace=0.35)

        # ---------populating the figure------------------------------
        axis_dirimg = plot_direct_image(full_hdu, axis_dirimg, args)
        axis_1dspec = plot_1d_spectra(od_hdu, axis_1dspec, args)
        if 'Ha' in args.available_lines: axes_sfr = plot_sfr_map(full_hdu, axes_sfr, args)
        if all([line in args.available_lines for line in ['Ha', 'Hb']]): axes_dust = plot_dust_map(full_hdu, axes_dust, args)
        if all([line in args.available_lines for line in ['OIII-4363', 'OII', 'Hb']]): axes_Z_R23 = plot_Z_R23_map(full_hdu, axes_Z_R23, args)
        if 'OIII-4363' in args.available_lines: axes_Z_Te = plot_Z_Te_map(full_hdu, axes_Z_Te, args)

        # ---------decorating the figure------------------------------


        print(f'Completed id {this_id} in {timedelta(seconds=(datetime.now() - start_time2).seconds)}, {len(args.id) - index - 1} to go!')

    os.chdir(args.code_dir)
    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
