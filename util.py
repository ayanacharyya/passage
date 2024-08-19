##!/usr/bin/env python3

"""

    Filename :   util.py
    Notes :      Contains various generic utility functions and classes used by the other scripts in PASSAGE, including a function to parse args
    Author :    Ayan
    Created: 11-06-24
    Last modified: 11-06-24

"""

from header import *

# --------------------------------------------------------------------------------------------------------------------
def parse_args():
    '''
    Parse command line arguments. Returns args object
    '''

    parser = argparse.ArgumentParser(description='Produces emission line maps for JWST-PASSAGE data.')

    # ---- common args used widely over the full codebase ------------
    parser.add_argument('--input_dir', metavar='input_dir', type=str, action='store', default=None, help='Where do the input files reside?')
    parser.add_argument('--output_dir', metavar='output_dir', type=str, action='store', default=None, help='Where do you want to store the outputs?')
    parser.add_argument('--system', metavar='system', type=str, action='store', default='hd', help='Which file system is the code being run on?')
    parser.add_argument('--code_dir', metavar='code_dir', type=str, action='store', default='/Users/acharyya/Work/astro/ayan_codes/passage/', help='Where is the source code?')
    parser.add_argument('--clobber', dest='clobber', action='store_true', default=False, help='Over-write existing plots? Default is no.')
    parser.add_argument('--silent', dest='silent', action='store_true', default=False, help='Suppress some generic print statements? Default is no.')
    parser.add_argument('--keep', dest='keep', action='store_true', default=False, help='Keep existing plot windows open? Default is no.')

    parser.add_argument('--field', metavar='field', type=str, action='store', default='Par008', help='Which passage field? Default is Par50')
    parser.add_argument('--id', metavar='id', type=str, action='store', default='100', help='Object ID. Default is 100')

    # ------- args added for make_spectra_from_beam.py ------------------------------
    parser.add_argument('--include_photometry', dest='include_photometry', action='store_true', default=False, help='Include photometry while computing fit parameters? Default is no.')
    parser.add_argument('--zmin', metavar='zmin', type=float, action='store', default=0.5, help='minimum of redshift range within which to search for lines; default is 0.5')
    parser.add_argument('--zmax', metavar='zmax', type=float, action='store', default=1.0, help='maximum of redshift range within which to search for lines; default is None')
    parser.add_argument('--line_list', metavar='line_list', type=str, action='store', default='all', help='Which emission lines to look for? Default is all') # OR set default to 'Lya,OII,Hb,OIII,Ha,Ha+NII,SII,SIII,PaB,He-1083,PaA'
    parser.add_argument('--pixscale', metavar='pixscale', type=float, action='store', default=0.04, help='Pixel scale (in arcsec/pixel) of the thumbnails produced; default is 0.04')

    # ------- args added for read_line_catalog.py ------------------------------
    parser.add_argument('--nbins', metavar='nbins', type=int, action='store', default=30, help='No. of bins for plotting the histogram. Default is 30')
    parser.add_argument('--fontsize', metavar='fontsize', type=int, action='store', default=10, help='fontsize of plot labels, etc.; default is 15')
    parser.add_argument('--mag_lim', metavar='mag_lim', type=float, action='store', default=None, help='magnitude limit above which to search for targets; default is None')

    # ------- args added for plot_footprints.py ------------------------------
    parser.add_argument('--bg_file', metavar='bg_file', type=str, action='store', default=None, help='Which file to be used for plotting the background image?')
    parser.add_argument('--plot_zcosmos', dest='plot_zcosmos', action='store_true', default=False, help='Overplot the (thousands of) zCOSMOS targets? Default is no.')

    # ------- args added for make_spectra_from_beams.py ------------------------------
    parser.add_argument('--skip_sep', dest='skip_sep', action='store_true', default=False, help='Skip the Source Extraction Pipeline? Default is no.')
    parser.add_argument('--do_all_obj', dest='do_all_obj', action='store_true', default=False, help='Reduce spectra and make beam files for ALL detected objects? Default is no.')
    parser.add_argument('--re_extract', dest='re_extract', action='store_true', default=False, help='Re-extract ALL objects to be re-extracted? Default is no.')

    # ------- args added for run_passagepipe.py ------------------------------
    parser.add_argument('--dry_run', dest='dry_run', action='store_true', default=False, help='Do a dry run i.e., only make the config file without actually running PASSAGEPipe? Default is no.')
    parser.add_argument('--start_step', metavar='start_step', type=int, action='store', default=1, help='Starting step for PASSAGEPipe. Default is first step')
    parser.add_argument('--do_download', dest='do_download', action='store_true', default=False, help='Do stage 1 of PASSAGEPipe? Default is no.')
    parser.add_argument('--do_prep', dest='do_prep', action='store_true', default=False, help='Do stage 2 of PASSAGEPipe? Default is no.')
    parser.add_argument('--do_image', dest='do_image', action='store_true', default=False, help='Do stage 3 of PASSAGEPipe? Default is no.')
    parser.add_argument('--do_grism', dest='do_grism', action='store_true', default=False, help='Do stage 4 of PASSAGEPipe? Default is no.')
    parser.add_argument('--do_extract', dest='do_extract', action='store_true', default=False, help='Do stage 5 of PASSAGEPipe? Default is no.')
    parser.add_argument('--do_post', dest='do_post', action='store_true', default=False, help='Do stage 6 of PASSAGEPipe? Default is no.')
    parser.add_argument('--do_upload', dest='do_upload', action='store_true', default=False, help='Do stage 7 of PASSAGEPipe? Default is no.')
    parser.add_argument('--download_from_mast', dest='download_from_mast', action='store_true', default=False, help='Download from MAST? Default is no.')
    parser.add_argument('--clobber_download', dest='clobber_download', action='store_true', default=False, help='Over-write existing files during downloading from MAST? Default is no.')
    parser.add_argument('--redo_level1', dest='redo_level1', action='store_true', default=False, help='Re-do Level1 of processing during stage 1 of PASSAGEPipe? Default is no.')
    parser.add_argument('--filters', metavar='filters', type=str, action='store', default=None, help='Which filters are included for this field? Default is None')
    parser.add_argument('--magmin', metavar='magmin', type=float, action='store', default=16, help='magnitude lower limit for refined magnitude search during PASSAGEPipe; default is 16')
    parser.add_argument('--start_id', metavar='start_id', type=int, action='store', default=0, help='Starting ID of the object whose spectra is to be extracted. Default is 0')
    parser.add_argument('--stop_id', metavar='stop_id', type=int, action='store', default=10000, help='Stopping ID of the object whose spectra is to be extracted. Default is all IDs till the end of the list')

    # ------- args added for make_diagnostic_maps.py ------------------------------
    parser.add_argument('--plot_target_frame', dest='plot_target_frame', action='store_true', default=False, help='Annotate plot axes in the object/target frame of reference? Default is no.')
    parser.add_argument('--arcsec_limit', metavar='arcsec_limit', type=float, action='store', default=1.0, help='Half box size (in arcsec) of the thumbnails to plot; default is 1.5')
    parser.add_argument('--vorbin', dest='vorbin', action='store_true', default=False, help='Voronoi bin the 2D emission line maps? Default is no.')
    parser.add_argument('--voronoi_snr', metavar='voronoi_snr', type=float, action='store', default=3, help='Target SNR to Voronoi bin the emission line maps to; default is 3')
    parser.add_argument('--voronoi_line', metavar='voronoi_line', type=str, action='store', default='Ha', help='Which emission line to be used for computing the Voronoi bins? Default is None i.e., the given emission line itself')
    parser.add_argument('--flam_max', metavar='flam_max', type=float, action='store', default=None, help='Maximum y-axis limit for f_lambda (in units of 1e-19 ergs/s/cm^2/A); default is None')
    parser.add_argument('--plot_radial_profiles', dest='plot_radial_profiles', action='store_true', default=False, help='Plot radial profiles corresponding to the 2D maps? Default is no.')
    parser.add_argument('--snr_cut', metavar='snr_cut', type=float, action='store', default=None, help='Impose an SNR cut on the emission line maps to; default is None')
    parser.add_argument('--only_seg', dest='only_seg', action='store_true', default=False, help='Cut out the emission line plots corresponding to the grizli segmentation map? Default is no.')
    parser.add_argument('--write_file', dest='write_file', action='store_true', default=False, help='Write the measured quantities to a master dataframe? Default is no.')
    parser.add_argument('--plot_mappings', dest='plot_mappings', action='store_true', default=False, help='Plot emission line locations as per MAPPINGS predictions (will lead to crowding of many lines)? Default is no.')
    parser.add_argument('--hide', dest='hide', action='store_true', default=False, help='Hide (do not display) the plots just made? Default is no.')

    # ------- args added for get_field_stats.py ------------------------------
    parser.add_argument('--EW_thresh', metavar='EW_thresh', type=float, action='store', default=300, help='EW threshold to consider good detection for emission line maps; default is 300')
    parser.add_argument('--do_all_fields', dest='do_all_fields', action='store_true', default=False, help='Include ALL available fields? Default is no.')
    parser.add_argument('--plot_venn', dest='plot_venn', action='store_true', default=False, help='Plot Venn diagrams? Default is no.')
    parser.add_argument('--merge_visual', dest='merge_visual', action='store_true', default=False, help='Include visually inspected dataframe for Venn diagrams? Default is no.')
    parser.add_argument('--plot_conditions', metavar='plot_conditions', type=str, action='store', default='detected', help='Which conditions are plotted in the Venn diagram? Default is None')

    # ------- wrap up and processing args ------------------------------
    args = parser.parse_args()
    if args.line_list != 'all': args.line_list = [item for item in args.line_list.split(',')]
    if 'Par' in args.field: args.field = f'Par{int(args.field.split("Par")[1]):03d}'
    args.id = [int(item) for item in args.id.split(',')]

    if args.system == 'hd' and not os.path.exists('/Volumes/Elements/'): args.system = 'local'
    if args.line_list == 'all': args.line_list = ['Lya', 'OII', 'Hb', 'OIII-4363', 'OIII', 'Ha', 'NII','Ha+NII', 'SII', 'SIII', 'PaD','PaG','PaB','HeI-1083','PaA']

    root_dir = '/Users/acharyya/Work/astro/passage' if args.system == 'local' else '/Volumes/Elements/acharyya_backup/Work/astro/passage' if 'hd' in args.system else '/Users/acharyya/Library/CloudStorage/GoogleDrive-ayan.acharyya@inaf.it/My Drive/passage' if 'gdrive' in args.system else ''
    args.root_dir = Path(root_dir)

    if args.input_dir is None:
        args.input_dir = args.root_dir / 'passage_data/'
    if args.output_dir is None:
        args.output_dir = args.root_dir / 'passage_output/'

    args.input_dir = Path(args.input_dir)
    args.output_dir = Path(args.output_dir)
    args.code_dir = Path(args.code_dir)

    if args.filters is None:
        if args.field in filter_dict: args.filters = filter_dict[args.field]
        else: args.filters = ['F115W'] # default
    else:
        args.filters = args.filters.split(',')

    args.plot_conditions = args.plot_conditions.split(',')

    return args

# -----------------------------------------------------------------------------
def fix_filter_names(filters):
    '''
    Fixes filter nomenclature, if needed, and sorts them by wavelength
    '''
    filters = np.atleast_1d(filters)
    for index, filter in enumerate(filters):
        if filter[0] != 'F': filter = 'F' + filter
        if filter[-1] not in ['N', 'M', 'W']: filter = filter + 'W'  # by default assume wide band filter
        filters[index] = filter
    filters = np.sort(filters)

    return filters

# ------------------------------------------------------------------------------------
def get_zranges_for_filters(line, filters=['F115W', 'F150W', 'F200W']):
    '''
    Computes the redshift range in which a certain emission line is available, for a given JWST filter,
    considering the lambda boundaries where sensitivity drops to 50% (by default corresponds to the set of filters used in PASSAGE)
    Input filter name/s must be from the known JWST filters
    NIRISS filter data taken from Table 1 in https://jwst-docs.stsci.edu/jwst-near-infrared-imager-and-slitless-spectrograph/niriss-instrumentation/niriss-filters#gsc.tab=0
    Returns min and max redshift/s
    '''
    filter_dict = {'F090W':[0.796, 1.005], 'F115W':[1.013, 1.283], 'F150W':[1.330, 1.671], 'F200W':[1.751, 2.226], 'F277W':[2.413, 3.143], 'F356W':[3.140, 4.068], 'F444W':[3.880, 5.023], \
                   'F140M':[1.331, 1.480], 'F158M':[1.488, 1.688], 'F380M':[3.726, 3.931], 'F430M':[4.182, 4.395], 'F480M':[4.668, 4.971]} # microns
    filters = fix_filter_names(filters)

    z_limits = []
    for index, filter in enumerate(filters):
        obs_wave_range = np.array(filter_dict[filter])
        z_min, z_max = get_zrange_for_line(line, obs_wave_range=obs_wave_range * 1e3)  # x1e3 to convert um to nm
        if index == 0:
            z_limits.extend([z_min, z_max])
        else:
            if z_min <= z_limits[-1]: z_limits[-1] = z_max
            else: z_limits.extend([z_min, z_max])

    z_limits = np.array(z_limits).reshape((int(len(z_limits)/2), 2))
    return z_limits, filters

# ------------------------------------------------------------------------------------
def print_zranges_for_filters(lines, filters=['F115W', 'F150W', 'F200W']):
    '''
    Prints the redshift range in which a certain emission line is available, for a given JWST filter,
    considering the lambda boundaries where sensitivity drops to 50% (by default corresponds to the set of filters used in PASSAGE)
    Input filter name/s must be from the known JWST filters
    NIRISS filter data taken from Table 1 in https://jwst-docs.stsci.edu/jwst-near-infrared-imager-and-slitless-spectrograph/niriss-instrumentation/niriss-filters#gsc.tab=0
    Prints min and max redshift/s
    '''
    lines = np.atleast_1d(lines)
    filters = fix_filter_names(filters)
    print(f'For given filters {filters}..')

    for line in lines:
        z_limits, filters = get_zranges_for_filters(line, filters=filters)
        print(f'{line} is captured between' + ' '.join(['z=[%.2f, %.2f]' %(zr[0], zr[1]) for zr in z_limits]))
        print('\n')

# ------------------------------------------------------------------------------------
def get_zrange_for_line(line, obs_wave_range=[800, 2200]):
    '''
    Computes the redshift range in which a certain emission line is available, for a given observed wavelength window (by default corresponds to NIRISS WFSS)
    Input wavelengths must be in nm
    Returns min and max redshift
    '''
    rest_wave_dict = {'Lya': 121.6, 'OII': 372.7, 'Hd': 434.0, 'OIII-4363': 436.3, 'Hb': 486.1, 'OIII': 500.7, 'Ha+NII': 655.5, 'Ha': 656.2, 'SII': 671.7,
                      'SIII': 953.0, 'PaD': 1004.6, 'PaG': 1093.5, 'PaB': 1281.4,
                      'PaA': 1874.5}  # approximate wavelengths in nm

    if type(line) in [float, int, np.float64, np.int64]: rest_wave = line
    elif type(line) in [str, np.str_]: rest_wave = rest_wave_dict[line]

    z_min = (obs_wave_range[0] / rest_wave) - 1
    z_max = (obs_wave_range[1] / rest_wave) - 1

    return max(0, z_min), z_max

# ------------------------------------------------------------------------------------
def print_zrange_for_lines(lines, obs_wave_range=[800, 2200]):
    '''
    Prints the redshift range in which a certain emission line is available, for a given observed wavelength window (by default corresponds to NIRISS WFSS)
    Input wavelengths must be in nm
    '''
    lines = np.atleast_1d(lines)

    for line in lines:
        z_min, z_max = get_zrange_for_line(line, obs_wave_range=obs_wave_range)
        if type(line) is float or type(line) is int: line = '%.2f nm' %line

        print(f'Line: {line}, z=[{z_min:.2f}, {z_max:.2f}]')

# ---------------------------------------------------------------------------
def get_kpc_from_arc_at_redshift(arcseconds, redshift):
    '''
    Function to convert arcseconds on sky to physical kpc, at a given redshift
    '''
    cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)
    d_A = cosmo.angular_diameter_distance(z=redshift)
    kpc = (d_A * arcseconds * u.arcsec).to(u.kpc, u.dimensionless_angles()).value # in kpc
    print('%.2f arcseconds corresponds to %.2F kpc at target redshift of %.2f' %(arcseconds, kpc, redshift))
    return kpc

# ------------------------------------------------------------------------------------
def get_passage_obslist(filename=None, max_delta=0.02):
    '''
    Returns the list of PASSAGE observations available on MAST, as pandas dataframe
    '''
    if filename is None: filename = '/Users/acharyya/Work/astro/passage/passage_data/Proposal_IDs_1571.csv'
    df = pd.read_csv(filename, skiprows=4, usecols=['obs_id', 's_ra', 's_dec', 'filters', 't_exptime'])

    df = df.rename(columns={'s_ra':'ra', 's_dec':'dec'})
    df = df.sort_values(by='ra')
    df[['mode', 'filter']] = df['filters'].str.split(';', expand=True)
    df = df.drop('filters', axis=1)

    df['ra_group'] = df.sort_values('ra')['ra'].diff().gt(max_delta).cumsum()
    df['dec_group'] = df.sort_values('dec')['dec'].diff().gt(max_delta).cumsum()

    ra_list, dec_list = np.round(df[['ra', 'dec', 'ra_group', 'dec_group']].groupby(['ra_group', 'dec_group']).mean().reset_index()[['ra', 'dec']], 2).values.transpose()
    print(f'{len(df)} total PASSAGE pointings available, with {len(ra_list)} unique fields')

    return df, ra_list, dec_list

# ------------------------------------------------------------------------------------
def extract_field_from_obslist(index, df, ra_list, dec_list, max_delta=0.02):
    '''
    Returns the list of PASSAGE observations for the same field, as pandas dataframe
    '''
    df_sub = df[(np.abs(df['ra'] - ra_list[index]) <= max_delta) & (np.abs(df['dec'] - dec_list[index]) <= max_delta)].drop(['ra_group', 'dec_group'], axis=1)
    return df_sub

# -------------------------------------------------------------------------------------------------------
def mast_query(request):
    '''
    Perform a MAST query.

        Parameters
        ----------
        request (dictionary): The MAST request json object

        Returns head,content where head is the response HTTP headers, and content is the returned data
    This function has been borrowed from the MAST tutorial at https://mast.stsci.edu/api/v0/MastApiTutorial.html
    This function is just a place holder for now, not actually used in this script
    '''

    # Base API url
    request_url = 'https://mast.stsci.edu/api/v0/invoke'

    # Grab Python Version
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent": "python-requests/" + version}

    # Encoding the request as a json string
    req_string = json.dumps(request)
    req_string = urlencode(req_string)

    # Perform the HTTP request
    resp = requests.post(request_url, data="request=" + req_string, headers=headers)

    # Pull out the headers and response content
    head = resp.headers
    content = resp.content.decode('utf-8')

    return head, content

# -------------------------------------------------------------------------------------------------------
def watson_function(img_hdul_orig, img_hdul_new):
    '''
    Borrowed from Peter Watson
    Not yet used in this script
    '''
    xpix, ypix = img_hdul_orig[1].data.shape
    xx = np.asarray([0, 0, xpix, xpix,0])
    yy = np.asarray([0, ypix, ypix, 0, 0])
    orig_celestial = pywcs.WCS(img_hdul_orig[0].header).celestial
    new_celestial = pywcs.WCS(img_hdul_new[0].header).celestial
    x_p, y_p = pywcs.utils.pixel_to_pixel(orig_celestial, new_celestial, yy, xx)

    # extent_in_sky_coords = pywcs.WCS(img_hdul_orig[0].header).calc_footprint

    return x_p, y_p

# -------------------------------------------------------------------------------------------------------
def copy_from_hd_to_local(files_to_move=['*.txt', '*.png']):
    '''
    To copy all less heavy files (all txt and png) of each existing field in the HD to corresponding location locally
    '''
    hd_path = Path('/Volumes/Elements/acharyya_backup/Work/astro/passage/passage_output')
    local_path = Path('/Users/acharyya/Work/astro/passage/passage_output')

    available_fields = [os.path.split(item)[-1] for item in glob.glob(str(hd_path / 'Par*'))]

    for index, field in enumerate(available_fields):
        print(f'Doing field {field} which is {index+1} out of {len(available_fields)}..')
        dest_dir =local_path / field
        dest_dir.mkdir(parents=True, exist_ok=True)
        for files in files_to_move:
            target_files = str(hd_path / field / files)
            command = f'cp {target_files} {str(dest_dir)}/.'
            print(command)
            try: dummy = subprocess.check_output([command], shell=True)
            except: print('No such files. Skipping..')

    print('All done')




