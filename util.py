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
    parser.add_argument('--input_dir', metavar='input_dir', type=str, action='store', default='/Users/acharyya/Work/astro/passage/passage_data/', help='Where do the input files reside?')
    parser.add_argument('--output_dir', metavar='output_dir', type=str, action='store', default='/Users/acharyya/Work/astro/passage/passage_output/', help='Where do you want to store the outputs?')
    parser.add_argument('--code_dir', metavar='code_dir', type=str, action='store', default='/Users/acharyya/Work/astro/ayan_codes/passage/', help='Where is the source code?')
    parser.add_argument('--clobber', dest='clobber', action='store_true', default=False, help='Over-write existing plots? Default is no.')
    parser.add_argument('--silent', dest='silent', action='store_true', default=False, help='Suppress some generic pritn statements? Default is no.')
    parser.add_argument('--keep', dest='keep', action='store_true', default=False, help='Keep existing plot windows open? Default is no.')

    parser.add_argument('--field', metavar='field', type=str, action='store', default='Par008', help='Which passage field? Default is Par50')
    parser.add_argument('--id', metavar='id', type=int, action='store', default=100, help='Object ID. Default is 100')

    # ------- args added for get_spectra_from_beam.py ------------------------------
    parser.add_argument('--include_photometry', dest='include_photometry', action='store_true', default=False, help='Include photometry while computing fit parameters? Default is no.')
    parser.add_argument('--zmin', metavar='zmin', type=float, action='store', default=0.5, help='minimum of redshift range within which to search for lines; default is 0.5')
    parser.add_argument('--zmax', metavar='zmax', type=float, action='store', default=1.0, help='maximum of redshift range within which to search for lines; default is None')
    parser.add_argument('--line_list', metavar='line_list', type=str, action='store', default='all', help='Which emission lines to look for? Default is all') # OR set default to 'Lya,OII,Hb,OIII,Ha,Ha+NII,SII,SIII,PaB,He-1083,PaA'

    # ------- args added for read_line_catalog.py ------------------------------
    parser.add_argument('--nbins', metavar='nbins', type=int, action='store', default=30, help='No. of bins for plotting the histogram. Default is 30')
    parser.add_argument('--fontsize', metavar='fontsize', type=int, action='store', default=15, help='fontsize of plot labels, etc.; default is 15')
    parser.add_argument('--mag_lim', metavar='mag_lim', type=float, action='store', default=None, help='magnitude limit above which to search for targets; default is None')

    # ------- args added for plot_footprints.py ------------------------------
    parser.add_argument('--bg_file', metavar='bg_file', type=str, action='store', default=None, help='Which file to be used for plotting the background image?')
    parser.add_argument('--plot_zcosmos', dest='plot_zcosmos', action='store_true', default=False, help='Overplot the (thousands of) zCOSMOS targets? Default is no.')

    # ------- wrap up and processing args ------------------------------
    args = parser.parse_args()
    if args.line_list is not 'all': args.line_list = [item for item in args.line_list.split(',')]
    args.field = f'Par{int(args.field.split("Par")[1]):03d}'
    args.input_dir = Path(args.input_dir)
    args.output_dir = Path(args.output_dir)
    args.code_dir = Path(args.code_dir)

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

    filters = fix_filter_names(filters)
    print(f'For given filters {filters}..')

    if len(np.atleast_1d(lines)) == 1: lines = [lines]
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
    rest_wave_dict = {'OII': 372.7, 'Hb': 434.0, 'OIII': 436.3, 'Ha+NII': 655.5, 'Ha': 656.2, 'SII': 671.7,
                      'SIII': 953.0, 'PaD': 1004.6, 'PaG': 1093.5, 'PaB': 1281.4,
                      'PaA': 1874.5}  # approximate wavelengths in nm

    if type(line) is float or type(line) is int: rest_wave = line
    elif type(line) is str: rest_wave = rest_wave_dict[line]

    z_min = (obs_wave_range[0] / rest_wave) - 1
    z_max = (obs_wave_range[1] / rest_wave) - 1

    return max(0, z_min), z_max

# ------------------------------------------------------------------------------------
def print_zrange_for_lines(lines, obs_wave_range=[800, 2200]):
    '''
    Prints the redshift range in which a certain emission line is available, for a given observed wavelength window (by default corresponds to NIRISS WFSS)
    Input wavelengths must be in nm
    '''
    if len(np.atleast_1d(lines)) == 1: lines = [lines]

    for line in lines:
        z_min, z_max = get_zrange_for_line(line, obs_wave_range=obs_wave_range)
        if type(line) is float or type(line) is int: line = '%.2f nm' %line

        print(f'Line: {line}, z=[{z_min:.2f}, {z_max:.2f}]')

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
    orig_celestial = wcs.WCS(img_hdul_orig[0].header).celestial
    new_celestial = wcs.WCS(img_hdul_new[0].header).celestial
    x_p, y_p = wcs.utils.pixel_to_pixel(orig_celestial, new_celestial, yy, xx)

    # extent_in_sky_coords = wcs.WCS(img_hdul_orig[0].header).calc_footprint

    return x_p, y_p

