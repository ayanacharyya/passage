'''
    Filename: plot_footprints.py
    Notes: Plots all PASSAGE field footprints and overlays with existing large survey footprints from MAST
           This script is still W.I.P and does not yet work (for now, the workaround is plotting the footprints manually using DS9)
    Author : Ayan
    Created: 18-06-24
    Last modified: 18-06-24
    Example: run read_line_catalog.py --input_dir /Users/acharyya/Work/astro/passage/passage_data/ --output_dir /Users/acharyya/Work/astro/passage/passage_output/ --field Par21
             run plot_footprints.py
'''

from header import *
from util import *
from make_region_files import get_passage_footprints

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
def overlay_footprint(img_hdul_orig, img_hdul_new):
    '''
    Borrowed from Peter Watson
    '''
    xpix, ypix = img_hdul_orig[1].data.shape
    xx = np.asarray([0, 0, xpix, xpix,0])
    yy = np.asarray([0, ypix, ypix, 0, 0])
    x_celestial = wcs.WCS(img_hdul_orig[0].header).celestial
    y_celestial = wcs.WCS(img_hdul_new[0].header).celestial
    x_p, y_p = wcs.utils.pixel_to_pixel(x_celestial, y_celestial, yy, xx)

    # extent_in_sky_coords = wcs.WCS(img_hdul_orig[0].header).calc_footprint

    return x_p, y_p

# -------------------------------------------------------------------------------------------------------
def mast_query(request):
    '''
    Perform a MAST query.

        Parameters
        ----------
        request (dictionary): The MAST request json object

        Returns head,content where head is the response HTTP headers, and content is the returned data
    This function has been borrowed from the MAST tutorial at https://mast.stsci.edu/api/v0/MastApiTutorial.html
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
def plot_footprints(df, args, fig=None):
    '''
    Plots the footprint for a given survey
    Returns fig handle
    '''
    wwt = WWTQtClient()
    wwt = WWTQtClient(block_until_ready=True)

    if fig is None: fig, ax = plt.subplots(figsize=(8, 10))
    else: ax = fig.gca()

    # -------do the plotting---------
    for index in len(df):
        footprint = df.loc[index]
        hdf = SkyCoord(footprint.ra, footprint.dec, unit=u.deg)
        wwt.background = wwt.imagery.visible.sdss
        wwt.center_on_coordinates(hdf, fov=2.8 * u.arcmin)

        fov = wwt.add_fov(wwt.instruments.jwst_niriss, center=hdf, rotate=footprint.pa * u.deg, line_color=footprint.color)

    # ------axes aesthetics-------------
    ax.set_xlabel('RA', fontsize=args.fontsize)
    ax.set_ylabel('Dec', fontsize=args.fontsize)

    figname = args.output_dir / 'footprints.png'
    fig.savefig(figname)
    print(f'Saved plot to {figname}')
    plt.show(block=False)

    return fig

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # -----reading the passage footprints---------------
    filenames = list(args.input_dir.glob('*catalog*'))
    if len(filenames) > 0: df_passage = get_passage_footprints(filenames[0])
    else: sys.exit('No survey file in %s' %args.input_dir)

    # ------querying MAST for common surveys to get their footprints---------
    surveys_dict = {}#{'GOODS-N':'red', 'GOODS-S':'crimson'}
    df_others = pd.DataFrame(columns={'name':[], 'ra':[], 'dec':[], 'size':[], 'pa':[], 'color':[]})

    for survey in surveys_dict.keys():
        resolver_request = {'service': 'Mast.Name.Lookup', 'params': {'input': survey, 'format': 'json'}}
        headers, resolved_object_string = mast_query(resolver_request)
        resolved_object = json.loads(resolved_object_string)
        pp.pprint(resolved_object)

        ra = resolved_object['resolvedCoordinate'][0]['ra']
        dec = resolved_object['resolvedCoordinate'][0]['decl']

        ra_width = 0
        dec_width = 0
        pa = 0

        df_others.iloc[len(df_others)] = [survey, ra, dec, ra_width, dec_width, pa, surveys_dict[survey]]

    # -----making the plot--------------
    df_master = df_passage.concat(df_others)
    fig = plot_footprints(df_master, args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
