'''
    Filename: plot_NB_results.py
    Notes: Makes plots of results from NebulaBayes for different PASSAGE and GLASS glaaxies used in Acharyya+26
    Author : Ayan
    Created: 17-03-26
    Example: run plot_NB_results.py
'''
from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_data_path(field, args):
    '''
    Returns Path object with full path to the data based on a given field name
    '''
    if 'Par' in field: survey = 'passage'
    elif 'glass' in field: survey = 'glass'

    return args.root_dir / f'{survey}_output/' / f'{args.version_dict[survey]}' / f'{field}'

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.version_dict = {'passage': args.drv, 'glass': 'orig'}

    nbins = 5
    msize = 50
    # --------getting object list----------
    Par28_objects = [300, 1303, 1849, 2867]
    glass_objects = [1721, 1983, 1991, 1333]

    passage_objlist = [['Par028', item] for item in Par28_objects]
    glass_objlist = [['glass-a2744', item] for item in glass_objects]
    objlist = passage_objlist + glass_objlist

    # ------------setting up figure------------
    params = ['12 + log O/H', 'log P/k', 'log U']
    fig, axes = plt.subplots(len(params), 1, figsize=(5, 7), layout='constrained', sharex=True)

    # -------looping through objects-----------
    for index, obj in enumerate(objlist):
        field = obj[0]
        objid = obj[1]
        print(f'\nDoing object {field}-{objid} which is {index + 1} of {len(objlist)} objects..')

        data_path = get_data_path(field, args)
        file_path = data_path / f'{objid:05d}_NB_results' / f'_ellbin_{nbins}bins' / 'best_model_catalogs'
        files = glob.glob(str(file_path) + '/*.csv', recursive=True)
        files.sort()


        # -----------read in individual estimate file-----------
        for index2, file in enumerate(files):
            print(f'\tReading file {Path(file).stem}..')
            df = pd.read_csv(file, index_col='Parameter')
            
            # ----------looping over parameters----------------------
            for index3, param in enumerate(params):
                quant = df.loc[param, 'Estimate']
                quant_low = df.loc[param, 'CI68_low']
                quant_high = df.loc[param, 'CI68_high']

                axes[index3].scatter(index, quant, s=msize, marker='s' if 'glass' in field else 'o', ec='k', lw=1, c=index2, cmap='rainbow', vmin=0, vmax=nbins)
                try: axes[index3].errorbar(index, quant, yerr=[[quant - quant_low], [quant_high - quant]], fmt='none', c='grey', lw=0.5)
                except: pass


    # ---------annotate figure-----------------
    for index3, param in enumerate(params):
        axes[index3].set_ylabel(param, fontsize=args.fontsize)
    
    axes[-1].set_xticks(np.arange(len(objlist)))
    axes[-1].set_xticklabels([item[1] for item in objlist], fontsize=args.fontsize, rotation=45)

    # ---------saving figure-----------
    save_fig(fig, args.root_dir / 'zgrad_paper_plots', 'NB_all_results.png', args)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

