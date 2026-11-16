'''
    Filename: compute_NB_1dstack.py
    Notes: Computes metallicity using NebulaBayes using integraed emission line fluxes for all available 1D stacks 
    Author : Ayan
    Created: 16-11-26
    Example: run compute_NB_1dstack.py --system ssd
             run compute_NB_1dstack.py --system ssd --nb_mass_prior
             run compute_NB_1dstack.py --system ssd --nb_prior
'''

from header import *
from util import *
from plot_stacked_maps import compute_Z_NB

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')
    args.dered_in_NB = True
    args.fits_dir = args.output_dir / 'catalogs' / 'stacked_spec_ms'
    args.NB_text = '_resetBD' if args.reset_bd else '_Zprior' if args.nb_prior else '_Z_massprior' if args.nb_mass_prior else '_Z_centralprior' if args.nb_central_prior else ''
    args.bin_text = f'stacked_integrated{args.NB_text}'
    args.line_list = ['OII', 'OIII', 'Hb']

    # -----dict for converting line label names to those acceptable to NB---------
    mylabel_to_stacklabel_dict = {'OII':'OII', 'OIII':'OIIIboth', 'Hb':'Hb', 'Ha':'HaNII', 'SII':'SII'}
    line_label_dict = {'OII':'OII3726_29', 'Hb':'Hbeta', 'OIII':'OIII5007', 'Ha':'Halpha', 'SII':'SII6716_31'}
    line_label_dict.update({k.upper(): v for k, v in line_label_dict.items()})

    line_label_array = [line_label_dict[line] for line in args.line_list]
    line_waves_array = [rest_wave_dict[line] * 10 for line in args.line_list] # factor of 10 to convert from nm to Angstroms
    line_int_array = []

    # --------finding files---------------
    files = glob.glob(str(args.fits_dir) + '/meas_logm*.dat')
    df_master = pd.DataFrame(columns=['log_mass_min', 'log_mass_max', 'log_mass'])

    # ------looping over stacks-----------
    for index, file in enumerate(files):
        df = pd.read_csv(file, sep=r'\s*\|\s*', engine='python').iloc[:, 1:-1]
        log_mass_min, log_mass_max = map(float, Path(file).stem.split('_')[-1].split('-'))
        log_mass_mean = np.mean(np.array([log_mass_min, log_mass_max]))

        print(f'\nComputing metallicity for stack {log_mass_min} <= log_mass <= {log_mass_max} which is {index + 1} of {len(files)} stacks..')
        this_line_int_array = []
        
        for line in args.line_list:
            line_in_df = mylabel_to_stacklabel_dict[line]
            obj = df[df['line'] == line_in_df]

            flux = obj['flux_norm'].values[0]# * obj['flux_scale'].values[0]
            flux_u = obj['flux_norm_err'].values[0]# * obj['flux_scale'].values[0]
            
            # ----------deblending flux--------------------
            if line_in_df == 'OIIIboth': # special treatment for OIII 5007 line
                ratio_5007_to_4959 = 2.98 # from grizli source code
                factor = ratio_5007_to_4959 / (1 + ratio_5007_to_4959)
                print(f'\tCorrecting OIII for 4959 component, by factor of {factor:.3f}')
                flux = flux * factor
                flux_u = flux_u * factor
                        
            # -------preparing flux list for NB--------------------
            this_line_int_array.append(ufloat(flux, flux_u))
        line_int_array.append(this_line_int_array)

        # ------appending to master df----------
        thisrow = {'log_mass_min': log_mass_min, 'log_mass_max': log_mass_max, 'log_mass': log_mass_mean}
        df_master.loc[len(df_master)] = thisrow

    # ----------calling NB----------------------
    line_int_array = np.transpose(line_int_array) # because NB needs nlines x nobj shape
    args.this_stack_mass = df_master['log_mass'].values
    logOH_int = compute_Z_NB(line_label_array, line_waves_array, line_int_array, args)
    plt.close('all') # close all NB corner plots

    # ------appending to master df----------
    df_master['logOH'] = unp.nominal_values(logOH_int.data)
    df_master['logOH_u'] = unp.std_devs(logOH_int.data)

    # ----------saving master dataframe-----------
    outfilename = args.output_dir / 'catalogs' / f'logOH_1dstacking{args.NB_text}.csv'
    if len(df_master) > 0:
        df_master.to_csv(outfilename, index=None)
        print(f'Saved results as {outfilename}')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')