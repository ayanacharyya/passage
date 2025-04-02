##!/usr/bin/env python3

"""

    Filename :   util.py
    Notes :      Contains various generic utility functions and classes used by the other scripts in PASSAGE, including a function to parse args
    Author :    Ayan
    Created: 11-06-24
    Last modified: 11-06-24

"""

from header import *

# -------------------------------------------------------------------------------------------------------
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

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
    parser.add_argument('--forpaper', dest='forpaper', action='store_true', default=False, help='Format plots to paper quality? Default is no.')
    parser.add_argument('--fortalk', dest='fortalk', action='store_true', default=False, help='Format plots suitable for putting in talks? Default is no.')
    parser.add_argument('--drv', metavar='drv', type=str, action='store', default='v0.1', help='Which data reduction version? Default v0.1')
    parser.add_argument('--fontsize', metavar='fontsize', type=int, action='store', default=10, help='fontsize of plot labels, etc.; default is 15')

    parser.add_argument('--field', metavar='field', type=str, action='store', default='Par3', help='Which passage field? Default is Par50')
    parser.add_argument('--do_only_fields', metavar='do_only_fields', type=str, action='store', default=None, help='Which passage field? Default is Par50')
    parser.add_argument('--id', metavar='id', type=str, action='store', default=None, help='Object ID. Default is None')

    # ------- args added for make_spectra_from_beam.py ------------------------------
    parser.add_argument('--include_photometry', dest='include_photometry', action='store_true', default=False, help='Include photometry while computing fit parameters? Default is no.')
    parser.add_argument('--zmin', metavar='zmin', type=float, action='store', default=0.5, help='minimum of redshift range within which to search for lines; default is 0.5')
    parser.add_argument('--zmax', metavar='zmax', type=float, action='store', default=1.0, help='maximum of redshift range within which to search for lines; default is None')
    parser.add_argument('--line_list', metavar='line_list', type=str, action='store', default='all', help='Which emission lines to look for? Default is all') # OR set default to 'Lya,OII,Hb,OIII,Ha,Ha+NII,SII,SIII,PaB,He-1083,PaA'
    parser.add_argument('--pixscale', metavar='pixscale', type=float, action='store', default=0.04, help='Pixel scale (in arcsec/pixel) of the thumbnails produced; default is 0.04')

    # ------- args added for read_line_catalog.py ------------------------------
    parser.add_argument('--nbins', metavar='nbins', type=int, action='store', default=30, help='No. of bins for plotting the histogram. Default is 30')
    parser.add_argument('--mag_lim', metavar='mag_lim', type=float, action='store', default=None, help='magnitude limit above which to search for targets; default is None')

    # ------- args added for make_region_files.py ------------------------------
    parser.add_argument('--survey', metavar='survey', type=str, action='store', default='passage', help='Which survey to be used for making the region files?')
    parser.add_argument('--split_regions_by_fields', dest='split_regions_by_fields', action='store_true', default=False, help='Split regions into different files based on if they are contained within each PASSAGE fields? Default is no.')

    # ------- args added for plot_footprints.py ------------------------------
    parser.add_argument('--bg', metavar='bg', type=str, action='store', default='COSMOS', help='Which survey to be used to plot background? Default is COSMOS')
    parser.add_argument('--bg_file', metavar='bg_file', type=str, action='store', default=None, help='Which file to be used for plotting the background image?')
    parser.add_argument('--bg_image_dir', metavar='bg_image_dir', type=str, action='store', default=None, help='Which folder to be used for looking for the background image?')
    parser.add_argument('--plot_zcosmos_objects', dest='plot_zcosmos_objects', action='store_true', default=False, help='Overplot the (thousands of) zCOSMOS targets? Default is no.')
    parser.add_argument('--plot_cosmos2020_objects', dest='plot_cosmos2020_objects', action='store_true', default=False, help='Overplot the (millions of) COSMOS2020 targets? Default is no.')
    parser.add_argument('--plot_cosmoswebb_objects', dest='plot_cosmoswebb_objects', action='store_true', default=False, help='Overplot the (millions of) COSMOSWebb targets? Default is no.')
    parser.add_argument('--only_passage_regions', dest='only_passage_regions', action='store_true', default=False, help='Overplot ONLY the PASSAGE Par regions? Default is no.')
    parser.add_argument('--fg_image_dir', metavar='fg_image_dir', type=str, action='store', default=None, help='Which folder to be used for looking for the foreground image?')
    parser.add_argument('--fg_file', metavar='fg_file', type=str, action='store', default=None, help='Which file to be used for plotting the foreground images?')
    parser.add_argument('--plot_fg_data', dest='plot_fg_data', action='store_true', default=False, help='Overplot the data of the foreground file? Default is no.')

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
    parser.add_argument('--remake_figures', dest='remake_figures', action='store_true', default=False, help='Re-make the figures made by PASSAGEPipe and tarball them? Default is no.')
    parser.add_argument('--fit_redshift', dest='fit_redshift', action='store_true', default=False, help='Fit redshifts for this field (will take significantly longer)? Default is no.')

    # ------- args added for make_diagnostic_maps.py ------------------------------
    parser.add_argument('--plot_target_frame', dest='plot_target_frame', action='store_true', default=False, help='Annotate plot axes in the object/target frame of reference? Default is no.')
    parser.add_argument('--arcsec_limit', metavar='arcsec_limit', type=float, action='store', default=1.0, help='Half box size (in arcsec) of the thumbnails to plot; default is 1.5')
    parser.add_argument('--vorbin', dest='vorbin', action='store_true', default=False, help='Voronoi bin the 2D emission line maps? Default is no.')
    parser.add_argument('--voronoi_snr', metavar='voronoi_snr', type=float, action='store', default=3, help='Target SNR to Voronoi bin the emission line maps to; default is 3')
    parser.add_argument('--voronoi_line', metavar='voronoi_line', type=str, action='store', default='Hb', help='Which emission line to be used for computing the Voronoi bins? Default is None i.e., the given emission line itself')
    parser.add_argument('--flam_max', metavar='flam_max', type=float, action='store', default=None, help='Maximum y-axis limit for f_lambda (in units of 1e-19 ergs/s/cm^2/A); default is None')
    parser.add_argument('--plot_radial_profiles', dest='plot_radial_profiles', action='store_true', default=False, help='Plot radial profiles corresponding to the 2D maps? Default is no.')
    parser.add_argument('--snr_cut', metavar='snr_cut', type=float, action='store', default=0., help='Impose an SNR cut on the emission line maps to; default is 0')
    parser.add_argument('--only_seg', dest='only_seg', action='store_true', default=False, help='Cut out the emission line plots corresponding to the grizli segmentation map? Default is no.')
    parser.add_argument('--write_file', dest='write_file', action='store_true', default=False, help='Write the measured quantities to a master dataframe? Default is no.')
    parser.add_argument('--plot_mappings', dest='plot_mappings', action='store_true', default=False, help='Plot emission line locations as per MAPPINGS predictions (will lead to crowding of many lines)? Default is no.')
    parser.add_argument('--hide', dest='hide', action='store_true', default=False, help='Hide (do not display) the plots just made? Default is no.')
    parser.add_argument('--trim_filter_by_wavelength_factor', metavar='trim_filter_by_wavelength_factor', type=float, action='store', default=0.05, help='Impose a trimming factor on wavelength for each filter, to be applied on both blue and red ends of a filter; default is 0.05 (i.e. 5%)')
    parser.add_argument('--plot_starburst', dest='plot_starburst', action='store_true', default=False, help='Plot the starbursty-ness map instead of the full diagnostic figure? Default is no.')
    parser.add_argument('--plot_slope_vs_mass', dest='plot_slope_vs_mass', action='store_true', default=False, help='Plot the slope of Ha/F115W profile vs stellar mass instead of the full diagnostic figure? Default is no.')
    parser.add_argument('--plot_vorbin', dest='plot_vorbin', action='store_true', default=False, help='Plot the voronoi bins? Default is no.')
    parser.add_argument('--plot_snr', dest='plot_snr', action='store_true', default=False, help='Plot the SNR map for a given 2D plot? Default is no.')
    parser.add_argument('--plot_metallicity', dest='plot_metallicity', action='store_true', default=False, help='Plot the metallicity map instead of the full diagnostic figure? Default is no.')
    parser.add_argument('--test_cutout', dest='test_cutout', action='store_true', default=False, help='Plot the cutout 2D clear image as a testing phase? Default is no.')
    parser.add_argument('--plot_direct_filters', dest='plot_direct_filters', action='store_true', default=False, help='Plot the direct filter images instead of the full diagnostic figure? Default is no.')
    parser.add_argument('--debug_vorbin', dest='debug_vorbin', action='store_true', default=False, help='Do extra plots and prints for debugging voronoi binning? Default is no.')
    parser.add_argument('--do_not_correct_flux', dest='do_not_correct_flux', action='store_true', default=False, help='Skip the step where it corrects for certain belnded line fluxes e.g., OIII5007, Ha, SII 6717? Default is no.')
    parser.add_argument('--plot_AGN_frac', dest='plot_AGN_frac', action='store_true', default=False, help='Plot AGN fraction 2D map (based on BPT diagram)? Default is no.')
    parser.add_argument('--diverging_cmap', metavar='diverging_cmap', type=str, action='store', default='vik', help='Which diverging colormap to use (out of managua, vanimo, lisbon, berlin)? Default is cork')
    parser.add_argument('--do_not_correct_pixel', dest='do_not_correct_pixel', action='store_true', default=False, help='Skip the step where it corrects for pixel offset in the emission lines compared to direct images? Default is no.')
    parser.add_argument('--Zbranch', metavar='Zbranch', type=str, action='store', default='low', help='Which R23 branch to be used (choose between high/low)? Default is low')
    parser.add_argument('--plot_ionisation_parameter', dest='plot_ionisation_parameter', action='store_true', default=False, help='Plot the plot_ionisation_parameter map along with metallicity? Default is no.')
    parser.add_argument('--ignore_combined_method', dest='ignore_combined_method', action='store_true', default=False, help='Ignore the combined method (S6 of KD02) while computing R23 metallicity and rely solely on R23? Default is no.')
    parser.add_argument('--Zdiag', metavar='Zdiag', type=str, action='store', default='R23', help='Which metallicity diagnostic to use (choose between R23,R3,O3S2,O3O2,Te,P25,NB? Default is R23')
    parser.add_argument('--debug_Zdiag', dest='debug_Zdiag', action='store_true', default=False, help='Make additional plots to debug the metallicity diagnostic implementation? Default is no.')
    parser.add_argument('--mask_agn', dest='mask_agn', action='store_true', default=False, help='Mask out the AGN-dominated pixels from all metallicity estimates? Default is no.')
    parser.add_argument('--plot_circle_at_arcsec', metavar='plot_circle_at_arcsec', type=float, action='store', default=None, help='Radius in arcseconds of a circle to be plotted on every 2D map; default is None')
    parser.add_argument('--plot_ratio_maps', dest='plot_ratio_maps', action='store_true', default=False, help='Plot the line ratio maps for a given 2D plot? Default is no.')
    parser.add_argument('--AGN_diag', metavar='AGN_diag', type=str, action='store', default='VO87', help='Which AGN-SF BPT-like diagnostic to use (choose between VO87,H21,O2O3,O2Hb,Ne3O2? Default is VO87')
    parser.add_argument('--use_variable_N2Ha', dest='use_variable_N2Ha', action='store_true', default=False, help='Use variable Ha/(NII + Ha) ratio across the face of the galaxy, to compute the Ha for the x-axis of BPT diagram, instead of constant 0.82? Default is no.')
    parser.add_argument('--no_text_on_plot', dest='no_text_on_plot', action='store_true', default=False, help='Skip putting text annotations on plot2D? Default is no.')
    parser.add_argument('--plot_models', dest='plot_models', action='store_true', default=False, help='Overplot MAPPINGS photoionisation models on the BPT diagram? Default is no.')
    parser.add_argument('--plot_DIG', dest='plot_DIG', action='store_true', default=False, help='Plot DIG diagnostics? Default is no.')
    parser.add_argument('--output_subdir', metavar='output_subdir', type=str, action='store', default=None, help='Any specific subdirectory (with output_dir) to put output files in; default is None (i.e. files will be put in output_dir)')
    parser.add_argument('--use_original_NB_grid', dest='use_original_NB_grid', action='store_true', default=False, help='Use the original, unmodified NebulaBayes grid? Default is no.')

    # ------- args added for get_field_stats.py ------------------------------
    parser.add_argument('--EW_thresh', metavar='EW_thresh', type=float, action='store', default=300.0, help='Rest-frame EW threshold to consider good detection for emission line maps; default is 300')
    parser.add_argument('--SNR_thresh', metavar='SNR_thresh', type=float, action='store', default=10.0, help='SNR threshold to consider good detection for emission line maps; default is 10')
    parser.add_argument('--log_SFR_thresh', metavar='log_SFR_thresh', type=float, action='store', default=0, help='SFR threshold (in log) to consider highly star-forming; default is 0')
    parser.add_argument('--log_sSFR_thresh', metavar='log_sSFR_thresh', type=float, action='store', default=-9, help='specific SFR threshold (in log) to consider highly star-forming; default is 0')
    parser.add_argument('--do_all_fields', dest='do_all_fields', action='store_true', default=False, help='Include ALL available fields? Default is no.')
    parser.add_argument('--merge_visual', dest='merge_visual', action='store_true', default=False, help='Include visually inspected dataframe for Venn diagrams? Default is no.')
    parser.add_argument('--plot_conditions', metavar='plot_conditions', type=str, action='store', default='detected', help='Which conditions are plotted in the Venn diagram? Default is None')
    parser.add_argument('--plot_EW_hist', dest='plot_EW_hist', action='store_true', default=False, help='Plot EW histograms for each line? Default is no.')
    parser.add_argument('--clobber_venn_df', dest='clobber_venn_df', action='store_true', default=False, help='Over-write existing dataframe which is a result of intersection of Venn diagram? Default is no.')
    parser.add_argument('--a_thresh', metavar='a_thresh', type=float, action='store', default=2.4, help='a_image (semi-major axis in pixels) threshold to consider good detection for emission line maps; default is 2.4')
    parser.add_argument('--plot_pie', dest='plot_pie', action='store_true', default=False, help='Plot pie chart instead of Venn diagram? Default is no.')
    parser.add_argument('--plot_sunburst', dest='plot_sunburst', action='store_true', default=False, help='Plot sunburst chart instead of Venn diagram? Default is no.')
    parser.add_argument('--plot_columns', dest='plot_columns', action='store_true', default=False, help='Plot (user input) columns vs each other for the entire cross-matched dataset? Default is no.')
    parser.add_argument('--clobber_cosmos', dest='clobber_cosmos', action='store_true', default=False, help='Over-write existing COSMOS dataframes (subset of columns)? Default is no.')

    # ------- args added for make_cosmos_plots.py ------------------------------
    parser.add_argument('--xcol', metavar='xcol', type=str, action='store', default='lp_mass_best', help='Column name in COSMOS2020 catalog to be used as the quantity on x-axis? Default is lp_mass_best')
    parser.add_argument('--ycol', metavar='ycol', type=str, action='store', default='lp_SFR_best', help='Column name in COSMOS2020 catalog to be used as the quantity on y-axis? Default is lp_SFR_best')
    parser.add_argument('--colorcol', metavar='colorcol', type=str, action='store', default='ez_z_phot', help='Column name in COSMOS2020 catalog to be used as the color axis? Default is ez_z_phot')
    parser.add_argument('--only_download', dest='only_download', action='store_true', default=False, help='Perform only the downloading (from GDrive) step for each field? Default is no.')

    # ------- args added for make_passage_plots.py ------------------------------
    parser.add_argument('--plot_BPT', dest='plot_BPT', action='store_true', default=False, help='Plot BPT? Default is no.')
    parser.add_argument('--plot_separately', dest='plot_separately', action='store_true', default=False, help='Plot BPTs separetly for each object? Default is no.')
    parser.add_argument('--plot_flux_vs_mag', dest='plot_flux_vs_mag', action='store_true', default=False, help='Plot line flux vs mag for each object? Default is no.')
    parser.add_argument('--foggie_comp', dest='foggie_comp', action='store_true', default=False, help='Plot Zgrad vs redshift to exact same limits as the FOGGIE plot, for comparison? Default is no.')
    parser.add_argument('--use_only_good', dest='use_only_good', action='store_true', default=False, help='Use only the pre-determined good galaxies? Default is no.')
    parser.add_argument('--plot_full_BPT', dest='plot_full_BPT', action='store_true', default=False, help='Plot BPT for the full speccat of a given PASSAGE field? Default is no.')
    parser.add_argument('--log_colorcol', dest='log_colorcol', action='store_true', default=False, help='Take log of the quantity before color coding? Default is no.')
    parser.add_argument('--plot_mass_excitation', dest='plot_mass_excitation', action='store_true', default=False, help='Plot mass-excitation diagram for the sample? Default is no.')
    parser.add_argument('--plot_full_mass_excitation', dest='plot_full_mass_excitation', action='store_true', default=False, help='Plot mass-excitation diagram for the full speccat of a given PASSAGE field? Default is no.')

    # ------- args added for compute_stellar_mass.py ------------------------------
    parser.add_argument('--plot_transmission', dest='plot_transmission', action='store_true', default=False, help='Plot transmission curves for all filters? Default is no.')
    parser.add_argument('--plot_SED', dest='plot_SED', action='store_true', default=False, help='Plot SED for all filters? Default is no.')
    parser.add_argument('--plot_niriss_direct', dest='plot_niriss_direct', action='store_true', default=False, help='Plot 2D NIRISS direct images? Default is no.')
    parser.add_argument('--plot_filter_cutouts', dest='plot_filter_cutouts', action='store_true', default=False, help='Plot 2D image cutouts? Default is no.')
    parser.add_argument('--plot_cutout_errors', dest='plot_cutout_errors', action='store_true', default=False, help='Plot 2D uncertainty maps of cutouts? Default is no.')
    parser.add_argument('--plot_all', dest='plot_all', action='store_true', default=False, help='Plot cutouts for ALL filters? Default is no.')
    parser.add_argument('--fit_sed', dest='fit_sed', action='store_true', default=False, help='Fit SEDs? Default is no.')
    parser.add_argument('--run', metavar='run', type=str, action='store', default='first_try', help='Which run label should be assigned to the SED fit (this decides subfolders where results are stored)? Default is "first_try"')
    parser.add_argument('--ncpus', metavar='ncpus', type=int, action='store', default=1, help='No. of processors to use for BAGPIPES. Default is all (4)')
    parser.add_argument('--clobber_sed_photcat', dest='clobber_sed_photcat', action='store_true', default=False, help='Over-write existing phot cat to be used for SED fitting with bagpipes? Default is no.')
    parser.add_argument('--use_only_bands', metavar='use_only_bands', type=str, action='store', default=None, help='Which bands to be used for SED fitting with bagpipes? Default is None, i.e., use all bands')
    parser.add_argument('--plot_restframe', dest='plot_restframe', action='store_true', default=False, help='Plot fitted SED in restframe? Default is no.')
    parser.add_argument('--log_x', dest='log_x', action='store_true', default=False, help='Plot x-axis of SED in log scale? Default is no.')
    parser.add_argument('--test_sed', metavar='test_sed', type=int, action='store', default=None, help='Fit and plot just one object fo the given ID SED as a test? Default is None')
    parser.add_argument('--include_cosmoswebb', dest='include_cosmoswebb', action='store_true', default=False, help='Include COSMOS Webb filters in the SED fitting? Default is no.')
    parser.add_argument('--do_field', metavar='do_field', type=str, action='store', default=None, help='To run SED fitting on all overlapping objects of a given PASSAGE field? Default is None, i.e., use all fields')

    # ------- args added for make_toy_model.py ------------------------------
    parser.add_argument('--res', metavar='res', type=str, action='store', default='0.2', help='Resolution element to convolve toy model by, as a fraction of the full radial extent; default is 0.2')
    parser.add_argument('--num_line', metavar='num_line', type=str, action='store', default='sii', help='Which line to be used as the numerator? Default SII')
    parser.add_argument('--den_line', metavar='den_line', type=str, action='store', default='ha', help='Which line to be used as the denominator? Default Ha')
    parser.add_argument('--num_snr', metavar='num_snr', type=float, action='store', default=1e5, help='SNR with which to add noise for the numerator line profile; default is 1e5 i.e., basically no noise')
    parser.add_argument('--den_snr', metavar='den_snr', type=float, action='store', default=1e5, help='SNR with which to add noise for the denominator line profile; default is 1e5 i.e., basically no noise')

    # ------- args added for make_mappings_grid.py ------------------------------
    parser.add_argument('--ynum_line', metavar='ynum_line', type=str, action='store', default='OIII', help='Which line to be used as the numerator on y-axis? Default OIII')
    parser.add_argument('--yden_line', metavar='yden_line', type=str, action='store', default='Hb', help='Which line to be used as the denominator on y-axis? Default Hb')
    parser.add_argument('--xnum_line', metavar='xnum_line', type=str, action='store', default='SII', help='Which line to be used as the numerator on x-axis? Default SII')
    parser.add_argument('--xden_line', metavar='xden_line', type=str, action='store', default='Ha', help='Which line to be used as the denominator on x-axis? Default Ha')
    parser.add_argument('--mappings_dir', metavar='mappings_dir', type=str, action='store', default='/Users/acharyya/Work/astro/Mappings', help='Where do MAPPINGS files reside?')
    parser.add_argument('--geometry', metavar='geometry', type=str, action='store', default='s', help='Models corresponding to which geometry to be used; choose between s (spherical) and p(plane parallel)? Default s')
    parser.add_argument('--iso', metavar='iso', type=str, action='store', default='P', help='Use isobaric (P) or isodensity (d) models? Default P')
    parser.add_argument('--quantity1', metavar='quantity1', type=str, action='store', default='Z', help='Column name for quantity1 in the ratios dataframe; Default is Z')
    parser.add_argument('--quantity2', metavar='quantity2', type=str, action='store', default='log(q)', help='Column name for quantity3 in the ratios dataframe; Default is log(q)')
    parser.add_argument('--quantity3', metavar='quantitt3', type=str, action='store', default='log(P/k)', help='Column name for quantity3 in the ratios dataframe; Default is log(P/k)')
    parser.add_argument('--slice_at_quantity1', metavar='slice_at_quantity1', type=str, action='store', default=None, help='Value of quantity1 to which the model will be curtailed to? Default is None, i.e., use all values')
    parser.add_argument('--slice_at_quantity2', metavar='slice_at_quantity2', type=str, action='store', default=None, help='Value of quantity2 to which the model will be curtailed to? Default is None, i.e., use all values')
    parser.add_argument('--slice_at_quantity3', metavar='slice_at_quantity3', type=str, action='store', default=None, help='Value of quantity3 to which the model will be curtailed to? Default is None, i.e., use all values')
    parser.add_argument('--plot_grid', dest='plot_grid', action='store_true', default=False, help='Plot grid of two line ratios? Default is no.')
    parser.add_argument('--plot_model', dest='plot_model', action='store_true', default=False, help='Plot a line ratio vs model parameters? Default is no.')
    parser.add_argument('--annotate', dest='annotate', action='store_true', default=False, help='Annotate the ratio grid plot with arrows? Default is no.')
    parser.add_argument('--fit_y_envelope', dest='fit_y_envelope', action='store_true', default=False, help='Fit an envelope for the maximum y envelope in a given model grid? Default is no.')

    # ------- args added for extract_spectrum_from_grism.py ------------------------------
    parser.add_argument('--extract_arcsec', metavar='extract_arcsec', type=float, action='store', default=0.5, help='Spatial (cross-dispersion) extent in arcseconds from within which grism 2D spectra will be extracted; default is 0.5')
    parser.add_argument('--debug_zero_order', dest='debug_zero_order', action='store_true', default=False, help='Debug mode to find the zero order location in the grism image? Default is no.')

    # ------- wrap up and processing args ------------------------------
    args = parser.parse_args()
    if 'v' not in args.drv: args.drv = 'v' + args.drv
    if args.line_list != 'all': args.line_list = [item for item in args.line_list.split(',')]

    args.field_arr = args.field.split(',')
    for index in range(len(args.field_arr)):
        if 'Par' in args.field_arr[index]: args.field_arr[index] = f'Par{int(args.field_arr[index].split("Par")[1]):03d}'
    args.field = args.field_arr[0]

    if args.id is not None: args.id = [int(item) for item in args.id.split(',')]
    if args.slice_at_quantity1 is not None: args.slice_at_quantity1 = [float(item) for item in args.slice_at_quantity1.split(',')]
    if args.slice_at_quantity2 is not None: args.slice_at_quantity2 = [float(item) for item in args.slice_at_quantity2.split(',')]
    if args.slice_at_quantity3 is not None: args.slice_at_quantity3 = [float(item) for item in args.slice_at_quantity3.split(',')]

    if args.system == 'hd' and not os.path.exists('/Volumes/Elements/'): args.system = 'local'
    if args.line_list == 'all': args.line_list = ['Lya', 'OII', 'NeIII-3867', 'Hb', 'OIII-4363', 'OIII', 'Ha', 'NII','Ha+NII', 'SII', 'ArIII-7138', 'SIII', 'PaD','PaG','PaB','HeI-1083','PaA']

    root_dir = '/Users/acharyya/Work/astro/passage' if 'local' in args.system else '/Volumes/Elements/acharyya_backup/Work/astro/passage' if 'hd' in args.system else '/Users/acharyya/Library/CloudStorage/GoogleDrive-ayan.acharyya@inaf.it/My Drive/passage' if 'gdrive' in args.system else ''
    args.root_dir = Path(root_dir)

    if 'glass' in args.field: survey_name = 'glass'
    else: survey_name = 'passage'
    if args.input_dir is None:
        args.input_dir = args.root_dir / f'{survey_name}_data/'
    if args.output_dir is None:
        args.output_dir = args.root_dir / f'{survey_name}_output/'
    if 'glass' not in args.field:
        args.input_dir = args.input_dir / args.drv
        args.output_dir = args.output_dir / args.drv

    args.input_dir = Path(args.input_dir)
    args.output_dir = Path(args.output_dir)
    args.code_dir = Path(args.code_dir)
    args.mappings_dir = Path(args.mappings_dir)

    if args.filters is None:
        if 'glass' in args.field:
            args.filters = ['F115W', 'F150W', 'F200W']
        else:
            if 'available_filters_for_field_dict' not in locals(): available_filters_for_field_dict = get_passage_filter_dict(args)
            args.filters = available_filters_for_field_dict[args.field]
    else:
        args.filters = args.filters.split(',')

    args.plot_conditions = args.plot_conditions.split(',')
    args.res = [float(item) for item in args.res.split(',')]

    if args.fortalk:
        print(f'Setting up plots for talks..')
        setup_plots_for_talks()

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
    Returns min and max redshift/s
    '''
    filters = fix_filter_names(filters)

    z_limits = []
    for index, filter in enumerate(filters):
        obs_wave_range = np.array(filter_waverange_dict[filter])
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

# -------------------------------------------------------------------------------------------------------
def get_fluxcols(args):
    '''
    Function to load or generate the list of filters and correpsonding flux columns in COSMOS2020 catalog
    Returns list of columns, and optionally the full cosmos2020 dataframe
    '''
    filepath = args.input_dir / 'COSMOS' / 'cosmos_fluxcols.npy'

    if os.path.exists(filepath):
        print(f'Reading flux columns from existing {filepath}')
        fluxcols = np.load(filepath)
        df_cosmos = None
    else:
        print(f'{filepath} does not exist, so preparing the list..')
        df_cosmos = read_COSMOS2020_catalog(args=args, filename=args.input_dir / 'COSMOS' / 'COSMOS2020_CLASSIC_R1_v2.2_p3.fits')

        all_flux_cols = [item for item in df_cosmos.columns if 'FLUX' in item and item != 'FLUX_RADIUS' and 'FLUXERR' not in item]
        filters = [item[:item.find('FLUX')] for item in all_flux_cols]
        fluxcols = [item + 'FLUX_AUTO' if item + 'FLUX_AUTO' in df_cosmos.columns else item + 'FLUX' for item in filters]
        fluxcols = list(dict.fromkeys(fluxcols)) # to remove duplicates
        np.save(filepath, fluxcols)

    return fluxcols, df_cosmos

# -------------------------------------------------------------------------------------------------------
def read_COSMOSWebb_catalog(args=None, filename=None, aperture=1.0):
    '''
    Reads in the zCOSMOS galaxy catalog
    Returns as pandas dataframe
    '''
    aperture_dict = {0.1:0, 0.25:1, 0.5:2, 1.0:3, 1.5:4}
    if filename is None:
        if args is None: input_dir = '/Users/acharyya/Work/astro/passage/passage_data'
        else: input_dir = args.input_dir
        filename = Path(input_dir) / 'COSMOS' / 'COSMOS_Web_for_Ayan_Dec23.fits'

    print(f'Reading in {filename}, might take a while..')
    start_time2 = datetime.now()

    data = fits.open(filename)
    table = Table(data[1].data)

    multi_index_columns = [item for item in table.colnames if len(table[item].shape) > 1]
    single_index_columns = list(set(table.columns) - set(multi_index_columns))

    df  = table[single_index_columns].to_pandas()
    for thiscol in multi_index_columns: df[thiscol] = table[thiscol][:, aperture_dict[aperture]]

    df = df.rename(columns={'ID_SE++':'id', 'RA_DETEC':'ra', 'DEC_DETEC':'dec'})
    print(f'Completed reading COSMOSWebb catalog in {timedelta(seconds=(datetime.now() - start_time2).seconds)}')

    return df

# -------------------------------------------------------------------------------------------------------
def read_COSMOS2020_catalog(args=None, filename=None):
    '''
    Reads in the zCOSMOS galaxy catalog
    Returns as pandas dataframe
    '''
    if filename is None:
        if args is None: input_dir = '/Users/acharyya/Work/astro/passage/passage_data'
        else: input_dir = args.input_dir
        filename = Path(input_dir) / 'COSMOS' / 'COSMOS2020_CLASSIC_R1_v2.2_p3_subsetcolumns.fits'

    if not os.path.exists(filename) or (args is not None and args.clobber_cosmos): make_COSMOS_subset_table(filename, args)

    print(f'Reading in {filename}, might take a while..')
    start_time2 = datetime.now()

    data = fits.open(filename)
    table = Table(data[1].data)
    df  = table.to_pandas()
    df = df.rename(columns={'ID':'id', 'ALPHA_J2000':'ra', 'DELTA_J2000':'dec'})
    print(f'Completed reading COSMOS2020 catalog in {timedelta(seconds=(datetime.now() - start_time2).seconds)}')

    return df

# -------------------------------------------------------------------------------------------------------
def read_zCOSMOS_catalog(args=None, filename=None):
    '''
    Reads in the zCOSMOS galaxy catalog
    Returns as pandas dataframe
    '''
    if filename is None:
        if args is None: input_dir = '/Users/acharyya/Work/astro/passage/passage_data'
        else: input_dir = args.input_dir
        filename = input_dir / 'COSMOS/zCOSMOS-DR3' / 'zCOSMOS_VIMOS_BRIGHT_DR3_CATALOGUE.fits'

    data = fits.open(filename)
    table = Table(data[1].data)
    df  = table.to_pandas()
    df = df.rename(columns={'OBJECT_ID':'id', 'RAJ2000':'ra', 'DEJ2000':'dec'})

    return df

# -------------------------------------------------------------------------------------------------------
def make_COSMOS_subset_table(filename, args):
    '''
    Reads in the massive COSMOS2020 catalog and makes a smaller table with subset of columns and saves it
    '''
    suffix = '_subsetcolumns'
    filename = str(filename)
    if suffix in filename: filename = filename[:filename.find(suffix)] + '.fits'

    # -------determining flux column other columns to extract from df_cosmos-------
    fluxcols, _ = get_fluxcols(args)
    lp_cols_suffix = ['med', 'med_min68', 'med_max68', 'best']
    lp_cols = np.ravel([f'lp_{item}_{suffix}' for item in ['mass', 'SFR', 'sSFR'] for suffix in lp_cols_suffix])
    ez_cols_suffix = ['', '_p160', '_p500', '_p840']
    ez_cols = np.ravel([f'ez_{item}{suffix}' for item in ['mass', 'sfr', 'ssfr'] for suffix in ez_cols_suffix])
    flux_and_err_cols = np.ravel([[item, item.replace('FLUX', 'FLUXERR')] for item in fluxcols])
    cols_to_extract = np.hstack((['ID', 'ALPHA_J2000', 'DELTA_J2000', 'ID_COSMOS2015', 'ez_z_phot', 'lp_MK', 'lp_zBEST'], lp_cols, ez_cols, flux_and_err_cols)).tolist()

    print(f'Trying to read in  {filename}; can take a while..')
    data = fits.open(filename)
    table = Table(data[1].data)
    table_sub = table[cols_to_extract]

    outfilename = str(filename).split('.fits')[0] + suffix + '.fits'
    table_sub.write(outfilename, overwrite=True)
    print(f'Saved subset table as {outfilename}')

# ----------------------------------------------------------------------------------------------------------
def get_crossmatch(df1, df2, sep_threshold=1., df1_idcol='id', df2_idcol='id'):
    '''
    Determines crossmatch between two dataframes df1 and df2
    df1 and df2 should have df1_idcol, and df2_idcol respectively, and they each have columns: ra, dec
    sep_threshold is in arcseconds
    Returns cross matched dataframe with IDs from both dataframes
    '''
    df1_coords = SkyCoord(df1['ra'], df1['dec'], unit='deg')
    df2_coords = SkyCoord(df2['ra'], df2['dec'], unit='deg')
    nearest_id_in_df2, sep_from_nearest_id_in_df2, _ = df1_coords.match_to_catalog_sky(df2_coords)

    df_crossmatch = pd.DataFrame({'df1_id': df1[df1_idcol].values, 'df2_id': df2[df2_idcol].iloc[nearest_id_in_df2].values, 'sep': sep_from_nearest_id_in_df2.arcsec})
    df_crossmatch = df_crossmatch[df_crossmatch['sep'] < sep_threshold]  # separation within XX arcsecond
    df_crossmatch = df_crossmatch.sort_values('sep').drop_duplicates(subset='df2_id', keep='first').reset_index(drop=True)  # to avoid multiple df1 objects being linked to the same df2 object

    return df_crossmatch

# -------------------------------------------------------------------------------------------------------
def split_COSMOS_subset_table_by_par(args):
    '''
    Reads in the subset of columns of COSMOS2020 catalog and splits it into smaller tables with only objects that are overlapping with individual PASSAGE fields
    '''
    # -------reading in the COSMOS2020 (sub)catalog------
    filename = Path(args.input_dir) / 'COSMOS' / 'COSMOS2020_CLASSIC_R1_v2.2_p3_subsetcolumns.fits'
    data = fits.open(filename)
    table_cosmos = Table(data[1].data)

    df_cosmos = table_cosmos.to_pandas()
    df_cosmos = df_cosmos.rename(columns={'ALPHA_J2000':'ra', 'DELTA_J2000':'dec'})

    field_list = [os.path.split(item[:-1])[1] for item in glob.glob(str(args.input_dir / args.drv / 'Par*') + '/')]
    field_list += [f'Par{item:03d}' for item in passage_fields_in_cosmos]
    field_list = list(np.unique(field_list))
    field_list.sort(key=natural_keys)

    for index, thisfield in enumerate(field_list):
        print(f'Starting {index+1} of {len(field_list)} fields..')
        # -------determining path to photometric catalog------
        product_dir = args.input_dir / args.drv / thisfield / 'Products'
        catalog_file = product_dir / f'{thisfield}_photcat.fits'

        if os.path.exists(catalog_file):
            # -------reading in photometric catalog------
            catalog = GTable.read(catalog_file)
            df = catalog['id', 'ra', 'dec'].to_pandas()
            df['passage_id'] = thisfield + '-' + df['id'].astype(str)  # making a unique combination of field and object id

            # -------cross-matching RA/DEC of both catalogs------
            df_crossmatch = get_crossmatch(df, df_cosmos, sep_threshold=1.0, df1_idcol='passage_id', df2_idcol='id')
            df_crossmatch = df_crossmatch.rename(columns={'df1_id': 'passage_id', 'df2_id': 'ID'})

            if len(df_crossmatch) > 0:
                table_crossmatch = Table.from_pandas(df_crossmatch)
                table_cosmos_thisfield = join(table_cosmos, table_crossmatch, keys='ID')

                outfilename = args.input_dir / 'COSMOS' /  args.drv / f'cosmos2020_objects_in_{thisfield}.fits'
                table_cosmos_thisfield.write(outfilename, overwrite=True)
                print(f'Saved subset table as {outfilename}')
            else:
                print(f'No overlapping objects found in field {thisfield}')
        else:
            print(f'{catalog_file} does not exist, so skipping {thisfield}.')

# -------------------------------------------------------------------------------------------------------
def split_COSMOSWebb_table_by_par(args, filename=None):
    '''
    Reads in the COSMOSWebb catalog and splits it into smaller tables with only objects that are overlapping with individual PASSAGE fields
    '''
    # -------reading in the COSMOS2020 (sub)catalog------
    if filename is None: filename = Path(args.input_dir) / 'COSMOS' / 'COSMOS_Web_for_Ayan_Dec24.fits'
    data = fits.open(filename)
    table_cosmos = Table(data[1].data)

    df_cosmos = table_cosmos[['ID_SE++', 'RA_DETEC', 'DEC_DETEC']].to_pandas() # using ID_SE++ instead of ID because it turns out that ID is not unique in the COSMOSWeb catalog
    df_cosmos = df_cosmos.rename(columns={'RA_DETEC':'ra', 'DEC_DETEC':'dec'})

    field_list = [os.path.split(item[:-1])[1] for item in glob.glob(str(args.input_dir / args.drv / 'Par*') + '/')]
    field_list += [f'Par{item:03d}' for item in passage_fields_in_cosmos]
    field_list = list(np.unique(field_list))
    field_list.sort(key=natural_keys)

    # The following method is slower because it involves cross matching all COSMOS objects with..
    # ..objects in ParXXX, but it leads to an accurate list of objects within ParXXX.
    for index, thisfield in enumerate(field_list):
        print(f'Starting {index+1} of {len(field_list)} fields..')
        # -------determining path to photometric catalog------
        product_dir = args.input_dir / args.drv / thisfield / 'Products'
        catalog_file = product_dir / f'{thisfield}_photcat.fits'

        if os.path.exists(catalog_file):
            # -------reading in photometric catalog------
            catalog = GTable.read(catalog_file)
            df = catalog['id', 'ra', 'dec'].to_pandas()
            df['passage_id'] = thisfield + '-' + df['id'].astype(str)  # making a unique combination of field and object id

            # -------cross-matching RA/DEC of both catalogs------
            df_crossmatch = get_crossmatch(df, df_cosmos, sep_threshold=1.0, df1_idcol='passage_id', df2_idcol='ID_SE++')
            df_crossmatch = df_crossmatch.rename(columns={'df1_id': 'passage_id', 'df2_id': 'ID_SE++'})

            if len(df_crossmatch) > 0:
                table_crossmatch = Table.from_pandas(df_crossmatch)
                table_cosmos_thisfield = join(table_cosmos, table_crossmatch, keys='ID_SE++')
                outfilename = args.input_dir / 'COSMOS' /  args.drv / f'cosmoswebb_objects_in_{thisfield}.fits'
                table_cosmos_thisfield.write(outfilename, overwrite=True)
                print(f'Saved subset table as {outfilename}')
            else:
                print(f'No overlapping objects found in field {thisfield}')
        else:
            print(f'{catalog_file} does not exist, so skipping {thisfield}.')

    '''
    # The following method is faster because it does not do cross matching, just takes all objects contained..
    # ..within the region file corresponding to ParXXX, but since the region files are approximate..
    # ..the resulting list of COSMOS objects "contained" within ParXXX is not accurate either.
    for index, thisfield in enumerate(field_list):
        print(f'Starting {index+1} of {len(field_list)} fields..')
        # -------determining path to direct images------
        product_dir = args.input_dir / args.drv / thisfield / 'Products'
        image_files = glob.glob(str(product_dir) + f'/{thisfield}*clear*sci*.fits')

        if len(image_files) > 0:
            image_file = Path(image_files[0])
            hdul = fits.open(image_file)
            source_wcs = pywcs.WCS(hdul[0].header)

            region_dir = image_file.parent.parent / 'Regions'
            region_dir.mkdir(parents=True, exist_ok=True)
            region_file = region_dir / Path(filename.stem + '.reg')
            source_wcs.footprint_to_file(region_file, color='cyan', width=1)

            sky_region = Regions.read(region_file, format='ds9')[0]
            contained_ids = sky_region.contains(SkyCoord(df_cosmos['ra'], df_cosmos['dec'], unit='deg'), pywcs.WCS(hdul[0].header))
            table_contained = table_cosmos[contained_ids]

            if len(table_contained) > 0:
                outfilename = args.input_dir / 'COSMOS' / args.drv / f'cosmoswebb_objects_in_{thisfield}.fits'
                table_contained.write(outfilename, overwrite=True)
                print(f'Saved subset table as {outfilename}')
            else:
                print(f'No overlapping objects found in field {thisfield}')
        else:
                print(f'No image file exists, so skipping {thisfield}.')
    '''

# -------------------------------------------------------------------------------------------------------
def get_passage_filter_dict(args=None, filename=None):
    '''
    Reads PASSAGE spreadsheet and and returns a dictionary with PASSAGE field names and which filters are present in them
    '''
    if filename is None:
        input_dir = Path('/Volumes/Elements/acharyya_backup/Work/astro/passage/passage_data') if args is None else args.input_dir
        filename = input_dir / 'JWST PASSAGE Cycle 1 - Cy1 Executed.csv'

    df = pd.read_csv(filename)
    df = df[['Par#', 'Obs Date', 'Filter']]
    df = df.fillna(method='ffill')
    df = df[~ df['Obs Date'].str.contains('SKIPPED')]

    dictionary = {item: np.unique(df[df['Par#'] == item]['Filter']).tolist() for item in np.unique(df['Par#'])}

    return dictionary

# ------------------------------------------------------------------------
def setup_plots_for_talks():
    '''
    Function to setup plto themes etc for talks
    '''
    plt.style.use('cyberpunk')
    background_for_talks = 'cyberpunk'  # 'dark_background' #'Solarize_Light2' #
    plt.style.use(background_for_talks)
    new_foreground_color = '#FFF1D0'
    #plt.rcParams['grid.color'] = new_foreground_color
    plt.rcParams['text.color'] = new_foreground_color
    plt.rcParams['xtick.color'] = new_foreground_color
    plt.rcParams['ytick.color'] = new_foreground_color
    plt.rcParams['xtick.color'] = new_foreground_color
    plt.rcParams['axes.titlecolor'] = new_foreground_color
    plt.rcParams['axes.labelcolor'] = new_foreground_color
    plt.rcParams['axes.edgecolor'] = new_foreground_color
    plt.rcParams['figure.edgecolor'] = new_foreground_color
    plt.rcParams['savefig.edgecolor'] = new_foreground_color
    plt.rcParams['axes.linewidth'] = 2

    new_background_color = '#120000'
    plt.rcParams['axes.facecolor'] = new_background_color
    plt.rcParams['figure.facecolor'] = new_background_color
    plt.rcParams['savefig.facecolor'] = new_background_color
    plt.rcParams['grid.alpha'] = 0.5
    plt.rcParams['grid.linewidth'] = 0.3

# ------------------------------------------------------------------------------------------------------
def get_files_in_url(url, ext='fits', auth=None):
    '''
    Function to get list of fits files in a given URL
    Returns list of urls of the files
    From https://stackoverflow.com/questions/11023530/python-to-list-http-files-and-directories
    '''
    if auth is None: page = requests.get(url).text
    else: page = requests.get(url, auth=auth).text
    soup = BeautifulSoup(page, 'html.parser')
    url_list = [url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]

    return url_list

# ------------------------------------------------------------------------------------------------------
def download_files_from_url(url, outdir, ext='fits', match_strings=[''], auth=None):
    '''
    Function to download all fits files that contain <amch_string> in the filename from a given URL
    Saves the downloaded files in outdir
    '''
    start_time = datetime.now()
    outdir = Path(outdir)
    url_list = get_files_in_url(url, ext=ext, auth=auth)
    len_orig = len(url_list)
    for match_string in match_strings: url_list = [item for item in url_list if match_string in item]
    print(f'Found {len(url_list)} matching files out of {len_orig} at URL {url}')

    for index, this_url in enumerate(url_list):
        target_file = os.path.split(this_url)[-1]
        if os.path.isfile(outdir / target_file) or os.path.isfile(outdir / Path(target_file).stem):
            print(f'Skipping file {target_file} ({index + 1} of {len(url_list)}) because it already exists at target location')
        else:
            start_time2 = datetime.now()
            print(f'Downloading file {target_file} which is {index + 1} of {len(url_list)}..')
            if auth is None: response = requests.get(this_url)
            else: response = requests.get(this_url, auth=auth)
            with open(outdir / target_file, 'wb') as file: file.write(response.content)
            if target_file.endswith('.gz'):
                print(f'Unzipping downloaded file..')
                unzip_and_delete(outdir / target_file, outdir)
            print(f'Completed this download in {timedelta(seconds=(datetime.now() - start_time2).seconds)}')


    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

# ------------------------------------------------------------------------------------------------------
def get_sky_region_from_fits_header(header, CDELT1='CD1_1', CDELT2='CD2_2', ORIENT='ORIENTAT'):
    '''
    Function to make an astro RectanguleSkyRegion from a given fits header
    Returns SkyRegion
    '''
    center_ra = header['CRVAL1'] + (header['NAXIS1']/2 - header['CRPIX1']) * header[CDELT1]
    center_dec = header['CRVAL2'] + (header['NAXIS2']/2 - header['CRPIX2']) * header[CDELT2]
    width = np.abs(header[CDELT1]) * header['NAXIS1']
    height = np.abs(header[CDELT2]) * header['NAXIS2']
    angle = header[ORIENT] if ORIENT in header else 0.

    sky_center = SkyCoord(center_ra, center_dec, unit='deg')
    sky_region = RectangleSkyRegion(center=sky_center, width=width * u.deg, height=height * u.deg, angle=angle * u.deg)

    return sky_region

# ------------------------------------------------------------------------------------------------------
def is_point_in_region(sky_coord, data, CDELT1='CD1_1', CDELT2='CD2_2', ORIENT='ORIENTAT'):
    '''
    Function to check if an input sky coordinate lies within the footprint of a given fits file
    Returns True/False
    '''
    if type(data) == str: # in case the input is the filename
        data = fits.open(data)
    header = data[0].header
    sky_region = get_sky_region_from_fits_header(header, CDELT1=CDELT1, CDELT2=CDELT2, ORIENT=ORIENT)

    contains = sky_region.contains(sky_coord, pywcs.WCS(header))

    return contains

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
def get_distance_from_line(xdata, ydata, func, method='K01'):
    '''
    Computes distance of each object in the given xdata and ydata (line ratios) arrays, from a given line func(x)
    Returns the distance as an array
    '''
    print(f'Computing distance form Kewley+2001 line on the BPT diagram..')
    x = np.linspace(-2, 1, 100)
    y = func(x, method)

    min_dist_arr = []
    for P in zip(xdata.flatten(), ydata.flatten()):
        min_idxs, distances = min_distance(x, y, P)
        if len(min_idxs) > 0: min_dist = distances[min_idxs[0]]
        else: min_dist = np.nan
        min_dist_arr.append(min_dist)

    return np.reshape(min_dist_arr, np.shape(xdata))

# -----------------------------------------------------------------
def rebin(array, dimensions=None, scale=None):
    """ Return the array ``array`` to the new ``dimensions`` conserving flux the flux in the bins
    The sum of the array will remain the same

    >>> ar = numpy.array([
        [0,1,2],
        [1,2,3],
        [2,3,4]
        ])
    >>> rebin(ar, (2,2))
    array([
        [1.5, 4.5]
        [4.5, 7.5]
        ])
    Raises
    ------

    AssertionError
        If the totals of the input and result array don't agree, raise an error because computation may have gone wrong

    Reference
    =========
    +-+-+-+
    |1|2|3|
    +-+-+-+
    |4|5|6|
    +-+-+-+
    |7|8|9|
    +-+-+-+
    """
    if dimensions is not None:
        if isinstance(dimensions, float):
            dimensions = [int(dimensions)] * len(array.shape)
        elif isinstance(dimensions, int):
            dimensions = [dimensions] * len(array.shape)
        elif len(dimensions) != len(array.shape):
            raise RuntimeError('')
    elif scale is not None:
        if isinstance(scale, float) or isinstance(scale, int):
            dimensions = map(int, map(round, map(lambda x: x * scale, array.shape)))
        elif len(scale) != len(array.shape):
            raise RuntimeError('')
    else:
        raise RuntimeError('Incorrect parameters to rebin.\n\trebin(array, dimensions=(x,y))\n\trebin(array, scale=a')
    if np.shape(array) == dimensions: return array  # no rebinning actually needed
    import itertools
    # dY, dX = map(divmod, map(float, array.shape), dimensions)

    result = np.zeros(dimensions)
    for j, i in itertools.product(*map(range, array.shape)):
        (J, dj), (I, di) = divmod(j * dimensions[0], array.shape[0]), divmod(i * dimensions[1], array.shape[1])
        (J1, dj1), (I1, di1) = divmod(j + 1, array.shape[0] / float(dimensions[0])), divmod(i + 1,
                                                                                            array.shape[1] / float(
                                                                                                dimensions[1]))

        # Moving to new bin
        # Is this a discrete bin?
        dx, dy = 0, 0
        if (I1 - I == 0) | ((I1 - I == 1) & (di1 == 0)):
            dx = 1
        else:
            dx = 1 - di1
        if (J1 - J == 0) | ((J1 - J == 1) & (dj1 == 0)):
            dy = 1
        else:
            dy = 1 - dj1
        # Prevent it from allocating outide the array
        I_ = np.min([dimensions[1] - 1, I + 1])
        J_ = np.min([dimensions[0] - 1, J + 1])
        result[J, I] += array[j, i] * dx * dy
        result[J_, I] += array[j, i] * (1 - dy) * dx
        result[J, I_] += array[j, i] * dy * (1 - dx)
        result[J_, I_] += array[j, i] * (1 - dx) * (1 - dy)
    allowError = 0.1
    if array.sum() > 0: assert (array.sum() < result.sum() * (1 + allowError)) & (array.sum() > result.sum() * (1 - allowError))
    return result

# --------------------------------------------------------------------------------------------------------------------
def get_custom_cmap(cmap_name, cmap_path=None):
    '''
    Computes custom colormaps, code borrowed from Eduardo Vitral

    SET COLORMAPS FROM:
    Crameri, F., G.E. Shephard, and P.J. Heron (2020)
    The misuse of colour in science communication,
    Nature Communications, 11, 5444.
    https://doi.org/10.1038/s41467-020-19160-7
    ---
    Crameri, F.: Geodynamic diagnostics, scientific visualisation
    and StagLab 3.0, Geosci. Model Dev. Discuss.,
    https://doi.org/10.5194/gmd-2017-328, 2018.

    Returns the colormap
    '''

    if cmap_path is None: cmap_path = HOME / 'Work/astro/ayan_codes/ScientificColourMaps8'
    cpal_quant = np.loadtxt(cmap_path / cmap_name / f'{cmap_name}.txt')
    mycmap_quant = mplcolors.ListedColormap(cpal_quant, name=cmap_name)

    return mycmap_quant

# --------------------------------------------------------------------------------------------------------------------
def get_combined_cmap(breaks, cmaps, new_name='my_colormap'):
    '''
    Combines an arbitrary number of matplotlib colormaps at given break points
    Returns a new colormap
    Adapted from https://stackoverflow.com/questions/31051488/combining-two-matplotlib-colormaps
    '''
    colors = []
    for index, cmap in enumerate(cmaps):
        ncolors = int(256 * (breaks[index + 1] - breaks[index]))
        this_colors = mplcolormaps[cmap](np.linspace(breaks[index], breaks[index + 1], ncolors))
        colors.append(this_colors)
    colors = np.vstack(colors)
    new_cmap = mplcolors.LinearSegmentedColormap.from_list(new_name, colors)

    return new_cmap

# --------------------------------------------------------------------------------------------------------------------
def unzip_and_delete(zip_file, destination_path):
    '''
    Unzips given zip file to given destination path and removes the zip file
    '''
    zip_file = Path(zip_file)
    print(f'Unpacking {zip_file} in to {destination_path}..')
    if zip_file.suffix == '.gz':
        with gzip.open(zip_file, 'rb') as gz_file:
            with open(zip_file.parent / zip_file.stem, 'wb') as out_file:
                shutil.copyfileobj(gz_file, out_file)
    else:
        shutil.unpack_archive(zip_file, destination_path)
    os.remove(zip_file)  # remove zipped files after unzipping

# --------------------------------------------------------------------------------------------------------------------
def move_field_after_download(field_arr, args=None):
    '''
    Unzips and moves the reduced data from Box from Downloads folder to appropriate location by making the right dorectories
    '''
    origin_dir = HOME / 'Downloads'
    field_arr = np.atleast_1d(field_arr)

    for index, field in enumerate(field_arr):
        start_time = datetime.now()
        if 'Par' in str(field): field = f'Par{int(field.split("Par")[1]):03d}'
        else: field = f'Par{field:03d}'
        print(f'Starting field {field} which is {index + 1} out of {len(field_arr)}..')

        try:
            if args is None: destination_dir = Path('/Volumes/Elements/acharyya_backup/Work/astro/passage/passage_data/v0.5/data')
            else: destination_dir = args.input_dir / 'v0.5/data'
            Path(destination_dir / 'linelist').mkdir(exist_ok=True, parents=True)

            unzip_and_delete(origin_dir / f'{field}.zip', destination_dir)

            Path(destination_dir / field / 'DATA/DIRECT_GRISM').mkdir(exist_ok=True, parents=True)
            unzip_and_delete(destination_dir / field / f'{field}_spec1D.tar.gz', destination_dir / field)
            unzip_and_delete(destination_dir / field / f'{field}_spec2D.tar.gz', destination_dir / field)

            cat_files = glob.glob(str(destination_dir / field / '*cat.fits'))
            for this_cat_file in cat_files:
                print(f'Moving {os.path.split(this_cat_file)[-1]}..')
                shutil.move(this_cat_file, destination_dir / field / 'DATA/DIRECT_GRISM')

            fits_files = glob.glob(str(destination_dir / field / '*.fits'))
            for this_fits_file in fits_files:
                print(f'Moving {os.path.split(this_fits_file)[-1]}..')
                shutil.move(this_fits_file, destination_dir / field / 'DATA')

            filename = f'{field}lines.dat'
            print(f'Moving {filename}..')
            shutil.move(destination_dir / field / filename, destination_dir / 'linelist')
        except Exception as e:
            print(f'Probably skipping {field} due to following error: {e}')
            continue

        print(f'Completed moving {field} in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

# --------------------------------------------------------------------------------------------------
def adjust_lightness(color, amount=0.5):
    '''
    Lightens (if amount > 1) or darkens (if amount < 1) matplitlib colors
    Adapted from https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
    '''
    c = colorsys.rgb_to_hls(*mplcolors.to_rgb(color))
    color = colorsys.hls_to_rgb(c[0], max(0.1, min(0.9, amount * c[1])), c[2])
    return color

# --------------------------------------------------------------------------------------------------
args = parse_args()
available_filters_for_field_dict = get_passage_filter_dict(args)