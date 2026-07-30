'''
    Filename: make_toy_model.py
    Notes: Makes a toy model of radial profile of emission line ratios and how they get affected by convolution and noise
    Author : Ayan
    Created: 07-03-25
    Example: run make_toy_model.py --drv 0.5 --num_line sii --den_line ha --num_snr 5 --den_snr 2 --res 0.2
             run make_toy_model.py --drv 0.5 --num_snr 2 --den_snr 5 --res 0.1,0.2,0.5
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ------model parameters-----------------------
    # assuming exponential profiles, with central value (intercept) and scale length (re) in arbitrary units
    line_model_dict = {'ha':{'name':'H alpha', 'intercept':1e-18, 're':0.50},
                       'nii':{'name':'N II', 'intercept':5e-19, 're':0.35},
                       'sii':{'name':'S II', 'intercept':1e-19, 're':0.30}} # dict for storing each line's name, slope (dex/length), intercept (dex)
    df_lines = pd.DataFrame(line_model_dict)
    num_color, den_color, ratio_color = 'salmon', 'cornflowerblue', 'sienna'
    xlim = 1
    npix_tot = 100 * xlim
    #np.random.seed(444)

    # ------declare the figure-----------------------
    num_snr_text = f'_SNR_{args.num_snr:.0f}' if args.num_snr is not None else ''
    den_snr_text = f'_SNR_{args.den_snr:.0f}' if args.den_snr is not None else ''
    fig, axes = plt.subplots(1,2, figsize=(12, 6))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.1, top=0.95, wspace=0.2)
    figname = args.output_dir / 'plots' / f'toy_model_ratio_{args.num_line.lower()}{num_snr_text}-{args.den_line.lower()}{den_snr_text}_PSF_{",". join(str(i) for i in args.res)}.png'
    xarr = np.linspace(0, xlim, npix_tot)

    # ------make the ideal model-----------------------
    num_profile = df_lines[args.num_line.lower()]['intercept'] * np.exp(- xarr / df_lines[args.num_line.lower()]['re'])
    den_profile = df_lines[args.den_line.lower()]['intercept'] * np.exp(- xarr / df_lines[args.den_line.lower()]['re'])
    ratio_profile = num_profile / den_profile

    # ------plot the ideal model-----------------------
    ls, lw, alpha = 'dashed', 2, 0.5
    axes[0].plot(xarr, num_profile, ls=ls, lw=lw, alpha=alpha, c=num_color, label=f'Input {df_lines[args.num_line.lower()]["name"]}{num_snr_text}')
    axes[0].plot(xarr, den_profile, ls=ls, lw=lw, alpha=alpha, c=den_color, label=f'Input {df_lines[args.den_line.lower()]["name"]}{den_snr_text}')    
    axes[1].plot(xarr, ratio_profile, ls=ls, lw=lw, alpha=alpha, c=ratio_color, label=f'Input {df_lines[args.num_line.lower()]["name"]}/{df_lines[args.den_line.lower()]["name"]}')

    # ------plot the noise limits-----------------------
    if args.num_snr is not None and args.den_snr is not None:
        if args.num_snr > 1: axes[0].fill_between(xarr, np.log10(num_profile * (1 - 1 / args.num_snr)), np.log10(num_profile * (1 + 1 / args.num_snr)), lw=0, alpha=alpha - 0.2, color=num_color)
        else: axes[0].fill_between(xarr, -30, np.log10(num_profile * (1 + 1 / args.num_snr)), lw=0, alpha=alpha - 0.2, color=num_color)
        
        if args.den_snr > 1: axes[0].fill_between(xarr, np.log10(den_profile * (1 - 1 / args.den_snr)), np.log10(den_profile * (1 + 1 / args.den_snr)), lw=0, alpha=alpha - 0.2, color=den_color)
        else: axes[0].fill_between(xarr, -30, np.log10(den_profile * (1 + 1 / args.den_snr)), lw=0, alpha=alpha - 0.2, color=den_color)

        ratio_snr = (args.num_snr * args.den_snr) / np.sqrt(args.num_snr ** 2 + args.den_snr ** 2) # for f=x/y, snr_f = (snr_x * snr_y) / sqrt(snr_x**2 + snr_y**2)
        if ratio_snr > 1: axes[1].fill_between(xarr, np.log10(ratio_profile * (1 - 1 / ratio_snr)), np.log10(ratio_profile * (1 + 1 / ratio_snr)), lw=0, alpha=alpha - 0.2, color=ratio_color)
        else: axes[1].fill_between(xarr, -3, np.log10(ratio_profile * (1 + 1 / ratio_snr)), lw=0, alpha=alpha - 0.2, color=ratio_color)

        # ------generate noisy gaussian data from the model-----------------------
        num_data = np.random.normal(loc=num_profile, scale=num_profile / args.num_snr) # treating 1-sigma std dev = signal/SNR
        den_data = np.random.normal(loc=den_profile, scale=den_profile / args.den_snr) # treating 1-sigma std dev = signal/SNR
        ratio_data = num_data / den_data
        
        # ------plot the noisy model-----------------------
        ls, lw, alpha = 'solid', 0.5, 0.8
        axes[0].plot(xarr, np.log10(num_data), ls=ls, lw=lw, alpha=alpha, c=num_color)
        axes[0].plot(xarr, np.log10(den_data), ls=ls, lw=lw, alpha=alpha, c=den_color)
        axes[1].plot(xarr, np.log10(ratio_data), ls=ls, lw=lw, alpha=alpha, c=ratio_color, label=f'Noisy data')

    # ------looping over different resolutions-----------------------
    for index, res in enumerate(args.res):
        #res = res * df_lines['ha']['re'] # res are in units of Halpha Re, converting them to distance units
        print(f'Doing res={res}, which is {index+1} of {len(args.res)}')
        # ---------------convolve the model-----------------------
        npix_kernel = int(npix_tot * res)
        if False:
            num_smoothed = convolve(num_data, Gaussian1DKernel(npix_kernel))
            den_smoothed = convolve(den_data, Gaussian1DKernel(npix_kernel))
        elif False:
            mode = 'nearest'
            num_smoothed = gaussian_filter1d(num_data, npix_kernel, mode=mode)
            den_smoothed = gaussian_filter1d(den_data, npix_kernel, mode=mode)
        else:
            mode = 'nearest'
            num_smoothed = gaussian_filter1d(num_profile, npix_kernel, mode=mode)
            den_smoothed = gaussian_filter1d(den_profile, npix_kernel, mode=mode)
        ratio_smoothed = num_smoothed / den_smoothed

        # -------------plot the convolved model-----------------------
        ls, lw, alpha = 'solid', 1, 1
        brightness_factor = (index + 1) / len(args.res)
        axes[0].plot(xarr, num_smoothed, ls=ls, lw=lw, alpha=alpha, c=adjust_lightness(num_color, amount=brightness_factor))
        axes[0].plot(xarr, den_smoothed, ls=ls, lw=lw, alpha=alpha, c=adjust_lightness(den_color, amount=brightness_factor), label=f'PSF {res:.1f}')
        axes[1].plot(xarr, ratio_smoothed, ls=ls, lw=lw, alpha=alpha, c=adjust_lightness(ratio_color, amount=brightness_factor), label=f'PSF {res:.1f}')

    # ------annotate figure-----------------------
    #axes[0].set_yscale('log')
    axes[1].set_yscale('log')

    axes[0].legend(fontsize=args.fontsize)
    axes[0].set_xlabel('Normalised distance', fontsize=args.fontsize)
    axes[0].set_ylabel('Line flux', fontsize=args.fontsize)
    axes[0].tick_params(axis='both', which='major', labelsize=args.fontsize)
    axes[0].set_xlim(0, 1)

    axes[1].legend(fontsize=args.fontsize)
    axes[1].set_xlabel('Normalised distance', fontsize=args.fontsize)
    axes[1].set_ylabel('Line ratio', fontsize=args.fontsize)
    axes[1].tick_params(axis='both', which='major', labelsize=args.fontsize)
    axes[1].set_xlim(0, 1)

    # --------for talk plots--------------
    if args.fortalk:
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    fig.savefig(figname, transparent=args.fortalk)
    print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
