'''
    Filename: make_toy_model.py
    Notes: Makes a toy model of radial profile of emission line ratios and how they get affected by convolution and noise
    Author : Ayan
    Created: 07-03-25
    Example: run make_toy_model.py --drv 0.5 --num_line sii --den_line ha --num_snr 5 --den_snr 2 --res 0.2
             run make_toy_model.py --drv 0.5 --num_snr 5 --den_snr 2 --res 0.1,0.2,0.5
'''

from header import *
from util import *

start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # ------model parameters-----------------------
    # assuming linear profile for lines in the linear scale, so essentially exponential profiles
    line_model_dict = {'nii':{'name':'N II', 'slope':0.3, 'intercept':1.0}, 'sii':{'name':'S II', 'slope':0.2, 'intercept':0.1}, 'ha':{'name':'H alpha', 'slope':1.0, 'intercept':0.8}} # dict for storing each line's name, slope (dex/length), intercept (dex)
    df_lines = pd.DataFrame(line_model_dict)
    num_color, den_color, ratio_color = 'salmon', 'cornflowerblue', 'sienna'
    npix_tot = 100
    np.random.seed(444)

    # ------declare the figure-----------------------
    fig, axes = plt.subplots(1,2, figsize=(10, 6))
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.1, top=0.95, wspace=0.2)
    figname = args.output_dir / 'plots' / f'toy_model_ratio_{args.num_line.lower()}_SNR_{args.num_snr:.0f}-{args.den_line.lower()}_SNR_{args.den_snr:.0f}_PSF_{",". join(str(i) for i in args.res)}.png'
    xarr = np.linspace(0, 1, npix_tot)

    # ------make the ideal model-----------------------
    log_num_profile = df_lines[args.num_line.lower()]['slope'] * xarr + df_lines[args.num_line.lower()]['intercept']
    log_den_profile = df_lines[args.den_line.lower()]['slope'] * xarr + df_lines[args.den_line.lower()]['intercept']
    log_ratio_profile = np.log10((10 ** log_num_profile) / (10 **log_den_profile)) # essentially log_num_profile - log_den_profile

    # ------plot the ideal model-----------------------
    ls, lw, alpha = 'solid', 0.5, 0.5
    axes[0].plot(xarr, log_num_profile, ls=ls, lw=lw, alpha=alpha, c=num_color, label=f'Input {df_lines[args.num_line.lower()]["name"]}; SNR={args.num_snr:.0f}')
    axes[0].plot(xarr, log_den_profile, ls=ls, lw=lw, alpha=alpha, c=den_color, label=f'Input {df_lines[args.den_line.lower()]["name"]}; SNR={args.den_snr:.0f}')
    axes[1].plot(xarr, log_ratio_profile, ls=ls, lw=lw, alpha=alpha, c=ratio_color, label=f'Input {df_lines[args.num_line.lower()]["name"]}/{df_lines[args.den_line.lower()]["name"]}')

    # ------plot the noise limits-----------------------
    axes[0].fill_between(xarr, log_num_profile * (1 - 1 / args.num_snr), log_num_profile * (1 + 1 / args.num_snr), lw=0, alpha=alpha - 0.2, color=num_color)
    axes[0].fill_between(xarr, log_den_profile * (1 - 1 / args.den_snr), log_den_profile * (1 + 1 / args.den_snr), lw=0, alpha=alpha - 0.2, color=den_color)

    ratio_snr = (args.num_snr * args.den_snr) / np.sqrt(args.num_snr ** 2 + args.den_snr ** 2) # for f=x/y, snr_f = (snr_x * snr_y) / sqrt(snr_x**2 + snr_y**2)
    axes[1].fill_between(xarr, log_ratio_profile * (1 - 1 / ratio_snr), log_ratio_profile * (1 + 1 / ratio_snr), lw=0, alpha=alpha - 0.2, color=ratio_color)

    # ------generate noisy gaussian data from the model-----------------------
    log_num_data = np.random.normal(loc=log_num_profile, scale=log_num_profile / args.num_snr) # treating 1-sigma std dev = signal/SNR
    log_den_data = np.random.normal(loc=log_den_profile, scale=log_den_profile / args.den_snr) # treating 1-sigma std dev = signal/SNR
    log_ratio_data = log_num_data - log_den_data

    # ------plot the noisy model-----------------------
    s, lw, alpha = 10, 0.5, 0.8
    axes[0].scatter(xarr, log_num_data, s=s, lw=lw, alpha=alpha, c=num_color)
    axes[0].scatter(xarr, log_den_data, s=s, lw=lw, alpha=alpha, c=den_color)
    axes[1].scatter(xarr, log_ratio_data, s=s, lw=lw, alpha=alpha, c=ratio_color, label=f'Noisy data')

    # ------looping over different resolutions-----------------------
    for index, res in enumerate(args.res):
        print(f'Doing res={res}, which is {index+1} of {len(args.res)}')
        # ------convolve the model-----------------------
        npix_kernel = int(npix_tot * res)
        kernel = Gaussian1DKernel(npix_kernel)
        # kernel = Box1DKernel(npix_kernel)
        log_num_data = convolve(log_num_data, kernel)
        log_den_data = convolve(log_den_data, kernel)
        log_ratio_data = log_num_data - log_den_data

        # ------plot the convolved model-----------------------
        ls, lw, alpha = 'solid', 1, 1
        brightness_factor = (index + 1) / len(args.res)
        axes[0].plot(xarr, log_num_data, ls=ls, lw=lw, alpha=alpha, c=adjust_lightness(num_color, amount=brightness_factor))
        axes[0].plot(xarr, log_den_data, ls=ls, lw=lw, alpha=alpha, c=adjust_lightness(den_color, amount=brightness_factor), label=f'PSF {res:.1f}')
        axes[1].plot(xarr, log_ratio_data, ls=ls, lw=lw, alpha=alpha, c=adjust_lightness(ratio_color, amount=brightness_factor), label=f'PSF {res:.1f}')

    # ------annotate figure-----------------------
    axes[0].legend(fontsize=args.fontsize)
    axes[0].set_xlabel('Normalised distance', fontsize=args.fontsize)
    axes[0].set_ylabel('Log line flux', fontsize=args.fontsize)
    axes[0].tick_params(axis='both', which='major', labelsize=args.fontsize)

    axes[1].legend(fontsize=args.fontsize)
    axes[1].set_xlabel('Normalised distance', fontsize=args.fontsize)
    axes[1].set_ylabel('Log line ratio', fontsize=args.fontsize)
    axes[1].tick_params(axis='both', which='major', labelsize=args.fontsize)

    # --------for talk plots--------------
    if args.fortalk:
        try: mplcyberpunk.make_scatter_glow()
        except: pass

    fig.savefig(figname, transparent=args.fortalk)
    print(f'\nSaved figure as {figname}')
    plt.show(block=False)

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
