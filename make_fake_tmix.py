'''
    Filename: make_fake_tmix.py
    Notes: To plot a cartoon/dummy plot for tmix vs Z-SFR slope, for JWST Cycle 5 proposal submission
    Author : Ayan
    Created: 15-10-25
    Example: run make_fake_tmix.py
'''

from matplotlib import pyplot as plt
import numpy as np

plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = True
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['ytick.major.size'] = 4.5
plt.rcParams['ytick.minor.size'] = 2
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['ytick.minor.width'] = 1.5

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['xtick.major.size'] = 4.5
plt.rcParams['xtick.minor.size'] = 2
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['xtick.minor.width'] = 1.5

plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 1.5

# -------------------------------------------------------------------------------------------------------
if __name__ == "__main__":

    # fake data
    logZ_logSFR_slope = np.linspace(-1.5, 1.5, 10)
    #log_t_dyn = 6 - 0.2 * (logZ_logSFR_slope - 1.7)**2
    log_t_dyn = 5.5 + 0.8 * logZ_logSFR_slope
    log_t_dyn_u = 1
    log_t_turb = 1.2 + 0.9 * log_t_dyn
    log_t_turb_u = 1

    color1, color2, fontsize = 'salmon', 'sienna', 25

    fig, ax = plt.subplots(1, figsize=(7, 7))
    fig.subplots_adjust(left=0.15, right=0.98, bottom=0.15, top=0.98)
                        
    ax.plot(logZ_logSFR_slope, log_t_dyn, color=color1, lw=2, ls='--', label=r'$t_{\mathrm{dyn}}$ calibration')
    ax.plot(logZ_logSFR_slope, log_t_turb, color=color2, lw=2,  ls='--', label=r'$t_{\mathrm{turb}}$ calibration')
    #ax.fill_between(logZ_logSFR_slope, log_t_dyn - log_t_dyn_u, log_t_dyn + log_t_dyn_u, color=color1, alpha=0.3)   

    ax.set_xlabel(r'Slope of $\log$ (O/H) + 12 vs $\log$ $\Sigma_{\rm SFR}$', fontsize=fontsize, fontdict={'color':'royalblue'})
    ax.set_ylabel(r'$\log$ $t_{dyn}$ & $t_{turb}$ (from H$\alpha$ kinematics)', fontsize=fontsize/1.2)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.legend(fontsize=fontsize, loc='upper left')

    ax.set_ylim(3.5, 7)

    fig.savefig('/Users/acharyya/Documents/writings/telescope_proposals/JWSTCy5/Figures/tmix_cartoon.png')
    plt.show(block=False)

