from astropy.io import fits
from astropy.table import Table, vstack, join, hstack
import numpy as np
from pathlib import Path
from os import PathLike
from astropy.wcs import WCS
import astropy.visualization as astrovis
import os
from matplotlib.ticker import MultipleLocator

# root_dir = Path(os.getenv("ROOT_DIR"))
# root_dir = Path("/media/watsonp/ArchivePJW/backup/data")
root_dir = Path("/Users/acharyya/Work/astro/passage/passage_data/vpjw")
field = "Par682"
date = "2026_02_18"

field_name = f"passage-{field.lower()}"

extractions_dir = Path(root_dir) / f"{field}"

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


aanda_columnwidth = 256.0748 / 72.27
aanda_textwidth = 523.5307 / 72.27


import cmcrameri.cm as cmc


def read_quant_maps_fits(filename):
    """
    Reads a multi-extension FITS file and reconstructs the quant_maps dictionary.
    Returns:
    logOH_map and sfr_map: 2D masked arrays
    extent: to be used for plt.imshow()
    """
    params_dict = {
        "logOH": {"data": "LOG_OH", "err": "LOG_OH_ERR", "mask": "LOG_OH_MASK"},
        "sfr": {"data": "SFR", "err": "SFR_ERR", "mask": "SFR_MASK"},
    }
    quant_maps = {"logOH": np.nan, "sfr": np.nan}

    print(f"Reading existing maps fits file from {filename}")
    hdul = fits.open(filename)

    for label in list(quant_maps.keys()):
        hdu_data = hdul[params_dict[label]["data"]]
        # hdu_err = hdul[params_dict[label]['err']]
        hdu_mask = hdul[params_dict[label]["mask"]]

        this_quant = np.ma.masked_where(hdu_mask.data, hdu_data.data)
        quant_maps.update({label: this_quant})

    logOH_map = quant_maps["logOH"]
    sfr_map = quant_maps["sfr"]

    left = (
        hdu_data.header["CRVAL1"]
        + (0.5 - hdu_data.header["CRPIX1"]) * hdu_data.header["CDELT1"]
    )
    right = (
        hdu_data.header["CRVAL1"]
        + (hdu_data.header["NAXIS1"] + 0.5 - hdu_data.header["CRPIX1"])
        * hdu_data.header["CDELT1"]
    )
    bottom = (
        hdu_data.header["CRVAL2"]
        + (0.5 - hdu_data.header["CRPIX2"]) * hdu_data.header["CDELT2"]
    )
    top = (
        hdu_data.header["CRVAL2"]
        + (hdu_data.header["NAXIS2"] + 0.5 - hdu_data.header["CRPIX2"])
        * hdu_data.header["CDELT2"]
    )
    extent = (left, right, top, bottom)

    return logOH_map, sfr_map, extent


def setup_aanda_style(dark: bool = False):
    """
    A helper function to setup the A&A style.

    Parameters
    ----------
    dark : bool, optional
        Use a dark plotting style, by default `False`.
    """

    import matplotlib as mpl

    rc_params = {
        "font.family": "serif",
        "font.size": 7,
        "figure.figsize": (aanda_columnwidth, 3),
        "text.usetex": True,
        "ytick.right": True,
        "ytick.direction": "in",
        "ytick.minor.visible": True,
        "ytick.labelsize": 8,
        "xtick.top": True,
        "xtick.direction": "in",
        "xtick.minor.visible": True,
        "xtick.labelsize": 8,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "legend.fontsize": 8,
        "lines.linewidth": 1.0,
        "lines.markersize": 3.5,
        "image.interpolation": "none",
        "text.latex.preamble": (
            r"""
        \usepackage{amsmath}
        \usepackage{txfonts}
        \usepackage{siunitx}
        %
        \DeclareMathAlphabet{\mathsc}{OT1}{cmr}{m}{sc}
        \def\testbx{bx}%
        \DeclareRobustCommand{\ion}[2]{%
        \relax\ifmmode
        \ifx\testbx\f@series
        {\mathbf{#1\,\mathsc{#2}}}\else
        {\mathrm{#1\,\mathsc{#2}}}\fi
        \else\textup{#1\,{\mdseries\textsc{#2}}}%
        \fi}
        %
        """
        ),
    }

    if dark:
        rc_fonts |= {
            "text.color": "white",
            "axes.facecolor": "0C1C23",  # axes background color
            "axes.edgecolor": "#e1e9ec",  # axes edge color
            "axes.labelcolor": "#e1e9ec",
            "grid.color": "#e1e9ec",
            "legend.edgecolor": "#e1e9ec",
            "legend.facecolor": "inherit",
            "legend.labelcolor": "#e1e9ec",
            "xtick.color": "#e1e9ec",
            "ytick.color": "#e1e9ec",
            "figure.facecolor": (0.0, 0.0, 0.0, 0.0),
            "savefig.facecolor": (0.0, 0.0, 0.0, 0.0),
        }

    plt.rcdefaults()
    mpl.rcParams.update(rc_params)

    return


setup_aanda_style()

obj_id = 3174

with fits.open(
    extractions_dir / "stack" / f"passage-par682_{obj_id:0>5}.stack.fits"
) as stack_hdul:
    # stack_hdul.info()

    plot_kwargs = dict(
        origin="lower",
        # cmap="binary",
        cmap=cmc.vik,
    )

    fig = plt.figure(constrained_layout=True, figsize=(0.6 * aanda_textwidth, 2.5))

    subfigs = fig.subfigures(2, 1, height_ratios=[1, 1.2])

    axs_maps = subfigs[1].subplots(1, 3, sharey=True)

    with fits.open(
        extractions_dir / "full" / f"passage-par682_{obj_id:0>5}.full.fits"
    ) as full_hdul:

        obj_z = full_hdul[0].header["REDSHIFT"]

        # axs_maps[1].imshow(np.log10(sfr))

        # axs_maps[2].imshow(Z_map)

        # quant_slice = (slice(0, None, None), slice(-1, None, -1))
        quant_slice = (slice(0, None, None), slice(-1, None, -1))

        logOH_map, sfr_map, extent = read_quant_maps_fits(
            #extractions_dir / "AA_quants" / f"passage-par682_{obj_id:0>5}_quants.fits"
            extractions_dir / "quants" / f"passage-par682_{obj_id:0>5}_quants.fits"
        )
        print (extent)
        extent = (extent[0], extent[1], -extent[2], -extent[3])
        sfr = axs_maps[1].imshow(
            sfr_map[quant_slice],
            extent=extent,
            cmap="plasma",
            origin="lower",
        )
        clb_sfr = plt.colorbar(sfr, ax=axs_maps[1])
        # clb_sfr.ax.set_title(r'log\,SFR')
        bbox_kwargs = dict(
            facecolor="white", edgecolor="k", boxstyle="round", alpha=0.9
        )
        # annotate_kwargs = dict(
        #     xy=(0.95, 0.95),
        #     xycoords="axes fraction",
        #     ha="right",
        #     va="top",
        #     bbox=bbox_kwargs,
        # )
        annotate_kwargs = dict(
            xy=(0.05, 0.05),
            xycoords="axes fraction",
            ha="left",
            va="bottom",
            bbox=bbox_kwargs,
        )
        axs_maps[1].annotate(r"log\,SFR", **annotate_kwargs)

        logOH = axs_maps[2].imshow(
            logOH_map[quant_slice],
            extent=extent,
            cmap="viridis",
            origin="lower",
        )
        plt.colorbar(logOH, ax=axs_maps[2])
        axs_maps[2].annotate(r"log\,(O/H)+12 [R23]", **annotate_kwargs)

        from astrocolour import ColourImage

        img_crop = (
            # slice(-1,None, -1),
            # slice(0,None, None),
            # slice(-1,None, -1)
            slice(17, -17),
            slice(-17, 17, -1),
        )

        axs_maps[0].imshow(
            ColourImage(
                [
                    full_hdul["DSCI", f"{f}-CLEAR"].data[img_crop]
                    for f in ["F115W", "F150W", "F200W"]
                ]
            ),
            origin="lower",
            extent=extent,
        )

        wcs = WCS(full_hdul["DSCI", "F115W-CLEAR"])

        from astropy.cosmology import LambdaCDM
        import astropy.units as u

        cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        arcsec_kpc_scale = cosmo.arcsec_per_kpc_proper(obj_z)

        # print(50 * (arcsec_kpc_scale / wcs.proj_plane_pixel_scales()[0]))
        print(wcs.proj_plane_pixel_scales()[0].to(u.arcsec))
        # print(50 * arcsec_kpc_scale * u.kpc)
        # print (1/arcsec_kpc_scale)
        print(arcsec_kpc_scale / wcs.proj_plane_pixel_scales()[0].to(u.arcsec))

        len_kpc = 1*u.kpc

        print (arcsec_kpc_scale / wcs.proj_plane_pixel_scales()[0].to(u.arcsec)*len_kpc)


        scale_bar_coords = np.array([-.3,0.3])
        axs_maps[-1].errorbar(
            *scale_bar_coords,
            c="k",
            fmt="none",
            # markersize=10,
            zorder=10,
            xerr=(arcsec_kpc_scale*len_kpc).value,
            capsize=2
        )


        axs_maps[-1].annotate(rf"{int(2*len_kpc.value)}\,kpc", (scale_bar_coords[0],scale_bar_coords[1]*1.1), ha="center", fontsize=8)
        
    axs_spec = subfigs[0].subplots(1, 2, width_ratios=[98, 138], sharey=True)

    crop_size = 20
    crop = (
        # slice(10,-10),
        slice(crop_size + 3, -crop_size),
        slice(0, None),
    )
    
    # Include padding for lines here
    label_dict = {
        r"\ion{H}{\alpha}": 0.655,
        r"\ion{H}{\beta}": 0.483,
        r"\ion{H}{\gamma}": 0.436,
        r"[\ion{O}{iii}]": 0.505,
        r"[\ion{S}{ii}]": 0.673
    }

    for ax, filt_name in zip(
        axs_spec,
        # [ax_f150w, ax_f200w],
        ["F150W", "F200W"],
    ):

        stack_data = stack_hdul["SCI", filt_name].data[crop]

        header = stack_hdul["SCI", filt_name].header

        extent = [header["WMIN"], header["WMAX"], 0, stack_data.shape[0]]
        # print(extent)
        ax.set_xlabel(r"$\lambda_{\rm{obs}}$ (\textmu m)")
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.set_yticklabels([])
        # ax.set_xticklabels([])

        ax.imshow(
            # stack_hdul["SCI", "F200W"].data[crop],
            stack_data,
            **plot_kwargs,
            norm=astrovis.ImageNormalize(
                stretch=astrovis.LinearStretch(),
                # vmin=-0.005,
                # vmax=np.nanpercentile(stack_data, 99.9),
                vmax=np.nanpercentile(stack_data, 99.9),
                vmin=-0.5 * np.nanpercentile(stack_data, 99.9),
            ),
            aspect="auto",
            extent=extent,
        )
        for lab, rf_lam in label_dict.items():
            if (rf_lam*(1+obj_z)>extent[0]) and (rf_lam*(1+obj_z)<extent[1]):
                print (lab)

                ax.annotate(rf"${lab}$", (rf_lam*(1+obj_z), extent[3]*1.05),annotation_clip=False, ha="center")


    plt.savefig(extractions_dir / "figs" / "2d_grism_spec_3714.pdf")


    plt.show()
    plt.close()
