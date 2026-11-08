#!/usr/bin/env python
"""Diagnose the FastSpecFit versus GSWLC-X2 SFR offset for star-forming BGS galaxies.

For SV3/bright BGS targets matched to GSWLC-X2 (Salim et al. 2016, 2018)
and classified as star-forming by the conservative rest-frame 0.1(g-r)
color cut (util.fit_gr_cut), FastSpecFit SFRs are higher than GSWLC-X2 by
~0.45 dex (Section 5.4, Figure compare-sfr-external). This script makes a
six-panel reference figure of the offset, Delta log SFR (FastSpecFit minus
GSWLC-X2, the sign convention used throughout Section 5.4), against:

  (a) the stellar-mass offset, Delta log M* (FastSpecFit minus GSWLC-X2);
  (b) the attenuation offset, Delta tau_V (FastSpecFit minus GSWLC-X2);
  (c) redshift;
  (d) FastSpecFit specific SFR (note that FastSpecFit's SFR enters both
      axes, which builds in part of the correlation);
  (e) FastSpecFit stellar mass;

and (f) shows the distributions of log SFR minus the dust-corrected,
aperture-corrected H-alpha SFR (util.halpha_sfr; S/N>3 in H-alpha and
H-beta) for FastSpecFit, GSWLC-X2, Zou et al. (CIGALE, Loa), and Siudek et
al. (CIGALE-AGN, Iron; mass and SFR quality flags applied).

The summary statistics behind the Section 5.4 narrative are also printed,
including the H-alpha comparison in bins of aperture correction. At the
time of writing: Delta log SFR = +0.45 +/- 0.24 dex, Delta log M* = +0.21
dex, Delta log sSFR = +0.20 dex; Spearman(Delta SFR, Delta tau_V) = 0.50;
FastSpecFit, GSWLC-X2, Zou, and Siudek SFRs exceed the H-alpha SFRs by
+0.57, +0.11, +0.41, and +0.34 dex (common sample), independent of
aperture correction. The Zou offset relative to GSWLC-X2 tracks the
attenuation difference even more tightly (Spearman 0.82; Zou median
tau_V = 0.64 versus 0.39 for FastSpecFit and 0.28 for GSWLC-X2). Our
interpretation is a dust-age degeneracy shared by all three optical+WISE
fits despite their different attenuation laws: without the rest-frame UV
(GALEX) constraints that GSWLC-X2 has, a fit can trade more dust for
younger stars, which inflates the dust-corrected SFR.

Inputs (all prepared by code/prepare-external.py; nothing new is needed):
    external/gswlcx2-sv3-bright.fits
    external/zouhu-loa-sv3-bright.fits  (FastSpecFit H-alpha/H-beta fluxes
                                         and Zou SFRs and tau_V, joined on TARGETID)
    external/cigaleagn-iron-sv3-bright.fits  (Siudek SFRs, joined on TARGETID)
plus the SV3 fastspec catalogs, which util.gr_cut_coeff reads to fit the
star-forming / quiescent color cut.

Usage (from the repo root):

    python code/gswlc-sfr-diagnostics.py [--verbose]

Output: tex/figures/gswlc-sfr-diagnostics.pdf (reference figure; not
included in the paper).

"""
import os, sys, argparse, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, join

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from util import (plot_style, nmad, running_binstat, halpha_sfr, gr_color,
                  gr_cut_value, gr_cut_coeff, good_redshift_by_class,
                  cigaleagn_quality)

REPODIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXTDIR  = os.path.join(REPODIR, 'external')
FIGDIR  = os.path.join(REPODIR, 'tex', 'figures')

BGS_COLOR = '#1B7837'  # matches TARGET_CLASS_COLORS['BGS'] in build-figures.py


def read_sample(survey='sv3', coeff=None, verbose=False):
    """Star-forming BGS galaxies matched to GSWLC-X2, with H-alpha SFRs.

    Parameters
    ----------
    survey : str
        DESI survey flavor of the prepared external files.
    coeff : array_like, optional
        Coefficients of the 0.1(g-r) star-forming / quiescent cut. If None,
        fit them with ``util.gr_cut_coeff``.
    verbose : bool
        Print progress.

    Returns
    -------
    astropy.table.Table
        Matched star-forming sample with derived columns LOGSFR, LOGSFR_GSWLC,
        LOGSFR_ZOU, LOGSFR_SIUDEK, TAUV_ZOU, LOGSFR_HA (NaN where H-alpha or H-beta fails the S/N
        cut), DLOGSFR, DLOGM, DTAUV, and LOGSSFR.

    """
    import fitsio

    if coeff is None:
        coeff = gr_cut_coeff(survey, verbose=verbose)

    gfile = os.path.join(EXTDIR, f'gswlcx2-{survey}-bright.fits')
    zfile = os.path.join(EXTDIR, f'zouhu-loa-{survey}-bright.fits')
    if verbose:
        print(f'Reading {gfile}')
        print(f'Reading {zfile}')
    gsw = Table(fitsio.read(gfile))
    sfile = os.path.join(EXTDIR, f'cigaleagn-iron-{survey}-bright.fits')
    if verbose:
        print(f'Reading {sfile}')
    zou = Table(fitsio.read(zfile, columns=[
        'TARGETID', 'HALPHA_FLUX', 'HALPHA_FLUX_IVAR', 'HBETA_FLUX',
        'HBETA_FLUX_IVAR', 'SFR_ZOUHU', 'TAUV_ZOUHU']))
    siu = Table(fitsio.read(sfile))
    siu = siu[cigaleagn_quality(siu)]['TARGETID', 'SFR_CIGALEAGN']
    cat = join(gsw, zou, keys='TARGETID', join_type='left')
    cat = join(cat, siu, keys='TARGETID', join_type='left')

    good = (good_redshift_by_class(cat, survey)['BGS'] &
            (cat['SFR'] > 0) & (cat['SFR_GSWLCX2'] > 0))
    gr = gr_color(cat)
    with np.errstate(invalid='ignore'):
        sf = good & np.isfinite(gr) & (gr < gr_cut_value(cat['LOGMSTAR'], cat['Z'], coeff))
    cat = cat[sf]

    with np.errstate(divide='ignore', invalid='ignore'):
        cat['LOGSFR']       = np.log10(np.asarray(cat['SFR'], float))
        cat['LOGSFR_GSWLC'] = np.log10(np.asarray(cat['SFR_GSWLCX2'], float))
        sfr_zou = np.ma.filled(cat['SFR_ZOUHU'].astype(float), np.nan)
        cat['LOGSFR_ZOU']   = np.where(sfr_zou > 0, np.log10(sfr_zou), np.nan)
        sfr_siu = np.ma.filled(cat['SFR_CIGALEAGN'].astype(float), np.nan)
        cat['LOGSFR_SIUDEK'] = np.where(sfr_siu > 0, np.log10(sfr_siu), np.nan)
    cat['TAUV_ZOU'] = np.ma.filled(cat['TAUV_ZOUHU'].astype(float), np.nan)
    log_ha, good_ha = halpha_sfr(cat, snr_cut=3.0)
    cat['LOGSFR_HA'] = np.where(good_ha, log_ha, np.nan)

    # FastSpecFit minus GSWLC-X2
    cat['DLOGSFR'] = cat['LOGSFR'] - cat['LOGSFR_GSWLC']
    cat['DLOGM']   = cat['LOGMSTAR'] - cat['LOGMSTAR_GSWLCX2']
    cat['DTAUV']   = cat['TAUV'] - cat['TAUV_GSWLCX2']
    cat['LOGSSFR'] = cat['LOGSFR'] - cat['LOGMSTAR']

    if verbose:
        print(f'Star-forming BGS galaxies matched to GSWLC-X2: {len(cat):,}')
    return cat


def print_summary(cat):
    """Print the summary statistics quoted in the module docstring."""
    from scipy.stats import spearmanr

    dS, dM = np.asarray(cat['DLOGSFR']), np.asarray(cat['DLOGM'])
    print(f'N = {len(cat):,}')
    print(f'Delta log SFR  = {np.median(dS):+.3f} +/- {nmad(dS):.3f} dex')
    print(f'Delta log M*   = {np.median(dM):+.3f} +/- {nmad(dM):.3f} dex')
    print(f'Delta log sSFR = {np.median(dS - dM):+.3f} +/- {nmad(dS - dM):.3f} dex')
    print(f'tau_V median: FastSpecFit {np.median(cat["TAUV"]):.2f}, '
          f'GSWLC-X2 {np.median(cat["TAUV_GSWLCX2"]):.2f}')
    print(f'Spearman(Delta SFR, Delta M*)    = {spearmanr(dS, dM)[0]:.2f}')
    print(f'Spearman(Delta SFR, Delta tau_V) = {spearmanr(dS, cat["DTAUV"])[0]:.2f}')

    ha = np.asarray(cat['LOGSFR_HA'])
    hgood = np.isfinite(ha)
    print(f'H-alpha comparison (S/N>3 in H-alpha and H-beta): N = {hgood.sum():,}')
    cols = (('FastSpecFit', 'LOGSFR'), ('GSWLC-X2', 'LOGSFR_GSWLC'),
            ('Zou et al.', 'LOGSFR_ZOU'), ('Siudek et al.', 'LOGSFR_SIUDEK'))
    common = hgood.copy()
    for _, col in cols:
        common &= np.isfinite(np.asarray(cat[col]))
    for label, col in cols:
        x = np.asarray(cat[col])
        k = hgood & np.isfinite(x)
        r, rc = x[k] - ha[k], x[common] - ha[common]
        print(f'  {label:13s} - H-alpha: {np.median(r):+.3f} +/- {nmad(r):.3f} dex (N={k.sum():,}); '
              f'common sample: {np.median(rc):+.3f} (N={common.sum():,})')

    tz = np.asarray(cat['TAUV_ZOU'])
    k = np.isfinite(tz) & np.isfinite(np.asarray(cat['LOGSFR_ZOU']))
    dZ = np.asarray(cat['LOGSFR_ZOU'] - cat['LOGSFR_GSWLC'])[k]
    dTZ = tz[k] - np.asarray(cat['TAUV_GSWLCX2'])[k]
    print(f'Zou tau_V median = {np.median(tz[k]):.2f}; Zou minus GSWLC-X2: '
          f'Delta log SFR = {np.median(dZ):+.3f}, '
          f'Spearman(Delta SFR, Delta tau_V) = {spearmanr(dZ, dTZ)[0]:.2f}')

    ap = np.asarray(cat['APERCORR'], float)
    print('  versus aperture correction (median FastSpecFit / GSWLC-X2 minus H-alpha):')
    for a, b in ((1, 2), (2, 3), (3, 5), (5, np.inf)):
        k = hgood & (ap >= a) & (ap < b)
        if k.sum() == 0:
            continue
        print(f'    APERCORR [{a}, {b}): N={k.sum():5,}  '
              f'{np.median(cat["LOGSFR"][k] - ha[k]):+.2f} / '
              f'{np.median(cat["LOGSFR_GSWLC"][k] - ha[k]):+.2f}')


def _trend_panel(ax, x, y, xrange, nbins, xlabel, ylim):
    """Scatter plot of y versus x with the running median and interquartile range."""
    ax.scatter(x, y, s=3, color=BGS_COLOR, alpha=0.25, linewidths=0,
               rasterized=True, zorder=1)
    k = np.isfinite(x) & np.isfinite(y)
    cen, med, lo, hi = running_binstat(x[k], y[k], bins=nbins, xrange=xrange,
                                       min_count=30)
    ax.fill_between(cen, lo, hi, color='k', alpha=0.15, lw=0, zorder=2)
    ax.plot(cen, med, color='k', lw=2, zorder=3)
    ax.axhline(0, color='k', ls='--', lw=1, zorder=0)
    ax.axhline(np.median(y[k]), color='C3', ls=':', lw=1.5, zorder=0)
    ax.set_xlim(xrange)
    ax.set_ylim(ylim)
    ax.set_xlabel(xlabel)


def make_figure(cat, outfile):
    """Six-panel diagnostic figure; see the module docstring."""
    from scipy.stats import spearmanr

    plot_style(talk=True, font_scale=0.8, palette='colorblind')
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.5))
    fig.subplots_adjust(hspace=0.35, wspace=0.12)
    ylim = [-1.2, 2]
    dS = np.asarray(cat['DLOGSFR'])
    dslabel = r'$\Delta\log_{10}\mathrm{SFR}$ [FSF $-$ GSWLC-X2]'

    panels = [
        (np.asarray(cat['DLOGM']), [-0.4, 0.8], 12,
         r'$\Delta\log_{10}\mathcal{M}_{*}$ [FSF $-$ GSWLC-X2]'),
        (np.asarray(cat['DTAUV']), [-0.75, 1.5], 12,
         r'$\Delta\tau_V$ [FSF $-$ GSWLC-X2]'),
        (np.asarray(cat['Z']), [0.0, 0.3], 10, 'Redshift'),
        (np.asarray(cat['LOGSSFR']), [-11.5, -8.5], 12,
         r'$\log_{10}(\mathrm{sSFR}/\mathrm{yr}^{-1})$ [FSF]'),
        (np.asarray(cat['LOGMSTAR']), [8.0, 11.75], 12,
         r'$\log_{10}(\mathcal{M}_{*}\,h^{-2}/\mathcal{M}_{\odot})$ [FSF]'),
    ]
    letters = 'abcdef'
    for ip, (x, xr, nb, xl) in enumerate(panels):
        ax = axes.flat[ip]
        _trend_panel(ax, x, dS, xr, nb, xl, ylim)
        k = np.isfinite(x)
        ax.text(0.04, 0.95, f'({letters[ip]})  '
                rf'$\rho_\mathrm{{S}}={spearmanr(x[k], dS[k])[0]:+.2f}$',
                transform=ax.transAxes, va='top', ha='left', fontsize='small',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.75, pad=2))
        if ip % 3 == 0:
            ax.set_ylabel(dslabel)
        else:
            ax.tick_params(labelleft=False)

    # (f) SFR minus H-alpha SFR
    ax = axes.flat[5]
    ha = np.asarray(cat['LOGSFR_HA'])
    bins = np.arange(-1.5, 2.0 + 0.05, 0.1)
    for (label, col), color in zip((('FastSpecFit', 'LOGSFR'),
                                    ('GSWLC-X2', 'LOGSFR_GSWLC'),
                                    ('Zou et al.', 'LOGSFR_ZOU'),
                                    ('Siudek et al.', 'LOGSFR_SIUDEK')),
                                   (BGS_COLOR, 'C0', 'C1', 'C4')):
        x = np.asarray(cat[col])
        k = np.isfinite(ha) & np.isfinite(x)
        r = x[k] - ha[k]
        ax.hist(r, bins=bins, histtype='step', lw=2, color=color,
                label=f'{label} ({np.median(r):+.2f})')
        ax.axvline(np.median(r), color=color, ls=':', lw=1.5)
    ax.axvline(0, color='k', ls='--', lw=1)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, 1.45 * ax.get_ylim()[1])
    ax.set_xlabel(r'$\log_{10}\mathrm{SFR} - \log_{10}\mathrm{SFR}(\mathrm{H}\alpha)$')
    ax.set_ylabel('Number')
    ax.yaxis.set_label_position('right')
    ax.yaxis.tick_right()
    ax.legend(loc='upper right', fontsize='x-small', frameon=False)
    ax.text(0.04, 0.95, f'(f)  $N={np.isfinite(ha).sum():,}$', transform=ax.transAxes,
            va='top', ha='left', fontsize='small')

    fig.suptitle(f'Star-forming BGS galaxies matched to GSWLC-X2 ($N={len(cat):,}$)',
                 y=0.97)
    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f'Wrote {outfile}')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__.split('\n')[0])
    parser.add_argument('--survey', default='sv3',
                        help='DESI survey flavor of the prepared external files.')
    parser.add_argument('--outfile', default=os.path.join(FIGDIR, 'gswlc-sfr-diagnostics.pdf'),
                        help='Output figure.')
    parser.add_argument('--verbose', action='store_true', help='Print progress.')
    args = parser.parse_args()

    warnings.filterwarnings('ignore', category=RuntimeWarning)
    cat = read_sample(args.survey, verbose=args.verbose)
    print_summary(cat)
    make_figure(cat, args.outfile)


if __name__ == '__main__':
    main()
