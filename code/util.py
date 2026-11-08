"""Utility functions for reading DESI/FastSpecFit VAC catalogs.

Each read_* function returns an astropy Table with a standardized data model
so that figure-making scripts do not need to know catalog-specific details.

Cosmology note: LOGMSTAR and SFR in the FastSpecFit VACs are stored with h=1
(i.e., M_star in units of h^{-2} M_sun). Conversion to H0=70 adds
+2*log10(1/0.7) ~ +0.31 dex; to Planck 2018 (h=0.674) adds ~+0.34 dex.
Comparisons across catalogs must use a consistent h; each read_* function
documents which convention it returns.

"""
import os
import numpy as np
from glob import glob
from astropy.table import Table, hstack, vstack
import fitsio

# At NERSC we read the official (read-only) VAC directories on CFS. Off NERSC
# (e.g., a laptop), we read from a local mirror under $DESI_ROOT, laid out to
# match the official VAC directory structure.
_ON_NERSC = 'NERSC_HOST' in os.environ

# Per-specprod catalog directory paths.  Add entries here as new specprods
# are released; set a value to None when the path is not yet known.
if _ON_NERSC:
    _DESI_VAC = '/dvs_ro/cfs/cdirs/desi/vac'
    _SPECPROD_CONFIG = {
        'loa': {
            'fastspec': f'{_DESI_VAC}/dr2/fastspecfit/loa/v2.0/catalogs',
            'fastphot': f'{_DESI_VAC}/dr2/fastphot/loa/v2.0/catalogs',
        },
        'iron': {
            'fastspec': '/dvs_ro/cfs/cdirs/desi/vac/dr1/fastspecfit/iron/v3.0/catalogs',
            'fastphot': None,
        },
        'fuji': {
            'fastspec': '/dvs_ro/cfs/cdirs/desi/vac/edr/fastspecfit/fuji/v3.2/catalogs',
            'fastphot': None,
        },
    }
else:
    if 'DESI_ROOT' not in os.environ:
        raise RuntimeError('DESI_ROOT must be set to locate VAC catalogs off of NERSC.')
    _DESI_VAC = os.path.join(os.environ['DESI_ROOT'], 'vac')
    _SPECPROD_CONFIG = {
        'loa': {
            'fastspec': f'{_DESI_VAC}/dr2/fastspecfit/loa/v2.0/catalogs',
            'fastphot': f'{_DESI_VAC}/dr2/fastphot/loa/v2.0/catalogs',
        },
        'iron': {
            'fastspec': f'{_DESI_VAC}/dr1/fastspecfit/iron/v4.0/catalogs',
            'fastphot': None,
        },
        'fuji': {
            'fastspec': None,
            'fastphot': None,
        },
    }

DEFAULT_SPECPROD = 'loa'

# Earlier public fastspec VACs, read only for the release-to-release comparison
# in Appendix B (see prepare-external.py --vacs). Keyed by (specprod, version);
# values are (catalog directory, FITS extensions to join). Catalogs made with
# fastspecfit<3.0 (data model 2.5.0) have no SPECPHOT extension.
_ARCHIVE_FASTSPEC = {
    ('iron', 'v3.0'): ('/dvs_ro/cfs/cdirs/desi/vac/dr1/fastspecfit/iron/v3.0/catalogs',
                       ['METADATA', 'SPECPHOT', 'FASTSPEC']),
    ('iron', 'v2.1'): ('/dvs_ro/cfs/cdirs/desi/vac/dr1/fastspecfit/iron/v2.1/catalogs',
                       ['METADATA', 'FASTSPEC']),
    ('fuji', 'v3.2'): ('/dvs_ro/cfs/cdirs/desi/vac/edr/fastspecfit/fuji/v3.2/catalogs',
                       ['METADATA', 'FASTSPEC']),
}

# survey+program combinations split into 12 nside=1 healpix files
_SPLIT_COMBOS = {('main', 'bright'), ('main', 'dark')}

# Columns always returned from METADATA regardless of the caller's ``columns``
# argument.  All survey-specific targeting bit names are listed here; columns
# absent from a given file are silently skipped by _read_extensions.
_DEFAULT_COLUMNS = frozenset({
    'TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX',
    'RA', 'DEC', 'Z', 'ZERR', 'ZWARN', 'DELTACHI2', 'SPECTYPE',
    'COADD_FIBERSTATUS',
    # OII doublet columns needed for the ELG redshift quality cut
    'OII_3726_FLUX', 'OII_3729_FLUX', 'OII_3726_FLUX_IVAR', 'OII_3729_FLUX_IVAR',
    # main-survey targeting bits
    'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET', 'SCND_TARGET', 'ETC_TARGET',
    # SV1
    'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'SV1_MWS_TARGET',
    'SV1_SCND_TARGET', 'SV1_ETC_TARGET',
    # SV2
    'SV2_DESI_TARGET', 'SV2_BGS_TARGET', 'SV2_MWS_TARGET',
    'SV2_SCND_TARGET', 'SV2_ETC_TARGET',
    # SV3
    'SV3_DESI_TARGET', 'SV3_BGS_TARGET', 'SV3_MWS_TARGET',
    'SV3_SCND_TARGET', 'SV3_ETC_TARGET',
    # CMX
    'CMX_TARGET',
})


def _catfiles(specprod, vactype, survey, program, vacdir=None):
    """Return sorted list of catalog paths for a specprod/survey/program combination.

    If `vacdir` is given, it overrides the directory from _SPECPROD_CONFIG.

    """
    if vacdir is None:
        cfg = _SPECPROD_CONFIG.get(specprod)
        if cfg is None:
            raise ValueError(f'Unknown specprod {specprod!r}. '
                             f'Known specprods: {list(_SPECPROD_CONFIG)}')
        vacdir = cfg[vactype]
        if vacdir is None:
            raise ValueError(f'VAC path not yet configured for specprod={specprod!r}, '
                             f'vactype={vactype!r}')

    basename = f'{vactype}-{specprod}-{survey}-{program}'
    if (survey, program) in _SPLIT_COMBOS:
        pattern = os.path.join(vacdir, f'{basename}-nside1-hp??.fits')
        files = sorted(glob(pattern))
        if not files:
            raise FileNotFoundError(f'No files found matching {pattern}')
        return files
    path = os.path.join(vacdir, f'{basename}.fits')
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    return [path]


def _read_extensions(filepath, extensions, columns=None):
    """Read and hstack FITS extensions from one file.

    The first extension (METADATA) is always read in full. For subsequent
    extensions, `columns` limits which columns are loaded; columns already
    present from METADATA are never duplicated.

    Parameters
    ----------
    filepath : str
    extensions : list of str
        Extension names to read, e.g. ['METADATA', 'SPECPHOT', 'FASTSPEC'].
    columns : list of str or None
        Columns to read from extensions[1:]. None means read all.

    Returns
    -------
    astropy.table.Table

    """
    parts = []
    seen = set()

    with fitsio.FITS(filepath) as fits:
        for i, ext in enumerate(extensions):
            ext_cols = set(c.upper() for c in fits[ext].get_colnames())

            #if i == 0:
            #    # METADATA: always full
            #    read_cols = None
            if columns is not None:
                want = (set(c.upper() for c in columns) & ext_cols) - seen
                if not want:
                    seen.update(ext_cols)
                    continue
                read_cols = sorted(want)
            else:
                read_cols = sorted(ext_cols - seen)
                if not read_cols:
                    continue

            data = fits[ext].read(columns=read_cols)
            t = Table(data)
            parts.append(t)
            seen.update(t.colnames)

    return hstack(parts) if len(parts) > 1 else parts[0]


def _read_vac(specprod, vactype, extensions, survey, program, columns, verbose,
              vacdir=None):
    """Internal driver: read and vstack all catalog files for a survey/program."""
    files = _catfiles(specprod, vactype, survey, program, vacdir=vacdir)
    chunks = []
    for f in files:
        if verbose:
            print(f'Reading {f}')
        chunks.append(_read_extensions(f, extensions, columns=columns))
    cat = vstack(chunks) if len(chunks) > 1 else chunks[0]
    if verbose:
        print(f'  {len(cat):,} rows')
    return cat


# ---------------------------------------------------------------------------
# Public read functions
# ---------------------------------------------------------------------------

def read_fastspec(survey='main', program='dark', specprod=DEFAULT_SPECPROD,
                  columns=None, version=None, verbose=False):
    """Read the fastspec VAC for a given specprod, survey, and program.

    Joins the METADATA, SPECPHOT, and FASTSPEC extensions into a single Table.
    For main-bright and main-dark the 12 nside=1 healpix files are stacked
    transparently.

    Parameters
    ----------
    survey : str
        Survey name: 'main', 'sv1', 'sv2', 'sv3', or 'special'.
    program : str
        Program name: 'bright', 'dark', or 'backup'.
    specprod : str
        Spectroscopic production name (default: 'loa'). Must be a key in
        _SPECPROD_CONFIG.
    columns : list of str, optional
        Columns to load from SPECPHOT and FASTSPEC. METADATA is always read
        in full. If None, all columns from all three extensions are returned.
    version : str, optional
        Read an earlier public VAC listed in _ARCHIVE_FASTSPEC (e.g., 'v3.0'
        with specprod='iron') instead of the default catalog for `specprod`.
    verbose : bool
        Print file names and row counts while reading.

    Returns
    -------
    astropy.table.Table
        Merged catalog. Key columns always present: TARGETID, SURVEY, PROGRAM,
        HEALPIX, RA, DEC, Z, ZERR, ZWARN, SPECTYPE, and all survey-specific
        targeting bit columns that exist in the file (DESI_TARGET, BGS_TARGET,
        MWS_TARGET, SV3_DESI_TARGET, etc.). Additional columns from SPECPHOT
        and FASTSPEC are controlled by ``columns``.

    """
    if columns is not None:
        columns = list(set(columns) | _DEFAULT_COLUMNS)
    vacdir = None
    extensions = ['METADATA', 'SPECPHOT', 'FASTSPEC']
    if version is not None:
        key = (specprod, version)
        if key not in _ARCHIVE_FASTSPEC:
            raise ValueError(f'Unknown VAC {specprod}/{version}. '
                             f'Known: {list(_ARCHIVE_FASTSPEC)}')
        vacdir, extensions = _ARCHIVE_FASTSPEC[key]
    return _read_vac(
        specprod, 'fastspec', extensions,
        survey, program, columns, verbose, vacdir=vacdir,
    )


def corner_plot(plotdata, labels, ranges, bins=50, truths=None, sigmas=None,
                titles=None, unity=False, diag_ylabel='N',
                contour_levels=None, contour_lw=1.5, smooth=1.0,
                cmap='Blues', show_residuals=True, groups=None,
                split_contours=False,
                figsize=None, suptitle='', subplots_adjust=None):
    """Corner-style N×N plot: histograms on the diagonal, Hess+contours on the lower triangle.

    Adapted from fastspecfit.qa._corner_plot for catalog-scale datasets.
    Off-diagonal panels show a log-stretched 2D histogram (Hess diagram) with
    smoothed density contours at specified cumulative enclosed fractions overlaid.

    Parameters
    ----------
    plotdata : array_like, shape (nsamples, ndim)
        Data array, one column per parameter. Ignored when ``groups`` is provided.
    labels : list of str
        Axis labels, one per parameter.
    ranges : list of (lo, hi) tuples
        Axis limits for each parameter.
    bins : int
        Number of bins along each axis for the 2D histogram and diagonal histograms.
    truths : list of float or None
        Reference values: vertical lines on diagonal, crosshairs on off-diagonal.
        Skipped if None. Not used in groups mode.
    sigmas : list of float or None
        1-sigma uncertainties drawn as dashed lines on the diagonal.
        Only used when truths is not None.
    titles : list of str or None
        Titles above each diagonal histogram. Skipped if None.
    unity : bool
        Draw a dashed 1:1 reference line on each off-diagonal panel.
    diag_ylabel : str
        Y-axis label on the leftmost diagonal histogram.
    contour_levels : list of float or None
        Cumulative enclosed fractions at which to draw contour lines, e.g.
        [0.5, 0.75, 0.95, 0.995]. Default: [0.5, 0.75, 0.95, 0.995].
    contour_lw : float
        Line width for contours.
    smooth : float
        Gaussian smoothing sigma (in bins) applied to the 2D histogram before
        computing contour levels. Set to 0 to disable.
    cmap : str
        Colormap for the Hess diagram (single-group mode only).
    show_residuals : bool
        If True, annotate each off-diagonal panel with Δ = col_y − col_x
        statistics.  In single-group mode: one-line text annotation.  In
        groups mode: a per-group legend with colored Line2D handles.
    groups : list of dict or None
        When provided, diagonal panels show per-group step histograms instead
        of a single gray histogram.  Each dict must have keys ``'data'``
        (array shape (n, ndim)), ``'color'`` (str), and ``'label'`` (str).
        Off-diagonal behavior is controlled by ``split_contours``.
    split_contours : bool
        Only meaningful when ``groups`` is provided.  If False (default),
        off-diagonal panels show the all-data Hess diagram (requires valid
        ``plotdata``).  If True, off-diagonal panels show per-group contours
        with no Hess background.
    figsize : tuple of (float, float) or None
        Figure size in inches as (width, height). Default is (3*ndim, 3*ndim),
        minimum 6×6.
    suptitle : str
        Figure suptitle.
    subplots_adjust : dict or None
        Forwarded to fig.subplots_adjust; if None, tight_layout is used.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from scipy.ndimage import gaussian_filter

    if contour_levels is None:
        contour_levels = [0.5, 0.75, 0.95, 0.995]

    use_groups = groups is not None and len(groups) > 0
    use_groups_offdiag = use_groups and split_contours
    ndim = len(labels)

    # plotdata is needed for off-diagonal Hess whenever split_contours is False
    if not use_groups or not split_contours:
        plotdata = np.asarray(plotdata)

    if figsize is None:
        _size = max(3 * ndim, 6)
        figsize = (_size, _size)
    fig, axes = plt.subplots(ndim, ndim, figsize=figsize)
    ax = np.array(axes).reshape((ndim, ndim))

    for yi in range(ndim):
        for xi in range(ndim):
            a = ax[yi, xi]
            if xi > yi:
                a.set_visible(False)
                continue

            lo_x, hi_x = ranges[xi]
            lo_y, hi_y = ranges[yi]

            if xi == yi:
                if use_groups:
                    for grp in groups:
                        a.hist(grp['data'][:, xi], bins=bins, range=(lo_x, hi_x),
                               color=grp['color'], histtype='step', lw=2)
                else:
                    a.hist(plotdata[:, xi], bins=bins, range=(lo_x, hi_x),
                           color='gray', alpha=0.75, edgecolor='k')
                    if truths is not None:
                        a.axvline(truths[xi], color='C0', lw=2, ls='-')
                        if sigmas is not None:
                            a.axvline(truths[xi] + sigmas[xi], color='C0', lw=1, ls='--')
                            a.axvline(truths[xi] - sigmas[xi], color='C0', lw=1, ls='--')
                if titles is not None:
                    a.set_title(titles[xi], fontsize=8)
                a.set_xlim(lo_x, hi_x)
                if xi == 0:
                    a.set_ylabel(diag_ylabel)
                else:
                    a.tick_params(labelleft=False)
            else:
                if use_groups_offdiag:
                    # Per-group contours (no Hess background)
                    for grp in groups:
                        H, xedges, yedges = np.histogram2d(
                            grp['data'][:, xi], grp['data'][:, yi],
                            bins=bins, range=[[lo_x, hi_x], [lo_y, hi_y]],
                        )
                        xc = 0.5 * (xedges[:-1] + xedges[1:])
                        yc = 0.5 * (yedges[:-1] + yedges[1:])
                        Hs = gaussian_filter(H, smooth) if smooth > 0 else H
                        flat = np.sort(Hs.flatten())[::-1]
                        cumsum = np.cumsum(flat)
                        total = cumsum[-1]
                        if total > 0 and contour_levels:
                            lvls = []
                            for frac in contour_levels:
                                idx = np.searchsorted(cumsum, frac * total)
                                lvls.append(flat[min(idx, len(flat) - 1)])
                            lvls = sorted(v for v in set(lvls) if v > 0)
                            if lvls:
                                a.contour(xc, yc, Hs.T, levels=lvls,
                                          colors=grp['color'], linewidths=contour_lw)
                    if show_residuals:
                        from matplotlib.lines import Line2D
                        handles = []
                        for grp in groups:
                            resid = grp['data'][:, yi] - grp['data'][:, xi]
                            med = np.median(resid)
                            mu  = np.mean(resid)
                            sig = np.std(resid)
                            lbl = (f"{grp['label']}: "
                                   f"$\\Delta={med:+.3f}\\,({mu:+.3f}\\pm{sig:.3f})$")
                            handles.append(Line2D([0], [0], color=grp['color'],
                                                  lw=2, label=lbl))
                        a.legend(handles=handles, loc='upper left',
                                 fontsize='x-small', framealpha=0.75,
                                 handlelength=1.5)
                else:
                    # Hess diagram: log-stretched 2D histogram
                    H, xedges, yedges = np.histogram2d(
                        plotdata[:, xi], plotdata[:, yi],
                        bins=bins, range=[[lo_x, hi_x], [lo_y, hi_y]],
                    )
                    xc = 0.5 * (xedges[:-1] + xedges[1:])
                    yc = 0.5 * (yedges[:-1] + yedges[1:])
                    a.pcolormesh(xedges, yedges, H.T, norm=LogNorm(vmin=1), cmap=cmap)
                    Hs = gaussian_filter(H, smooth) if smooth > 0 else H
                    flat = np.sort(Hs.flatten())[::-1]
                    cumsum = np.cumsum(flat)
                    total = cumsum[-1]
                    if total > 0 and contour_levels:
                        lvls = []
                        for frac in contour_levels:
                            idx = np.searchsorted(cumsum, frac * total)
                            lvls.append(flat[min(idx, len(flat) - 1)])
                        lvls = sorted(v for v in set(lvls) if v > 0)
                        if lvls:
                            a.contour(xc, yc, Hs.T, levels=lvls,
                                      colors='k', linewidths=contour_lw)
                    if truths is not None:
                        a.axvline(truths[xi], color='C0', lw=1, ls='-', alpha=0.75)
                        a.axhline(truths[yi], color='C0', lw=1, ls='-', alpha=0.75)
                    if show_residuals:
                        resid = plotdata[:, yi] - plotdata[:, xi]
                        med = np.median(resid)
                        mu  = np.mean(resid)
                        sig = np.std(resid)
                        txt = f'$\\Delta={med:+.3f}\\,({mu:+.3f}\\pm{sig:.3f})$'
                        a.text(0.04, 0.96, txt, transform=a.transAxes,
                               fontsize='x-small', va='top', ha='left',
                               bbox=dict(facecolor='white', edgecolor='none',
                                         alpha=0.75, pad=2))

                if unity:
                    lo = max(lo_x, lo_y)
                    hi = min(hi_x, hi_y)
                    a.plot([lo, hi], [lo, hi], color='k', lw=1, ls='-')

                a.set_xlim(lo_x, hi_x)
                a.set_ylim(lo_y, hi_y)
                if xi == 0:
                    a.set_ylabel(labels[yi])
                else:
                    a.tick_params(labelleft=False)

            if yi == ndim - 1:
                a.set_xlabel(labels[xi])
            else:
                a.tick_params(labelbottom=False)

    if suptitle:
        fig.suptitle(suptitle)
    if subplots_adjust:
        fig.subplots_adjust(**subplots_adjust)
    else:
        fig.tight_layout()

    return fig


def hess_contours(ax, x, y, xrange, yrange, bins=50, smooth=1.0,
                  contour_levels=None, cmap='Blues', contour_lw=1.5,
                  contour_color='k', outlier_ms=2, background=True):
    """Hess diagram with smoothed cumulative contours, matching corner_plot style.

    Points that fall outside the outermost contour are drawn as individual small
    markers (size ``outlier_ms``); pass ``outlier_ms=0`` to suppress them.
    The pcolormesh background is rasterized so PDF output stays compact.
    Pass ``background=False`` to draw contours and outlier scatter only (no
    pcolormesh), which is useful when overlaying multiple samples on one axes.
    """
    from matplotlib.colors import LogNorm
    from scipy.ndimage import gaussian_filter

    if contour_levels is None:
        contour_levels = [0.5, 0.75, 0.95, 0.995]

    H, xedges, yedges = np.histogram2d(x, y, bins=bins,
                                        range=[xrange, yrange])
    xc = 0.5 * (xedges[:-1] + xedges[1:])
    yc = 0.5 * (yedges[:-1] + yedges[1:])
    if background:
        ax.pcolormesh(xedges, yedges, H.T, norm=LogNorm(vmin=1), cmap=cmap,
                      rasterized=True, zorder=1)

    Hs = gaussian_filter(H, smooth) if smooth > 0 else H
    flat = np.sort(Hs.flatten())[::-1]
    cumsum = np.cumsum(flat)
    total = cumsum[-1]
    lvls = []
    if total > 0 and contour_levels:
        for frac in contour_levels:
            idx = np.searchsorted(cumsum, frac * total)
            lvls.append(flat[min(idx, len(flat) - 1)])
        lvls = sorted(v for v in set(lvls) if v > 0)
        if lvls:
            ax.contour(xc, yc, Hs.T, levels=lvls,
                       colors=contour_color, linewidths=contour_lw, zorder=3)

    # individual points outside the outermost contour, drawn above the Hess
    # background (zorder=2) but below contours (zorder=3)
    if lvls and outlier_ms > 0:
        in_range = ((x >= xrange[0]) & (x <= xrange[1]) &
                    (y >= yrange[0]) & (y <= yrange[1]))
        xi = np.clip(np.digitize(x[in_range], xedges) - 1, 0, len(xc) - 1)
        yi = np.clip(np.digitize(y[in_range], yedges) - 1, 0, len(yc) - 1)
        outside = Hs[xi, yi] < lvls[0]
        ax.scatter(x[in_range][outside], y[in_range][outside],
                   s=outlier_ms, c=[contour_color], alpha=0.8,
                   linewidths=0, zorder=2, rasterized=True)


def nmad(x):
    """Normalized median absolute deviation: 1.4826 * median(|x - median(x)|)."""
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def running_binstat(x, y, bins=20, xrange=None, q_lo=25, q_hi=75, min_count=5):
    """Bin y on x; return bin centers, running median, and percentile bounds.

    Parameters
    ----------
    x, y : array_like
        Input arrays (must be 1-D and the same length; should be finite).
    bins : int
        Number of equal-width bins.
    xrange : (float, float) or None
        Range of x. Passed to binned_statistic; defaults to (x.min(), x.max()).
    q_lo, q_hi : float
        Lower and upper percentiles for the shaded band (default: 25, 75 → IQR).
    min_count : int
        Bins with fewer than this many points are set to NaN.

    Returns
    -------
    centers : ndarray, shape (bins,)
        Bin midpoints.
    median : ndarray, shape (bins,)
        Median of y in each bin.
    p_lo, p_hi : ndarray, shape (bins,)
        q_lo-th and q_hi-th percentiles of y in each bin.
        Bins with count < min_count contain NaN.
    """
    from scipy.stats import binned_statistic

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    kw = dict(bins=bins, range=xrange)

    cnt, edges, _ = binned_statistic(x, y, statistic='count', **kw)
    med, _, _     = binned_statistic(x, y, statistic='median', **kw)
    lo,  _, _     = binned_statistic(x, y,
                                     statistic=lambda v: np.percentile(v, q_lo), **kw)
    hi,  _, _     = binned_statistic(x, y,
                                     statistic=lambda v: np.percentile(v, q_hi), **kw)

    centers = 0.5 * (edges[:-1] + edges[1:])
    sparse = cnt < min_count
    for arr in (med, lo, hi):
        arr[sparse] = np.nan

    return centers, med, lo, hi


def wilson_binomial_ci(k, n, z=1.0):
    """Wilson binomial confidence interval (asymmetric, suitable for plt.errorbar).

    Parameters
    ----------
    k : array_like
        Number of successes in each bin.
    n : array_like
        Total number of trials in each bin.
    z : float
        Number of standard deviations (default 1.0 → 68% CI).

    Returns
    -------
    lo_err, hi_err : ndarray
        Distances below and above the fraction p = k/n to the CI bounds.
        Pass as ``yerr=[lo_err, hi_err]`` to ``plt.errorbar``.
    """
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    p = k / n
    z2 = z * z
    n_plus = n + z2
    center = (k + z2 / 2.0) / n_plus
    margin = (z / n_plus) * np.sqrt(k * (n - k) / n + z2 / 4.0)
    lo = np.clip(center - margin, 0.0, 1.0)
    hi = np.clip(center + margin, 0.0, 1.0)
    return p - lo, hi - p


def halpha_sfr(cat, snr_cut=3.0):
    """Dust-corrected Hα SFR for objects with significant Hα and Hβ detections.

    Applies an aperture correction (APERCORR), derives E(B-V) from the Balmer
    decrement using a Cardelli+O'Donnell (R_V=3.1) extinction curve, and converts
    to SFR via the Hao et al. (2011, Table 2) Kroupa-IMF calibration, converted
    to a Chabrier IMF with the Madau & Dickinson (2014) factors (0.63/0.67).

    Parameters
    ----------
    cat : astropy.table.Table
        Must include HALPHA_FLUX, HALPHA_FLUX_IVAR, HBETA_FLUX, HBETA_FLUX_IVAR,
        APERCORR, Z.  Fluxes must be in units of 1e-17 erg/s/cm².
    snr_cut : float
        Minimum S/N required on both Hα and Hβ.

    Returns
    -------
    log_sfr : numpy.ndarray
        log10 SFR [M_sun/yr], NaN for objects that do not pass the S/N cut.
    good : numpy.ndarray of bool
        True where both lines exceed snr_cut.
    """
    from fastspecfit.cosmo import TabulatedDESI

    ha_flux  = np.array(cat['HALPHA_FLUX'],      dtype=float)
    hb_flux  = np.array(cat['HBETA_FLUX'],       dtype=float)
    ha_ivar  = np.array(cat['HALPHA_FLUX_IVAR'], dtype=float)
    hb_ivar  = np.array(cat['HBETA_FLUX_IVAR'], dtype=float)
    apercorr = np.array(cat['APERCORR'],         dtype=float)
    z        = np.array(cat['Z'],                dtype=float)

    with np.errstate(invalid='ignore'):
        good = (
            (ha_flux * np.sqrt(ha_ivar) > snr_cut) &
            (hb_flux * np.sqrt(hb_ivar) > snr_cut) &
            (ha_flux > 0) & (hb_flux > 0) & (apercorr > 0)
        )

    # Cardelli+O'Donnell R_V=3.1: k(Hα)=2.53 at 6563 Å, k(Hβ)=3.61 at 4861 Å
    k_ha, k_hb = 2.53, 3.61
    with np.errstate(divide='ignore', invalid='ignore'):
        ebv = np.where(
            good,
            np.maximum(2.5 / (k_hb - k_ha) * np.log10(ha_flux / hb_flux / 2.86), 0.0),
            0.0,
        )

    # Aperture-corrected, dust-corrected Hα flux [erg/s/cm²]
    f_ha_corr = ha_flux * apercorr * 10.0 ** (0.4 * k_ha * ebv) * 1e-17

    # Luminosity distance [Mpc/h] → cm; h=1 convention matches LOGMSTAR in VAC
    MPC_TO_CM = 3.0857e24
    dlum_cm = TabulatedDESI().luminosity_distance(z) * MPC_TO_CM

    l_ha = 4.0 * np.pi * dlum_cm ** 2 * f_ha_corr  # [erg/s]

    # Hao+11 Table 2: log C(Hα) = -41.257 for Kroupa (0.1-100 Msun, 100 Myr);
    # Kroupa -> Chabrier: log10(0.63/0.67) ≈ -0.027 dex (Madau & Dickinson 2014).
    log_c_ha = -41.257 + np.log10(0.63 / 0.67)

    with np.errstate(divide='ignore', invalid='ignore'):
        log_sfr = np.where(good, np.log10(l_ha) + log_c_ha, np.nan)

    return log_sfr, good


GR_CUT_ZBINS     = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2]
GR_CUT_MSTARBINS = np.arange(8.5, 12.0 + 0.01, 0.25)
GR_CUT_COLORLIM  = [-0.2, 1.4]


def gr_color(cat):
    """Rest-frame 0.1(g-r) color from the FastSpecFit ABSMAG01 columns.

    Returns NaN where either magnitude is missing (zero or non-finite).

    """
    g = np.asarray(cat['ABSMAG01_SDSS_G'], dtype=float)
    r = np.asarray(cat['ABSMAG01_SDSS_R'], dtype=float)
    bad = ~np.isfinite(g) | ~np.isfinite(r) | (g == 0) | (r == 0)
    return np.where(bad, np.nan, g - r)


def gr_cut_value(logmstar, z, coeff):
    """Dividing 0.1(g-r) color at a given stellar mass (h=1) and redshift.

    color_cut = a0 + a1 * (logmstar - 10) + a2 * (z - 0.1)

    """
    a0, a1, a2 = coeff
    return a0 + a1 * (np.asarray(logmstar) - 10.) + a2 * (np.asarray(z) - 0.1)


GR_CUT_DELTA = 0.20


def _fit_plane(cells, ycol, wcol, nclip=3, sigclip=3.0, floor=0.02):
    """Weighted, sigma-clipped fit of y = a0 + a1 (logM - 10) + a2 (z - 0.1).

    Returns the coefficients and a boolean mask of the cells kept in the fit.

    """
    A = np.vstack([np.ones(len(cells)), cells['MSTARCEN'] - 10.,
                   cells['ZCEN'] - 0.1]).T
    y = np.asarray(cells[ycol], dtype=float)
    w = np.sqrt(np.asarray(cells[wcol], dtype=float))

    def _solve(used):
        coeff, *_ = np.linalg.lstsq(A[used] * w[used, None], y[used] * w[used],
                                    rcond=None)
        return coeff

    used = np.ones(len(cells), bool)
    for _ in range(nclip):
        coeff = _solve(used)
        resid = y - A @ coeff
        newused = np.abs(resid) < sigclip * max(nmad(resid[used]), floor)
        if np.array_equal(newused, used) or newused.sum() < 4:
            break
        used = newused
    return _solve(used), used


def fit_gr_cut(logmstar, z, color, delta=GR_CUT_DELTA, zbins=GR_CUT_ZBINS,
               mstarbins=GR_CUT_MSTARBINS, colorlim=GR_CUT_COLORLIM, dcolor=0.02,
               smooth=2.0, nmin=200, nred_min=100, min_sep=0.15,
               min_prominence=0.05):
    """Fit a conservative, redshift-dependent star-forming / quiescent color cut.

    The cut runs parallel to the red sequence, ``delta`` mag blueward of its
    ridge line. We derive it in two passes over cells in (redshift, stellar
    mass), using a smoothed 0.1(g-r) histogram in each cell:

    1. Valley. In cells with two prominent peaks separated by at least
       ``min_sep`` mag, we take the color at the minimum between them. A plane
       fit to these valley colors gives a provisional red/blue boundary. The
       valley is only measurable in a minority of cells (mostly at low mass
       and low redshift), so we use it only to separate the red sequence from
       the blue cloud in the next step.
    2. Ridge. In every cell with at least ``nred_min`` galaxies redward of the
       provisional boundary, we take the most populous histogram peak redward
       of that boundary as the red-sequence ridge color, and fit a plane to
       these ridge colors.

    The returned cut is the ridge plane shifted blueward by ``delta``.
    Galaxies bluer than the cut are classified as star-forming. With the
    default ``delta=0.20``, the cut lies about 0.07 mag blueward of the
    median valley, which keeps green-valley galaxies out of the star-forming
    sample.

    Each plane has the form a0 + a1 (logM - 10) + a2 (z - 0.1), and is fit
    with each cell weighted by sqrt(N), iteratively rejecting 3-sigma
    outliers.

    Parameters
    ----------
    logmstar, z, color : array_like
        Stellar mass (log10, h=1), redshift, and rest-frame 0.1(g-r) color.
        Non-finite entries are ignored.
    delta : float
        Blueward offset of the cut from the red-sequence ridge [mag].
    zbins, mstarbins : array_like
        Bin edges in redshift and stellar mass.
    colorlim : (float, float)
        Color range of the histogram.
    dcolor : float
        Histogram bin width [mag].
    smooth : float
        Gaussian smoothing kernel in units of histogram bins.
    nmin : int
        Minimum number of galaxies in a cell.
    nred_min : int
        Minimum number of galaxies redward of the provisional boundary for a
        cell to contribute to the ridge fit.
    min_sep : float
        Minimum color separation between the two peaks in the valley pass [mag].
    min_prominence : float
        Minimum peak prominence as a fraction of the smoothed histogram maximum.

    Returns
    -------
    coeff : ndarray, shape (3,)
        Coefficients (a0, a1, a2) of the star-forming / quiescent cut.
    fit : dict
        ``valley_coeff`` and ``ridge_coeff`` (plane coefficients), and the
        per-cell tables ``valley`` (ZCEN, MSTARCEN, N, VALLEY, PEAK_BLUE,
        PEAK_RED, USED) and ``ridge`` (ZCEN, MSTARCEN, N, NRED, RIDGE, USED),
        where USED is False for cells rejected by the sigma clipping.

    """
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks
    from astropy.table import Table

    logmstar = np.asarray(logmstar, dtype=float)
    z        = np.asarray(z, dtype=float)
    color    = np.asarray(color, dtype=float)
    ok = np.isfinite(logmstar) & np.isfinite(z) & np.isfinite(color)

    edges   = np.arange(colorlim[0], colorlim[1] + dcolor / 2, dcolor)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # smoothed histogram and peaks for every populated cell
    hcells = []
    for z1, z2 in zip(zbins[:-1], zbins[1:]):
        for m1, m2 in zip(mstarbins[:-1], mstarbins[1:]):
            I = ok & (z >= z1) & (z < z2) & (logmstar >= m1) & (logmstar < m2)
            if I.sum() < nmin:
                continue
            H, _ = np.histogram(color[I], bins=edges)
            Hs = gaussian_filter1d(H.astype(float), smooth)
            peaks, props = find_peaks(Hs, prominence=min_prominence * Hs.max())
            hcells.append(dict(zcen=0.5 * (z1 + z2), mcen=np.median(logmstar[I]),
                               I=I, Hs=Hs, peaks=peaks,
                               prom=props['prominences']))

    # pass 1: valley between the two most prominent peaks
    vrows = []
    for c in hcells:
        if len(c['peaks']) < 2:
            continue
        top = np.sort(c['peaks'][np.argsort(c['prom'])[-2:]])
        if centers[top[1]] - centers[top[0]] < min_sep:
            continue
        valley = top[0] + np.argmin(c['Hs'][top[0]:top[1] + 1])
        vrows.append((c['zcen'], c['mcen'], c['I'].sum(), centers[valley],
                      centers[top[0]], centers[top[1]]))
    if len(vrows) < 4:
        raise ValueError(f'Only {len(vrows)} cells show a bimodal color '
                         'distribution; cannot fit the color cut.')
    vcells = Table(rows=vrows, names=['ZCEN', 'MSTARCEN', 'N', 'VALLEY',
                                      'PEAK_BLUE', 'PEAK_RED'])
    vcoeff, vcells['USED'] = _fit_plane(vcells, 'VALLEY', 'N', floor=dcolor)

    # pass 2: red-sequence ridge redward of the provisional boundary
    rrows = []
    for c in hcells:
        vline = vcoeff[0] + vcoeff[1] * (c['mcen'] - 10.) + vcoeff[2] * (c['zcen'] - 0.1)
        nred = np.sum(color[c['I']] > vline)
        red = c['peaks'][centers[c['peaks']] > vline]
        if nred < nred_min or len(red) == 0:
            continue
        ridge = red[np.argmax(c['Hs'][red])]
        rrows.append((c['zcen'], c['mcen'], c['I'].sum(), nred, centers[ridge]))
    if len(rrows) < 4:
        raise ValueError(f'Only {len(rrows)} cells have a measurable red '
                         'sequence; cannot fit the color cut.')
    rcells = Table(rows=rrows, names=['ZCEN', 'MSTARCEN', 'N', 'NRED', 'RIDGE'])
    rcoeff, rcells['USED'] = _fit_plane(rcells, 'RIDGE', 'NRED', floor=dcolor)

    coeff = rcoeff - np.array([delta, 0., 0.])
    fit = dict(valley_coeff=vcoeff, ridge_coeff=rcoeff, valley=vcells,
               ridge=rcells)
    return coeff, fit


def gr_cut_sample(survey='sv3', verbose=False):
    """BGS, LRG, and ELG galaxies used to fit the 0.1(g-r) color cut.

    Reads the bright and dark fastspec catalogs for ``survey`` and keeps
    objects that pass their own class's redshift-quality cut
    (``good_redshift_by_class``) and have a stellar mass and finite color.

    Returns
    -------
    astropy.table.Table
        Columns include Z, LOGMSTAR, and GR (rest-frame 0.1(g-r)).

    """
    cols = ['LOGMSTAR', 'ABSMAG01_SDSS_G', 'ABSMAG01_SDSS_R',
            'OII_3726_FLUX', 'OII_3729_FLUX',
            'OII_3726_FLUX_IVAR', 'OII_3729_FLUX_IVAR']
    chunks = []
    for program in ('bright', 'dark'):
        cat = read_fastspec(survey, program, specprod=DEFAULT_SPECPROD,
                            columns=cols, verbose=verbose)
        gz = good_redshift_by_class(cat, survey)
        keep = (gz['BGS'] | gz['LRG'] | gz['ELG']) & (cat['LOGMSTAR'] > 0)
        cat = cat[keep]
        cat['GR'] = gr_color(cat)
        chunks.append(cat[np.isfinite(cat['GR'])])
    cat = vstack(chunks)
    if verbose:
        print(f'Color-cut sample: {len(cat):,} galaxies')
    return cat


def gr_cut_coeff(survey='sv3', verbose=False, return_all=False):
    """Fit the star-forming / quiescent 0.1(g-r) cut (see ``fit_gr_cut``).

    Returns the cut coefficients (a0, a1, a2); with ``return_all=True`` also
    returns the ``fit`` dictionary from ``fit_gr_cut`` and the fitting sample.

    """
    cat = gr_cut_sample(survey, verbose=verbose)
    coeff, fit = fit_gr_cut(cat['LOGMSTAR'], cat['Z'], cat['GR'])
    fmt = '{:.3f} {:+.3f} (logM - 10) {:+.3f} (z - 0.1)'
    for label, key in (('valley', 'valley'), ('ridge', 'ridge')):
        print(f'0.1(g-r) {label:6s}: ' + fmt.format(*fit[f'{key}_coeff']) +
              f"  [{fit[key]['USED'].sum()} of {len(fit[key])} cells used]")
    print('0.1(g-r) cut   : ' + fmt.format(*coeff))
    if return_all:
        return coeff, fit, cat
    return coeff


def cigaleagn_quality(d, sfr=True):
    """Siudek et al. (CIGALE-AGN) quality mask.

    Always applies the recommended stellar-mass cut, 0.2 <= FLAG_MASSPDF <= 5.
    With ``sfr=True`` (for the SFR comparison) also applies the SFR cut,
    0.2 <= FLAG_SFRPDF < 5. The mass cut matters for the SFR comparison too:
    it removes ~100k objects whose SFR is pinned at exactly 1 Msun/yr
    (h=0.704), which pass the SFR cut but have FLAG_MASSPDF > 1e7. We do not
    apply the SFR cut to the mass comparison, because it preferentially
    removes quiescent galaxies (~70% of LRGs) and biases that sample.

    """
    fmass = np.asarray(d['FLAG_LOGMSTAR_CIGALEAGN'], dtype=float)
    good = (fmass >= 0.2) & (fmass <= 5.0)
    if sfr:
        fsfr = np.asarray(d['FLAG_LOGSFR_CIGALEAGN'], dtype=float)
        good &= (fsfr >= 0.2) & (fsfr < 5.0)
    return good


def jiyan_p1p3(log_nii_ha, log_sii_ha, log_oiii_hb):
    """Ji & Yan (2020) P1 and P3 line-ratio projections for AGN/SF separation.

    Parameters
    ----------
    log_nii_ha  : array_like
        log10([NII] 6584 / Hα)
    log_sii_ha  : array_like
        log10(([SII] 6716 + [SII] 6731) / Hα)
    log_oiii_hb : array_like
        log10([OIII] 5007 / Hβ)

    Returns
    -------
    P1, P3 : ndarray, ndarray
        P1 =  0.63 * log([NII]/Hα) + 0.51 * log([SII]/Hα) + 0.59 * log([OIII]/Hβ)
        P3 = −0.46 * log([NII]/Hα) − 0.37 * log([SII]/Hα) + 0.81 * log([OIII]/Hβ)
    """
    n = np.asarray(log_nii_ha,  dtype=float)
    s = np.asarray(log_sii_ha,  dtype=float)
    o = np.asarray(log_oiii_hb, dtype=float)
    P1 =  0.63 * n + 0.51 * s + 0.59 * o
    P3 = -0.46 * n - 0.37 * s + 0.81 * o
    return P1, P3


def _good_fiberstatus(cat):
    """Fiber status mask: allow only RESTRICTED and VARIABLE bits."""
    from desispec.maskbits import fibermask
    okmask = fibermask.mask('RESTRICTED|VARIABLE')
    return (cat['COADD_FIBERSTATUS'] & okmask) == cat['COADD_FIBERSTATUS']


def good_redshift_by_class(cat, survey, fiberstatus_cut=True, ignore_emline=False):
    """Per-class redshift-quality-and-membership masks.

    Matches the LSS redshift-quality cuts in desispec.validredshifts (as of
    2026-08), substituting FastSpecFit's OII_3726/3729 fluxes for the
    emlinefit OII_FLUX used there.

    Applies the standard DESI per-class cuts:
      BGS  — ZWARN==0 & DELTACHI2>40 & Z<0.8
      LRG  — ZWARN==0 & Z<1.5 & DELTACHI2>15 (LRG or LGE targeting bit)
      ELG  — OII S/N cut: log10(OII_FLUX * sqrt(OII_FLUX_IVAR)) > 0.9 - 0.2*log10(DELTACHI2)
              (falls back to ZWARN==0 & Z<1.6 & DELTACHI2>15 if OII columns are
              absent or ``ignore_emline=True``; no redshift-range cut is applied
              when the OII S/N cut is used)
      Other — ZWARN==0 & Z>0.001
    BGS, LRG, and ELG additionally require SPECTYPE!='STAR' & Z>0.001 to
    reject stellar contaminants in the target bits.

    Each returned mask already combines that class's targeting-bit
    membership with its own quality cut above. This matters because DESI
    targeting bits are not mutually exclusive (an object can carry, e.g.,
    both a BGS and an LRG bit, especially in SV3): naively taking a
    combined ``good_redshift`` mask (the OR of all four classes) and then
    slicing it by targeting bit alone can leak an object into the wrong
    class's sample if it passed a *different* class's cut. Use this
    function directly whenever building a per-class sample; use
    ``good_redshift`` only when a single combined "is this redshift good
    for its own class" mask is wanted.

    Parameters
    ----------
    cat : astropy.table.Table
        Must include ZWARN, DELTACHI2, Z, SPECTYPE, all survey-appropriate
        targeting bit columns, and COADD_FIBERSTATUS (if fiberstatus_cut=True).
    survey : str
        DESI survey flavor: 'sv1', 'sv3', 'main', or 'special'.
    fiberstatus_cut : bool
        Apply fiber status mask (RESTRICTED and VARIABLE bits allowed).
    ignore_emline : bool
        Skip the OII emission-line cut for ELGs and use DELTACHI2>15 instead.

    Returns
    -------
    dict of numpy.ndarray of bool, keyed by 'BGS', 'LRG', 'ELG', 'Other'
    """
    if survey == 'sv3':
        from desitarget.sv3.sv3_targetmask import desi_mask
        desi_col, bgs_col = 'SV3_DESI_TARGET', 'SV3_BGS_TARGET'
    elif survey == 'sv1':
        from desitarget.sv1.sv1_targetmask import desi_mask
        desi_col, bgs_col = 'SV1_DESI_TARGET', 'SV1_BGS_TARGET'
    elif survey in ('main', 'special'):
        from desitarget.targets import desi_mask
        desi_col, bgs_col = 'DESI_TARGET', 'BGS_TARGET'
    else:
        raise ValueError(f'Unknown survey {survey!r}; expected sv1, sv3, main, or special.')

    is_bgs = cat[bgs_col] != 0
    is_lrg = (cat[desi_col] & int(desi_mask['LRG'])) != 0
    if 'LGE' in desi_mask.names():
        is_lrg |= (cat[desi_col] & int(desi_mask['LGE'])) != 0
    is_elg = (cat[desi_col] & int(desi_mask['ELG'])) != 0

    good_fs = _good_fiberstatus(cat) if fiberstatus_cut else np.ones(len(cat), bool)

    # reject stellar contaminants in the BGS/LRG/ELG target bits
    nonstar = (cat['SPECTYPE'] != 'STAR') & (cat['Z'] > 0.001)

    good_bgs = (cat['ZWARN'] == 0) & (cat['DELTACHI2'] > 40) & (cat['Z'] < 0.8) & good_fs & nonstar
    good_lrg = (cat['ZWARN'] == 0) & (cat['Z'] < 1.5) & (cat['DELTACHI2'] > 15) & good_fs & nonstar

    _OII_COLS = ('OII_3726_FLUX', 'OII_3729_FLUX',
                 'OII_3726_FLUX_IVAR', 'OII_3729_FLUX_IVAR')
    has_oii = all(c in cat.colnames for c in _OII_COLS)
    if not ignore_emline and has_oii and is_elg.any():
        oii_flux = cat['OII_3726_FLUX'] + cat['OII_3729_FLUX']
        # cast to float64 before multiplying to avoid float32 overflow
        ivar_a   = np.array(cat['OII_3726_FLUX_IVAR'], dtype=np.float64)
        ivar_b   = np.array(cat['OII_3729_FLUX_IVAR'], dtype=np.float64)
        ivar_sum = ivar_a + ivar_b
        with np.errstate(divide='ignore', invalid='ignore'):
            # harmonic sum of IVARs; zero wherever both (or either) are zero
            oii_ivar = np.where(ivar_sum > 0, ivar_a * ivar_b / ivar_sum, 0.0)
            oii_snr = np.log10(oii_flux * np.sqrt(oii_ivar))
            dc2     = np.log10(np.maximum(cat['DELTACHI2'], 1e-10))
        good_elg = ((oii_flux > 0) & (oii_ivar > 0) &
                    (oii_snr > 0.9 - 0.2 * dc2) &
                    good_fs & nonstar)
    else:
        good_elg = ((cat['ZWARN'] == 0) & (cat['Z'] < 1.6) &
                    (cat['DELTACHI2'] > 15) & good_fs & nonstar)

    good_other = (cat['ZWARN'] == 0) & (cat['Z'] > 0.001)

    return {
        'BGS':   is_bgs & good_bgs,
        'LRG':   is_lrg & good_lrg,
        'ELG':   is_elg & good_elg,
        'Other': (~(is_bgs | is_lrg | is_elg)) & good_other,
    }


def good_redshift(cat, survey, fiberstatus_cut=True, ignore_emline=False):
    """Combined redshift-quality mask: any class's cut, OR'd together.

    Equivalent to ``functools.reduce(operator.or_,
    good_redshift_by_class(...).values())``. Useful for a general-purpose
    "does this object have a good redshift for its own target class"
    filter, but do *not* use this to build a per-class sample by slicing
    on a targeting bit afterward --- targeting bits are not mutually
    exclusive, so an object can survive this mask via one class's cut
    while carrying a different class's bit. Use
    ``good_redshift_by_class`` directly for that.

    See ``good_redshift_by_class`` for the cut definitions and parameters.

    Returns
    -------
    numpy.ndarray of bool
    """
    masks = good_redshift_by_class(cat, survey, fiberstatus_cut=fiberstatus_cut,
                                   ignore_emline=ignore_emline)
    return masks['BGS'] | masks['LRG'] | masks['ELG'] | masks['Other']


def good_galaxies(cat, survey=None, fiberstatus_cut=True):
    """Boolean mask: good redshift quality and successful stellar mass fit.

    When ``survey`` is provided, applies per-class DELTACHI2 and fiber status
    cuts via ``good_redshift``; otherwise falls back to ZWARN==0 & Z>0.001.
    Always requires LOGMSTAR > 0.
    """
    if survey is not None:
        good_z = good_redshift(cat, survey, fiberstatus_cut=fiberstatus_cut)
    else:
        good_z = (cat['ZWARN'] == 0) & (cat['DELTACHI2'] > 15.)
    if 'LOGMSTAR' in cat.colnames:
        good_z = good_z & (cat['LOGMSTAR'] > 0)
    return good_z


def make_class_cmap(color, lighten=0.3):
    """White-to-lightened-color colormap for Hess diagram backgrounds.

    The endpoint is mixed ``lighten`` fraction toward white so that the
    densest bins stay visually lighter than the full class color, keeping
    the contours (drawn in the full color) clearly visible.

    """
    from matplotlib.colors import LinearSegmentedColormap, to_rgb
    rgb = np.array(to_rgb(color))
    light_end = tuple((1.0 - lighten) * rgb + lighten * np.ones(3))
    return LinearSegmentedColormap.from_list('', ['white', light_end])


def plot_style(talk=True, font_scale=1.0, palette=None, style='white'):
    """Set seaborn plot style and return (sns, color_palette).

    Parameters
    ----------
    talk : bool
        Use 'talk' context (larger fonts) if True, else 'paper'.
    font_scale : float
        Additional font scaling factor on top of the context default.
    palette : str or list or None
        Color palette passed to seaborn.set_palette().  None keeps the
        seaborn default.  Pass 'colorblind' for the Okabe-Ito palette,
        which is safe for deuteranopia, protanopia, and tritanopia.
    style : str
        Seaborn axes style. Default 'whitegrid' (this module's historical
        default); pass 'ticks' to match fastspecfit.qa.qa_fastspec's own
        style (gridlines off, plain tick axes).

    Returns
    -------
    (sns, colors) : (module, list)
        The seaborn module (so callers need only one import) and the
        current color palette as a list.
    """
    import seaborn as sns
    sns.set(context='talk' if talk else 'paper',
            style=style, font_scale=font_scale)
    if palette is not None:
        sns.set_palette(palette)
    return sns, sns.color_palette()


def target_label(objmeta, coadd_type, redshift):
    """Compact, publication-style target-ID and redshift label.

    Unlike fastspecfit.qa._target_label (which is tuned for the production
    QA figure's legend), this is meant to be read at a glance in a paper
    figure: survey/program/healpix (or tile/night/fiber), TARGETID, and
    redshift, one per line.

    Parameters
    ----------
    objmeta : astropy.table.Row
        Object metadata; must contain SURVEY/PROGRAM/HEALPIX (coadd_type
        'healpix', 'custom', or 'uniqpix') or TILEID/NIGHT/FIBER (coadd_type
        'cumulative', 'pernight', or 'perexp'), plus TARGETID in all cases.
    coadd_type : str
        Coadd type, as returned by fastspecfit.io.DESISpectra.coadd_type.
    redshift : float
        Redshift, appended as a final "z=" line.

    Returns
    -------
    str
        Newline-joined label, ready to pass directly to ax.text().

    """
    if coadd_type in ('healpix', 'custom', 'uniqpix'):
        id_line = f"{objmeta['SURVEY']}/{objmeta['PROGRAM']}, healpix {objmeta['HEALPIX']}"
    elif coadd_type in ('cumulative', 'pernight', 'perexp'):
        id_line = f"tile {objmeta['TILEID']}, night {objmeta['NIGHT']}, fiber {objmeta['FIBER']}"
    else:
        errmsg = f'Unrecognized coadd_type {coadd_type}!'
        raise ValueError(errmsg)

    return '\n'.join([
        #id_line,
        f"TARGETID {objmeta['TARGETID']}",
        f"$z={redshift:.4f}$",
    ])


def read_fastphot(survey='main', program='dark', specprod=DEFAULT_SPECPROD,
                  columns=None, verbose=False):
    """Read the fastphot VAC for a given specprod, survey, and program.

    Joins the METADATA and SPECPHOT extensions into a single Table.
    For main-bright and main-dark the 12 nside=1 healpix files are stacked
    transparently.

    Parameters
    ----------
    survey : str
        Survey name: 'main', 'sv1', 'sv2', 'sv3', or 'special'.
    program : str
        Program name: 'bright', 'dark', or 'backup'.
    specprod : str
        Spectroscopic production name (default: 'loa'). Must be a key in
        _SPECPROD_CONFIG.
    columns : list of str, optional
        Columns to load from SPECPHOT. METADATA is always read in full.
        If None, all columns are returned.
    verbose : bool
        Print file names and row counts while reading.

    Returns
    -------
    astropy.table.Table
        Merged catalog. Key columns always present: TARGETID, SURVEY, PROGRAM,
        HEALPIX, RA, DEC, Z, ZERR, ZWARN, SPECTYPE, and all survey-specific
        targeting bit columns that exist in the file (DESI_TARGET, BGS_TARGET,
        MWS_TARGET, SV3_DESI_TARGET, etc.). Additional columns from SPECPHOT
        are controlled by ``columns``.

    """
    if columns is not None:
        columns = list(set(columns) | _DEFAULT_COLUMNS)
    return _read_vac(
        specprod, 'fastphot',
        ['METADATA', 'SPECPHOT'],
        survey, program, columns, verbose,
    )


def target_class(objmeta, survey='main'):
    """Best-effort one-word target-class label for a single object.

    Checked in DESI dark/bright program priority order (QSO > LRG > ELG >
    BGS > MWS); an object can carry more than one targeting bit (e.g. a BGS
    target that also passes the LRG cut), so this is a display label for
    figure titles, not a full membership classification -- use
    good_redshift_by_class for that.

    Parameters
    ----------
    objmeta : astropy.table.Row
        Metadata row with DESI_TARGET/BGS_TARGET/MWS_TARGET (or the
        survey-prefixed equivalents for sv1/sv3).
    survey : str
        DESI survey flavor: 'sv1', 'sv3', 'main', or 'special'.

    Returns
    -------
    str
        'QSO', 'LRG', 'ELG', 'BGS', 'MWS', or 'Target' if no bit matches.
    """
    if survey == 'sv3':
        from desitarget.sv3.sv3_targetmask import desi_mask
        desi_col, bgs_col, mws_col = 'SV3_DESI_TARGET', 'SV3_BGS_TARGET', 'SV3_MWS_TARGET'
    elif survey == 'sv1':
        from desitarget.sv1.sv1_targetmask import desi_mask
        desi_col, bgs_col, mws_col = 'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'SV1_MWS_TARGET'
    elif survey in ('main', 'special'):
        from desitarget.targets import desi_mask
        desi_col, bgs_col, mws_col = 'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET'
    else:
        raise ValueError(f'Unknown survey {survey!r}; expected sv1, sv3, main, or special.')

    desi_target = int(objmeta[desi_col])
    is_qso = (desi_target & int(desi_mask['QSO'])) != 0
    is_lrg = (desi_target & int(desi_mask['LRG'])) != 0
    if 'LGE' in desi_mask.names():
        is_lrg |= (desi_target & int(desi_mask['LGE'])) != 0
    is_elg = (desi_target & int(desi_mask['ELG'])) != 0
    is_bgs = int(objmeta[bgs_col]) != 0
    is_mws = int(objmeta[mws_col]) != 0

    for is_class, label in ((is_qso, 'QSO'), (is_lrg, 'LRG'), (is_elg, 'ELG'),
                            (is_bgs, 'BGS'), (is_mws, 'MWS')):
        if is_class:
            return label
    return 'Target'


def fetch_cutout_cached(objmeta, cutoutdir, layer='ls-dr9', pixscale=0.262, timeout=15):
    """Fetch a Legacy Survey grz cutout JPEG and keep it on disk permanently.

    fastspecfit.qa._fetch_cutout downloads to a temp JPEG and always deletes
    it after loading (it's built for one-shot batch QA, not for reuse), so
    scripts that call it directly re-download on every run. This instead
    caches to ``cutoutdir/cutout-<targetid>.jpeg`` and only fetches if that
    file isn't already there.

    Parameters
    ----------
    objmeta : astropy.table.Row
        Metadata row with RA, DEC, TARGETID.
    cutoutdir : str
        Directory for the cached JPEG (must already exist).
    layer : str
        Legacy Survey viewer layer (e.g. 'ls-dr9').
    pixscale : float
        Pixel scale in arcsec/pixel.
    timeout : float
        Network timeout in seconds for the single fetch attempt.

    Returns
    -------
    numpy.ndarray
        RGB image array, shape (height, width, 3).
    """
    import urllib.request
    import matplotlib.image as mpimg

    width = int(30 / pixscale)
    height = int(width / 1.3)
    cutoutjpeg = os.path.join(cutoutdir, f"cutout-{objmeta['TARGETID']}.jpeg")

    if not os.path.isfile(cutoutjpeg):
        url = ('https://www.legacysurvey.org/viewer/jpeg-cutout?ra=' +
              f'{objmeta["RA"]}&dec={objmeta["DEC"]}&width={width}&height={height}&layer={layer}')
        with urllib.request.urlopen(url, timeout=timeout) as response:
            with open(cutoutjpeg, 'wb') as f:
                f.write(response.read())

    return mpimg.imread(cutoutjpeg)
