# Figure and Table Generation

Paper figures are written to `tex/figures/` and committed for inclusion in the paper. Reference and diagnostic figures that are not in the paper are written to `qa/` (see `qa/README.md`).

Most scripts run on NERSC (Perlmutter), where the DESI software environment and VAC paths are available. Run everything from the repo root. Each script's module docstring documents its exact invocation and options.

```bash
conda env create -f environment.yml   # once
conda activate fastspec
```

## build-figures.py

The primary figure-generation script; each flag builds one figure (`--all` rebuilds every figure in the paper). Stellar masses and SFRs are stored with **h=1** and a **Chabrier IMF**. Each figure function takes its own `survey` argument (default `'sv3'`), so switching one figure to the main survey is a one-line change to that function's default.

| Flag | Output |
|---|---|
| `--sps-models` | `tex/figures/sps-models.pdf` |
| `--compare-mstar` | `tex/figures/compare-mstar-fastphot.pdf` |
| `--compare-mstar-external` | `tex/figures/compare-mstar-external.pdf` |
| `--compare-emlines-external [--pull-denom {quadrature,ref,ext}]` | `tex/figures/compare-emlines-external.pdf`; `qa/compare-emlines-external-{ref,ext}.pdf` for the non-default pull denominators |
| `--compare-sfr-external` | `tex/figures/compare-sfr-external.pdf` |
| `--mstar-redshift` | `tex/figures/mstar-redshift.pdf` |
| `--bpt-agn` | `tex/figures/bpt-agn.pdf` |
| `--compare-vdisp [--vdisp-two-panel]` | `tex/figures/compare-vdisp.pdf` |
| `--vac-compare-table` | `tex/tables/table-vac-compare.tex` (Appendix B; requires `prepare-external.py --vacs`) |
| `--gr-mstar-cut-qa` | `qa/gr-mstar-cut-qa.pdf` |
| `--cosmos2020` | `qa/compare-cosmos2020.pdf` |
| `--ewoii-dn4000` | `qa/ewoii-dn4000.pdf` |
| `--sfr-mstar-bgs [--flag-agn]` | `qa/sfr-mstar-bgs.pdf` |
| `--broadhalpha-fraction-bgs [--binwidth DEX]` | `qa/broadhalpha-fraction-bgs.pdf` |

Common options: `--specprod` (default `loa`) and `--verbose`.

## Other scripts

| Script | Purpose |
|---|---|
| `example-qa.py` | Section 4 illustrative-example figures (`tex/figures/example-{bgs,elg,qso}.pdf`); fits a single target standalone from the files in `data/` |
| `linemasker-qa.py` | Line-masking / patch-fitting figure (`tex/figures/linemasker.pdf`, Section 3.3) |
| `pipeline-diagnostics.py` | Section 5.5 diagnostics: `--wavecal` and `--fluxcal` build samples in `external/`; `--wavecal-figure` writes `tex/figures/wave-residuals.pdf` |
| `prepare-external.py` | Cross-matches external catalogs (`--zouhu`, `--cigaleagn`, `--gswlcx2`, `--cosmos2020`, `--fpcatalog`, `--emlinefit`) and earlier VACs (`--vacs`) to the FastSpecFit VAC; see `external/README.md` |
| `gswlc-sfr-diagnostics.py` | Reference figure on the FastSpecFit vs. GSWLC-X2 SFR offset (`qa/`, Section 5.4) |
| `vdisp-qa.py` | Draft velocity-dispersion QA figure (`qa/`); inputs in `data/` |
| `qso-redshift-systematics.py` | QSO emission-line redshift systematics, for a separate proposal (`qa/`) |
| `util.py` | Shared catalog readers and plotting helpers |
