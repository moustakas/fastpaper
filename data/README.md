# Small Reference Data

Small per-target files needed to regenerate specific paper figures without a
full NERSC catalog: DESI coadd + Redrock spectra (`coadd-*.fits` /
`redrock-*.fits`, named `{survey}-{program}-{healpix}`), one Legacy Surveys
Tractor photometry tree (`north/`, `south/`), and a handful of files produced
by in-progress QA scripts. FITS files are gitignored (see repo root
`.gitignore`); this README is the manifest of what each committed pair is
for, since the filenames themselves only carry SURVEY-PROGRAM-HEALPIX, not
TARGETID or which figure they feed.

**Only files actually wired into a current figure/script belong directly in
`data/`.** Anything not currently used, but not yet safe to delete, goes in
[`extras/`](#extras) instead of being removed outright.

## Manifest

| File pair | TARGETID | z | Figure / script |
|---|---|---|---|
| `{coadd,redrock}-main-bright-17366` | 39627670102744387 | 0.3292 | Fig. 2, `linemasker.pdf` (§3.3) — `code/linemasker-qa.py` |
| `{coadd,redrock}-main-bright-27851` | 39628087574398347 | 0.0813 | Fig. 3, `example-bgs.pdf` (§4.1) — `code/example-qa.py` |
| `{coadd,redrock}-main-dark-29973` | 39627663165362386 | 0.9449 | Fig. 4, `example-elg.pdf` (§4.2) — `code/example-qa.py` |
| `{coadd,redrock}-main-dark-29972` | 39627657142340435 | 2.0731 | Fig. 5, `example-qso.pdf` (§4.3) — `code/example-qa.py` |
| `redrock-main-dark-17289`, `f-main-dark-17289-39627496647296130`, `qa-vdisp-chi2scan-39627496647296130.ecsv`, `cutout-39627496647296130.jpeg` | 39627496647296130 | — | Draft vdisp-section figure, not yet in `ms.tex` — `code/vdisp-qa.py` |

The exact regeneration command for each is given inline in `CLAUDE.md`
("Running the Analysis Scripts").

When you swap a target for one of these figures (as happened for Fig. 3,
BGS: 31152 → 40945 → 27851 in September 2026), update the TARGETID/z/HEALPix
here, in `CLAUDE.md`, and in the figure's caption in `ms.tex` together, and
move the old file pair into `extras/` below.

## `extras/`

Coadd/Redrock pairs downloaded at some point (candidate targets, earlier
choices for one of the Fig. 3-5 examples, etc.) that no figure or script
currently reads. Kept around rather than deleted in case they're useful
again; like everything else in `data/`, these are gitignored, so deleting
them (if space becomes an issue) has no git-history consequences.
