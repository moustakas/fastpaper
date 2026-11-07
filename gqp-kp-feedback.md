# GQP DR2 Key Paper: feedback and follow-ups

Cross-comparison of `gqp_dr2_kp.tex` ("DESI DR2 Results XX: The Physical Properties of the DESI Cosmological Tracer Populations") against `tex/ms.tex`, 2026-11-06. This file replaces `provabgs-kp-notes.md` (Nov 1). The KP draft is byte-identical to the one read then.

Contents:

1. Wiki comments for the KP authors (ready to paste)
2. ToDo list for `ms.tex`
3. PROVABGS: status and plan (carried over from the Nov 1 notes)

KP line numbers below are `.tex` source lines, not the PDF line numbers, so the wiki comments also give the section and a short quote.

---

## 1. Wiki comments for the KP authors

```text
=== John Moustakas ===

Thanks for a terrific paper, and for putting the FastSpecFit catalogs to such good use. My comments mostly concern how FastSpecFit is described, plus a few other things I noticed along the way.

 * '''Naming, throughout.''' FastSpecFit is a single code with two fitting modes, `fastspec` (spectrum plus photometry) and `fastphot` (photometry only), and the corresponding VACs are the `fastspec` and `fastphot` catalogs. The paper instead treats "!FastSpecFit" and "!FastPhot" as two separate codes (abstract, Table 2 columns, Sec. 3.2.3, Sec. 6), and also uses `\textsc{FastPhot}`, "FSF", and "FP" (e.g., Sec. 6.1.3, "3,168,984 for the FSF catalogue, and 2,892,212 for the FP catalogue"). Could you use "!FastSpecFit (`fastspec`)" and "!FastSpecFit (`fastphot`)" at first use, and `fastspec`/`fastphot` (or "!FastSpecFit spectrophotometric/photometric") thereafter, and drop FSF and FP? This matches the convention in the !FastSpecFit paper. Separately, the `\fsf` macro ends in `~`, which produces a double space before "Loa-v1.0" in the Appendix C title and text (tex L1578, L1585).
 * '''Table 2 (SED-fitting summary), stellar library.''' "C3K`$_a$`" should be `C3K_R3K`, a custom version of the C3K theoretical library at solar metallicity with no `$\alpha$`-enhancement (see https://github.com/moustakas/c3k).
 * '''Cosmology of the !FastSpecFit masses.''' `LOGMSTAR` and `SFR` in both VACs are computed with `$h=1$` (`$H_0=100$` km/s/Mpc), i.e., in units of `$h^{-2}\,M_\odot$` and `$h^{-2}\,M_\odot$`/yr. For a Planck cosmology this is a `$-2\log_{10}h \approx +0.34$` dex correction. The sign and size of your !FastSpecFit-minus-CIGALE offsets suggest you did apply it, but the paper never says so, and the stellar mass function, the HOD stellar-mass bins, and the BPT/KEx mass bins all depend on it. Could you confirm and state the correction explicitly (e.g., in Sec. 2.3.2)?
 * '''Sec. 2.3.2 and Appendix C, VAC versions.''' Please state which !FastSpecFit VAC is used in each analysis, by version and software tag: `Loa/v1.0` (!FastSpecFit 3.1.4, templates `chabrier-2.0.0`) or `Loa/v2.0` (!FastSpecFit 3.6.1, templates `chabrier-2.2.0`). Only `v2.0` is documented in the !FastSpecFit paper.
 * '''Sec. 2.3.3 and Sec. 5, AGN/Galaxy Classification VAC.''' Which !FastSpecFit version are the emission-line diagnostics in this VAC built on? If it is `Loa/v1.0`, note that the emission-line flux and EW inverse variances in that catalog are known to be wrong (fixed in !FastSpecFit 3.2.0), which would directly affect the S/N>3 cuts behind the BPT and KEx AGN fractions.
 * '''Sec. 3.2.3, !FastSpecFit paragraph (tex L520).''' Two small corrections. (1) !FastSpecFit does not use the spectroperfectionism method itself; it forward-models each spectrum with the resolution matrix that spectroperfectionism produces (Bolton & Schlegel 2010; Guy et al. 2023). (2) The `fastphot` results are independent of the DESI spectrum except for its redshift.
 * '''Sec. 3.2.4, aperture-correction footnote (tex L617).''' The column is `FLUX_SYNTH_{G,R,Z}`, not `SYNTH_FLUX`. Also, for the H`$\alpha$` SFR aperture correction (tex L625), note that the `fastspec` catalog already provides this quantity as `APERCORR_R` (and the `$grz$` median as `APERCORR`), which might be simpler and more consistent than recomputing it.
 * '''Sec. 3.2.4, BGS stellar masses (tex L627).''' "Slightly lower" undersells a 0.25-0.35 dex offset; maybe "systematically lower"?
 * '''Sec. 6.2.1 (tex L1207).''' It is true that the `fastspec` and `fastphot` stellar-mass uncertainties reflect only the propagated measurement noise, but the reason is not really "the limited number of templates." The fit is maximum-likelihood over a fixed model (five non-parametric SFH bins, solar metallicity, a fixed attenuation slope), and the Monte Carlo resampling only captures the noise for that model, not the systematic uncertainty from the model assumptions themselves. Something like "do not include systematic uncertainties from the SPS, SFH, and dust assumptions" would be more accurate.
 * '''Sec. 7.2, LRG galaxy-halo connection (tex L1296).''' "The CIGALE results are systematically offset towards higher stellar mass by about 0.2 dex" has the wrong sign. Sec. 3.2.4 (tex L627, L683) and Fig. (LRG stellar-mass comparison) both find CIGALE masses ~0.2 dex ''lower'' than !FastSpecFit.
 * '''Appendix C.''' A few points.
   * The appendix mixes up the LSS catalog versions (v1 versus v2, which differ in the BGS and LRG redshift ranges) with the !FastSpecFit VAC versions (`Loa/v1.0` versus `Loa/v2.0`, which differ in code, templates, and SFR definition). The LSS redshift extension has no bearing on the !FastSpecFit measurements themselves, since !FastSpecFit fits every extragalactic target independent of the LSS catalogs. Please separate the two and give the !FastSpecFit versions as above.
   * "Different integration windows for SFHs" should be "different SFR averaging timescales."
   * The SFR difference is a change in definition, so "!FastSpecFit Loa v1 underestimates the SFR" is not quite right. `Loa/v1.0` reports the instantaneous SFR of the best-fitting model, to which only the youngest (0-30 Myr) age bin contributes; `Loa/v2.0` reports the mean SFR over the last 100 Myr. The !FastSpecFit paper (Appendix B) quantifies the same effect relative to the DR1 and EDR VACs and recommends against combining SFRs across these versions.
 * '''PROVABGS stellar-mass offset.''' The PROVABGS-minus-!FastSpecFit offset in the BGS corner figure (+0.208 dex, NMAD 0.182) differs by ~0.37 dex from the commented-out pairwise table (-0.161, scatter 0.346, "case 0"), while the other !FastSpecFit pairs agree between the two. Could you check which one is right? Possible culprits are the fiber-correction case, the `$h$` normalization (PROVABGS uses `$H_0=67.8$`), or formed versus surviving stellar mass (the !FastSpecFit masses are surviving masses, including remnants).
 * '''Other things I noticed.'''
   * tex L603: `\textttt` is a typo (compile error).
   * Sec. 3.2.1 (tex L492): "{\tt BGS\_FAINT} or {\tt BGS\_FAINT}" should presumably include `BGS_BRIGHT`.
   * Table 2, CIGALE_MS AGN component: lists Dale et al. (2014), but the DR1 VAC page gives Fritz et al. (2006) for the AGN, with Dale et al. (2014) for dust emission.
   * Sec. 7.3 (tex L1306): "increasing by ~20% over the low-redshift \msun" is garbled.
   * Sec. 7.4 (tex L1319): refers to `fig:GHC_ELG` but should be `fig:GHC_QSO`.
   * The stellar mass function spans "5 decades" in the abstract and Sec. 6 but "six decades" in the introduction and Sec. 8.
   * The ELG redshift success rate is 72.15% in Sec. 8 but 72.24% in Table 1; the `BGS_BRIGHT` success rate is 98.81% in Sec. 6.2.3 but 98.83% in Table 1.
   * The median bolometric luminosity is `$10^{45.3}$` erg/s in Sec. 8 but `$10^{45.41}$` in Sec. 4.2.
```

---

## 2. ToDo list for `ms.tex`

One at a time, each to be reviewed before it lands. Line numbers as of today.

1. **Fill the PROVABGS `\todo` (Sec. 5.4, ~L3318).** Point to the KP for the broader comparison. Possible wording: "For a broader comparison of the stellar masses and SFRs from several SED-fitting codes applied to DESI, including PROVABGS, see DESI Collaboration et~al. (in preparation)." Closes the Saulder CWR item.
2. **Cite the KP and fix a typo in the introduction (~L617).** "luminosity and stellar functions" should be "luminosity and stellar mass functions", and add "DESI Collaboration et~al., in preparation" to that citation.
3. **Add an SFR caveat (Sec. 3.6, after eq. `sfr`, where the commented-out `\todo` sits).** The paper has no SFR validation; the KP is the main independent test. Possible wording: "Because our photometry does not reach the rest-frame UV or far-IR, these SED-based SFRs are only weakly constrained, especially for quiescent galaxies. DESI Collaboration et~al. (in preparation) find a scatter of at least 0.3~dex between the SFRs inferred by different codes, even for star-forming galaxies, so we recommend that users vet them carefully."
4. **Add an uncertainty caveat (Sec. 3.7).** The KP criticizes our mass uncertainties. Possible wording: "These uncertainties reflect only the propagated measurement noise for our fixed template set; they do not include systematic uncertainties from the SPS, SFH, and dust assumptions (\S\ref{ssec:templates})."
5. **Repeat the h=1 convention where users will see it (Sec. 5.1 or 5.2).** It's stated only once, in the introduction (~L649). Possible wording: "As noted in \S\ref{sec:intro}, \texttt{LOGMSTAR} and \texttt{SFR} assume $h=1$; for $H_0=67.4~\kms~\mathrm{Mpc}^{-1}$, add $-2\log_{10}h\approx0.34$~dex."
6. **Note the high-mass consequence of the `fastspec`/`fastphot` offset (Appendix A, after the "does not exceed 0.05~dex" sentence).** Possible wording: "Even an offset this small matters where the stellar mass function falls off exponentially; DESI Collaboration et~al. (in preparation) find that the \fastspec{} and \fastphot{} stellar mass functions differ by $\approx0.2$~dex in number density at $\log_{10}(\Mstar/\msun)\approx11.5$."
7. **Zou AGN component (Sec. 5.4, ~L3296).** The KP says the DR2 Zou catalog includes a Fritz+06 AGN component. Possible wording: "\citet{siudek24a} and the DR2 Zou et~al. VAC also include an AGN component \citep{fritz06a}." Wait for Hu to confirm.
8. **ELG mass scatter (Sec. 5.4, ~L3310).** Currently blamed on low-S/N photometry alone; the KP's `fastphot`-versus-CIGALE comparison (photometry on both sides) still shows 0.43--0.60~dex. Possible wording: "likely reflecting both the lower-S/N broadband photometry ... and differences in the model assumptions."
9. **Where the nebular templates are used (Sec. 3.1, ~L974).** The text says only K-corrections, but `continuum.py` also uses the line-including templates in the continuum fit itself (`fastspec` full fit at ~L2133, `fastphot` at ~L1618), where masked line pixels leave them affecting only the synthetic photometry. The Buhler item 6 response in `comments-26nov06.txt` has the same problem; the earlier Nov 4 response had it right. Possible wording: "We use the models with emission lines in the continuum fit, where masking line-contaminated pixels limits their effect to the synthetic photometry, and when computing K-corrections."
10. **Camera wavelength ranges (Sec. 2.1, ~L679--686).** B and R come from the instrument paper, but the Z upper edge (9824~\AA) comes from the pipeline extraction grid (the instrument paper gives 9800~\AA). Either use the instrument ranges throughout (B 3600--5930, R 5600--7720, Z 7470--9800~\AA) or the Guy+23 extraction ranges (B 3600--5800, R 5760--7620, Z 7520--9824~\AA). The extraction ranges match the 3600--9824~\AA{} used elsewhere in the paper.
11. **Document the Loa/v1.0 bugs (Appendix B, commented-out paragraph at ~L4033).** The VAC history table lists `Loa/v1.0` as a 2027 public release, and the KP used it, so the note about its emission-line inverse-variance and clipping bugs probably needs to come back out of the comment.
12. **Optional: cite the KP for AGN fractions (Sec. 6.1, BPT paragraph, ~L3537).** Its BPT, KEx, and WISE AGN fractions are built on our line fluxes.

---

## 3. PROVABGS: status and plan

Adding PROVABGS to Sec. 5.4 and Fig. 8 (`compare-mstar-external`) is on hold until I have access to the DR2/Loa catalogs. ToDo item 1 above is the stopgap.

Why it would still be worth adding: PROVABGS shares FastSpecFit's FSPS+MIST+Chabrier setup and has a flexible SFH, so it directly tests the Sec. 5.4 explanation that SFH priors drive the 0.15--0.30 dex offset from the CIGALE catalogs. If the true offset doesn't fit that explanation, the conclusion of that paragraph changes too.

### What's known

| Item | Value | Source |
|---|---|---|
| Model | FSPS, MIST, MILES+BaSeL, Chabrier; 4 SFH bases from NMF of Illustris histories plus a burst; 2 metallicity-history bases | Hahn+23 |
| Dust and lines | Charlot & Fall two-component dust, no dust emission; emission lines masked | `hahn24a` |
| Data | Joint spectrum + *grz* fit; single-band (*r*) fiber correction as a free parameter with a prior ("C1") | KP Appendix D |
| Cosmology | H0 = 67.8 in the paper (Planck13 in the code), so h=1 means -0.338 dex | `hahn24a`, `models.py` |
| Coverage | BGS only: no PROVABGS in the KP LRG or ELG panels; emulator covers z < 0.6 | KP figures, GitHub README |

**Public EDR catalog (fallback):** v1.0, SV3 BGS only, HDF5, matched on `TARGETID`.

- File: `/global/cfs/cdirs/desi/public/edr/vac/edr/provabgs/v1.0/BGS_ANY_full.provabgs.sv3.v0.hdf5`
- Redshift column: `Z_HP`
- Mass columns: `provabgs_logMstar` (100 posterior samples), `provabgs_logMstar_bf` (best fit)
- `PROVABGS_MCMC` includes `log(M_form)`, which the VAC page says differs from log stellar mass.
- Needs `h5py`, which isn't in `environment.yml`.

### Open problems

1. **Sample mismatch.** The only Loa run the KP describes is a test catalog of about 50,000 BGS galaxies, apparently drawn from the main-survey LSS catalogs (Loa LSS v1). Fig. 8 uses SV3, so the overlap may be small.
2. **Formed vs. surviving mass.** The `log M*` parameter normalizes the SFH, so it looks like formed mass (my inference). FastSpecFit reports surviving mass including remnants. Using the wrong column would bias PROVABGS high by roughly 0.2 dex, assuming about 40% mass return.
3. **Unexplained KP offset.** Raised in the wiki comments above. Don't quote the KP "largest masses" claim until it's settled.

### Questions for the catalog owner

- Which file and version is the Loa catalog? Does it cover SV3 or the full main survey, or only the test sample?
- Which column holds surviving mass, and does it include remnants?
- Which fiber-correction case does it use?
- What H0 was used, and are the masses already rescaled?
- Which SFR timescale applies? The KP table says 100 Myr, but the code default and both papers use 1 Gyr. This only matters if SFRs are compared later.

### Plan once the catalog is available

- **Catalog choice:** use Loa if it covers SV3 or the full main survey; otherwise use public EDR, which matches Fig. 8's SV3 sample (the Iron row already mixes data releases). Avoid the test sample unless it will be released.
- **`code/prepare-external.py`:** add `prepare_provabgs()` with per-specprod paths, like `prepare_zouhu()`, writing `LOGMSTAR_PROVABGS` and `LOGMSTAR_ERR_PROVABGS` from the posterior median and 16th--84th percentiles. Document it in `external/README.md`.
- **`code/build-figures.py`:** add a BGS-only fourth row to `compare_mstar_external()`; the 3x3 grid is hard-coded.
- **`tex/ms.tex`:** rewrite the "three CIGALE-based" and "all three... parametric SFH" text and the Fig. 8 caption.
- **Citations (BibDesk):** `hahn24a` is the stellar mass function paper, which fits the EDR catalog but not the method. The method paper (ApJ 945, 16; arXiv:2202.01809) isn't in `refs.bib`, and `hahn23a` is already BGS target selection, so probably `hahn23b`. Also add Kwon+23 for the emulator.

### Don't quote KP numbers

Its BGS offsets (-0.22 for Zou, -0.36 for Siudek) are close to ours, but its ELG Siudek offset is about 0.09 dex against our 0.28. Its samples and catalog runs also differ: its figure names suggest the 5-band Zou fit, while we use the 15-band one.

### Sources

- [Hahn+23, PROVABGS mock challenge](https://arxiv.org/abs/2202.01809)
- [Hahn+24 (`hahn24a`)](https://arxiv.org/abs/2306.06318)
- [EDR PROVABGS VAC page](https://data.desi.lbl.gov/doc/releases/edr/vac/provabgs/)
- [DR1 CIGALE VAC page](https://data.desi.lbl.gov/doc/releases/dr1/vac/cigale/)
- [PROVABGS GitHub repository](https://github.com/changhoonhahn/provabgs)
