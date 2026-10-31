# AD1 (visible methane and ammonia absorption of giant planets): sources check

Filed by the TODO thread on 2026-10-17 at about 1:58 PM EDT (as this thread's clock reads). The data folder is outside the repository: the build fetches it (gap plan section 4.11). Summary and the item notes are in docs/TODO.md (GEN.531, GEN.422).

Date: 2026-10-17, checked between 1:40 PM and 1:50 PM EDT as the lane clock reported (lane clocks drift). Read-only research, no repo change. Boss's pasted note (AI-written, 1:38 PM EDT) was treated as a list of unverified leads. Downloaded data and its README: `/mnt/project-files/data/ad1-methane/` (README in the same folder). AD1 as the plan states it: `docs/design/research/gap-filling-plan.md` row AD1, "band strengths of CH4, NH3, H2O in giants, Karkoschka tables not extracted", feeds GEN.422 and atmosphere colour; stand-in is neighbouring-band strengths with 10 percent.

Labels: **READ** = the file or paper text was opened here; **2ND-HAND** = read in a paper or code that cites the original; **NOT REACHED** = no reachable page gave it. Verdicts on the note's claims: verified, partly verified, contradicted, not found.

## 1. Verdict

1. **AD1 for methane can be filled from free data, and the data downloaded.** Two independent routes, both READ and consistent with each other:
   - **The K&T 2010 band model itself**, as k-tables by Patrick Irwin (Oxford), 14 May 2010: 3689 wavelengths, 400 nm to 5.0 micron, κ at 100, 198 and 296 K plus the line-spacing ratio, and correlated-k values on a 10 pressure by 11 temperature grid (1e-8 to 10 atm, 50 to 300 K). The count matches the paper's own abstract ("for each of ~3700 wavelengths, we present a temperature dependence").
   - **Karkoschka 1994 and 1998 as NASA PDS tables** (300 to 1050 nm, one effective cold methane spectrum, plus full-disk albedos of Jupiter, Saturn, Uranus, Neptune and Titan).
2. **The note's catalogue id is wrong.** The Karkoschka data is PDS volume **gbat_0001**, not go_0016 (go_0016 to go_0022 are Galileo orbiter images). It holds both the 1993 data (1994 paper) and the 1995 data (1998 paper); the NMSU web pages cite only the 1994 paper.
3. **Not covered by any free source found: ammonia and water coefficients in the visible.** K&T 2010 is methane only. Karkoschka 1998 gives ammonia only as band depths in Jupiter's albedo (equivalent widths 0.05, 0.32, 0.62, 0.81 nm for the 552, 645, 790 and 930 nm bands, read from Fig. 5); no coefficient table. Water: the 1994 suggestion of a 930 nm water band was withdrawn in 1998 (the band is probably ammonia). AD1 stays partly open for NH3 and H2O.
4. **Gaps in the methane data:** 5.0 to 5.5 micron (the paper says its data reach 5.5; the Oxford tables stop at 2000 cm-1 = 5.0 micron); below 450 nm every coefficient sits at a 1e-5 floor in the tables (65 wavenumbers), so the true weak absorption there is not resolved; the Oxford page says its K&T tables are "not yet fully verified"; Irwin 2016 reports K&T 2010 models 1.4 to 1.8 micron less well than the WKMC-80K line list (outside the visible).
5. **Several of the note's details are wrong** (section 4): Neptune paper citation, the eleventh data set, PMC5105178 is not a different paper, one pressure level in the filter table, the Uranus methane range.

## 2. What could be reached from here (measured in this lane)

| Host | Result |
|---|---|
| `atmos.nmsu.edu`, `pds-atmospheres.nmsu.edu` | reachable; full PDS tree downloadable |
| `users.ox.ac.uk/~atmp0035` (Irwin's page) | reachable; all k-table files downloaded (4 files, 10 to 19 MB each) |
| `github.com` through `git clone` / `git ls-remote` | reachable (NEMESIS cloned, 1.1 GB, not kept); the github.com web pages and `api.github.com` return 403; `raw.githubusercontent.com` reachable |
| `arxiv.org` (abs and PDF) | reachable; `export.arxiv.org/api/query` returned an empty answer |
| `api.crossref.org`, `doi.org` | reachable; doi.org redirects to ScienceDirect |
| `eutils.ncbi.nlm.nih.gov` | reachable and returns the **full text XML of PMC articles**; the `pmc.ncbi.nlm.nih.gov` web pages return a reCAPTCHA wall |
| `www.ebi.ac.uk/europepmc/webservices/rest/` | reachable (full text XML, supplementary file bundles); `europepmc.org/articles/...` returns 403 |
| `api.semanticscholar.org` | reachable, rate limited (429 after a few calls); abstracts of Elsevier papers elided |
| `experts.arizona.edu` | reachable; carries the **full abstracts** of K&T 2009 and K&T 2010 |
| `hitran.org`, `exomol.com`, `archive.stsci.edu` | front pages reachable (not used); `www.mast.stsci.edu` no connection; Zenodo web search reachable, `api.zenodo.org` no connection |
| refused | ScienceDirect (403), ADS (405), MDPI (403, so Wong 2023 not read), science.org (403), NASA PSG (`psg.gsfc.nasa.gov` answers "unsupported command / browser"; its data not reachable), `lpl.arizona.edu/~erich` (404) |

Note for other lanes: this list differs from the atmosphere check of 10-17 (that lane found github.com, PMC, NIST and Crossref refused); the proxy rules seem to have changed between lanes. The PMC route that works is `eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=PMCnnnn&rettype=xml` (full text, CC-licensed articles).

## 3. What was downloaded and what it covers

All under `/mnt/project-files/data/ad1-methane/`; the README there has source URLs, sizes, sha256 and terms.

### 3.1 PDS volume gbat_0001 (READ)

| File | Rows | Range | Resolution | Source paper |
|---|---|---|---|---|
| `data/1993.tab` | 1750 | 300.4 to 1000.0 nm | 0.4 nm sampling, 1 nm | Karkoschka 1994 (ESO, 23 to 26 July 1993) |
| `data/1995low.tab` | 1875 | 300.4 to 1050.0 nm | 0.4 nm sampling, 1 nm | Karkoschka 1998 (ESO, 6 to 10 July 1995) |
| `data/1995high.tab` | 4750 | 520.1 to 995.0 nm | 0.1 nm sampling, 0.4 nm | Karkoschka 1998 (Jupiter, Saturn, Uranus) |

Columns: vacuum wavelength (nm), air wavelength (nm), methane absorption coefficient (**1/(km-amagat)**), then full-disk albedos (Jupiter at 9.8 deg phase in 1993 and 6.8 deg in 1995; Saturn; Uranus and Neptune geometric; Titan). Methane k is zero below 405.6 nm (1993) and 401.2 nm (1995 low), and reaches 44 (1993), 47 (1995 low) and 51.6 (1995 high) 1/km-am at the strongest bands. It is **one effective cold-temperature spectrum** ("the corresponding temperature may be near 100 K", 1994 paper), so no pressure or temperature dependence.

Sample rows (copied data, `1995low.tab`: vac nm, air nm, k, Jupiter, Saturn, Uranus, Neptune, Titan):
```
 300.4  300.31 .0000 .2128 .2532 .5306 .6933 .0431
 675.2  675.00 .0169 .5456 .5903 .3694 .2720 .2930
1049.6 1049.29 .2599 .4007 .5402 .0710 .0441 .1779
```
The papers themselves (1994 and 1998, ASCII) are in the same volume's `document/` folder, so the numbers in section 4 for those two were checked against the full text.

### 3.2 Oxford K&T 2010 k-tables (READ)

Header: "k-table for Karkoschka and Tomasko (2010) band model, Patrick Irwin, AOPP, University of Oxford 14th May 2010". Files `ch4_karkoschka_vis.par` (19300 to 25000 cm-1, 400 to 518 nm, 25 cm-1 step, 50 cm-1 FWHM, 229 wavenumbers) and `ch4_karkoschka_IR.par` (2000 to 19295 cm-1, 518 nm to 5.0 micron, 5 cm-1 step, 10 cm-1 FWHM, 3460 wavenumbers). Together 3689 wavelengths. Units "(km-amagat)^-1".

- **Grid:** 10 g-ordinates; pressures 1e-8, 1e-7 ... 10 atm (10 values); temperatures 50 to 300 K in 25 K steps (11 values); each (wavenumber, temperature) block carries the band-model line and 10 pressure rows of 10 k values plus the transmission fitting error.
- **Band-model line:** wavenumber, wavelength (nm), κ(100 K), κ(198 K), κ(296 K), δ/α_D0. This is exactly the input the K&T Voigt-Goody model needs (open-access equations 1 to 5, READ; NEMESIS `tkark.f` implements it, READ). Sample (copied data, extracted CSV; wavenumber cm-1, nm, κ100, κ198, κ296, δ/α_D0):
```
19300.00, 517.99, 0.90E-03, 0.19E-02, 0.41E-02, 8.0
16100.00, 620.95, 0.561, 0.532, 0.505, 8.0
 2000.00, 4998.64, 0.0476, 0.0471, 0.0383, 17741
```
  (`ch4_karkoschka_tomasko2010_band_parameters.csv` has all 3689 rows.)
- **Strength across the range:** median κ(100 K) 2.2e-2 (500 to 700 nm), 0.54 (700 to 1000 nm), 7.5 (1.0 to 1.5 micron), 34 (1.5 to 2.5), 7.4 (2.5 to 5.0); largest 2.25e5 1/km-am near 3.27 micron.
- **Cross-check against the PDS tables:** κ(100 K) against the 1995 low-resolution PDS k at the same wavelengths, where both exceed 0.001: correlation of the logarithms 0.97 (500 to 700 nm), 0.995 (700 to 900 nm), 0.998 (900 to 1050 nm); median ratio within 0.04 dex (about 10 percent) in all three ranges. At the band centres: 619 nm 0.719 against 0.720, 727 nm 4.75 against 3.93, 790 nm 0.569 against 0.571, 890 nm 20.3 against 24.3. So the Oxford tables and Karkoschka 1998 agree (the 2010 abstract says its coefficients below 1000 nm are "generally close to Karkoschka 1998"), which supports reading the Oxford files as the real K&T 2010 parameters; this is a consistency test, not proof of provenance.

Whether the Oxford tables are the paper's own supplementary values: the paper's abstract speaks of about 3700 wavelengths and Collins et al. 2018 thank Karkoschka and Tomasko for "machine-readable tabulations ... in the Supplementary Materials" of K&T 2010 (READ, PMC6157968); the Oxford file has 3689 wavelengths with three-digit values. The supplement itself is on ScienceDirect and was **not reached**.

### 3.3 Other free items found

- **NEMESIS repository** (cloned with git): `raddata/karkoschka.dat` = the 1993 PDS methane column to 1000 nm (identical except the last row), extended to 1060 nm; `tkark.f` = K&T transmission; no K&T 2010 k-table inside (README: tables "under embargo ... available on a common server at some point"). NEMESIS README also says the visible ammonia bands (552, 645 nm) were originally modelled from Lutz and Owen (1980) and that this is "obsolete by newer ammonia line data" (those line data embargoed); Lutz and Owen was not looked up.
- **Oxford ammonia tables** on the same page: Irwin et al. 1999 (JQSRT 62, 193), 4000 to 11000 cm-1 (0.91 to 2.5 micron), self-broadened; and two tables at 25025 to 34975 cm-1 (286 to 400 nm). Neither reaches the 552, 645, 790 nm bands. Not downloaded.
- **Sithajan et al. 2025** (arXiv:2509.06496, READ abstract): 6195 methane features in Titan ESPRESSO spectra (R about 190000), detection limit about 0.02 km-am^-1, comparison with experimental data at similar temperature; and **Rianço-Silva et al. 2026** (arXiv:2607.10250): an empirical high-resolution optical methane line list (RRS-2026) from the same Titan data. Not needed for band strengths; useful as an independent check at 380 to 790 nm.
- **Sromovsky, Fry and Kim 2011** (Icarus 215, 292; arXiv:1503.02476, READ): independent Uranus methane retrieval, see section 4.

## 4. The note's claims, one by one

### 4.1 Bibliographic facts (Crossref, READ)

| Claim in the note | Result |
|---|---|
| K&T 2010: Icarus 205(2), 674 to 694, DOI 10.1016/j.icarus.2009.07.044, Feb 2010, Karkoschka and Tomasko, title as given | **Verified** (Crossref; Semantic Scholar record agrees, CorpusId 120940761 and id cd062ae0 verified; 148 citations there, 141 in Crossref, 140 Scopus on the Arizona page; the note's "148 to 176" upper end not found) |
| Karkoschka 1994: Icarus 111(1), 174 to 192, DOI 10.1006/icar.1994.1139, sole author, title as given | **Verified** (Crossref; 297 citations there, the note's "over 360" not found) |
| K&T 2009 Uranus: Icarus 202(1), 287 to 309, DOI 10.1016/j.icarus.2009.02.010, July 2009, Karkoschka and Tomasko | **Verified** (Crossref) |
| Karkoschka 1998: Icarus 133(1), 134 to 146, DOI 10.1006/icar.1998.5913, May 1998 | **Verified** (Crossref and the paper's own header) |
| K&T 2011 Neptune: "Icarus 212(2), 762 to 778, April 2011" | **Contradicted.** Crossref: DOI 10.1016/j.icarus.2010.08.013, Icarus **211**(1), **780 to 797**, January 2011. (The title "The haze and methane distributions on Neptune from HST-STIS spectroscopy" is right.) |
| "Accepted July 2009" for K&T 2010 | Not found (the DOI string carries 2009 but no acceptance date was read) |
| Corpus/S2 ids for 1994 and 1998 | Not checked |
| Irwin et al. 2022: JGR Planets 127(6), e2022JE007189, PMC9286428 | **Verified** (Crossref and PMC XML; CC BY 4.0) |
| Wong et al. 2023: Remote Sensing 15(3), 702 | **Verified** as a citation (Crossref: title "Deep Clouds on Jupiter", published 25 January 2023). Its content (the note says it details the K&T k-table parameterization) was **not read**: MDPI refuses |
| Collins et al. 2018 Science Advances, PMC6157968 | **Verified** (Crossref: Sci. Adv. 4(9), eaas9593, 7 September 2018; the DOI the note's text implies, "aas9930", is not it; CC BY-NC) |
| arXiv:1611.03257 (Irwin 2016/2017) | **Verified** (READ): Irwin, Wong, Simon, Orton, Toledo, "HST/WFC3 observations of Uranus' 2014 storm clouds and comparison with VLT/SINFONI and IRTF/SpeX", submitted to Icarus |
| arXiv:2607.10250 (Rianço-Silva 2026) and 2509.06496 (Sithajan 2025) | **Verified** as papers (READ abstracts). "6195 lines" is Sithajan's figure; "378.2 to 788.7 nm" not found in the abstracts |

### 4.2 Data hosts and files

| Claim | Result |
|---|---|
| PDS dataset **go_0016** = "Spectra of methane, jovian planets and Titan" | **Contradicted.** go_0016 to go_0022 are Galileo orbiter images. The volume is **gbat_0001** (catalogue text: "Spectrophotometry of the Jovian Planets and Titan"), data set id `ESO-J/S/N/U-SPECTROPHOTOMETER-4-V2.0`, whose catalogue summary says July 1995 data (the 1998 paper); it also contains the 1993 table |
| `atmos.nmsu.edu/planetary_datasets/jupiter_infrared.html` serves the file | **Verified**; but that page and the saturn, uranus, neptune and titan pages each serve **one planet** (four columns, 1993 data, cited as Karkoschka 1994); the combined 8-column tables are in the PDS volume |
| "Tab-delimited", wavelength starting 300.31 to 300.4 nm | Partly: fixed-width, space separated; vacuum 300.4 nm and air 300.31 nm both appear in the first row |
| Methane k in (km-am)^-1 with albedos of five bodies | **Verified** for the PDS volume |
| "Also cross-indexed in HITRAN-2024" | Not checked (not needed) |
| NEMESIS includes the K&T 2010 k-tables | **Partly verified.** The NEMESIS repository (git) does **not** contain them. The tables are on Irwin's Oxford page (`ch4_karkoschka_*.par.gz`), downloaded here. NEMESIS carries `karkoschka.dat` (the 1993 PDS column) and the `tkark.f` band-model routine |
| NASA PSG incorporates Karkoschka cross-sections | **Not found** (PSG host answers only "unsupported command / browser") |
| MAST holds the 19 August 2002 STIS cube of Uranus | Not checked (the 2009 abstract confirms the observation date; MAST host not reachable: `www.mast.stsci.edu` no connection, `archive.stsci.edu` front page only) |

### 4.3 K&T 2010 (AD1) content

Abstract read in full (University of Arizona repository, [experts.arizona.edu](https://experts.arizona.edu/en/publications/methane-absorption-coefficients-for-the-jovian-planets-from-labor/)).

| Claim | Result |
|---|---|
| "We use 11 data sets of methane transmission measurements within 0.4 to 5.5 micron" | **Verified** (word for word) |
| Eight laboratory data sets; two Huygens spectrometer sets (Tomasko et al. 2008b, PSS 56, 624) | **Verified** |
| "The eleventh data set incorporates planetary reflectance and transmission constraints from HST and earlier ground-based spectrophotometry" | **Contradicted.** The abstract: "The last data set is a set of **new Jupiter images by the Hubble Space Telescope to measure atmospheric transmission with Ganymede as the light source**" |
| Below 1000 nm results are close to Karkoschka 1998 and add temperature and pressure dependence; 882 to 902 nm inconsistency with laboratory data | **Verified** (abstract). The note omits the inconsistency; it matters for AD1 near the 890 nm band (the PDS and Oxford values differ by 20 percent there, section 3.2) |
| Parameterised with a correlated-k distribution and a Voigt-Goody band model | **Partly verified.** The paper provides band-model coefficients: for each of about 3700 wavelengths a temperature dependence (abstract); the Voigt-Goody formula with κ(ν,100), κ(ν,198), κ(ν,296) is printed in open-access Collins 2018 (eq. 1 to 5). The **correlated-k tables were computed by Irwin** (14 May 2010), not published in the paper itself |
| Huygens data cover 0.48 to 1.6 micron; Titan 70 to 94 K | Not found (abstract says "similar temperatures as many observations of Saturn, Uranus, Neptune and Titan") |
| Huygens data confirm Irwin et al. 2006 at laboratory temperatures; Irwin's extrapolation to Titan's low temperatures predicts coefficients up to 100 times too low | Not in the note; stated in the abstract and relevant: **do not extrapolate Irwin's 2006 temperature dependence to cold giants** |

### 4.4 Karkoschka 1994 and 1998 (READ in the PDS ASCII copies)

| Claim | Result |
|---|---|
| 1994: full-disk albedo spectra of the jovian planets and Titan, ESO July 1993, 300 to 1000 nm, 1 nm | **Verified** (plus S/N about 1000, 2 percent relative and 4 percent absolute accuracy, 40 to 60 Raman features per planet, three new weak bands on Uranus and Neptune) |
| 1998: ESO July 1995, 300 to 1050 nm, 0.4 nm between 520 and 995 nm and 1 nm elsewhere, 4 percent calibration, UBV, methane and ammonia bands, temperature | **Verified** (abstract, word for word). New 930 nm band "probably due to ammonia" |
| Both abstracts open with the same phrase | **Verified** (both start "Full-disk albedo spectra of the jovian planets and Titan were derived from observations at the European Southern Observatory in July ...") |
| Methane coefficients are in km-amagat^-1, derived indirectly | **Verified** (PDS label: "the indirect determination ... is described by Karkoschka and Tomasko (1992) ICARUS 97, 161-181", a paper not checked here) |
| 1994 strong Jupiter band "possibly water" | **Verified** (abstract); **withdrawn in 1998** (930 nm band probably ammonia; most of its lines are not water lines) |

### 4.5 Uranus numbers and filter coefficients

| Claim | Result |
|---|---|
| Uranus CH4 mole fraction at low and equatorial latitudes 3.2 to 4.0 percent ("0.04 plus or minus 0.01") | **Partly verified.** K&T 2009 abstract (READ): "The methane mixing ratio varied between **0.014 and 0.032 from high to low southern latitudes**". 3.2 percent is K&T's low-latitude figure; **4 percent is not K&T 2009's** but the deep value of Voyager occultation Model F (Lindal et al. 1987) and of Sromovsky, Fry and Kim 2011 ("a methane mixing ratio near 4% at 5 deg S", READ). "0.04 plus or minus 0.01" not found. Under the widest-range rule the usable bracket is 1.4 to 4 percent, with sources named |
| About 1.4 percent over the southern polar region | **Verified** (0.014, K&T 2009 abstract) |
| "Roughly double the about 2 percent from Voyager 2 radio occultation" | **Partly verified.** Lindal 1987's preferred Model D had 2.3 percent but Model F had 4 percent (READ in Sromovsky 2011). "Voyager gave 2 percent and K&T found double" is not what the sources say |
| Poleward of 30 to 45 degrees depletion | Not found in the abstract (the abstract gives a trend from low to high southern latitudes and an aerosol enhancement at 45 deg S); Sromovsky 2011 finds a model fit 30 deg S to 20 deg N with one profile and an increasing upper-tropospheric depletion closer to the pole |
| Aerosol opacity jumps by about 100 per bar near 1.2 to 2.0 bar | **Partly verified.** Abstract: "Below the 1.2 bar level, the aerosol opacity increased some 100 fold". "Per bar" and "2.0 bar" not found |
| Effective coefficients per HST filter (the note's table) | The numbers exist in **Table 1 of Sánchez-Lavega et al. 2016** (PMC5105178, READ), a Saturn paper that cites the K&T and Karkoschka coefficients (refs 33 to 35). **Verified:** F410M 410.8 nm k_eff 7.2e-5 1/km-am, P(tau=1) 1.210 bar; FQ937N 937.7 nm k_eff 0.0397, 24.866 bar; F225W 237.8 nm k_eff below 1e-6; FQ727N 727.7 nm k_eff 3.2302, 0.993 bar. **Contradicted:** the note puts 0.515 bar on F225W; the table has 0.112 bar for F225W and 0.515 bar for F336W. Full table: F225W 237.8 nm <1e-6, 0.112 bar; F336W 335.9 <1e-6, 0.515; F410M 410.8, 7.2e-5, 1.210; F502N 501.0, 0.0017, 2.760; F547M 545.2, 0.0221, 3.815; F689M 688.0, 0.2456, 5.772; FQ727N 727.7, 3.2302, 0.993; FQ750N 750.2, 0.0211, 13.183; FQ889N 889.4, 23.733, 0.138; FQ937N 937.7, 0.0397, 24.866. These are **effective coefficients for Saturn's atmosphere** (Rayleigh plus CH4 optical depth one), not new K&T data; use as a check only |
| Collins 2018: Jovian methane absorption inside MODTRAN, forcing 0.026 W/m2 global, +0.25 W/m2 localised (10 times), 2.75 times over 88 percent of Earth | **Verified** (READ, PMC6157968). Detail: MODTRAN version 5.3.0.0 with DISORT; the two line-by-line calculations are HITRAN-based and the K&T Voigt-Goody one; the equation set 1 to 5 is the open description of how to use κ at 100, 198 and 296 K |
| Irwin 2022 Fig. 1 re-plots K&T 2009 STIS Uranus spectra (0.3 to 1.0 micron) | **Partly verified.** Figure 1 shows the visible and near-infrared spectra of both planets; the text says the 2002 STIS Uranus data are those of K&T 2009 (READ). The re-plot not inspected |

## 5. What this means for AD1 (for the TODO thread and the lanes)

- **Methane (CH4):** the data exists in two free forms. For the atmosphere-colour work (reflected visible light of a giant, 0.4 to 1.0 micron), `κ(ν, T)` from the three-temperature band parameters with the Voigt-Goody formula (Collins 2018 eq. 1 to 5) gives absorption as a function of path, pressure, temperature and methane fraction, replacing the "neighbouring band 10 percent" stand-in; the PDS single-spectrum table (`gbat_0001`) is the simpler form and a test value. Source rules apply: where Oxford and PDS differ (890 nm band, up to 20 percent; weak bands below 700 nm, up to about 0.15 dex) both are verified, K&T 2010 is the newer, so it wins, and the PDS spread is the width for the draw.
- **Build-time download (Boss's decision of 2026-10-17):** both hosts answered here (`atmos.nmsu.edu`, `users.ox.ac.uk`); no licence text on either page, the PDS volume is public. The Oxford IR table is 19 MB compressed (68 MB text); the band-parameter CSV is 0.24 MB and holds everything the Voigt-Goody route needs. The Oxford page says "not yet fully verified" for these tables, and its host is a personal academic page (a stability risk): a build-time fetch should fall back to the labelled stand-in. For offline tests a few rows of the CSV are enough.
- **Ammonia and water:** AD1 remains open for the visible NH3 bands (552, 645, 790, 930, 1030 nm) and H2O; the only free numbers are the four equivalent widths above (read from a figure, labelled). The Irwin 1999 NH3 tables start at 910 nm (partly overlapping the 930 and 1030 nm bands, a lead, not read). The stand-in for ammonia stays.
- **Not usable under NO FIXED CONTENT:** fitting a given planet's measured albedo (the PDS albedo columns are Jupiter, Saturn, Uranus, Neptune and Titan data). The methane coefficients are physical data (laboratory plus Huygens), the albedo columns are not to be copied into the generator.
- **Citation fixes for the plan files:** the AD1 source is Karkoschka and Tomasko (2010), Icarus 205(2), 674 to 694; the Neptune paper is Icarus 211(1), 780 to 797 (2011), DOI 10.1016/j.icarus.2010.08.013; the Uranus paper is Icarus 202(1), 287 to 309 (2009); the dataset is PDS gbat_0001, not go_0016; the filter-table source is Sánchez-Lavega et al. 2016 (Nature Communications 7, 13262).
- **Free routes found for other gaps:** the Collins 2018 text (CC BY-NC) prints the full Voigt-Goody temperature interpolation; Sithajan 2025 and Rianço-Silva 2026 give an empirical Titan methane line list at 380 to 790 nm; `eutils.ncbi.nlm.nih.gov` and the Europe PMC REST host both return full text of open-access papers.

## 6. Not done

- K&T 2010 supplement itself (ScienceDirect refused) and Wong 2023 (MDPI refused) were not read; the Oxford file is the nearest copy found.
- Raw STIS cubes at MAST, PSG data, and Karkoschka and Tomasko 1992 (Icarus 97, 161) were not looked up.
- The Oxford tables were cross-checked only against the PDS tables, not against laboratory data or Huygens.
