# G4: ices, brines and licences (check of Boss's second note, claims G4.1 to G4.7)

Scope: body-related claims only. Labels follow `_brief.md` (READ, 2ND-HAND, COMPUTED, SNIPPET, RECALLED, NOT REACHED). Scripts and extracted tables: `parts/G-files/G4_*` (listed at the end). Reused from Part F: the transcribed Fray and Schmitt Table 2 (`F-files/fray_schmitt_table2_transcribed_and_checked.py`, NTRS 20210019502).

Key finding in one paragraph: the cold-end vapour-pressure papers exist, but the attribution is off (Blakley et al. is 2024, Planetary and Space Science 244, 105863, and covers only beta-O2 and NH3; N2, CO and CH4 are Grundy et al. 2024, Icarus 410, 115767). A 2026 preprint by Grundy's group applies a further +19 percent correction (ccold = 1.19) to all five species, so the "lower by 1.6" and "lower by 1.5" figures are already out of date. The "0 K to 800 K, 12 GPa" range belongs to Trusler 2011, not to Jager and Span 2012. The CO2 hydrate quadruple point matches one paper (Semenov et al. 2023). Two licence claims fail: NIST Standard Reference Data (WebBook, XCOM, ESTAR, PSTAR) is copyrighted by the Secretary of Commerce, not public domain, and the NTRS record for JPL Publication 19-5 says it is copyrighted with US-government-purpose use rights, not public domain.

## Verdicts

| id | claim as written (short) | verdict | the right figure (unit) | source | label |
|---|---|---|---|---|---|
| G4.1a | Grundy et al. 2024 and Blakley et al. 2025 replace Fray and Schmitt (F&S) extrapolations below 50 K for alpha-N2, alpha-CO, CH4 I, NH3, beta-O2 with QCM experiments | PARTLY | Two papers exist. Grundy et al. 2024, Icarus 410, 115767 (N2, CO, CH4). Blakley et al. **2024** (not 2025), Planet. Space Sci. 244, 105863 (NH3, beta-O2; arXiv 2403.09917 is this paper). Measured ranges: N2 23.3 to 31.1 K, CO 26.2 to 34.1 K, CH4 31.1 to 39.1 K (Grundy fig03 files); beta-O2 lower bound near 30 K; NH3 lower bound near 90 K, fitted over about 80 to 120 K (Blakley sec 2.2, 4.1). NH3 is therefore not "below 50 K". They are lab data in a narrow window plus a fitted extrapolation law (Lobo and Ferreira form), not a replacement below 50 K for all five. | Grundy 2024 abstract (NAU repository eprint 6246, Crossref DOI 10.1016/j.icarus.2023.115767); Blakley 2024 (arXiv 2403.09917 PDF, par.nsf.gov/servlets/purl/10610695) | READ |
| G4.1b | beta-O2 vapour pressures lower than legacy by an average factor 1.6 | CONFIRMED | "lower vapor pressures than previous work by an average factor of 1.6" (Blakley sec 4.2, Fig. 7; compared to F&S/Brown and Ziegler polynomial). NH3: factor 1.5 (sec 4.1). Both are the 2024 values; the 2026 correction (G4.1e) raises them by 4 to 39 percent depending on T (COMPUTED below) | Blakley 2024 sec 4.1, 4.2 | READ |
| G4.1c | sublimation enthalpy 9.2 percent higher | CONFIRMED | "Over the range of 20-50 K, our calculated values are higher by an average of 9.2% than the NIST accepted constant value of 9.26 kJ/mol". NH3: 2.74 percent higher than NIST 31.2 kJ/mol over 80 to 120 K. | Blakley sec 4.2, Fig. 8; sec 4.1 | READ |
| G4.1d | "E_a about 9.2 kJ/mol" (note's second wording) | CONTRADICTED | The two 9.2 values are unrelated. 9.26 kJ/mol is the old NIST constant for O2; the paper's own result is 9.2 percent above it, about 10.1 kJ/mol on average over 20 to 50 K (9.26 x 1.092 = 10.11, COMPUTED), and it varies with T. The paper gives L(T) as a fourth-degree function, not one constant. | Blakley sec 4.2; Table 1 | READ + COMPUTED |
| G4.1e | N2 "chamber background pressure correction factor Phi about 1.1875" | PARTLY | Phi in Grundy 2024 is a per-run polynomial scale on the gauge pressure (eq. 6: p_vap = Phi p_gauge sqrt(T_QCM/T_room) minus the dQ/dt term), fitted from the stepped deposition, not one constant, and its N2 value is not printed in anything reached. The only constant near 1.19 is a different correction: **ccold = 1.19 +/- 0.03**, the geometric back-condensation from cold parts of the QCM mount, in the 2026 Grundy et al. preprint (Appendix A), applied to all five species. 1.1875 itself is not found. | Grundy 2026 preprint (openknowledge.nau.edu/6360, "in press in Icarus 2026"; title Vapor Pressures of Low-Temperature C2Hx Hydrocarbon Ices), eq. 2 to 3, Appendix A; Blakley sec 3.1 | READ |
| G4.1f | CH4 reduced vapour pressure below 36 K | PARTLY | CH4 I is reduced over the whole measured and compared range, not only below 36 K. Original Grundy 2024 fit over F&S: 0.52 at 30 K, 0.62 at 35 K, 0.70 at 40 K, 0.83 at 50 K, 0.93 at 60 K, 1.01 at 70 K, 1.07 at 80 K (COMPUTED from the data package). With ccold = 1.19: 0.62, 0.74, 0.83, 0.99, 1.1 (30, 35, 40, 50, 60 K). Nothing special happens at 36 K for CH4. The 35.61 K transition is alpha to beta N2, a different molecule. | Grundy 2024 data package fig06 files (openknowledge.nau.edu/6246) | READ + COMPUTED |
| G4.1g | CO "largest relative error among super-volatiles" | CONFIRMED | Abstract: "carbon monoxide being the most discrepant of the three species, almost an order of magnitude lower than had been thought." Ratio to F&S (original fit): 0.067 at 25 K (15 times lower), 0.146 at 30 K, 0.25 at 35 K, 0.39 at 40 K, 0.69 at 50 K, 1.02 at 60 K. N2 ratio: 1.24 at 20 K, 0.97 at 22 K, 0.52 at 30 K, 0.42 at 34 K. | Grundy 2024 abstract; fig04 to fig06 files | READ + COMPUTED |
| G4.2a | Jager and Span 2012: validity "0 K to 800 K and up to 12 GPa" | CONTRADICTED | That range is the title and abstract of a different paper: Trusler, J. Phys. Chem. Ref. Data 40, 043105 (2011), "Equation of State for Solid Phase I of Carbon Dioxide Valid for Temperatures up to 800 K and Pressures up to 12 GPa" (abstract: "valid at temperatures from 0 to 800 K and at pressures from the solid-fluid phase boundary to 12 GPa"). Jager and Span's abstract states no numbers: "behaves physically reasonably within its range of validity". Their full text was not reachable; the range stays unknown. A 2025 JCED paper (acs.jced.5c00260) builds a new solid-I EoS valid "up to 400 K and 1.2 GPa" and lists Jager and Span and Trusler as comparison models. | Crossref abstract of 10.1063/1.3664915 (Trusler); OpenAlex/TU Dresden abstract of 10.1021/je2011677; OpenAlex abstract of 10.1021/acs.jced.5c00260 | READ (abstracts) |
| G4.2b | coupling with Span and Wagner 1996 | CONFIRMED | Abstract: sublimation and melting pressures are "predicted accurately using Gibbs energies calculated from the new dry ice equation and from the fundamental equation by Span and Wagner for the fluid phase"; they agree with the S&W reference equations "mostly within the uncertainty"; "very sensitive to variations of the Gibbs energy"; a different fluid EoS "has a significant impact". J. Chem. Eng. Data 57 (2), 590 to 597, published 6 January 2012. | Crossref/OpenAlex abstract | READ (abstract) |
| G4.3a | H2 critical temperature 33.19 K (check 33.145 K) | PARTLY | NIST WebBook phase-change table, H2: **Tc 33.18 K** (+/- 0.2 K), Pc 13.00 bar. The modern reference-EoS value is **33.145 K** (EoS as stored in CoolProp `Hydrogen.json`, Pc 1.2964 MPa). Equilibrium-H2 measured value 32.994 K (old temperature scale). 33.19 is not in any source read; use 33.145 K (EoS) or 33.18 K (WebBook), difference 0.1 percent. | webbook.nist.gov/cgi/cbook.cgi?ID=C1333740 ; CoolProp dev/fluids/Hydrogen.json ; J. Res. NBS 47, 75 | READ / 2ND-HAND |
| G4.3b | He critical temperature 5.20 K (5.195) | CONFIRMED | WebBook He: Tc 5.2 K, Pc 2.274 bar. EoS value 5.1953 K, Pc 0.2283 MPa (CoolProp `Helium.json`). | WebBook C7440597; CoolProp Helium.json | READ / 2ND-HAND |
| G4.3c | HD critical temperature about 35.9 K | CONFIRMED | **35.908 K**, Pc 14.645 atm, Vc 62.8 cm3/mol; e-H2 32.994 K, 12.770 atm; e-D2 38.262 K, 16.282 atm. WebBook lists HD only with a triple point 16.6 K; Hoge and Arnold 1951 Table 9: HD triple 16.60 K at 92.8 mmHg. | J. Res. NBS 47, 75 (RP2229) PDF; J. Res. NBS 47, 63 (RP2228) Table 9 | READ |
| G4.3d | H2, HD, He "non-condensable above 10 K" | PARTLY | Wrong as stated for H2 and HD. H2 solid exists below the triple point **13.95 to 13.96 K** (WebBook; EoS 13.957 K) at **0.0721 bar** (7.2 kPa; EoS 7358 Pa); H2 liquid exists up to Tc; H2 boils at 20.27 K at 1 atm (Hoge and Arnold Table 9, e-H2). HD: triple 16.60 K, boils 22.14 K. Equilibrium-D2 triple 18.69 K. So H2 and HD are condensible below about 14 and 17 K (solid) and 20 and 22 K (1 atm liquid). He has no solid at its own vapour pressure: triple (lambda) point of the liquid 2.1768 K, 5039 Pa (CoolProp), liquid below Tc 5.195 K at 1 atm; solid He needs about 25 bar (RECALLED, not checked). What the generator needs: **no H2 ice on any body whose surface is warmer than the H2 triple point (about 14 K at the triple pressure); H2 vapour pressure at 10 K is of order 100 Pa**, so even a 10 K surface cannot hold H2 ice against vacuum (ESTIMATE: Clausius-Clapeyron from the triple point, RECALLED sublimation enthalpy about 1 kJ/mol; not read). | WebBook C1333740; Hoge and Arnold 1951; CoolProp | READ |
| G4.4a | CO2 hydrate quadruple point V-Lw-I-H: T = 271.60 K, P = 1.044 MPa | CONFIRMED | Matches Semenov et al., Int. J. Mol. Sci. 24, 9321 (2023) (PMC10253777): direct measurement, "T = 271.60 K, P = 1.044 MPa", combined standard uncertainty +/- 0.023 K and +/- 0.021 MPa. Earlier indirect value Nema et al. 2017 (Fluid Phase Equilib. 441, 49): 271.6 +/- 0.2 K, 1.04 +/- 0.02 MPa (SNIPPET). Not found in Sloan and Koh, Mohammadi and Richon, Hashimoto, Chen or Nakano texts; Nagashima et al. are cited in the paper for V-I-H data. A value I recalled for Sloan and Koh (273.1 K, 1.256 MPa) was not found in any source: do not use it (RECALLED). | PMC10253777 full text | READ |
| G4.4b | uncertainty +/- 0.02 K, +/- 0.002 MPa | PARTLY | Paper: +/- 0.023 K and +/- 0.021 MPa (10 times larger in pressure). The tighter +/- 0.01 K and +/- 0.002 MPa appear only in a Mendeley data record by the same group (SNIPPET). Upper quadruple point Q2 (V-Lw-L(CO2)-H): **283.13 K, 4.494 MPa** (intersection of fits through the authors' earlier data and refs 57, 58 of that paper, Fig. in sec 1/2); roughly 283 K, 4.5 MPa in simulation papers (SNIPPET). | PMC10253777 figure caption | READ |
| G4.4c | the claim includes the salts MgCl2, RbCl, NaCl | UNVERIFIABLE | Tried: the Semenov paper (no salt data), web search for Q1 in NaCl, MgCl2 and RbCl systems. Found only inhibition studies (Holzammer and Braeuer 2019 Raman; Dholabhai, Sami and Cordeiro et al. for NaCl and MgCl2 hydrate equilibria above Q2 or without Q1). Nothing on RbCl. Missing: a primary paper that gives Q1 or the eutectic shift for CO2 hydrate in salt solutions. Not needed for v1 (hydrate in brines is later work, Part F item F4). | n/a | NOT REACHED |
| G4.5a | Mg(ClO4)2 eutectic about 205 K (literature 206 K at 44 wt percent) | CONFIRMED | **206 +/- 1 K at 44.0 wt% Mg(ClO4)2** (Chevrier, Hanley, GRL 36, L10202, theoretical eutectic). 205 K with eutectic water activity aw = 0.55 (Gough et al. 2011, as quoted by Martinez et al. arXiv 1811.00862); 206 K, aw = 0.56 (Rivera-Valentin et al. arXiv 2012.00100 ref 6, i.e. Chevrier 2009); 212 K, aw 0.53 in another group's table (Pal and Kereszturi arXiv 2112.08287). Spread 205 to 212 K between authors. | Crossref abstract 10.1029/2009gl037497; arXiv preprints | READ (abstract), 2ND-HAND |
| G4.5b | Ca(ClO4)2 eutectic about 197 K | PARTLY | **198 K**, aw = 0.52 (Pestova et al. 2005; Marion et al. 2010, as quoted by Martinez et al. 1811.00862 and Rivera-Valentin et al. (arXiv 2012.00100)). 199 K, aw 0.51 in Pal and Kereszturi (arXiv 2112.08287, Table 1). 197 K is not found. Eutectic composition not reached (52.4 wt% appears in Martinez et al. as the brine composition at about 205 K, not as a stated eutectic). | arXiv 1811.00862, 2012.00100, 2112.08287 | 2ND-HAND |
| G4.5c | NaClO4 eutectic about 236 K | CONFIRMED | **236 +/- 1 K at 52 wt% NaClO4** (Chevrier et al. 2009 abstract). | Crossref abstract 10.1029/2009gl037497 | READ |
| G4.5d | ln DRH = (dH_sol/R)(1/T - 1/T0) + ln DRH(T0) | UNVERIFIABLE | Form is the integrated van't Hoff relation; no text I read writes it for DRH. What is read: DRH of Mg(ClO4)2 hydrate rises as T falls (about 42 percent near 273 K to about 64 percent near 223 K, Gough 2011, SNIPPET); Gu et al. 2017 (RSC Adv.) find 42.8 to 40.5 percent from 278 to 303 K (SNIPPET). Fitting the Gough pair to the form gives an effective dH/R = 513 K (about 4.3 kJ/mol, COMPUTED from SNIPPET inputs, not usable). Nuding et al. 2014 (Icarus 243, 420) give Ca(ClO4)2 DRH from 223 to 273 K: 5 to 55 percent, highest at the lowest T (SNIPPET); Martinez et al. use their fits as DRH(T) and ERH(T). Needs Gough 2011 (EPSL 312, 371) and Nuding 2014 full texts. | see left | SNIPPET / NOT REACHED |
| G4.6 | Venus H2SO4-H2O clouds 70 to 95 wt percent (typical 75 to 85); activities by Tabazadeh 1997, Zeleznik 1991, Ayers 1980; cloud base evaporation 48 to 52 km | PARTLY | Wt range: published retrievals 73 to 87 wt% in the upper clouds and 85 to 96 wt% in the lower clouds; best-fit droplets about 75 wt% (SNIPPET from a 2025 paper and several Venus papers; the Titov 2018 review not read). Cloud-opacity drop at **48 km** where the acid evaporates (Grassi review arXiv 1804.08340 sec 4.1.1, 2ND-HAND); "48 to 52" is not supported. Source mapping: The activity formulation is valid for **185 to 260 K** and RH above 1 percent: stratospheric aerosol, **not Venus cloud temperatures** (cloud base is warm, above 300 K). Zeleznik 1991 (JPCRD 20, 157) covers aqueous H2SO4 "to 350 K" and all compositions: usable. The fit is the H2SO4 vapour pressure over 98.01 mass-% acid, 338 to 445 K: relevant to the base of the deck, not an activity-coefficient model. | Crossref abstracts 10.1029/97gl01879, 10.1063/1.555899, 10.1029/gl007i006p00433; arXiv 1804.08340 | READ (abstracts) / SNIPPET |
| G4.7a | Mainz Spectral Atlas = CC BY 4.0, cite Keller-Rudek 2013 | CONFIRMED (frozen copy) | The Zenodo record 6951 (frozen snapshot of 24 June 2013, 126 MB zip) has licence `cc-by-4.0`. The live site pages read state no licence; they ask to cite Keller-Rudek, Moortgat, Sander, Sorensen, ESSD 5, 365 (2013), doi 10.5194/essd-5-365-2013. The ESSD paper itself is CC BY **3.0**. So: CC BY 4.0 is stated only for the 2013 snapshot. | zenodo.org/api/records/6951; uv-vis-spectral-atlas-mainz.org/uvvis/ ; essd.copernicus.org/articles/5/365/2013 | READ |
| G4.7b | Leiden photodissociation database = "academic citation request / open access", cite Heays 2017 | PARTLY | Citation request confirmed: "When using material from this site, please cite one of the overview papers. In cases of specific molecules, please cite the original source of the cross sections." Overview list: Hrodmarsson and van Dishoeck 2023 (A&A), Heays, Bosman, van Dishoeck 2017 (A&A). No licence or "open access" statement on the pages read (index, downloads, acknowledgements, references). Data are compiled from many third-party sources, so redistribution status is unknown. | home.strw.leidenuniv.nl/~ewine/photo/ (index, downloads, references) | READ |
| G4.7c | JPL Publication 19-5 = "public domain (US Government work)", free redistribution | CONTRADICTED | NTRS record 20210006316 (Burkholder et al., Evaluation No. 19): distribution PUBLIC, but the copyright block says `belongsToUsGov: false`, `belongsToPublisher: true` and the approver comment "Items in TRS are protected by copyright, but are furnished with U.S. government purpose use rights." Full text is not disseminated by NTRS (metadata only). jpldataeval.jpl.nasa.gov is refused by this environment's proxy (502) and jpl.nasa.gov/copyrights returned 403, so the site's own copyright FAQ was not read. Not public domain on the evidence; redistribution of derived tables unknown. | api NTRS 20210006316 | READ (metadata) |
| G4.7d | "Bremen molecular spectra files" = proprietary | UNVERIFIABLE | The IUP Bremen database pages (reference spectra, DLR database, Databases index) state no licence or terms in the text read; the CH4 data page needs a login ("Please log on"); O3 2011 page asks to check for updates and contact V. Gorshelev. "Proprietary" is an inference, not a stated term. | iup.uni-bremen.de/gruppen/molspec/databases/ | READ (absence) |
| G4.7e | "NMSU Karkoschka band files" = academic attribution | PARTLY | Host is right: atmos.nmsu.edu/planetary_datasets/ serves Karkoschka (1994, Icarus 111, 174) methane absorption coefficient and Uranus albedo text files (300 to 1000 nm), source line only, no licence. Karkoschka and Tomasko (2010, Icarus 205, 674) coefficients are available as derived k-tables from P. Irwin's Oxford page (users.ox.ac.uk/~atmp0035/methane.html), flagged "not yet fully verified", no licence. The Arizona page of E. Karkoschka (pirlwww.lpl.arizona.edu/~erich/) is a 1998 placeholder with no data. "Academic attribution" is not a stated term. | as left | READ (absence) |

## Figures for the generator

### F1 cold-end vapour pressure (answers G4.1; needed by `p_sat(species, T)`)

What changes in item F1. Part F's transcribed F&S coefficients are now independently confirmed: the F&S polynomials tabulated in the Grundy 2024 data package (fig04 to fig06) agree with the Part F transcription of N2-1, CO-1 and CH4-1 to 0.5 percent (20 to 60 K; `G4_qcm_check.out.txt`). So the transcription can be relabelled READ (F&S curves) plus COMPUTED (match).

Recommended structure of `p_sat` for the super-volatiles (PARAMETRISATION of which source to use, not a physical constant):

| species | below the QCM window | QCM window (lab data) | above the window |
|---|---|---|---|
| alpha-N2 | refit needed; do not use F&S below 22 K without a flag (original fit crosses F&S at about 21.5 K) | 23.3 to 31.1 K | F&S set 1 or Brown and Ziegler up to the alpha-beta transition 35.61 K (Grundy fig07 gives L(T), 7.18 kJ/mol near 20 K) |
| alpha-CO | same | 26.2 to 34.1 K | transition from corrected-QCM to F&S near 55 to 60 K (ratio crosses 1.0 at 55 K with ccold, 60 K without) |
| CH4 I | same | 31.1 to 39.1 K | crossover to F&S near 50 to 60 K |
| beta-O2 | Blakley fit | about 30 K upward | Blakley fit joined with Aoyama and Kanda / Hoge data (37.6 to 53.1 K) |
| NH3 | Blakley fit | about 90 K upward (fit uses 80 to 120 K) | F&S set 1 above 162 K (lab data 162 to 195 K) |

Table of corrected values (READ + COMPUTED): `G-files/G4_qcm_corrected_vapor_pressure.csv` (114 rows: T, F&S, Grundy 2024 original Lobo and Ferreira fit, same times 1.19, ratio). Spot values (bar):

| species | T (K) | F&S 2009 | Grundy 2024 fit | times ccold 1.19 | ratio corrected / F&S |
|---|---|---|---|---|---|
| alpha-N2 | 25 | 8.42e-11 | 6.14e-11 | 7.31e-11 | 0.87 |
| alpha-N2 | 30 | 3.90e-8 | 2.01e-8 | 2.39e-8 | 0.61 |
| alpha-CO | 30 | 3.61e-9 | 5.27e-10 | 6.27e-10 | 0.17 |
| alpha-CO | 40 | 1.01e-5 | 3.90e-6 | 4.64e-6 | 0.46 |
| CH4 I | 36 | 3.01e-10 | 1.91e-10 | 2.27e-10 | 0.76 |
| CH4 I | 50 | 2.98e-6 | 2.47e-6 | 2.94e-6 | 0.99 |

Direct QCM measurements (Grundy 2024 fig03, not corrected for ccold; T in K, P in bar): N2 23.29 3.25e-12, 25.24 7.82e-11, 27.17 9.36e-10, 29.12 8.89e-9, 31.08 6.30e-8. CO 26.20 4.57e-12, 28.14 4.74e-11, 30.09 5.99e-10, 32.06 4.81e-9, 34.06 3.45e-8. CH4 31.07 1.21e-12, 33.05 9.12e-12, 35.05 6.51e-11, 37.07 4.06e-10, 39.09 2.41e-9.

Functional form (used by both papers): ln p[bar] = A - B/T + C ln T + D2 T + D3 T^2 + D4 T^3. Printed coefficients:

| species | A | B (K) | C | D2 (1/K) | D3 (1/K^2) | D4 (1/K^3) | source |
|---|---|---|---|---|---|---|---|
| beta-O2, 2024 | 15.29 +/- 0.29 | 1166.2 +/- 8.6 | -0.75587 | 0.14188 | -1.8665e-3 | 7.582e-6 | Blakley Table 1 (READ) |
| NH3, 2024 | -5.55 +/- 0.45 | 3605 +/- 31 | 4.82792 | -0.024895 | 2.1669e-5 | -2.3575e-8 | Blakley Table 1 (READ) |
| alpha-N2, with ccold | -3.35 +/- 0.19 | 817 +/- 5 | not printed | not printed | not printed | not printed | Grundy 2026 Table 3 (READ, A and B only) |
| alpha-CO, with ccold | -6.150 +/- 0.036 | 971.6 +/- 2.1 | not printed | | | | same |
| CH4 I, with ccold | -2.666 +/- 0.027 | 1127.8 +/- 2.0 | not printed | | | | same |
| beta-O2, with ccold | 14.85 +/- 0.23 | 1147 +/- 9 | take Blakley's | | | | same |
| NH3, with ccold | -5.73 +/- 0.57 | 3570 +/- 27 | take Blakley's | | | | same |

C, D2, D3, D4 for N2, CO, CH4 are printed in the Grundy 2024 paper (Table; not reachable), not in the 2026 preprint. Reconstruction for use until the paper is read (COMPUTED, `G4_fit_updated_lf.py`, fitting the tabulated 2024 curve times 1.19 with the 2026 A and B fixed; max error 0.3 to 0.6 percent in ln p inside the fitted ranges; **valid only inside** 20 to 35 K, 25 to 61 K, 30 to 90 K; these are curve-fit constants, not the authors' physical C..D4):

| species | C | D2 | D3 | D4 | range (K) |
|---|---|---|---|---|---|
| alpha-N2 | 5.4606 | -0.27791 | 0.0042387 | -3.7668e-5 | 20 to 35 |
| alpha-CO | 5.9227 | -0.0961 | 7.746e-5 | 5.0999e-7 | 25 to 61 |
| CH4 I | 3.5975 | -0.025762 | -1.7561e-4 | 1.0998e-6 | 30 to 90 |

Effect of the 2026 correction on O2 and NH3 (COMPUTED, `G4_blakley_update_ratio.py`, assumes the 2024 C..D4 are unchanged): beta-O2 new/old = 1.39 at 25 K, 1.22 at 30 K, 1.12 at 35 K, 1.04 at 40 K, 0.99 at 45 K. NH3 new/old = 1.29 at 80 K, 1.23 at 90 K, 1.19 at 100 K, 1.12 at 120 K, 1.07 at 140 K. So the "factor 1.6" for O2 is about 1.3 to 1.5 at 30 K after the correction, and the corrected NH3 factor at 90 K is 1.5/1.23 = 1.2 (COMPUTED). The 2026 paper is a preprint ("in press in Icarus 2026"); the Elsevier version was not read. Honest status: ccold is a geometry factor measured from photographs for one apparatus (1.19 +/- 0.03), so the corrected figures are the authors' latest best estimate, not yet an independent confirmation. Use the corrected set as default and keep the F&S set as a switchable alternative; the generator's tests should accept either within a factor 2.

Enthalpy of sublimation: alpha-N2 7.18 to 7.21 kJ/mol over 20 to 35.6 K; beta-N2 about 7.09 kJ/mol at 35.6 K (Grundy fig07). beta-O2 about 10.1 kJ/mol average, NH3 about 32.0 kJ/mol average over 80 to 120 K (31.2 x 1.0274, COMPUTED).

Cross-check (COMPUTED, `G4_qcm_check.py`): Lobo and Ferreira NH3 from Blakley Table 1 divided by F&S set 1 (Part F coefficients) is 0.57 at 80 K, 0.63 at 90 K, 0.68 at 100 K, 0.88 at 140 K, consistent with "lower by an average factor 1.5" over 80 to 120 K. The beta-O2 curve from Table 1 gives 7.07e-11 bar at 30 K, 1.40e-6 bar at 40 K (no beta-O2 F&S coefficients in Part F to compare).

### CO2 sublimation (answers G4.2c; test target for p_sat and for CO2 frost)

Span and Wagner 1996, J. Phys. Chem. Ref. Data 25, 1509, eq. 3.12 (READ from the NIST SRD reprint, p. 1521): ln(p_sub/p_t) = (T_t/T) (a1 th + a2 th^1.9 + a3 th^2.9), th = 1 - T/T_t; T_t = 216.592 K, p_t = 0.51795 MPa, a1 = -14.740846, a2 = 2.4327015, a3 = -5.3061778. Stated uncertainty (OCR of the reprint; the first bound is garbled and read as 185 K): +/- 2.5 Pa for 185 K <= T <= T_t, +/- 100 Pa for 170 to 185 K, +/- 50 Pa below 170 K. The equation was fitted to data above 154 K; low-temperature data are matched to 0.01 Pa for T > 85 K and 0.0001 Pa below. Normal sublimation temperature 194.6857 +/- 0.0030 K (Bedford 1984); eq. 3.12 returns 194.6855 K (also COMPUTED here: 194.68546 K).

COMPUTED values (`G4_co2_sublimation.py`):

| T (K) | p_sub (Pa) | note |
|---|---|---|
| 80 | 2.169e-6 | |
| 100 | 1.3408e-2 | |
| 120 | 3.736 | |
| **150** | **842.65** | test target |
| 170 | 9960.4 | |
| 185 | 43935 | |
| **195** | **103971** | test target (above 1 atm: the normal sublimation point is 194.6857 K, 101327 Pa) |
| 210 | 327089 | |
| 216.592 | 517950 | triple point |

Cross-check of Part F item F1: F&S CO2 set 1 (Part F transcription) over S&W is 0.999 at 150 K, 0.996 at 170 K, 0.995 at 195 K, 0.976 at 210 K, 1.07 at 125 K and 1.41 at 100 K; set 2 over S&W is 1.38 at 150 K, 1.09 at 170 K, 1.006 at 195 K, 1.0045 at 210 K. So **set 1 is the one for 150 to 200 K, set 2 is only valid above about 190 K, and below about 130 K both F&S sets are 7 to 40 percent (set 1) or more (set 2) above Span and Wagner**. Recommendation: use S&W eq. 3.12 for CO2 (three constants, one equation, accurate to 0.1 percent at 150 to 200 K) and keep F&S for the others. This answers the open question in Part F "selection rule between the two sets for CO2".

### Light gases, critical and triple points (answers G4.3; Part E materials-table row for H2 already matches)

| gas | Tc (K) | Pc | T triple (K) | P triple | normal boiling (K) | source | label |
|---|---|---|---|---|---|---|---|
| H2 (normal, WebBook) | 33.18 (+/- 0.2) | 13.00 bar | 13.95 / 13.96 | 0.0721 bar | 20.27 (e-H2, Hoge and Arnold) | NIST WebBook C1333740 | READ |
| H2 (reference EoS) | 33.145 | 1.2964 MPa | 13.957 | 7357.8 Pa | n/a | CoolProp Hydrogen.json (Leachman) | 2ND-HAND |
| HD | 35.908 | 14.645 atm | 16.60 | 92.8 mmHg | 22.143 | Hoge and Arnold 1951 | READ |
| D2 (equilibrium) | 38.262 (EoS 38.34) | 16.282 atm | 18.69 (WebBook 18.65 to 18.73) | 128.5 mmHg | 23.63 | same; WebBook C7782390 | READ |
| He-4 | 5.2 (WebBook) / 5.1953 (EoS) | 2.274 bar / 0.2283 MPa | lambda point 2.1768 | 5039 Pa | not read | WebBook C7440597; CoolProp | READ / 2ND-HAND |

Part E's row `gas_H2 ... 13.95, 20.39, 13.95, 0.0721, 33.18, 13.0` matches the WebBook values for Tt, Pt, Tc and Pc (so no change); its Tb 20.39 K differs from e-H2 20.27 K by 0.6 percent (normal-H2 boils 20.39 K per other tables, RECALLED, not checked here).

### CO2 hydrate points (answers G4.4)

| point | T (K) | P (MPa) | uncertainty | source | label |
|---|---|---|---|---|---|
| Q1 (V-Lw-I-H) | 271.60 | 1.044 | +/- 0.023 K, +/- 0.021 MPa | Semenov et al. 2023 | READ |
| Q2 (V-Lw-L-H) | 283.13 | 4.494 | from fit intersection, not quoted | same, figure caption | READ |

### Brine freezing-point lookup (answers G4.5; for hydrosphere and Mars-type brines)

Pure-salt eutectics. The generator should draw within a spread because authors differ by 1 to 7 K (PARAMETRISATION: normal about the stated value with the stated spread; Solar System / Mars use is calibration only).

| salt | eutectic T (K) | eutectic composition | eutectic aw | source | label |
|---|---|---|---|---|---|
| NaClO4 | 236 +/- 1 | 52 wt% | not read | Chevrier et al. 2009 | READ |
| Mg(ClO4)2 | 206 +/- 1 (author range 205 to 212) | 44.0 wt% | 0.55 to 0.56 | Chevrier 2009; Gough 2011 (205 K) via arXiv preprints; Pal and Kereszturi (212 K, aw 0.53) | READ / 2ND-HAND |
| Ca(ClO4)2 | 198 (author range 198 to 199) | not reached (52.4 wt% quoted at about 205 K in Martinez et al.) | 0.51 to 0.52 | Pestova 2005, Marion 2010 via arXiv 1811.00862 and 2012.00100; arXiv 2112.08287 (199 K) | 2ND-HAND |
| H2SO4 (Venus) | see G4.6; Zeleznik 1991 freezing-point curve, T to 350 K and 0 to 100 wt% | not read | n/a | Crossref abstract | READ (abstract) |
| NaCl, MgCl2, other chlorides | not in this task | | | | NOT REACHED |

Toner and Catling 2016 (Geochim. Cosmochim. Acta 181, 164; PDF read, faculty.washington.edu): water activity of NaClO4, Ca(ClO4)2 and Mg(ClO4)2 solutions from heat capacities to 178 K; aw in NaClO4 rises as T falls (0.6 to 0.7 for up to 15 molal below 200 K), in Ca and Mg perchlorate it falls. Brines reach about 195 K for Mg(ClO4)2 (SNIPPET). Useful for the aw(T) of a brine, not read in detail.

### Venus acid clouds (answers G4.6)

Keep: H2SO4-H2O solution, 75 to 85 wt% in the upper cloud as a draw range (PARAMETRISATION, Solar System calibration only); cloud-deck lower edge where the acid evaporates, near 48 km on Venus (a result of temperature and acid vapour pressure, not a constant: derive from T(z) and the Ayers/Zeleznik vapour pressure). Do not use Tabazadeh 1997 for cloud decks; use Zeleznik 1991 (valid to 350 K).

## Licence table (G4.7 and the data this design will load)

"Can we ship derived tables in a GPL/open repo?" is answered plainly: YES, NO, or UNKNOWN. I am not a lawyer; where the source states nothing the answer is UNKNOWN and the action is to ask the maintainer.

| source | licence as stated | can we ship derived tables in a GPL/open repo? | evidence |
|---|---|---|---|
| MPI-Mainz UV/VIS Spectral Atlas | Zenodo snapshot 6951 (2013-06-24): `cc-by-4.0`. Live site: no licence stated on home, contact or cross-section index pages; asks to cite Keller-Rudek et al. 2013 (ESSD paper is CC BY 3.0) | YES for the 2013 snapshot with attribution (CC BY 4.0 permits copying, derived works and redistribution if credit and a change notice are kept). Live-site data added after 2013: UNKNOWN | zenodo.org/api/records/6951 ("license": {"id": "cc-by-4.0"}); uv-vis-spectral-atlas-mainz.org/uvvis/ ("please cite our publication") |
| Leiden Photodissociation / photoionisation database (van Dishoeck, Heays) | no licence stated; "When using material from this site, please cite one of the overview papers" | UNKNOWN (citation request only; many entries repackage third-party cross sections) | home.strw.leidenuniv.nl/~ewine/photo/index.html, downloads.html |
| JPL Publication 19-5 | NTRS 20210006316: "Items in TRS are protected by copyright, but are furnished with U.S. government purpose use rights"; `belongsToUsGov` false, `belongsToPublisher` true. Site copyright FAQ not reachable | UNKNOWN to NO as public domain; do not assume. Individual facts (a rate constant, a cross section) with citation are conventional use; bulk copy of tables: ask JPL | ntrs.nasa.gov/api/citations/20210006316 |
| Bremen IUP molecular spectra | no licence or terms found; one page needs login | UNKNOWN | iup.uni-bremen.de/gruppen/molspec/databases/ (pages read) |
| Karkoschka methane files (NMSU host, Oxford k-tables) | none stated, source line only | UNKNOWN | atmos.nmsu.edu/planetary_datasets/uranus_infrared.html; users.ox.ac.uk/~atmp0035/methane.html |
| NIST WebBook (SRD 69) | "Copyright for NIST Standard Reference Data is governed by the Standard Reference Data Act"; "Standard Reference Data (SRD) are copyrighted by the U.S. Secretary of Commerce on behalf of the United States of America. All rights reserved. None of our SRD may be reproduced, stored in a retrieval system or transmitted ... without prior permission." (nist.gov/srd/public-law). The nist.gov/open/license page adds that NIST data that is not SRD, created by NIST employees, is not under US copyright and carries a worldwide royalty-free right to copy, derive and distribute with acknowledgement | **NOT as public domain.** Using a few individual values with citation is normal practice; bulk reproduction of tables is UNKNOWN and needs written permission. Prefer computing from published primary sources (EoS, Span and Wagner etc.) and citing NIST for spot checks | nist.gov/srd/public-law ; nist.gov/open/license ; webbook.nist.gov/chemistry/ ("Note: This site is covered by copyright") |
| NIST XCOM (SRD 8), ESTAR / PSTAR (SRD 124) | pages read carry the Berger et al. introductions and the NIST header; no separate licence text found. NIST lists these among SRD | UNKNOWN; treat as SRD (permission for bulk copy). The underlying cross-section tables are also in published reports (Berger and Hubbell) that can be cited | physics.nist.gov/PhysRefData/Xcom/Text/intro.html, Star/Text/intro.html |
| PubChem | NLM/NCBI: "NCBI itself places no restrictions on the use or distribution of the data" (molecular databases policy, quoted in the NCBI policies page text I read) but records "incorporate material contributed or licensed by individuals, companies, or organizations that may be protected" and users "are expected to adhere to the terms ... asserted by the copyright holder"; US-government-created NLM web content is public domain, with NLM acknowledgement requested. The PubChem docs pages themselves were behind a captcha | UNKNOWN per record; basic facts (formula, mass, melting point as a fact) are low risk; check the "data source" licence of each property; do not bulk-copy contributor text | ncbi.nlm.nih.gov/home/about/policies/ (read); PubChem docs (captcha) |
| USGS | "USGS-authored information is in the public domain; credit requested; third-party images and some datasets are copyrighted or carry other terms" | YES for USGS-authored data with credit; NO assumption for items marked as third-party. Official page returned 403 here | SNIPPET from usgs.gov FAQ "are usgs reports copyrighted"; not READ |
| Hypatia Catalog | no licence or terms found; acknowledgement text asks to cite Hinkel et al. 2014 (AJ 148, 54) and the individual catalogs used (the catalog merges more than 80 literature sources) | UNKNOWN | hypatiacatalog.com (home, acknowledgements section) ; SNIPPET for the cite request |
| Zenodo records | each record has its own `license` field (e.g. `cc-by-4.0` for record 6951); Zenodo does not impose one | per record: read the field via `zenodo.org/api/records/<id>` and store it in the loader | READ for 6951 |
| IAPWS releases (e.g. R14-08, sublimation pressure of ice) | "Publication in whole or in part is allowed in all countries provided that attribution is given to the International Association for the Properties of Water and Steam" | YES with attribution (formulas and coefficients) | R14-08(2011) PDF, first page (iapws.org/documents/release/MeltSub.download) |
| NASA Glenn thermodynamic polynomials (NASA/TP-2002-211556, McBride, Zehe, Gordon) | "Unclassified - Unlimited"; NTRS distribution PUBLIC; written by NASA Glenn employees (US government work) | YES (US-government work); cite the report. The NTRS metadata flag `belongsToUsGov` is false for this record, so that flag is not reliable | NTRS 20020085330 (text), report block 12a |
| Burcat database (garfield.chem.elte.hu) | host does not resolve from this environment | UNKNOWN (NOT REACHED; RECALLED: free with acknowledgement) | n/a |
| Fray and Schmitt 2009 coefficients | published in an Elsevier paper (Planet. Space Sci. 57, 2053); also reproduced in a NASA/Icarus paper supplement (NTRS 20210019502; `belongsToPublisher: true`, so the NTRS copy is not a government work) | numerical coefficients of a fit are facts, low risk with citation; do not copy the paper's tables or text. A better option is to refit our own polynomial to computed values (Part F did the check) | NTRS record metadata |
| Grundy 2024, Blakley 2024, Grundy 2026 data and Lobo-Ferreira coefficients | data archived at openknowledge.nau.edu (eprints 6246, 6254, 6360); no licence statement on the pages read | coefficients, as numbers with citation: low risk. Archive redistribution: UNKNOWN | openknowledge.nau.edu/id/eprint/6246 |

Change needed to Part F item F5 ("NIST, PubChem, USGS public domain"): NIST WebBook / XCOM / ESTAR / PSTAR are copyrighted SRD, PubChem is per record, USGS public domain is only a SNIPPET. The materials-table loader should store a `licence` and a `licence_evidence_url` per row and treat `UNKNOWN` as "do not ship the table, ship a loader that fetches at build time or recomputes from a cited formula".

## Still unreachable (Boss must research)

1. Grundy et al. 2024 (Icarus 410, 115767) main text and Table of Lobo-Ferreira C, D2, D3, D4 for N2, CO, CH4, plus the Phi values (journal page refused; the accepted manuscript is not in the NAU repository, only data and README). Gives: the exact coefficients to replace my reconstruction, and the Phi per run.
2. Jager and Span 2012 full text: validity range of the solid CO2 EoS. Gives: whether it can be used above 195 K or in the deep ice layer of a CO2-rich body. Trusler 2011 (0 to 800 K, up to 12 GPa) is the alternative with a stated range.
3. Gough et al. 2011 (EPSL 312, 371) and Nuding et al. 2014 (Icarus 243, 420): DRH(T) fits and the Ca(ClO4)2 eutectic composition. Gives the DRH equation and its constants. Pestova et al. 2005 and Marion et al. 2010 for the Ca(ClO4)2 and NaClO4 phase diagrams (Marion 2010 is Icarus 207, 675).
4. JPL data-evaluation site copyright FAQ (jpldataeval.jpl.nasa.gov/ and jpl.nasa.gov/copyrights are blocked from here), and NIST written permission for bulk SRD reuse (srd@nist.gov).
5. A primary source for Q1 of CO2 hydrate in NaCl, MgCl2 (and RbCl) brines (only inhibition work found).
6. Titov et al. 2018 (Space Sci. Rev. 214, 126) for the Venus H2SO4 wt% versus altitude.
7. Hypatia and Leiden: ask the maintainers whether redistribution of derived tables is allowed.

## Counts

26 claim rows: CONFIRMED 10 (G4.1b, G4.1c, G4.1g, G4.2b, G4.3b, G4.3c, G4.4a, G4.5a, G4.5c, G4.7a), PARTLY 10 (G4.1a, G4.1e, G4.1f, G4.3a, G4.3d, G4.4b, G4.5b, G4.6, G4.7b, G4.7e), CONTRADICTED 3 (G4.1d, G4.2a, G4.7c), UNVERIFIABLE 3 (G4.4c, G4.5d, G4.7d).

Verification labels used in this file: READ about 38, 2ND-HAND about 8, COMPUTED 9, SNIPPET 9, RECALLED 4 (listed above), NOT REACHED 4.

## Files (all under `/mnt/project-files/design-drafts/parts/G-files/`)

- `G4_qcm_check.py`, `.out.txt`: Part F transcription of F&S N2/CO/CH4 against the data package; Blakley Table 1 against F&S for NH3.
- `G4_qcm_table_csv.py`, `G4_qcm_corrected_vapor_pressure.csv`: F&S, Grundy 2024 fit and the 1.19-corrected values per kelvin.
- `G4_fit_updated_lf.py`, `.out.txt`: reconstruction of C..D4 (flagged COMPUTED).
- `G4_blakley_update_ratio.py`, `.out.txt`: effect of the 2026 update on O2 and NH3.
- `G4_co2_sublimation.py`, `.out.txt`: Span and Wagner eq. 3.12 values and the F&S CO2 cross-check.
- `G4_grundy2024_datapackage_README.txt`: copy of the NAU data README (file list).
- Inputs downloaded outside the repo folder (not saved): Grundy 2024 data package `data_package.zip` (openknowledge.nau.edu/6246/2), Blakley PDF, Grundy 2026 preprint PDF, NIST reprints RP2228 and RP2229, Span and Wagner reprint.
