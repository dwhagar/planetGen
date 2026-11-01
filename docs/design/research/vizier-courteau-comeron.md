# VizieR search for [R0908] and Comerón thick-disc tables (GR5 and GR3 rank 2)

Status: research result. Written 11:40 AM EDT on 2026-10-17 as the VizieR search thread's clock reported (the thread "VizieR table search (2026-10-17)", started 11:29 AM EDT and finished 11:39 AM EDT as the coordinator reported; read-only) and filed by the TODO thread on 2026-10-17. It answers Boss's decision 6 of the gap-filling plan ([gap-filling-plan.md](gap-filling-plan.md) section 8): the plan's GR5 fill stays, and GR3 rank 2 becomes the [R0865] fit. The scripts and downloaded copies sit in `/mnt/project-files/reports/vizier-courteau-comeron-2026-10-17-files/`. The replacement of the stand-ins in GEN.401 and GEN.402 is a new stand-in version after GEN.401 merges (Core 1's gate ends about 12:20 PM EDT).

Written 11:40 AM EDT on 2026-10-17 for the coordinator and TODO. Serves GEN.401 and GEN.402. Nothing here changes code, TODO.md or docs/plan. The result replaces the plan's fills only after GEN.401 merges (Core 1's gate runs until about 12:20 PM EDT).

Scripts, the downloaded copies and the raw outputs are in `/mnt/project-files/reports/vizier-courteau-comeron-2026-10-17-files/` (listed at the end).

## 1. Answer in short

1. **[R0908]: the 64-galaxy table is not on VizieR and not reachable anywhere I could open.** The paper in the plan, astro-ph/0610721, is the Letter ApJ 655, L21 (not ApJ 671, 203), and its own text uses 154 spirals, 24 ellipticals, 54 lenticulars and 11 dwarfs in Fig. 1 and concentrations for 77 galaxies in Fig. 2. The number 64 does not appear in it. The authors' file (www.astro.queensu.ca/~courteau/data/VSigmaC28.txt) is refused by the proxy (connection reset, HTTP 503). The residual scatter about Eq. 5 (which needs the SDSS concentration C28 of each galaxy) therefore **cannot be computed**.
2. **What I could compute instead is the raw scatter of log10(Vc/σ0) for the spirals that feed Courteau's compilation**, from two of its source samples that are reachable: [R3619] (VizieR J/A+A/366/68, table1; 69 usable spirals) and [R3540] (Table 1 of astro-ph/0503649, transcribed from the PDF; 40 high-surface-brightness plus 8 low-surface-brightness spirals). The standard deviation of log10(Vc/σ0) is **0.114 dex (PMS01, N = 69), 0.119 dex (P05 HSB, N = 40), 0.123 dex for the 113 de-duplicated spirals, 0.118 dex without the LSB galaxies (N = 105)**. About Hubble-type means (a weak stand-in for concentration) it is **0.104 dex** (PMS01, five type groups).
3. **Against the plan's GR5 fill (log-normal 0.09 dex, bracket 0.06 to 0.12): the raw data sit at the top edge of the bracket and above the centre.** The values 0.104 to 0.123 are inside or within 0.003 dex of the 0.12 edge, so the plan's rule ("the table gives an RMS outside 0.06 to 0.12") is **not triggered**, but 0.09 is below every number I could measure. These numbers include measurement error and the unexplained concentration dependence, so the true scatter about Eq. 5 is smaller and unknown (see section 3.4).
4. **Comerón thick discs: the usable table is [R0865], on VizieR as J/A+A/610/A5** (tableb.dat 124 galaxies with at least two distinct discs, tablec.dat 17 without; per-galaxy circular velocity, gas, thin, thick and central-mass masses, scale heights and scale lengths). It supersedes the 2011 sample (47 galaxies) and the 2012 sample (70) by re-fitting them with a better point spread function and more than doubling the sample to 141.
5. **Fit for GR3 rank 2 (table B, N = 124): log10(M_thick/M_thin) = 2.113 − 1.086 log10(Vc), scatter 0.268 dex** (slope error ± 0.178, correlation coefficient r = −0.48). Against the plan's 0.2 dex scatter the data give **0.27 dex, outside the plan's value** (0.32 dex if the 17 single-disc galaxies of table C are added).
6. **The plan's hot-component bound does not hold galaxy by galaxy.** Built from [R0860: Eq. 1] on each galaxy's own stellar mass, the bound is exceeded by **48 of 124 galaxies (39 percent)**; the median galaxy sits 0.10 dex below it, the 84th percentile 0.18 dex above it, the maximum 0.57 dex above. The bound is a statement about the mean hot-component mass, which has 0.13 dex scatter itself. The plan's "any seed breaks the bound (a test, not a tolerance)" would reject 39 percent of real galaxies, so the test needs a margin or should be applied to the mean relation only (section 4.4).

## 2. What I checked and what is not reachable

| Item | Result |
|---|---|
| VizieR web service (`vizier.cds.unistra.fr/viz-bin/votable`, `/VizieR?-source=`) | reachable; catalogue list queries work with `-source=` wildcards and keywords |
| VizieR TAP (`tapvizier.cds.unistra.fr/TAPVizieR/tap/sync`) | HTTP 503 "No server is available" at about 11:30 AM EDT |
| CDS archive (`cdsarc.cds.unistra.fr/ftp/J/...` ReadMe and data files) | reachable |
| arXiv abstracts, PDFs and the arXiv API | reachable |
| Authors' page (www.astro.queensu.ca/~courteau/data/VSigmaC28.txt) | not reachable: connection reset (https) and HTTP 503 (http) |
| J/ApJ/671/203 | not found in VizieR (the relay's journal reference is for a different paper; the Letter is ApJ 655, L21) |
| J/ApJ/655/L21 (the actual Letter) | not found in VizieR |
| Keyword and author searches for "Courteau" | only [R0906] (J/AJ/114/2402, 304 Sb-Sc galaxies, no σ0) and unrelated catalogues |
| J/ApJ/631/785 [R3540], J/ApJ/578/90 [R1357], J/MNRAS/341/L44 [R0227], J/MNRAS/371/1912 [R0311], J/A+AS/144/53 [R2494] | none in VizieR (all 404 on the archive) |
| J/A+A/366/68 [R3619] | found: table1 (93 rows), table3 (155 rows) |
| [R0861] and 2012 (ApJ 759, 98; [R0863]) | no VizieR table (J/ApJ/741/28 and J/ApJ/759/98 not found, archive 404); the tables exist only in the papers and are superseded by 2018 |
| [R0864] as named in the brief (A&A 562, A121) | that id is **ARRAKIS**, the atlas of resonance rings (J/A+A/562/A121, 1343 features), not thick discs. The thick-disc paper of 2014 is A&A 571, A58 [R0864]; it has no VizieR table (archive 404) |
| [R0865] | **found**: J/A+A/610/A5 (Comerón, Salo, Knapen, "The reports of thick discs' deaths are greatly exaggerated"), [R0865] |

The relay named "Courteau, Dutton, van [R4492]". That is a different paper ("Scaling relations of spiral galaxies"). The Vc/σ0 work in the plan is the Letter by Courteau, McDonald, Widrow and Holtzman, ApJ 655, L21 (2007), [R0908]. I did not look for a VizieR table of the other paper because the plan's equation and citation are the Letter's.

## 3. GR5: scatter of log10(Vc/σ0)

### 3.1 Where the Letter's data come from (read from astro-ph/0610721 section 2)

Courteau's compilation is built on the P05 list (40 HSB and 8 LSB spirals, 20 ellipticals with dynamical-model Vc and 5 ellipticals with HI Vc), extended by the 66 spirals of [R3619] (PMS01), plus [R0311] and the Southern galaxies of [R0227]. Fig. 1 shows 154 spirals, 24 ellipticals, 54 lenticulars and 11 dwarf irregulars. Concentrations C28 were measured from SDSS i-band images for 32 spirals, 12 ellipticals and 33 lenticulars. The Letter prints no scatter, no per-galaxy table and no "64".

### 3.2 Data used

- **PMS01, VizieR J/A+A/366/68 table1**, columns `Vmrot` (bytes 43 to 45, km/s, "maximum velocity of rotation, from PGC96") and `Sig0` (bytes 47 to 49, km/s, "central velocity dispersion from Hypercat"), plus `MType`. 93 rows, of which 69 have both values positive (the Letter says 66 spirals; the three extra galaxies may have been dropped by Courteau's cuts, which are not stated).
- **P05, Table 1 of [R3540]** [R3540], columns σc and Vc for the 40 HSB and the first 8 LSB rows. Parsed from pdftotext output; the 48 rows are in `parse_p05.py` output form and the excerpts are saved.
- Four galaxies appear in both samples (NGC 2844, NGC 5055, NGC 7331, NGC 7606); the union takes the P05 row for them (113 galaxies).

### 3.3 Numbers (log10 of Vc over σ0, standard deviation with N − 1)

| Sample | N | mean | standard deviation (dex) |
|---|---|---|---|
| PMS01 spirals (Vmrot and Sig0) | 69 | +0.206 (ratio 1.61) | **0.114** |
| P05 HSB spirals | 40 | +0.246 (1.76) | **0.119** |
| P05 LSB spirals | 8 | +0.360 (2.29) | 0.107 |
| P05 HSB + LSB | 48 | +0.265 (1.84) | 0.124 |
| Union PMS01 + P05, de-duplicated | 113 | +0.230 (1.70) | **0.123** |
| Union without LSB | 105 | +0.220 (1.66) | **0.118** |

Further PMS01 results (`analyse_pms01.out.txt`): the 16th to 84th percentile of log10(Vc/σ0) is +0.131 to +0.298; the mean by morphological group rises from +0.139 (Sa) through +0.177 (Sab, Sb) and +0.267 (Sbc) to +0.281 (Sc), which is the dependence on concentration the Letter describes; **the RMS about the five type means is 0.104 dex**; a fit of log(Vc/σ0) against log σ0 leaves 0.094 dex.

### 3.4 How to read these against the plan

- The brief asked for the RMS of log(Vc/σ0). The Letter's relation is log(Vc/σ0) = 0.63 − 0.11 C28 (Eq. 5), so the generator's scatter is the scatter **about that line at fixed C28** (or fixed B/T). The raw standard deviations above are **larger than that**, because they include the whole morphology trend (the mean changes from 1.4 to 1.9 in Vc/σ0 across Hubble types) and the measurement error of an HI line width and an aperture-dependent σ0. Removing even a crude proxy (five type groups) lowers 0.114 to 0.104. I cannot estimate the final residual without C28.
- **Compared with the plan (centre 0.09 dex, bracket 0.06 to 0.12):** the measurable scatters (0.094 to 0.123) lie inside the bracket or touch its upper edge (0.123 for the 113-galaxy union). The plan's replacement rule needs a value outside the bracket, so the fill stays. The centre 0.09 is slightly below the data; 0.10 to 0.11 would match the type-corrected PMS01 value and P05.
- **A consistency check the Letter allows:** the relation Vc/σ0 = 2 (1 − B/T) for spirals. At the sample mean ratio of 1.6 to 1.8 this means B/T of about 0.1 to 0.2, which is a plausible spiral mix and agrees with the plan's use of the row 6 B/T law.
- **Selection bias to state:** HI width and stellar-dispersion samples favour galaxies with well-defined flat curves and measurable absorption lines. PMS01 is mostly the Heraudeau sample of bright spirals. The LSB galaxies (mean ratio 2.29) sit at the high end, as the Letter says.

### 3.5 What would settle it

The file `VSigmaC28.txt` from the authors' page (blocked here) or the 77-galaxy concentration list; or SDSS i-band r20 and r80 for the 113 spirals above. Neither is reachable from this session.

## 4. GR3 rank 2: thick-to-thin ratios from [R0865]

### 4.1 Table and columns (J/A+A/610/A5, ReadMe in `c18/ReadMe`)

Two data files with identical columns: `tableb.dat` (124 galaxies with at least two distinct discs) and `tablec.dat` (17 without). Columns used: `vc` (km/s, maximum circular velocity, references to HyperLeda, Springob, Courtois, Bedregal, ATLAS3D and others, listed in note 3), `GasM`, `ThinM`, `ThickM`, `CMCM` (10^9 solar masses), `sc-thin`, `sc-thick` (scale heights in arcsec), and for each of total, thin and thick disc the scale lengths `h1` to `h4` of up to four exponential sections with their break radii. Distances, absolute 3.6 µm magnitudes and the fitting types are also given. Source: 3.6 µm images of edge-on galaxies from S4G and its early-type extension, 141 galaxies in all. Masses assume (Υ_thick/Υ_thin) at 3.6 µm = 1.2, which the authors call the smallest and most conservative value; M_T/M_t scales roughly with that ratio, so every mass ratio below scales with it.

### 4.2 Fits (table B only unless stated; `analyse_c18.out.txt`)

| Relation | N | Fit | Scatter about the fit | r |
|---|---|---|---|---|
| log(M_thick/M_thin) against log Vc | 124 | 2.113 − 1.086 log Vc (slope ± 0.178) | **0.268 dex** | −0.48 |
| same, tables B and C | 141 | 2.382 − 1.190 log Vc (± 0.193) | 0.318 dex | −0.46 |
| log(M_thick/M_thin) against log M* (M* = thin + thick + central mass) | 124 | 2.849 − 0.302 log M* (± 0.053) | 0.272 dex | −0.46 |
| log(height ratio) against log Vc | 124 | 0.240 + 0.189 log Vc (± 0.086) | 0.129 dex | 0.20 |
| log(length ratio) against log Vc, outermost exponential section of each disc | 104 | −1.638 + 0.817 log Vc (± 0.233) | 0.316 dex | 0.33 |

Medians: thick-to-thin scale-height ratio **4.24** (mean 4.61, scatter 0.13 dex, nearly independent of Vc); thick-to-thin scale-length ratio **1.46** using each disc's outermost section (1.57 for the six galaxies where both discs are a single exponential, Type I, too few to fit). The scale-length ratio depends on which section is chosen; the first sections (`h1`) are often hundreds of arcseconds for a flat inner region and are not usable.

Binned ratio of thick to thin mass (median, 16th to 84th percentile):

| Vc (km/s) | N | median | 16 to 84 percent |
|---|---|---|---|
| below 100 | 16 | 1.13 | 0.81 to 2.02 |
| 100 to 150 | 55 | 0.73 | 0.42 to 1.35 |
| 150 to 200 | 44 | 0.45 | 0.23 to 0.77 |
| above 200 | 9 | 0.42 | 0.14 to 0.94 |

These agree with the statements printed in the 2018 paper (comparable masses below 120 km/s, 0.2 to 1 above, up to about 3 at the low-mass end).

### 4.3 Comparison with the plan's fills (gap-filling-plan section 4.4 and the GR3 row)

Values at the three masses of the plan's table (plan Vc about 89, 181, 255 km/s). The "upper bound" is the plan's M_h/(M* − M_h) with Eq. 1 of [R0860], recomputed here (0.98, 0.64, 0.53; the plan prints 0.96, 0.64, 0.54).

| M* | Plan upper bound (no classical bulge) | Plan, Yoachim and Dalcanton formula | **[R0865] fit** | Comerón fit, 1-sigma range (× or ÷ 1.86) |
|---|---|---|---|---|
| 3e9 (Vc 89) | 0.98 | 0.69 | **0.99** | 0.53 to 1.84 |
| 5e10 (Vc 181) | 0.64 | 0.14 | **0.46** | 0.25 to 0.86 |
| 2e11 (Vc 255) | 0.53 | 0.06 | **0.32** | 0.17 to 0.60 |

- **The fit lies close to the bound with B/T about 0.1** (plan values 0.76, 0.48, 0.38 at the three masses), not with the Yoachim and Dalcanton formula at high mass (0.14 and 0.06 are 3 to 5 times lower than the measured 0.46 and 0.32). That resolves the conflict the plan describes in favour of [R0865] on the data side, and the plan's morphology explanation (high B/T early-type spirals have small thick discs) is consistent with the fact that Comerón's central mass concentration is a separate component from the thick disc (it carries the bulge mass), so a large bulge leaves less for the thick disc in the bound.
- **Scatter: 0.27 dex measured, 0.2 dex in the plan.** The plan's fill is therefore too narrow; 0.27 dex (0.32 with the single-disc galaxies) is the data value. This falls outside the plan's value of 0.2 dex, so the plan's own rule ("a published fit equation with scatter replaces the scatter only") applies and this is that fit, with the caveats in 4.4 and the Υ ratio of 1.2.
- **Length ratio:** the plan uses 1.15 (Vc/100)^−0.2 from Yoachim and Dalcanton; the [R0865] outermost-section median is 1.46 with 0.32 dex scatter, rising with Vc (slope 0.82 ± 0.23). These differ in direction of the Vc trend (Yoachim falls, Comerón rises). The Comerón numbers depend on the section choice, so I do not recommend replacing the Yoachim length ratio without a decision on that choice.
- **Height ratio:** the plan uses 2.5 ± 30 percent (Yoachim and Dalcanton). [R0865] gives a median of 4.24 (scatter 0.13 dex, which is 35 percent). These differ by a factor 1.7; the Comerón heights are from edge-on 3.6 µm decompositions with their own definition of the scale height, so the two are not the same measurement.

### 4.4 The hot-component bound against the same galaxies

- The thick disc plus the central mass concentration (the hot component in [R0860]'s sense) of these 124 galaxies matches [R0860: Eq. 1] on M* on average (mean deviation 0.000 dex, standard deviation 0.132 dex, median +0.016). That is **not an independent check**, since the 2021 relation was built from the same family of fits; it only shows that the plan transcribed Eq. 1 correctly.
- Taking the bound **galaxy by galaxy** with that mean relation, 48 of 124 galaxies (39 percent) have M_thick/M_thin above the bound, 22 percent by more than the 0.13 dex scatter of the hot component itself, and the worst case is 0.57 dex above. The median is 0.10 dex below. So the plan's test text "any seed breaks the bound" will fail for roughly four draws in ten if the generator reproduces the observed scatter. Options for the plan (not my decision): apply the bound to the mean ratio only; add a 0.13 dex margin to the bound (then 22 percent of real galaxies still exceed it); or bound the hot-component fraction instead of the ratio. The bound as a mean statement is supported.

## 5. Corrections to the plan's references

| In the plan or brief | What I found |
|---|---|
| [R0908] = ApJ 671, 203, "64 galaxies" | [R0908] is ApJ 655, L21; no "64" in its text; samples are 154/24/54/11 (Fig. 1) and 32/12/33 with concentrations (Fig. 2). `docs/design/galaxy/galaxy-laws-research.md` row "(64 galaxies)" and R5 should drop the number |
| [R0864] = A&A 562, A121 | ARRAKIS ring atlas; the thick-disc paper is A&A 571, A58 [R0864] |
| [R0861], 2012 tables on VizieR | not on VizieR; superseded by J/A+A/610/A5 (2018, 141 galaxies) |
| [R0865] "prints ranges only" | the numbers are in the VizieR table J/A+A/610/A5 with per-galaxy masses, so a fit with scatter is possible |

## 6. Limits

- The Courteau residual scatter about Eq. 5 is not computed; only raw scatter of the source samples is. The sources are 113 of Fig. 1's 154 spirals (about 73 percent) with Vc and σ0 from several heterogeneous papers; Courteau notes a 24 km/s systematic offset in Vc and 9 km/s in σ0 between samples.
- P05 numbers were read from pdftotext output of the arXiv PDF (not a machine-readable table); I checked the count (40 plus 8) and spot-read the first and last rows, but did not re-transcribe by hand.
- The PMS01 type groups are crude (five groups, 11 to 22 galaxies each).
- The Comerón masses carry a systematic from the assumed Υ ratio of 1.2; the fit is for that choice. The fits treat vc as exact; its errors are not in the table, so the scatter includes them.
- Nothing here was compared with a Milky Way value, and no value was copied from one.

## 7. Files (`vizier-courteau-comeron-2026-10-17-files/`)

| File | What |
|---|---|
| `c18/ReadMe`, `c18/tableb.dat`, `c18/tablec.dat` | [R0865], J/A+A/610/A5, copied from cdsarc.cds.unistra.fr |
| `c18_meta.xml`, `c14_meta.xml` | VizieR metadata queries (the second shows 562/A121 is ARRAKIS) |
| `pms01/ReadMe`, `pms01/table1.dat`, `pms01/table3.dat` | [R3619], J/A+A/366/68 |
| `p05/pizzella2005-astro-ph-0503649.txt`, `p05/table1-hsb-excerpt.txt`, `p05/table1-lsb-excerpt.txt` | text of the arXiv PDF and the two Table 1 excerpts |
| `parse_c18.py`, `analyse_c18.py`, `analyse_c18.out.txt` | parsing and fits for GR3 |
| `analyse_pms01.py`, `analyse_pms01.out.txt` | PMS01 scatter, type groups |
| `parse_p05.py`, `analyse_courteau_samples.py`, `analyse_courteau_samples.out.txt` | P05 parsing and the combined scatter table |

Run each script with `python3 -I`; the command is in the file's docstring.
