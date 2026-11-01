# Galaxy types: one family per known galaxy type

**Status:** checked design, for Boss to read and approve before DB.66 and GEN.379 are built. Boss asked for the feature on 2026-10-16 at 8:48 PM EDT for version 9.4. GEN.378 (this document, written on 2026-10-21 in the Blocker 2 lane) turns the plan into the checked design: the family list, the draw order, one table row per parameter with its law, its source checked against the paper, its range and scatter, and the component it feeds, the mass limits, the reference galaxies, the columns DB.66 adds, and the research requests. Nothing here is built except the spiral family's first version (GEN.400, GEN.401, GEN.532, GEN.533), which is the first family of this design.

Boss (2026-10-16, 8:48 PM EDT): "For v9.4 or so, I want to add a feature that we have different galactic models for generation based on known galaxy types, as a chain of required items counted back from 9.4."

## 1. What it means under the standing rules

- **No fixed content** (Boss, 2026-10-16, 4:11 and 4:14 PM EDT). A galaxy type is a parameter family, never a copy of a named galaxy: every number a family needs is calculated from the galaxy's stellar mass and a few seeded draws through published scaling relations, with a draw only where the science gives a distribution. The real galaxies appear only as reference points the plausibility tests (TEST.204) compare against (section 9).
- **No home position** (Boss, 2026-10-16, 9:43 PM EDT: "we will choose viewing based on where sectors are generated"). No family stores a Sun or home position and the design has no home row; GEN.388 is dropped. Per-position calculations replace every Sun-radius use.
- **Source rules** (Boss): a verified source wins over an unverified one; two verified sources that disagree are both kept and the draw is widened over both (the gap plan's policy, `docs/design/research/gap-filling-plan.md` section 8); a value no paper gives stays a labelled placeholder (kind P) with a research request.
- **One rule set.** The density model (`galaxy/density.py`), the potential (GEN.338, 8.3), the halo (GEN.234), the satellites (GEN.369), the clusters (GEN.297), the scatters (GEN.301) and the sector profile already take a galaxy shape. A family changes the shape's contents, so the chain below makes the shape carry its components and leaves the readers alone.
- **Multiple galaxies stays the existing "later" note.** A database still holds one galaxy; the family is a property of that galaxy (multiple-galaxies.md).
- **Every scatter** is a log-normal about the median law, drawn through the one random wrapper (`planetgen.util.draw.Stream`) from `unit_seed(galaxy_seed, "<family>", "<link>")`, cut at 3 sigma, so one seed gives one galaxy and two seeds give two (the rule of GEN.400 and GEN.401).

## 2. How the sources were checked

The network was open to arXiv on 2026-10-21, so every paper below was read in this lane from its arXiv text or LaTeX source (tables read from the source, not from a rendered page). Tags in the last column of the tables: **V** the paper says what the row says; **C** the number was computed here from V values (arithmetic reproducible); **P** a labelled placeholder, a value no paper read gives; **N** not found in any paper read, which is a research request.

| Paper | arXiv | What it gave |
|---|---|---|
| [R2310], GAMA stellar mass functions by Hubble type | 1407.7555 | Tables 1 and 3: the mass function of each morphological class, which gives the family fractions |
| [R2558], GAMA bulges, discs and spheroids | 1607.01096 | Table of fits `R_e = a (M*/1e10)^b` by Hubble type and by component |
| [R4048], SDSS sizes | astro-ph/0301527 | Early-type fit and scatter, and a probable misprint in its Table 1 (section 6.3) |
| [R1656], bulge-to-disc ratios | 0805.3565 | Tables 3 and 4: S0 and spiral bulge and disc parameters by type |
| [R1267], ATLAS3D III | 1102.4444 | Fast and slow rotators: fractions and masses |
| [R4654], ATLAS3D XXIV | 1408.1099 | Intrinsic axis ratios of fast and slow rotators |
| [R1653], review of galaxy structure | 1108.0997 | Eq. 7: Sérsic index against luminosity for ellipticals |
| [R0655], ELVES dwarf structure | 2105.03435 | Table of mass-size fits and table of intrinsic shapes |
| ELVES-Field, isolated dwarf quenched fractions | 2602.16778 | Quenched fraction of isolated dwarfs against stellar mass; isolated dwarfs about 20 percent smaller |
| [R2561], barred S0 galaxies in Coma | 1401.3775 | Bar fraction of S0 galaxies |
| Greene, Strader and Ho 2020, intermediate-mass black holes | 1911.09678 | Occupation fraction of black holes in low-mass galaxies |
| [R1831], globular cluster catalogue | 1306.2247 | Specific frequency is about 1 for every galaxy type at intermediate mass |

Already checked for GEN.368 and cited from its report (`/mnt/project-files/reports/gen368-paper-check-2026-10-17.md`): [R0238, R1200, R2484, R3181], Kormendy and Ho 2013, [R1499, R2368, R3273, R4423]. Also used: `docs/design/galaxy/galaxy-laws-research.md` (GEN.406), `docs/design/research/gap-filling-plan.md`, `/mnt/project-files/reports/galaxy-size-sources-2026-10-21.md`.

## 3. The families

| Key | Family | What it has | What it leaves out |
|---|---|---|---|
| `spiral` | Spiral, unbarred or barred (`bar` is a draw inside the family, not a key), Hubble stage Sa to Sm | thin and thick disk, bulge, arms, bar by a draw, stellar halo, dark halo | |
| `lenticular` | Lenticular (S0) | thin and thick disk, a larger bulge, bar by a draw, old stars, almost no gas | arms, young stars |
| `elliptical` | Elliptical | one Sérsic spheroid with axis ratios (round and triaxial, or flat and oblate), old stars, dark halo | disk, arms, bar, young stars |
| `irregular` | Dwarf irregular | a clumpy exponential or low-n Sérsic body with young stars and gas | arms, bar, bulge, thick disk |
| `dwarf_spheroid` | Dwarf spheroid | a low-n Sérsic or Plummer body, old and metal-poor | disk, arms, bar, young stars |

Six outlines are drawn for them (UX.135): spiral, barred spiral, lenticular, elliptical, irregular and dwarf spheroid. A barred spiral is the `spiral` family with `bar` set.

## 4. The draw order

Each step reads only the steps before it. Steps 1 and 3 are GEN.400's `draw_mass_chain` (no second draw of either).

1. **Stellar mass** `M*` from the double Schechter function of [R0238] (inverse CDF in log mass), inside the mass limits of section 8.
2. **Family** from the fractions of section 5, given `M*`.
3. **Halo mass** `Mh` from `M*` by the stellar-to-halo relation of [R3181] with the approximation Boss chose (10:11 PM EDT on 2026-10-16), then `r200` and `c200`. The relation is fitted to all galaxies, so it is used for every family (research request R7).
4. **The family's shape**: the rows of section 6, each a median law and a scatter.
5. **Shared links** of section 7: black hole, nuclear cluster, globular clusters, satellites, metallicity, populations. Their laws are the spiral family's (GEN.368 rows 7 to 10, 13 and 14), with the family change stated per row.

A pinned family (`--family`, UX.134) fixes step 2 only. Every size and count is still drawn.

## 5. The family fractions (step 2)

**The law.** [R2310] give a single Schechter function for each of four morphological classes (Table 3), fitted to 2,711 GAMA galaxies with `M*` above 1e9 solar. The probability of a class at a stellar mass is its mass function divided by the sum of the four, so the fractions are calculated from fitted mass functions, not read from a binned table. With `x = M / M*_class`, the density per dex of a class is `phi*_class x^(alpha_class + 1) exp(-x)` (the common factor `ln 10` cancels in the ratio).

| Class | log M* | alpha | phi* (1e-3 per dex per Mpc^3) | Tag |
|---|---|---|---|---|
| E | 10.94 | -0.79 | 0.85 | V (Table 3) |
| S0-Sa | 10.25 | +0.87 | 2.38 | V |
| Sab-Scd | 10.09 | -0.01 | 3.57 | V |
| Sd-Irr | 9.57 | -1.36 | 3.40 | V |

**Check.** Integrating these four functions above a mass bound reproduces the paper's own Table 1 (number fractions above the bound) within 0.06: above 1e9 the functions give E 0.18, S0-Sa 0.19, Sab-Scd 0.28, Sd-Irr 0.35 against the table's 0.19, 0.18, 0.28, 0.31; above 1e10 they give 0.25, 0.41, 0.33, 0.01 against 0.24, 0.42, 0.32, 0.02; above 1e11, 0.78, 0.21, 0.00, 0.00 against 0.72, 0.22, 0.06, 0.00. So the density convention is right. [C]

**The fractions at a mass** (computed from the table above). [C]

| log M* | E | S0-Sa | Sab-Scd | Sd-Irr |
|---|---|---|---|---|
| 9.0 | 0.07 | 0.00 | 0.06 | 0.87 |
| 9.5 | 0.15 | 0.03 | 0.26 | 0.56 |
| 10.0 | 0.20 | 0.19 | 0.54 | 0.07 |
| 10.5 | 0.20 | 0.50 | 0.30 | 0.00 |
| 11.0 | 0.55 | 0.43 | 0.02 | 0.00 |
| 11.5 | 1.00 | 0.00 | 0.00 | 0.00 |

**From classes to families.**

- `E` becomes `elliptical`.
- `Sab-Scd` and `Sd-Irr` become `spiral` at or above 3e9 (the lower edge of the first family, Boss, 9:40 PM EDT on 2026-10-16) and `irregular` below it. A disk-dominated galaxy under 3e9 has no arm law to draw from (section 6.1).
- `S0-Sa` splits into `lenticular` and `spiral` (the Sa part). Kelvin et al. merge the two, so the split is derived: **0.55 lenticular and 0.45 Sa spiral, constant in mass** (range 0.50 to 0.62; [R0626, R2310] and ATLAS3D, no catalogue separates them by mass; `reports/s0-sa-galaxies-2026-10-21.md`; Decided by Boss, card tapped 8:19 PM EDT on 2026-10-21) [D], research request R1 answered in part (the mass dependence stays open). Below 3e9 the Sa part becomes `irregular`.
- **Below 1e9** (under Kelvin's limit) the family is `irregular` or `dwarf_spheroid`. The dwarf spheroid share is the isolated-dwarf quenched fraction of ELVES-Field (2602.16778 abstract): essentially zero above 1e7 solar and about 30 percent below it [V]. As a step: `p(dwarf_spheroid) = 0.30` below 1e7 and `0` from 1e7 to 1e9 [P for the shape of the step; the two levels are V]. Between 1e9 and 3e9 the Kelvin fractions apply with the rules above.

The quenched fraction rests on a few dozen isolated dwarfs ("the statistics on the quenched fraction are rather poor", the paper), so R2 asks for a better one.

**What the shape stand-in in GEN.401 becomes.** `shape_draw.py` today reads Kelvin's Table 1 (a cumulative fraction above a bound) as a local fraction to choose the spiral's stage (early or late). GEN.382 replaces it with this law: the spiral's Hubble stage `T` is drawn within its class (Sa from the S0-Sa share, Sb to Sc from Sab-Scd, Sd to Sm from Sd-Irr), uniformly over the stages of the class [P], and the bulge-to-disc ratio is that stage's median (section 6.1).

## 6. The shape table

*Feeds* names the component of the density model (`galaxy/density.py`, GEN.379) or the stored value that takes the number.

### 6.1 Spiral family (built in GEN.401, extended by GEN.382)

The rows of GEN.368's table that belong to this family are rows 4 to 6 (disk length and height, bulge, bar, arms). Rows 1 to 3 (the mass chain), 7 to 10 and 12 to 14 are shared by every family (section 7). Row 11 is removed.

| Parameter | Law | Range and scatter | Source, check | Feeds |
|---|---|---|---|---|
| Disk scale length | late-type size-mass relation `R50 = gamma M^alpha (1 + M/M0)^(beta - alpha)` (gamma 0.10, alpha 0.14, beta 0.39, M0 3.98e10), `Rd = R50 / 1.678` | 3e9 to 2e11; scatter 0.20 dex falling to 0.15 dex | [R4048], as in GEN.368 row 4 (Boss chose Shen, 10:09 PM EDT on 2026-10-16) [V]. Cross-check: [R2558] final z = 0 discs `R_e = 5.141 (M/1e10)^0.274` kpc [V], 8.0 kpc at 5e10 against Shen's 3.9: the two differ by 2, as in GEN.368 Q5 (Lange fits r-band elliptical-aperture sizes, Shen z-band circular ones) | thin disk length |
| Disk scale height | `hz = Rd / q`, `q = 7.3` | 0.13 dex | [R2484] [V] | thin disk height |
| Thick disk | mass, length and height ratios to the thin disk from the circular speed | see `structure_draw.py` | [R0865, R4824] (both kept, widened) [V] | thick disk |
| Hubble stage `T` | drawn within the class of section 5 | uniform in the class [P] | stage groups of [R1656: Table 3] | bulge, bar |
| Bulge-to-disc `log10(B/D)` | median of the stage with the table's half widths: Sa -0.34 (+0.40 -0.32), Sab -0.54 (+0.53 -0.41), Sb -0.60 (+0.28 -0.49), Sbc -0.82 (+0.28 -0.42), Sc -1.06 (+0.43 -0.34), Scd -1.23 (+0.75 -0.28), Sd -1.06 (+0.16 -0.50), Sdm -1.49, Sm -1.57, Irr -1.39 (the last three from 4, 2 and 1 galaxies) | `B/T = 1 / (1 + 1/(B/D))`; K-band flux ratios | [R1656: Table 4], read from the LaTeX source [V]. GEN.401 uses the two class medians, T = 1 to 3 (-0.51) and T = 6 to 9 (-1.19), also in the table [V] | bulge amplitude |
| Bulge size | `Re / h` 0.24 (+0.22 -0.10) for T = 1 to 3 and 0.22 (+0.08 -0.09) for T = 6 to 9; bulge Sérsic `n` 2.18 (+1.82 -0.77) and 1.15 (+0.95 -0.55) | | [R1656: Tables 3 and 4] [V] | bulge radius |
| Bar | present with the fraction of `structure_draw.bar_fraction`; length in disk scale lengths 1.29 +/- 0.54 (S0 to Sab), 0.57 +/- 0.44 (Sc to Sd) | 0.2 to 3 | [R2934] [V]; [R1289] [V]; `docs/design/galaxy/gen533-bar-sources.md` | bar |
| Arm count and pitch | count shares 3.9, 64.5, 17.8, 6.5, 7.3 percent for 1 to 5 arms; pitch falls with mass from 26 to 13.8 degrees over log M* 9.6 to 11.1 | 10 to 40 degrees; 5.5 degrees | [R1836]; Yu and Ho 2019 (they disagree on the sign against bulge share, kept as the 50/50 mixture of the gap plan, GX2) [V] | arms |

### 6.2 Lenticular family (GEN.384)

A smooth disk with a large bulge and no arms.

| Parameter | Law | Range and scatter | Source, check | Feeds |
|---|---|---|---|---|
| Total half-light radius (a check, not a draw) | `R_e = 2.574 (M*/1e10)^0.326` kpc (S(B)0-S(B)a) | 1.22 kpc at 1e9, 2.57 at 1e10, 6.84 at 2e11 | [R2558] [V]. The test compares the drawn galaxy's half-light radius with this | tests |
| Disk length | the early-type disc component relation `R_e,disc = 4.55 (M_disc/1e10)^0.247` kpc with `M_disc = (1 - B/T) M*`, `Rd = R_e,disc / 1.678` | scatter 0.15 dex [P: the late-type disk scatter, no S0 value read] | [R2556, R2558] [V]. Cross-check: Graham and Worley's S0 disc length `h = 2.83 (+3.10 -1.10)` kpc, 2.60 for the class T <= 0, against 2.5 kpc from the relation at `B/T = 0.30` and `M* = 1e10` [C] | thin disk length |
| Bulge-to-disc | `log10(B/D) = -0.37 (+0.40 -0.35)` for T <= 0 (S0 and S0/a alone: -0.54 and -0.31); `B/T` about 0.30 | | [R1656: Table 4] [V] | bulge amplitude |
| Bulge size and index | `Re / h = 0.17 (+0.24 -0.08)`, bulge `n = 2.08 (+1.76 -0.74)` | | [R1656: Tables 3 and 4] [V] | bulge radius |
| Disk height | as the spiral's, `q = 7.3`, 0.13 dex | | no S0 value found [P], R3 | thin disk height |
| Bar | present with probability 0.48 to 0.72, both kept | | [R2561] (the central S0 galaxies of Coma: 72 percent by decomposition, 48 percent by ellipse fitting) [V]; a cluster sample, so the isolated value is [P], R3 | bar |
| Arms | none, amplitude 0 | | | arms |
| Populations | old; gas fraction at the floor | | GEN.386 | profile |

### 6.3 Elliptical family (GEN.383, GEN.381)

| Parameter | Law | Range and scatter | Source, check | Feeds |
|---|---|---|---|---|
| Effective radius | `R_e = 2.114 (M*/1e10)^0.329` kpc (all ellipticals) | 1.0 kpc at 1e9, 2.1 at 1e10, 5.7 at 2e11; scatter 0.15 dex | [R2556, R2558] [V]. The same paper's high-mass fit (`M* >= 1e10`) is `1.382 (M*/1e10)^0.643`, 9.5 kpc at 2e11, a factor 1.7 above the single fit; its low-mass fit (`M* < 1e10`) is `1.978 (M*/1e10)^0.265`. The two are not continuous at 1e10, so both are kept: the draw's median is the single fit and the scatter is widened by the gap (R4). Scatter: [R4048] say the early-type scatter is "similar to that of bright late-type galaxies" and give 0.34 in ln R at high mass (0.15 dex) [V] | spheroid radius |
| Sérsic index `n` | `M_B = -9.4 log10 n - 14.3`, i.e. `n = 6.6, 4.0, 2.5, 1.5` at `M_B = -22, -20, -18, -16` [C] | `n` kept between 1 and 8 [P] | [R1653: Eq. 7], matched to the compilation of dwarf and bright ellipticals [V]. It is in luminosity: the step from `M*` to `M_B` needs a stellar mass-to-light ratio no paper read here gives [N], R5; the scatter is not printed | spheroid index |
| Axis ratios | two modes. Fast rotators (flat): intrinsic `q` Gaussian with mean 0.25 to 0.33 and sigma 0.11 to 0.14, `p = 0.99` (oblate). Slow rotators (round): `q` Gaussian with mean 0.61 to 0.66 and sigma 0.08 to 0.09, `p` from `ln(1 - p)` between -5 and -3 (0.95 to 0.99, mildly triaxial) | `q` cut to 0.05 to 1 | [R4654]: direct inversion 0.25 (sigma 0.14) and 0.63 (0.09), the 2-D triaxial fit 0.33 (0.11) and 0.66 (0.08), after removing four counter-rotating discs 0.61 (0.09); the fits are kept as a bracket [V]. `p` for slow rotators is only bracketed ("very strong indications that they are (mildly) triaxial"; the best fit sits at the grid edge, `mu_Y = -5`) [V] | spheroid axis ratios |
| Slow-rotator share | 14 percent of the early-type sample overall; 4 of the 36 slow rotators are low-mass counter-rotating discs (`Mdyn < 10^10.5`) and the other 32 have `Mdyn > 10^10.5`; about one fourth of early-type galaxies above `10^11.5` are fast | no value below 10^10 | [R1267] abstract [V]. The law against mass is a figure in the paper: placeholder logistic in `log M*`, zero at 1e10, 0.75 at 10^11.5 [P], R6 | spheroid axis ratios |
| Orientation | uniform on the sphere | | geometry | spheroid |
| Populations | old; gas fraction at the floor | | GEN.386 | profile |

Fast rotators are 86 percent of ATLAS3D early types, 20 percent of them are classified E, and 66 percent of all classified E are fast rotators [V, [R1267]]. That is why the elliptical family has a flat mode: a galaxy GAMA calls elliptical is more often flat than round.

**A misprint to watch.** [R4048: Table 1] prints the early-type fit `R50 = b M^a` with `a = 0.56` and `b = 3.47e-5`, which gives 34 kpc at 5e10, far from the paper's own Figure 11 and from [R2556, R2558]. With `b = 3.47e-6` the fit gives 3.4 kpc at 5e10, within 15 percent of Lange's high-mass fit (3.9 kpc), so the printed `b` is probably a typo. It is not used here; Lange et al. is the early-type size source. [C]

### 6.4 Irregular family (GEN.385)

| Parameter | Law | Range and scatter | Source, check | Feeds |
|---|---|---|---|---|
| Size | `log10(r_e / pc) = 0.667 + 0.296 log10(M*/Msun)` (the field fit, isolated dwarfs) | 5.5 < log M* < 8.5; scatter 0.215 dex; 277 pc at 1e6, 548 pc at 1e7, 1.08 kpc at 1e8, 2.1 kpc at 1e9 [C] | [R0655] table of fits [V]; the satellite fit is `1.071 + 0.247` with 0.181 dex. ELVES-Field says isolated dwarfs are about 20 percent smaller than satellites, which the field fit already shows (13 percent at 1e7) [V]. Above 1e8.5 the Shen late-type law continues it (1.8 kpc at 1e9 against 2.1) | body radius |
| Profile | Sérsic `n = 0.7` | `n <= 1` | [R0655]: early and late types below 1e8.5 have "essentially the same average index of about 0.7" [V] | body index |
| Intrinsic shape | `C/A = 0.35`, `B/A = 0.78` with the table's spreads 0.10 and 0.15 (late-type dwarfs: flattened, mildly triaxial) | | [R0655] shape table [V] | body axis ratios |
| Clumps (young share) | a clumped young population holding a share of the stellar mass | share [N] | GEN.386; no paper read gives the clump share (R8) | young population |
| Bulge, arms, bar, thick disk | none | | | |
| Metallicity | `<[Fe/H]> = -1.69 + 0.30 log10(M*/1e6)`, 0.17 dex rms (0.09 for irregulars) | | [R2368] [V] (GEN.368 row 10) | profile |

### 6.5 Dwarf spheroid family (GEN.385)

| Parameter | Law | Range and scatter | Source, check | Feeds |
|---|---|---|---|---|
| Size | `log10(r_e / pc) = 1.179 + 0.231 log10(M*/Msun)` (early-type satellites), times 0.8 for an isolated galaxy | 5.5 < log M* < 8.5; scatter 0.183 dex; 294 pc at 1e6, 500 pc at 1e7, 851 pc at 1e8 (isolated) [C] | [R0655] [V]; the factor 0.8 is ELVES-Field's "about 20 percent" smaller isolated dwarfs [V] | body radius |
| Profile | Sérsic `n = 0.7` below 1e8.5, rising above it | | [R0655] [V]; no law above 1e8.5 for dwarfs found [N], R9 | body index |
| Intrinsic shape | `C/A = 0.55`, `B/A = 0.85`, table spreads 0.03 and 0.12 (oblate) | | [R0655] shape table [V] | body axis ratios |
| Populations | old, metal-poor; same Kirby law | | GEN.386; [R2368] [V] | profile |

The ultra-faint end follows the floor the Galaxy dials lane set on 2026-10-21 (Segue 1 as the example, Boss 7:31 PM EDT: stellar mass 281 solar, scale length 7.1 pc), which belongs to the Scale dial, not to this table.

## 7. Links shared by every family

| Link | Law | How a family changes it | Source, check |
|---|---|---|---|
| Central black hole | `M_BH / 1e9 = 0.49 (M_bulge / 1e11)^1.16`, 0.29 dex | the spheroid mass is `B/T x M*` for a spiral or lenticular and `M*` for an elliptical. A disk-dominated galaxy with no classical bulge takes the Reines and Volonteri relation (GEN.368 row 7). Occupation fraction: 1 above 1e10, at least 0.5 from 1e9 to 1e10, and a draw from 0.1 to 0.5 below 1e9 [P below 1e9] | Kormendy and Ho 2013 [V]; the occupation fractions are a planetGen stand-in with no source. The hole stays exactly on the origin (Boss, 7:27 PM EDT on 2026-10-16) |
| Nuclear cluster | `log M_NSC = 0.48 log(M*/1e9) + 6.51`, 0.6 dex; nucleated share 80 percent or more near 1e9, falling both sides | none (GEN.368 row 8) | [R3273]; the Georgiev slope conflict is gap GX3 |
| Globular clusters | system mass `eta x Mh` (GEN.403) | **family factor 1 for every family.** [R1831] (422 galaxies: 248 ellipticals, 93 S0, 81 spirals and irregulars) find that in the mid range the numbers "fall along a well defined baseline value of `S_N` about 1 ... similar among all galaxy types", and that cluster formation "appears to be most nearly proportional to the galaxy halo mass" | [R1831] abstract [V]. The U-shape (more clusters per star in dwarfs and supergiants) follows from the halo-mass law and the stellar-to-halo relation, so no per-family multiplier is needed. GEN.387's family factor stays an input set to 1 (the input GEN.403 already takes), to change only if a test shows otherwise |
| Satellites | subhalo count from `Mh` (GEN.369) | none | [R0508] and the rest of GEN.368 row 10 |
| Metallicity and rates | mass-metallicity relations, star formation main sequence (GEN.405, GEN.386) | the young share is the family's (section 6) | GEN.368 rows 13 and 14 |

## 8. Mass limits and the 262 kpc reach

The limits are settings; these are the defaults. The grid reaches 262 kpc (decided 2026-10-15). The rule that sets a family's upper limit: **the radius holding 95 percent of the family's light, at the median size plus two sigma of its scatter, stays inside the reach.** That radius is `4.74 Rd` for an exponential disk and `R_e` times 4.5, 8.6, 14.4 and 22.3 for a Sérsic profile of index 2, 4, 6 and 8 (computed from the incomplete gamma function). [C]

| Family | Lower limit | Upper limit | Basis |
|---|---|---|---|
| spiral | 3e9 | 2e11 | Boss (9:40 PM EDT on 2026-10-16). The rule would allow far more (95 percent radius 34 kpc at plus two sigma at 2e11), so the Schechter cut-off, not the reach, ends it |
| lenticular | 1e9 | 2e11 | Kelvin's lower limit; the disk rule is far inside the reach (49 kpc at 2e11 at plus two sigma) |
| elliptical | 1e9 | 2e11 | the rule: with the high-mass size fit (9.5 kpc at 2e11) and `n = 6`, the 95 percent radius is 137 kpc at the median and 273 kpc at plus two sigma; the rule's limit is 1.9e11, so 2e11 sits 6 percent past it, inside the slack of the assumed index [C] |
| irregular | the Scale dial's floor | 3e9 | the spiral family's lower limit |
| dwarf spheroid | the Scale dial's floor | 3e9 | the same |

So the common top of 2e11 solar masses is not arbitrary: the elliptical family reaches the grid edge at about that mass. The dark halo is a profile with no cut-off, so its `r200` past the reach (GEN.368 section 5.13.1) is already handled. The Scale dial's mass range (`reports/galaxy-size-2026-10-21.md`) is applied on top of these defaults: it runs from 281 to 1.55e11 solar masses (Boss: the floor is half of Segue 1, the ceiling 1.5 times Andromeda), so the dial cuts the 2e11 top of the spiral, lenticular and elliptical ranges to 1.55e11 and the dwarf floors to 281. Below 1e10 solar masses of halo mass M200 (4.55e6 solar masses of stars) no halo law runs (Boss, 7:33 PM EDT on 2026-10-21), and the dials lane's hard limits apply to every family: scale height at least 4 pc, bulge radius at least 2 pc, outline never empty (`reports/galaxy-size-2026-10-21.md`). The Shape row of the Generate form is an enum whose only key today is `spiral` (flag `--shape`); the dials lane adds the five family keys exactly as section 3 spells them when this design is built (GEN.608, UX.134), and a shape pins the family only, never a size.

## 9. Reference galaxies for the plausibility tests (TEST.204)

These are the points each family's draw must be able to produce, never presets, and nothing else is compared. The values are read from the sources named when TEST.204 is built; this table fixes which galaxies and which quantities.

| Family | Reference galaxies | Quantities compared | Source to read the values from |
|---|---|---|---|
| spiral | the Milky Way; Andromeda | stellar mass, disk length, disk height, bulge-to-total, bar, arm count, halo mass | [R0422] (1602.07702); the Andromeda values in `reports/galaxy-size-sources-2026-10-21.md` [R0909, R0979, R1145, R4084] [V] |
| lenticular | the S0 galaxies of the ATLAS3D volume-limited sample | effective radius, bulge-to-total, bar | the ATLAS3D papers (per [R4654]) |
| elliptical | the round and flat early-type galaxies of the ATLAS3D sample; one giant elliptical with a measured black hole (M87) | effective radius, Sérsic index, axis ratios, black hole mass | the ATLAS3D papers; Kormendy and Ho 2013 table |
| irregular | the Large and Small Magellanic Clouds | stellar mass, size, metallicity | [R2972] (1204.1562) |
| dwarf spheroid | Fornax and Sculptor | stellar mass, size, metallicity | [R2972, R3400] (Local Volume Database, 2411.07424) |

Over many seeds each family's drawn parameters fall inside the ranges of sections 6 and 7 and show the stated slopes: mass-size, mass-metallicity, stellar-to-halo, black hole against bulge; a failure lists the parameter and the mass bin (like TEST.196). The family fractions of section 5 are tested against Kelvin's Table 1 within sampling error.

## 10. The columns DB.66 adds

`galaxy_shape` today carries the spiral's numbers (`disk_scale_length_pc`, `disk_scale_height_pc`, `bulge_scale_radius_pc`, `bulge_amplitude`, `arm_count`, `pitch_angle_rad`, `arm_amplitude`, `arm_level`, `nuclear_amplitude`, the spiral reference, `k_norm`). The bar is drawn in 8.1 by GEN.533 with no stored column (Boss, 12:56 PM EDT, card "Draw it in 8.1"). DB.66 adds, as one migration only:

| Column | Type | Meaning | Default for an existing spiral |
|---|---|---|---|
| `family` | VARCHAR(16) NOT NULL | one of `spiral`, `lenticular`, `elliptical`, `irregular`, `dwarf_spheroid` | `spiral` |
| `bar_present` | TINYINT NOT NULL | the bar draw of GEN.533 | the drawn value |
| `bar_axis_ratio_b`, `bar_axis_ratio_c`, `bar_angle_rad` | DOUBLE | the bar's axis ratios and angle | GEN.533's values |
| `thick_disk_amplitude` | DOUBLE | the thick disk's density amplitude | the structure draw's value |
| `spheroid_amplitude` | DOUBLE | the Sérsic spheroid's amplitude (0 removes it) | 0 |
| `spheroid_sersic_n`, `spheroid_effective_radius_pc` | DOUBLE | index and effective radius | 0 (unused) |
| `spheroid_axis_ratio_b`, `spheroid_axis_ratio_c` | DOUBLE | triaxial axes | 1, 1 |
| `spheroid_orientation_1_rad`, `spheroid_orientation_2_rad`, `spheroid_orientation_3_rad` | DOUBLE | the three orientation angles | 0 |
| `body_clump_share` | DOUBLE | the irregular family's young clump share | 0 |

The migration gives every existing galaxy the `spiral` family and these defaults, so no stored density changes. The component masses (`component_masses`, GEN.379) are calculated, not stored. A family that has no thick disk, bar or spheroid stores amplitude 0; an amplitude of 0 removes the component exactly.

As built (GEN.379): `galaxy/density.py` `density_components(shape)` returns the density model as a list of `DensityComponent`s in a fixed order (`thin_disk`, `thick_disk`, `bulge`, `halo`, `nuclear_cluster`, `nuclear_disc`), each with its population, amplitude, mass (the model's units, in closed form) and `mass_solar` (the share of `structure.stellar_mass_solar`; the six sum to it), and the scale lengths, heights, bar angle and halo slopes that shape it. `component_masses(shape)` returns the six `mass_solar` values by name, which is where the potential (GEN.338) reads its bulge, disk and halo masses from; the numerical integral the earlier `component_masses` was is now `integrated_masses`, kept for the math check and the tests. `component_values(position, shape)` and `component_arrays(x, y, z, shape)` give each component's value at a point or at arrays of points, and `relative_density`, `population_densities`, `population_densities_array` and `total_density_array` sum and split those, so there is one place that knows each component's formula; a component with amplitude 0 is exactly 0. Today's spiral gives today's density to the last digit (`src/tests/fixtures/density_identity.json` holds the density, the population split and the array forms of six shapes at 24 points, made before the rewrite) and the array forms run as fast as before. Closed-form masses agree with a quadrature of each profile to 1e-6 (the nuclear cluster's profile is held flat inside its softening radius, which keeps `softened_mass_fraction` of its nominal mass, 78 percent for the reference cluster: about one part in ten thousand of the galaxy's stars). `get_galaxy_shape` reads the DB.66 columns into `GalaxyShape` (`family`, `bar_present`, `bar_axis_ratio_b`, `bar_axis_ratio_c`, `bar_angle_rad`, `thick_disk_amplitude`, the spheroid numbers and `body_clump_share`) and `save_galaxy_shape` writes them. A NULL bar or thick disk number stays `None` on the shape and the structure draw's value stands; a number pins it, in the shape and in its `structure`. The Sersic spheroid (GEN.381) adds its own component to the list; no family draws one yet, so a stored spheroid amplitude above 0 is refused with a message until GEN.381 builds the component.

## 11. Research requests

Relations no paper read here gives. Each stays a labelled placeholder until answered; none blocks the build.

| Id | Request | Used by |
|---|---|---|
| R1 | The split of Kelvin's S0-Sa class into S0 and Sa as a number fraction against stellar mass (a classification that separates S0 from Sa) | section 5 |
| R2 | The quenched fraction of isolated dwarfs against stellar mass from a sample larger than a few dozen | section 5 |
| R3 | S0 disk height against length, and the bar fraction of field (not cluster) S0 galaxies | section 6.2 |
| R4 | The early-type size-mass relation above 1e10 as one continuous fit with its scatter against mass | section 6.3 |
| R5 | The Sérsic index against stellar mass (not luminosity) for ellipticals with its scatter, or a stated early-type stellar mass-to-light ratio in the B band | section 6.3 |
| R6 | The slow-rotator fraction against stellar mass (the ATLAS3D paper III figure data) | section 6.3 |
| R7 | The stellar-to-halo relation for early types and late types separately (Moster fits all galaxies together) | section 4 |
| R8 | The clump share of the stellar mass in dwarf irregular galaxies | section 6.4 |
| R9 | The Sérsic index of dwarfs above 1e8.5 | section 6.5 |
| R10 | A source for the stellar mass function below 1e8 (Baldry's fit is valid from 1e8) | the Scale dial's dwarf range |

## 12. Questions for Boss (defaults stand if unanswered)

- **Q1.** Decided (Boss, card tapped 8:19 PM EDT on 2026-10-21): the S0-Sa class splits 0.55 lenticular and 0.45 Sa spiral, constant in mass.
- **Q2.** The dwarf spheroid share of galaxies under 1e9: use the isolated-dwarf step (30 percent under 1e7, none from 1e7 to 1e9)? Decided: yes (Boss, card tapped 8:10 PM EDT on 2026-10-21).
- **Q3.** Is the common upper limit of 2e11 solar masses for every family right, given that the elliptical rule gives 1.9e11? Decided: yes (Boss, card tapped 8:11 PM EDT on 2026-10-21).

## 13. Defaults standing as decided

- The family is drawn from the seed; the Generate page may pin it ("auto", spiral, barred spiral, lenticular, elliptical, irregular, dwarf spheroid). Pinning fixes the family, never a size or a count: the laws still draw those from the seed.
- Barred or unbarred is a draw inside the spiral and lenticular families, not a family of its own.
- A satellite galaxy is made of the dwarf families (GEN.369 already draws profile by mass); a Plummer dwarf spheroid and an exponential dwarf irregular are the same objects the satellite laws draw.
- The root draw is the stellar mass; the halo mass follows from it by the stellar-to-halo relation, and GEN.368's table is read from that root.
- Each family has a simple hand-drawn static SVG outline in the icon library (UX.135, Boss, 2026-10-16, 9:00 PM EDT). Static, kept in the repo, no generation step, works on the offline page.
- The galaxy stores no home or Sun position. USR.11, a user's own home field, is a different thing and stays.
- Interacting and merging galaxies are not part of it.

## 14. The chain, counted back from 9.4

Each row is an item in TODO.md; a row sits in the earliest version that can hold it and be useful on its own. The build order rule applies: a need is earlier in its thread, in an earlier phase or in an earlier version.

| Step | Item | Version | Why there |
|---|---|---|---|
| 1 | GEN.378 the checked design (this document) | 8.1 phase 5 | the shape-draw build that follows reads it, so it is not built twice |
| 2 | DB.66 the shape row carries the family key and the component columns of section 10 (migration only) | 8.3 phase 1 | one migration beside the other galactic-gravity work; Core 1 owns the numbering |
| 3 | GEN.379 the density model as a list of components, each with an amplitude that may be 0, and `component_masses` | 8.3 phase 1 | the potential (GEN.338) reads each component's mass from it, so it comes first; today's spiral reproduces today's density |
| 4 | GEN.380 the type draw (section 5); GEN.381 the Sérsic spheroid component; UX.135 one hand-drawn outline per family | 9.4 phase 1 | nothing calls them before the families exist |
| 5 | GEN.382 spiral (section 6.1), GEN.383 elliptical (6.3), GEN.385 irregular and dwarf spheroid (6.4, 6.5) | 9.4 phase 2 | three independent threads |
| 6 | GEN.384 lenticular (6.2); GEN.386 populations; GEN.387 clusters, nuclear cluster and black hole (section 7) | 9.4 phase 3 | they read the families |
| 7 | GEN.389 creation draws the family and stores it | 9.4 phase 4 | needs every family |
| 8 | UX.134 the Generate page, MAP.242 the Galaxy Map, TEST.204 plausibility tests (section 9), DOC.59 docs | 9.4 phase 5 | four independent threads |

GEN.338 gains one need (GEN.379) and one test: a component of zero mass is legal, so an elliptical (no disk) needs no change in the potential.
