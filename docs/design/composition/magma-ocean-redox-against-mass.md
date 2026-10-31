# Magma-ocean redox against mass and the later mantle (GEN.482, AR3)

Status: research note, written on the night of 2026-10-24 (as the lane reported) by the Citations lane. Rows a1, a2 and a3 of the atmosphere possibility table (`scripts/atmosphere_recipes_data.py`, part 09 of `docs/design/composition/atmosphere-composition-parts/`) are checked against first-hand sources. GEN.482 (the check) changed nothing the generator draws; GEN.482b (section 7, done the same night) wired the findings; no epoch, because nothing in the generator draws the mixture yet (the golden battery is unchanged).

## 1. The answer in short

- **a1, the magma-ocean surface redox against mass.** The present mean, `clamp(-2 + 1.25 (log10 M + 2), -2, 0.5)`, is inside the published ranges for all three Solar System anchors (table 3.1: the Moon, Mars and Earth), within 0.2 log units of the middle of each. The row is therefore a checked fit to three sourced points, not a free guess. What the sources do not support: the width (0.5 flat is too narrow for Earth, where the published span is 3.0 log units wide, so a standard deviation of about 0.9), and the flat clamp at +0.5 above 1 Earth mass (the physics keeps rising up to the pressure where the calculations stop, about 140 GPa at the magma-ocean base). Both are now wired (section 7).
- **a2, the later upper mantle.** `3.5 + 2 log10 M` with width 1.0 passes through the three sourced anchors within one standard deviation (Earth IW+3.5, Mars about IW to IW+1, Moon IW to IW-2). It is a three-point fit and nothing more: no source gives the law above 1 Earth mass, and at 10 Earth masses the code gives +5.5, which nothing supports.
- **a1 and a2 are drawn independently, but the sources say they are linked**: the later mantle is never reduced below its own magma-ocean surface and the gap grows with mass (section 5). GEN.482b takes the larger of the two.
- **a3, small-body building blocks.** The check bracket in the row (`-1.4 to 0`) only covers the reduced half of Schaefer and Fegley's Table 2. The same table puts carbonaceous material (CI, CV) at IW+1.6 to IW+2.6 at 1 bar and 1500 K. Row a3's draw (normal, centre -1.25, width 0.22, from white-dwarf debris) is for differentiated rocky debris and does not describe C-complex bodies (section 6).
- **Na and K condensation temperatures** (the primary source the gap note asked for): Lodders 2003 Table 8, 50 percent condensation at 10^-4 bar: Na 958 K, K 1006 K, both in feldspar (section 8).
- **Not found:** the degassed 40Ar fraction against tectonic regime (no source with a number other than Earth's about 43 percent in part 07), Hirschmann 2012 (paywalled, Boss cannot fetch it), Gaillard and Scaillet 2014 coefficients (the earlier paper note says equations without coefficients).
- On the row states: a1 and a2 move from DESIGN to DERIVED in the table (state "partly verified", with the limits written in the row). a3 stays as it is (section 6). This does not change a number the generator draws.

## 2. Sources and what was read

All opened and read in full text for the numbers taken.

| Source | Access | Used for |
|---|---|---|
| Hirschmann 2022, Geochimica et Cosmochimica Acta 328, 221-241, "Magma oceans, iron and chromium redox, and the origin of comparatively oxidized planetary mantles" | public copy on the NSF Public Access Repository (par.nsf.gov, purl 10335656) | Table 1 (core-forming conditions), the Fe3+/FeT ranges of Earth and Mars, section 6.3 (surface redox of the magma ocean, Eq. 21 and the corrected Sossi relation Eq. 34), the redox of the present mantles |
| Deng, Du, Karki, Ghosh and Lee 2020, Nature Communications 11, 2007 | PMC7181735 | the pressure dependence (the molar volume difference stays positive to 140 GPa), the base conditions and the surface redox of the Moon, Mars and Earth |
| Sossi, Burnham, Badro, Lanzirotti, Science Advances 6, eabd1387 | uploaded by Boss | the 0.037 Fe3+/FeT anchor and the Eq. 33 relation (read earlier, GEN.481 and the citation pass) |
| Hirschmann 2021, GCA 313, 74 | uploaded by Boss | the iron-wuestite scale the repo already uses |
| Gaillard, Bernadou, Roskosz, Bouhifd, Marrocchi, Iacono-Marziano, Moreira, Earth and Planetary Science Letters 577, 117255 | HAL copy insu-03383209 | the three published views of the Earth's surface magma ocean (section 3.2) |
| Schaefer and Fegley 2017, The Astrophysical Journal 843, 120 | public copy on par.nsf.gov (purl 10038075) | Table 2 (oxygen fugacity fits of meteoritic material), the text on the end members |
| Lodders 2003, The Astrophysical Journal 591, 1220 | uploaded by Boss | Table 8 (condensation temperatures) |

One reading note. The PDF text of Hirschmann 2022 loses its minus signs (for example "IW1.2" for IW-1.2). The signs were restored by reading the page image of section 6.3 (page 235) and, for the Earth values of Eq. 34, by arithmetic with Eq. 34 itself: Fe3+/FeT 0.034 and 0.103 give delta IW -0.2 and +1.8, which are the printed values.

## 3. The surface redox of the magma ocean against mass (row a1)

### 3.1 The anchors

delta IW is log fO2 minus the iron-wuestite buffer (the repo uses the Hirschmann 2021 buffer, the scale of Hirschmann 2022).

| Body | Base of the ocean (P, delta IW) | Surface of the ocean (this note's reading) | Source | Row a1 now |
|---|---|---|---|---|
| Moon (0.0123 Earth masses) | 3.5 to 4.8 GPa; core-forming delta IW -2.2 to -1.5 (Hirschmann 2022 Table 1 and text); Deng base 5 GPa, -2 | about -2.2 to -1.2: little Fe3+ is stabilised at these pressures, so the surface stays near the core-forming value (Deng: lower end of the lunar basalts, IW-2; Hirschmann Table 1 and text) | Hirschmann 2022; Deng 2020 | -1.9 |
| Mars (0.107) | 14 GPa; core-forming -1.5 to -1.1 | Eq. 21: -1.4 to -0.7; corrected Sossi relation Eq. 34: -0.4 to +0.2 (Fe3+/FeT 0.026 to 0.038, section 6.3 and Fig. 9 inset); Deng: about IW | Hirschmann 2022; Deng 2020 | -0.7 |
| Earth (1) | 28 to 53 GPa, 2250 to 4100 K; core-forming -2.5 to -1.4 (hotter is more reduced) | Eq. 21: -1.2 to +1.4; Eq. 34: -0.2 to +1.8 (Fe3+/FeT 0.034 to 0.103); Sossi 2020: about +0.5 on its own scale, 0 on the corrected scale (Hirschmann 2022 section 6.3); Deng: at the lower bound of the present mantle (IW+3.5 plus or minus 1.0), so about +2.5 | Hirschmann 2022; Sossi 2020; Deng 2020 | +0.5 |

The widest honest span for Earth is the union of the first-hand model results, -1.2 to +1.8, middle +0.3. Deng's +2.5 comes from a cold temperature profile and from the older equation of state treatment; Hirschmann 2022 uses Deng's equation of state with a temperature function and the chromium correction, and is the newer work. By the project's rule (two verified sources in doubt: the newer one) the Hirschmann range is the centre and Deng's value is the upper tail, kept in the width. Hirschmann himself states that the uncertainties "remain too great to establish with certainty the redox state" of the atmospheres above Earth's and Mars's magma oceans.

### 3.2 The published views of Earth

Gaillard et al. 2022 lists three, and the table above sits between the first two. (i) One model: oxidised surface close to the FMQ buffer (Fe3+/Fe-total 0.03 to 0.1). (ii) Sossi et al. 2020: IW+0.5 plus or minus 0.7. (iii) Keppler and Golabek 2019: graphite flotation holding the shallow ocean down to IW-5 at low pressure. View (iii) is a different mechanism, not a spread of the same one: it is the only way the sources get far below -1.2. It stays out of the draw and is flagged as an alternative (limit of the row).

### 3.3 The check of the row

| Body | Published middle | Published half-range | Row a1 mean | Difference |
|---|---|---|---|---|
| Moon | -1.7 | 0.5 | -1.9 | -0.2 |
| Mars | -0.6 | 0.8 | -0.7 | -0.1 |
| Earth | +0.3 | 1.5 | +0.5 | +0.2 |

Every difference is inside half the published half-range, so the mean is consistent with all three. A least-squares line through the three middles has a slope of 1.0 per dex of mass against the row's 1.25; the data do not tell them apart (the Earth middle depends on an assumed core-forming temperature, 2300 to 4100 K, which Hirschmann names as the main control).

### 3.4 What the sources do not give

- **The width.** The row uses 0.5 (Deng's one-sigma). The published half-ranges are 0.5 (Moon), 0.8 (Mars) and 1.5 (Earth), which as uniform distributions are standard deviations of 0.29, 0.46 and 0.87. A width that grows with mass from 0.5 at 0.1 Earth masses to 0.9 at 1 Earth mass fits the three: `0.5 + 0.4 clamp(log10 M + 1, 0, 1)` (an L2 tweak: tuned to the published ranges, not to a statistic). The floor 0.22 the row carries from the Fe3+/FeT measurement stays.
- **Above 1 Earth mass.** The row clamps at +0.5. The physics (Deng: the volume difference of FeO1.5 and FeO stays positive to 140 GPa, so the surface of a deeper ocean stays more oxidised) says the surface keeps rising with the pressure at the ocean base; none of the sources gives a super-Earth number, and both Hirschmann and Deng stop at Earth. The only first-hand slope is the Mars-to-Earth one, 0.9 to 1.25 per dex. A pressure-only slope from the three Deng points is about 2 per dex. The two bracket the extension, 0 (the clamp) to 2 per dex, centre about 1. Where the base pressure passes 140 GPa the calculation does not apply; that is at about 5 Earth masses if the base pressure scales as M to the 0.55 (an assumption of this note, not a number read from a source). This is the part that stays a stand-in.
- **The graphite tail** (section 3.2) and the Cr and Ti effects on the crystallised mantle (Hirschmann 2022 section 6.4 and 6.5) are not mass laws and stay out.

## 4. The row states

- a1: state "partly verified"; the part 09 status moves from DESIGN to DERIVED. Limit: three Solar System anchors, the width and the extension above 1 Earth mass are not sourced, the graphite-flotation alternative is not drawn.
- a2: same, see section 5.

## 5. The later mantle (row a2)

Anchors (all first-hand):

| Body | Later upper mantle | Source | Row a2 now |
|---|---|---|---|
| Earth | IW+3.5, constant within about 1 log unit since the early Archean; above IW+3 in spinel peridotite (Hirschmann 2022 section 4) | Deng 2020 (refs 37, 6, 7); Hirschmann 2022 | +3.5 |
| Mars | the most reduced basalts record IW (Schmidt 2013; Nicklas 2021), many up to IW+5; the primitive upper mantle IW to IW+1 | Hirschmann 2022 introduction and section 4; Deng 2020 | +1.6 |
| Moon | lunar basalts IW to IW-2 | Deng 2020 (refs 8 to 10) | -0.3 |

With the width 1.0, Mars is +0.6 above the middle of its range and the Moon +0.7 above its middle; both are inside one standard deviation. The law is a straight-line fit to three points and says so.

Two structural facts from the sources that the row does not use:

1. Deng and Hirschmann find that for Mars and the Moon the later mantle is close to the magma-ocean surface (no large later oxidation), and for Earth later processes add at most 1 to 3 log units (bridgmanite disproportionation, late accretion, Hirschmann section 6.4). The gap (a2 minus a1) in the sources is about 0 to 2 at the Moon, 0 to 1 at Mars and 1 to 3 at Earth. In the code it is 1.6, 2.3 and 3.0 as means, and the two rows are drawn independently, so a draw can put the later mantle below its own magma-ocean surface.
2. Above 1 Earth mass the row is an unsupported extrapolation (5.5 at 10 Earth masses, clipped at 6).

## 6. Small-body building blocks (row a3)

The row, as it ships, is a normal distribution centred at -1.25 with width 0.22, from 22 polluted white dwarf values (Doyle 2020 and 2023, already in the repo), with the Schaefer and Fegley bracket `-1.4 to 0` as a check. Schaefer and Fegley 2017 Table 2 gives the fits of log10 fO2 against temperature for outgassing at 1 bar. Against the Frost 1991 iron-wuestite line (constants from PetroBuffer, `-27489/T + 6.702`; the Hirschmann 2021 scale the repo stores was not recomputed here and differs by a few tenths of a log unit), the polynomials give delta IW:

| T (K) | CI | CV | H | EH | eucrite |
|---|---|---|---|---|---|
| 1000 | +2.84 | +2.12 | -0.90 | -1.43 | -1.38 |
| 1200 | +2.48 | +1.85 | -0.67 | -1.19 | -1.11 |
| 1500 | +1.98 | +1.65 | -0.53 | -1.03 | -0.91 |
| 2000 | +1.22 | +1.57 | -0.51 | -0.98 | -0.80 |

The authors' text agrees: the most oxidised mixture is CI at "about QFM-1" (about IW+2.6 at 1500 K by the Frost 1991 buffers; the polynomials above give +2.0, so the text is read off a figure) and the most reduced is EH at "about IW-1", with all mixtures between the two. The unsupported part of row a3 is therefore not its reduced end (which the white-dwarf data and Schaefer and Fegley both give, -1.0 to -0.5 at 1500 K) but its use for every small body: carbonaceous (CI, CV, C-complex) material outgasses at IW+1.6 to +2.6. The two scales also differ (white-dwarf values are core-formation redox of differentiated parents; Schaefer and Fegley are outgassing of chondritic material), so no single number covers both. The row is not changed: nothing reads a3 for a body with an air today (section 7, point 5).

## 7. GEN.482b: what was wired (no epoch: the generator does not draw the mixture yet)

1. a1 width: `0.5 + 0.4 clamp(log10 M + 1, 0, 1)` (0.5 at 0.1 Earth masses and below, 0.9 at 1 and above).
2. a1 mean above 1 Earth mass: plus `clamp(log10 M, 0, 0.7324)`, that is 1.0 log unit per dex of mass up to 5.4 Earth masses (+1.23) and held there, where the base pressure passes the 140 GPa of Deng's calculations (assuming the base pressure scales as M to the 0.55). The bracket 0 to 2 per dex stays the stated uncertainty; this is the part that stays a stand-in.
3. a2 mean: `3.5 + 2 clamp(log10 M, -3, 0)`: the Earth value is held above 1 Earth mass (no source for a super-Earth's later mantle), and the 6 clip of the row is no longer reached.
4. `draw_redox` takes the larger of a2 and a1 for the later mantle, so a body never moves to a later mantle more reduced than its own magma-ocean surface (Deng 2020 and Hirschmann 2022: the later processes only oxidise).
5. a3: not changed. No code reads a3 for a body with an air today (an asteroid is airless and takes no mixture; the rocky path takes a3 only when a1 is absent), so a split by building-block type has no consumer. It waits for one (the first body that outgasses from a carbonaceous parent), with the carbonaceous range IW+1.6 to +2.6 of section 6 as its source.

## 8. Na and K condensation temperatures

Lodders 2003 Table 8 (equilibrium condensation temperatures for a solar-system composition gas at 10^-4 bar total pressure; 50 percent condensation): Na 958 K (host: feldspar, the initial phase albite NaAlSi3O8 dissolving), K 1006 K (feldspar, the initial phase KAlSi3O8 dissolving). For reference in the same table: Fe 1334 K (alloy), S 664 K (troilite), Zn 726 K, Mn 1158 K, Mg 1336 K, Si 1310 K. This answers the gap note's request for a primary source and ties "K depletion follows the nebular temperature at the formation distance" to numbers; it does not by itself give the planet's K abundance (that needs the temperature at the formation distance, the disc model's, not read here).

## 9. Limits

- Hirschmann 2022 states that the model is an extrapolation in temperature, pressure and melt composition (section 6.2) and that the two Fe3+ models differ by 0.7 to 1.3 log units in delta IW at the same Fe3+/FeT; the ranges above carry both.
- The pressure-scaling exponent 0.55 and the "about 5 Earth masses" are this note's assumptions.
- Gaillard and Scaillet 2014 and Hirschmann 2012 were not used: the first was read earlier (equations without coefficients), the second is paywalled and Boss cannot reach it.
- The 40Ar degassed fraction against tectonic regime was not found in any open source read; part 07's Earth value of about 43 percent stays the only anchor.

## 10. Licence

Hirschmann 2022 and Schaefer and Fegley 2017 are copyright of their publishers (Elsevier; the American Astronomical Society), read from the public access copies; only formulas, numbers and short quotations of fact were taken. Deng et al. 2020 is CC BY 4.0. Gaillard et al. 2022 is an open archive copy of an Elsevier paper; only numbers were taken. Lodders 2003 is AAS copyright; one table of values. See `DATA-LICENSE.md`.
