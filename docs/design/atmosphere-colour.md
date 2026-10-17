# Atmosphere colour: one consistent atmosphere, then its colour

**Status (2026-10-16):** plan; nothing is built. Boss (11:24 PM EDT, 11:25 PM EDT) made this the seminal feature of 9.0: the colour of the air of any world, worked out from the atmosphere it really has. GEN.408 checked the sources (done) and GEN.409 turns this note into the checked design. Boss (11:58 PM EDT) made composition the joint headline of 9.0: the master note is [composition-model.md](composition-model.md), and the research, design and table items moved down to 8.1. The base is [sky-colour-analysis.md](sky-colour-analysis.md).

## 1. What Boss decided

- Store the albedo and the greenhouse values used in generation; the migration works the formula back from the stored temperature and pressure, or draws from the planet class where that fails.
- Break the atmosphere's molar density into per-gas partial pressures and calculate the molar density from them, so both match and the light meets a known mixture.
- Build ozone into the system from real research; if the research is missing, build a generator from the oxygen and the star's UV; Boss researches on request.
- Sky colour and atmosphere colour are the same; the wording is "atmosphere colour"; it replaces the planet colour where that is what it meant.
- The colour is absolute, not adapted to the star's light (card, earlier that evening), and gas giants are in the first version.
- Haze, dust and similar scatterers get reasonable draws or derivations.
- One table of every gas and compound that can occur, including every compound the descriptions mention, so that atmosphere text, composition, density and colour agree: "all of the values consistent across the board".
- Nothing is fixed per class and nothing is copied from the Solar System; every input is a stored physical quantity or a seeded draw.

## 2. The chain

| Where | Item | What |
|---|---|---|
| 8.4 phase 1 | GEN.408 | Verify the sources and say what is missing (done 2026-10-17, [atmosphere-sources-check.md](atmosphere-sources-check.md)) |
| 8.1 phase 1 | GEN.409 | The checked design (this note) |
| 8.1 phase 2 | GEN.410 | The gas and compound table |
| 8.4 phase 2 | GEN.411 | Molar density from the partial pressures |
| 8.4 phase 1 | GEN.412 | Albedo and greenhouse multiplier recorded, with an inverse |
| 8.4 phase 2 | GEN.422 | Stop the recalculation redraw; order the steps so pressure, mixture and temperature agree after the freeze-out |
| 8.4 phase 2 | GEN.423 | Atmosphere density, scale height and the greenhouse multiplier follow from the mixture |
| 8.4 phase 2 | GEN.424 | The balance gas is a named gas with its own molar mass; habitability reads the stored mean molar mass |
| 8.4 phase 2 | DB.67 | Columns and migration for albedo, greenhouse multiplier, molar density |
| 8.4 phase 3 | GEN.413 | Ozone from the oxygen, the star's UV and the Chapman balance |
| 8.4 phase 3 | GEN.414 | Haze, dust and cloud optical depths |
| 8.4 phase 4 | DB.68 | Columns for ozone and the aerosol depths |
| 8.4 phase 4 | GEN.415 | Description text written from the stored mixture |
| 8.4 phase 4 | TEST.209 | One consistency sweep over every class |
| 8.4 phase 4 | DOC.63 | The data model docs |
| 9.0 phase 1 | GEN.416 | Blackbody colour of a star (first part of the photometry library) |
| 9.0 phase 1 | GEN.417 | Rayleigh optical depth of the mixture |
| 9.0 phase 1 | GEN.418 | Atmosphere colour of a solid world |
| 9.0 phase 1 | GEN.419 | Atmosphere colour of gas and ice giants at a reference level |
| 9.0 phase 2 | GEN.420 | The description row and swatch in every output |
| 9.0 phase 2 | MAP.244 | The maps use the atmosphere colour in place of the planet colour |
| 9.0 phase 2 | DOC.64 | Docs of the colour |
| 9.0 phase 2 | GEN.421 | The headline item |

The foundations sit in 8.4 (five phases, each blocking the next; the bulk composition chain of section 5 comes first) because the planet classes are reworked there (GEN.90) and the new classes should be born on the one table and the one mixture; the migrations take the revisions after the 0102 to 0105 stack of 8.1. VIEW.11 extends GEN.416 and stays in 9.1 with the view from a planet.

The read-only trace of 2026-10-17 (`/mnt/project-files/reports/composition-density-link-2026-10-17.md`) found that only the pressure flows from the stored mixture's side: the drawn molar density, atmosphere density, scale height and greenhouse multiplier never see the mix (molar density off by more than 10 percent for classes A, B, E, G, H, L and Q), the recalculation after orbit validation redraws the albedo, greenhouse multiplier, redox and mix, the freeze-out lowers the pressure after the temperature is set, and habitability uses 28 g/mol for every balance gas. GEN.411, GEN.422, GEN.423 and GEN.424 carry those fixes.

## 3. Open points for the design

1. **Giants:** the reference level for a gas or ice giant (cloud top or a fixed bar) and how the cloud depth sets it.
2. **Thick air:** the multiple-scattering step above tau of about 0.5 (a two-stream pass or the saturated limit).
3. **Inverse of the temperature formula:** whether the albedo and the greenhouse multiplier can both be worked back from the stored temperature and pressure (the temperature is one relation between them; whether the pressure fixes the multiplier is for GEN.412 to check in the code). Where only one can be solved, the migration draws the other from the class range, seeded.
4. **Ozone:** Earth's column near 300 Dobson units is the test target, not a preset; the oxygen and UV set it.
5. **Aerosols:** dust has no data today; its rule is a flagged parametrisation with a physical range.
6. **Stars:** the blackbody is the light; M dwarfs and cooler get a model-atmosphere table later.

## 4. Research: what the sources check found

The study of 2026-10-17 (GEN.408, [atmosphere-sources-check.md](atmosphere-sources-check.md)) went through 19 sources: 3 read in full (Betremieux and Kaltenegger, which is a 2015 paper, Sudarsky, Burrows and Pinto 2000, Meadows et al. 2018), 10 verified second-hand, 4 not reached and 2 citation only. Only arXiv papers (through `export.arxiv.org`), open-source code and PyPI package data could be read; journal sites were refused.

What it changed in the plan:

- **Rayleigh table (GEN.410, GEN.417):** each gas stores the temperature its refractivity formula was fitted at (288.15 K for N2, Ar, CO2, CO and CH4; 273.15 K for O2). The ratios to N2 are CO2 2.5, H2 0.20, He 0.013. The coded CO formula is half the measured value and is replaced. H2S and SO2 are not found. Earth air gives tau 0.0969 at 550 nm (computed, 2 percent tolerance).
- **Ozone (GEN.413):** use the Kozakis 2022 polynomial fit of the column against the oxygen by host star, not a power law; seeded spread about 30 percent; ozone is a Chappuis absorber along the sun path and matters at twilight, not at noon; the Earth figure is a test range of 279 to 300 Dobson units.
- **Aerosols (GEN.414):** dust is grey, haze goes as lambda^-1 to lambda^-2.3, cloud depth follows Venus (20 to 30) and Earth (about 10).
- **Giants (GEN.419):** Sudarsky classes by T_eff (III above about 350 K, IV 900 to 1500 K), colour from the Cahoy, Marley and Fortney table.
- **Meadows 2018** is the colour of the disk from space, a precedent for star plus Rayleigh only.

**Open, Boss's to research** (ranked by how much the chain depends on it):

1. Single scattering albedo of Mars dust (Lemmon 2004, Tomasko 1999) and the refractive index of iron-oxide dust.
2. Titan haze absolute optical depth and spectrum (Tomasko 2008, Doose 2016).
3. Methane, ammonia and water absorption strengths in the visible for giants (Karkoschka 1994; Karkoschka and Tomasko 2009, 2010).
4. Ozone cross sections with a clear licence (Serdyuchenko 2014 at IUP Bremen, Burkholder and Talukdar 1994 at MPI Mainz).
5. Refractivity and King factors of H2S, SO2 (Bogaard 1978), NH3 and CO.
6. Condensation curves of NH3, H2O and CH4 (Lodders and Fegley 2002, Lodders 2003 and 2010, Weidenschilling and Lewis 1973).
7. Bucholtz 1995 Table 3, and the originals of Peck and Huang 1977, Kasting and Donahue 1980, Segura 2003, Hulburt 1953, an Angstrom-exponent review.
8. Ice sublimation curves for N2, CO, CH4 and CO2, which the vapour-pressure relations reached by the mixing-ratio research (part 5) do not cover.
9. Per-comet abundance tables (volatile ratios by comet, not the few well-known ones).
10. Magma-ocean redox papers: the mantle redox of a young rocky planet and how it drifts (the code's mean is 3.5 above the iron-wuestite buffer, the building blocks sit near 0 to 1 below; GEN.451).
11. Widening the proxy (Project settings, environment network): doi.org, ADS, Science, Copernicus (AMT, ACP), Optica, ScienceDirect, NTRS and the Mainz and Bremen data hosts.

Until an item is researched, its design entry names the seeded range used meanwhile (rule: nothing fixed per class). Boss widened the network at about 12:02 AM EDT on 2026-10-17: arxiv.org itself is reachable now and the Atmosphere research lane (GEN.431) re-reads the four sources the check could not reach.

## 5. Bulk composition: the cause of mass, density and gravity

Boss (2026-10-16, 11:49 PM EDT) asked whether the system stores a body's composition and whether it drives mass, density and gravity. It does not: `composition` is a fixed text per class (`physics/planets.py` 481), the bulk density is a uniform draw per class (500) and the mass and gravity follow from it (see `/mnt/project-files/reports/composition-density-link-2026-10-17.md`). The chain below makes the composition the cause, and the atmosphere reads the gravity it derives (GEN.411 needs GEN.429).

| Where | Item | What |
|---|---|---|
| done | GEN.425 | Verify the mass-radius and equation-of-state sources ([bulk-composition-sources-check.md](bulk-composition-sources-check.md)) |
| 8.1 phase 1 | GEN.426 | The checked design (docs/design/bulk-composition.md) |
| 8.1 phase 2 | GEN.427 | The bulk materials table with equation-of-state constants |
| 8.4 phase 1 | GEN.428 | Composition fractions drawn, replacing the fixed text (generator epoch) |
| 8.4 phase 1 | GEN.429 | Density, mass, radius and gravity derived from the fractions |
| 8.4 phase 1 | DB.69 | Columns and migration for the fractions |
| 8.4 phase 1 | GEN.430 | Descriptions state the stored composition |
| 8.4 phase 4 | TEST.210 | Bulk composition sweep over every class |
| 8.4 phase 4 | DOC.65 | Docs |

Research needed (GEN.425 checked them, done 2026-10-17, see [bulk-composition-sources-check.md](bulk-composition-sources-check.md); only the sources it cannot reach may need Boss): mass-radius relations for solid planets (Seager 2007, Zeng and Sasselov 2013, Zeng, Sasselov and Jacobsen 2016, Zeng 2019), envelopes (Lopez and Fortney 2014, Fortney, Marley and Barnes 2007), equations of state (Vinet or Birch-Murnaghan constants, Saumon, Chabrier and van Horn 1995 for hydrogen and helium), the link from the star's Fe, Mg and Si to the core fraction (Dorn 2015, Thiabaud 2015, Santos 2015) and the water fraction against the snow line. All are recalled, not yet verified.
