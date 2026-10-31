# Body composition: research findings and the composition design (draft for the 9.0 headline)

**Status (2026-10-17):** research result of GEN.432 and GEN.433 (the Body composition research lane), merged here by the TODO thread so that every lane has it in the repository. It feeds the designs GEN.426 (bulk), GEN.435 (stars), GEN.436 (radiation) and GEN.457 (small bodies and moons), the tables GEN.427, GEN.438 and GEN.439, and the generation chain of [composition-model.md](composition-model.md). Times in it are as the lane's clock reported. Section 16 says where the TODO thread placed each finding and cluster. The five decisions of section 5 are taken (Boss's card taps of 1:21 to 1:24 AM EDT, as the coordinator relayed them).

Draft of 2026-10-17 (written between 12:00 AM and 1:10 AM EDT as the lane clock reported; lane clocks drift). Research only; the TODO thread merged it into the repository. Requested by Boss (2026-10-16, 11:58 PM EDT): composition of planets, moons, stars and everything else becomes a joint headline of the atmosphere-colour release; every body gets composition figures based on actual science, not made up where it can be avoided; drawing density and other bulk figures first and filling in a fitting composition is fine, and so is building up from the most basic.

This is the **synthesis**. Five detailed parts (A to E) and a sixth check part (F) carry the tables, formulas, verified code and per-number labels; read this file for decisions and structure, the parts for numbers.

| Part | File | Covers |
|---|---|---|
| A | [A-stars.md](body-composition-parts/A-stars.md) | Stars of every kind and phase, abundance sets, [X/Fe] patterns, X/Y/Z, white dwarfs, neutron stars, black holes, brown dwarfs |
| B | [B-rocky.md](body-composition-parts/B-rocky.md) | Rocky and iron planets: EOS constants, mass-radius, core mass fraction from the star, mantle mineralogy, verified solver |
| C | [C-volatile-and-giants.md](body-composition-parts/C-volatile-and-giants.md) | Water and ice worlds, ice phases, sub-Neptunes, ice and gas giants, brown dwarfs, rogue planets |
| D | [D-small-bodies-moons.md](body-composition-parts/D-small-bodies-moons.md) | Asteroids, comets, Kuiper belt objects, dwarf planets, moons of every origin |
| E | [E-ground-radiogenic-shielding.md](body-composition-parts/E-ground-radiogenic-shielding.md) | Crust and ground, regolith, albedo by band, radiogenic heat, surface radiation and shielding, one materials table (CSV) |
| F | [F-boss-note-check.md](body-composition-parts/F-boss-note-check.md) | Check of Boss's additional research note (1:22 AM EDT): vapour pressures of ices, magma-ocean redox, ICRP weights, comet ratios; new items F1 to F5 |
| H | [H-gen465-second-pass-tables.md](body-composition-parts/H-gen465-second-pass-tables.md) | The GEN.465 second extraction pass (UNSCEAR 2008, Brewer and Fischer 2016, Zenodo 10527170, Lodders 2003, Weisberg 2006) and the GEN.468 water-share result (merged 7:06 AM EDT 2026-10-23) |

Built on the GEN.425 source check (`/mnt/project-files/reports/bulk-composition-sources-check-2026-10-17.md`, called G425 below), which already verified Seager 2007, Zeng and Sasselov 2013, Zeng 2016 and 2019, Fortney 2007, Lopez and Fortney 2014, Thiabaud 2015, Santos 2015 and Dorn 2015; the parts cite its sections and do not repeat that work. All five parts agree with G425 where they overlap (differences are listed in section 7).

Reach: the first hours ran on a narrow proxy (arXiv, PyPI, raw GitHub); Boss widened it at about 12:02 AM EDT, after which the parts redid their NOT REACHED items. Still refused at the end: ADS, Science, Oxford Academic, ScienceDirect, A&A, Annual Reviews, Europe PMC, and the IOP and Elsevier full texts that sit behind bot checks. Labels used everywhere: READ (primary text read), 2ND-HAND (read in a citing paper or open-source code), COMPUTED (re-derived with Python from stated inputs), SNIPPET, RECALLED, NOT REACHED.





## 1. Result in one page

1. **The science exists for all of it, and most of it is now in tables with sources.** About 1,100 labelled rows in the five parts carry READ (about 550), 2ND-HAND (about 170) or COMPUTED (about 390), against about 60 RECALLED and about 45 NOT REACHED (counts in section 9). The weakest ground is listed in section 8, not hidden.
2. **One principle across all five parts: a body's composition is the cause; mass, density, radius, gravity, heat and radiation follow from it.** Every part recommends the same pipeline: the **star's abundances** (A) set the refractory makeup (Fe/Mg/Si, C/O, Th/U/K); the **snow line and formation distance** (C, D) set the ice and water share; the **structure solver** (B for rock and iron, C for ice, water and envelopes) turns composition and mass into radius, density, moment of inertia and core pressure; the **ground and radiogenic rules** (E) read the crust and mantle that B hands over.
3. **Both of Boss's routes work and the parts give both.** Route 1 (keep the drawn density/mass/radius and fill a composition that fits) is solvable with honest uncertainty: mass and radius fix the total iron share to about 0.01 but the core mass fraction only to about plus or minus 0.1 (B 2.3). Route 2 (draw the composition from the star, snow line and mass, derive the bulk figures) keeps everything mutually consistent by construction and is recommended for everything newly generated. Route 1 stays as the **inverse** (edit-everything, imports, validation).
4. **The existing code disagrees with the science in several places** (section 6): bulk-silicate-Earth radiogenic heat 37 to 50 percent too high, the ground-radiation draw not tied to the heat calculation, class D (ice) density drawn as rock, the water-layer depth 40 percent too deep at 10 percent water, a 10 km neutron-star radius floor with no support, a 10 Gyr main-sequence lifetime against 8 Gyr from the Sun's central hydrogen.
5. **The five decisions are taken** (section 5). Everything else carries a recommended default so work can proceed.
6. **What would still need Boss's own research** is down to publisher-walled papers (section 8), ranked; a second extraction pass can also take newly reachable sources (UNSCEAR 2008, VizieR, Zenodo data records).

## 2. What each body type is made of and what the generator must store

The cells give the composition model and the single source of the numbers. "S" is stored per body, "D" is derived on demand.

| Body | Composition model (what is drawn) | Drawn from | Derived | Part and section |
|---|---|---|---|---|
| Star, main sequence | X, Y, Z, plus Fe, Mg, Si, O, C, (Eu, Th, U, K) abundances | stored [Fe/H], population, age; branch (thin, high-alpha, halo, CEMP) and four [X/Fe] offsets drawn | X, Y, Z (Y = Y_P + 1.68 Z, Y_P = 0.2467), all element ratios, central H fraction | A 2.4 to 2.8, 3.2 |
| Star, evolved | surface C/N shift after first dredge-up, carbon stars, Wolf-Rayet N/C types, Pop III | star type, mass, phase | surface X, Y, C/O | A 2.10, 3.4 |
| White dwarf | DA/DB/DC/DQ/DZ, core C/O or O/Ne, pollution | initial-final mass relation from progenitor | mass-radius (Chandrasekhar integration, Sirius B to 0.1 percent), layers | A 2.11 |
| Neutron star | mass mixture, crust Fe/Ni then neutron-drip, core | stored mass | radius 11.6 to 13.1 km band, crust sequence | A 2.12 |
| Black hole | none beyond mass, spin, charge (no-hair); say so | stored mass, spin | r_s, ISCO, Hawking T, efficiency | A 2.13 |
| Brown dwarf | Y, Z of the natal gas; burning limits (D 13.55 MJ, Li 0.06 Msun, H 0.07 to 0.074 Msun, shifted by Z and Y) | host gas | radius (Baraffe table), limits | A 2.14, C 2.7 |
| Rocky/iron planet | core mass fraction, mantle Mg#, light elements in core, Mg/Si mineralogy, C/O flag | star's Fe/Si, Mg/Si, C/O (1:1 with scatter; low-probability super-Mercury branch) | R, g, v_esc, moment of inertia, core-mantle pressure, mantle oxides | B 2.1 to 2.8, 3 |
| Water/ice world, sub-Neptune | water share W of solids from snow-line position, formation distance and migrant probability; envelope fraction | snow line (per system), mass, formation distance | radius (closed form rms 1.2 percent, max 3.6), ocean/ice structure, base pressure | C 2.1 to 2.5, 3.3 |
| Ice and gas giant | heavy-element mass M_Z = 57.9 (M/MJ)^0.61 Earth masses, scatter 0.26 dex; Y_env from the natal gas | mass, host star | Z_env, radius from Fortney table plus inflation | C 2.6, 3.3 F |
| Rogue planet | the composition it was formed with; irradiation set to zero | formed branch | heat flow, surface state | C 2.8 |
| Asteroid belt | taxon mass shares (Bus-DeMeo classes with meteorite analogue, grain density, porosity, albedo) | position relative to snow line (kernel), per-system snow-line scatter | bulk density, water, long-lived and Al-26 heating | D 2.1 to 2.5 |
| Comet, Kuiper belt object, dwarf planet | refractory-to-ice ratio, volatile percentages relative to water, porosity | formation distance, size | bulk density, rock fraction, compaction | D 2.6, 2.7 |
| Moon | origin type (disk, impact, captured, ring, binary), ice fraction from the circumplanetary-disk ice line, iron share | parent mass, age, distance | density, ocean/ice shell | D 2.8 to 2.11 |
| Any ground | crust type (ultramafic, anorthositic, basaltic, felsic, ice, sulfur, carbon, evaporite, clay, oxidised, metal), grain state | mantle oxides (B), ice inventory (C, D), heat flux, water state | band albedo 0.4 to 0.7 um, regolith thermal inertia, K/Th/U, ground gamma | E 2.2 to 2.9, 3 |

## 3. The chain, in dependency order

```
 A star abundances ──► B core mass fraction, mantle oxides, mineralogy ──► E crust type, K/Th/U, ground
        │                      │                                              │
        │                      ▼                                              ▼
        │              structure solver ◄── C ice EOS, water share ──► radiogenic heat, ground gamma
        │                      │                  ▲                           │
        ▼                      ▼                  │                           ▼
   snow line per system ──► D small bodies, moons                    surface radiation (with atmosphere column)
        │                      │
        └───────────────► C / D water and ice shares
```

Layers (a worker cannot start a layer before the one above is stored):

1. **Research closed and decisions taken** (section 5). Pure documents.
2. **Tables as data files**: solar abundance set and [X/Fe] model (A), EOS constants and mantle density table (B), ice and water phase table and volatile ices (C), meteorite, taxon, comet-volatile and radiogenic tables (D), `materials.csv` with crust, oxide, albedo and radiogenic rows (E). Merge into the one materials table GEN.427 already plans (section 4).
3. **Star-side draws** (A): branch, [Mg/Fe], [Si/Fe], [O/Fe], [C/Fe], birth Th/U/K; plus the per-system snow line (C, D).
4. **Body-side draws** (B, C, D, E): core mass fraction and oxides, water share, ice species, taxon shares, moon origin, ground type, K/Th/U and Al-26.
5. **Derivations**: structure solver and closed-form radius, giant and brown-dwarf radius lookup, compaction and porosity, radiogenic heat, ground gamma, shielding column, albedo by band.
6. **Storage and migration**, **tests**, **docs**.

## 4. Where this overlaps the plan already in the repo

`docs/design/composition/atmosphere-colour.md` section 5 already holds GEN.425 (done), GEN.426 (this design), GEN.427 (materials table), GEN.428 (fractions drawn), GEN.429 (density, mass, radius derived), DB.69, GEN.430, TEST.210, DOC.65, plus GEN.411 (molar density from partial pressures; needs gravity from GEN.429). The research widens three of them:

| Existing item | What the research adds |
|---|---|
| GEN.426 design | Use this file and the five parts as its content; split the design note into the four notes in section 10 |
| GEN.427 materials table | E's 44-row CSV (`parts/E-...md` 2.10) is the base: formula, molar mass, density, specific heat, conductivity, melting/sublimation, latent heat, albedo, attenuation, tags. Add B's EOS rows (Seager Table 1, SLB via burnman), C's ice phases (Ih to X) and volatile ices, D's grain densities |
| GEN.428 fractions drawn | Becomes the five draws of layer 3 and 4 above, not one |
| GEN.429 derived bulk | One structure solver (B) with C's ice layers and Lopez and Fortney envelopes; closed forms as fast path and as tests (section 7 resolves which) |
| DB.69 | One migration set for A, B, C, D and E fields (section 10 cluster 1), with deterministic back-fill only if Boss wants old seeds kept |
| TEST.210 | The invariants in section 11 |

New groundwork that belongs with the 9.0 headline because the atmosphere is known by then: **atmospheric radiation** and **surface radiation** (E 2.9: ground gamma per ppm, air-equivalent column rule, existing `physics/radiation.py` proton range law matches PSTAR to 5 percent up to 400 MeV), and **radiogenic heat** replacing the single `a_rad` draw.

## 5. Decisions (taken by Boss, 1:24 AM EDT on 10-17, as the coordinator relayed them; card taps are not visible to this lane)

| # | Decision | Taken |
|---|---|---|
| 1 | Route | **Decided: composition first (route 2)** for everything newly generated; route 1 stays as the inverse for editing and imports |
| 2 | Water share of solids beyond the snow line | **Decided: keep the switch** (Marboeuf-derived 0.14 to 0.39 or the 0.5 convention). The default setting was not stated; the draft carries the Marboeuf-derived sequence as default. The Lodders 2003 full text was read on 2026-10-23 (GEN.468, section 5.1): the default stays, and the switch names its denominator |
| 3 | Radius reference level | **Decided: 1 bar** (scale-height correction for envelope tables) |
| 4 | Reference abundance scale | **Decided: Asplund 2021 proto-solar** log eps as the scale for ratios (not a preset); spread among A21, Lodders 2009 and Lodders 2019 recorded as the uncertainty on ratios |
| 5 | The two stellar figures | **Decided: replace both**: neutron-star radius from the 11.6 to 13.1 km band (no 10 km floor), main-sequence lifetime from the star's own mass and composition instead of a fixed 10 Gyr (A 2.12, 3.4) |

### 5.1 The water share beyond the snow line, researched (GEN.468, 2026-10-23)

The two values do not disagree about the chemistry. They count different things, and they start from different total ice to rock ratios.

**Read (Lodders 2003, ApJ 591, 1220, Table 11, equilibrium condensation, solar-system composition, mass percent of the total nebular mass):** silicates and oxides 0.303, metal and FeS 0.186 (rock 0.489); H2O ice 0.571, CH4 ice 0.330, NH3 ice 0.097 (ices 0.998); all condensate 1.487. Water ice over rock is 1.17 (2.09 for the Anders and Grevesse 1989 set, which the paper revised to about unity). Water ice condenses at 182 K at 10^-4 bar with 22.75 % of the oxygen already in rock. The photospheric composition gives the same ratios.

**Computed from those numbers:** the water share of the solids is 0.539 outside the water line (rock and water only), 0.494 once ammonia hydrate has condensed at 131 K, and 0.384 once methane has condensed too (all ices counted).

**Read (the Bern model, Emsenhuber et al. 2021, arXiv 2007.05561, section 3.3.3; Venturini et al. 2024 data, Zenodo 10527170):** the model uses the Marboeuf 2014 condensation model, notes that water is about 60 % of the ices by mass, and maps all ice species to water ice in the planet structure. The water mass fraction of its water-vapour planets tops out at 0.51 to 0.53 (median 0.509 at 10 to 30 AU, 95th percentile 0.529). So the convention of about 0.5 is the total ice share labelled as water; its water proper is about 0.3 of the solids, inside the Marboeuf-derived 0.14 to 0.39.

**What is explained.** (1) The convention counts water against rock plus water (or all ices lumped as water); the Marboeuf-derived range counts water against all solids. (2) What remains is a real factor of about two in the total ice to rock ratio: 2.04 for equilibrium chemistry with the carbon as CH4 (Lodders), 1 +/- 0.5 calibrated on comets (Marboeuf, second-hand). No source read says which one planetesimals follow.

**Decision for the default.** The switch default stays the Marboeuf-derived sequence (Boss, 1:21 AM EDT on 10-17). The switch is defined as the water share of all solids, with the option "0.5 convention" read as the water share of rock plus water (0.49 to 0.54 from Lodders 2003), so the two settings no longer pretend to be the same quantity. The one overlap is 0.38 to 0.39, where both say the same thing for a disc with every ice condensed. The 170 K snow-line centre with 0.2 dex of scatter (Boss, 10-18) is unchanged. Extracted tables and arithmetic: `design-drafts/gen465-second-pass/README-gen465-gen468-findings.md` in the shared project folder.

Lane-decidable, no card needed: neon following O or Mg and the O plateau at low [Fe/H] (A-R2); the super-Mercury branch probability (B, 0.05 to 0.23, a parametrisation); whether ages and metallicities are linked (A-D3, optional).

## 6. Findings against the existing code (to TODO as bugs or design items)

All read from the repo main of 2026-10-17 by the parts. Each is a failure the research lane does not own.

| # | Where | Finding | Size | Part |
|---|---|---|---|---|
| F1 | `physics/radiation.py` / radiogenic constants | Bulk-silicate-Earth radiogenic heat from the repo's values is 29.8 TW against 21.7 TW (Palme and O'Neill) or 19.8 TW (McDonough and Sun): 37 to 50 percent too high. The Urey ratio 0.67 against 0.47 to 0.52 | large | E 2.8 |
| F2 | `a_rad` draw | The ground-radiation factor (factor-of-4 draw) is independent of the heat calculation; real surface gamma spans about 650 from ultramafic to felsic | large | E 2.9a |
| F3 | class D density | "Frozen hydrocarbons and ice" draws 3.9 to 5.4 g/cm3; real icy bodies are about 1 to 2.5 | large | D 2.1, trace item 8 |
| F4 | `physics/hydrosphere.py` | Incompressible water-layer depth is about 40 percent too deep at a water fraction of 0.1; density constants are within 1 to 10 percent | medium | C 2.4 |
| F5 | stellar remnants | 10 km neutron-star radius floor has no support in the cited constraints; 10 Gyr main-sequence lifetime against 8 Gyr (decided, section 5) | small | A 3.4 |
| F6 | `tuning.py` | Disc metallicity gradient comment says 0.06 dex/kpc, Hayden 2015 gives 0.053; the `solids_factor` uses 10^[Fe/H] where Z/Z0 reaches 2 to 3 at [Fe/H] of -1 or lower | small | A 3.6 |
| F7 | giant radius | The repo's Neptunian exponent 0.519 against 0.589 in Chen and Kipping; scatter 10 and 7 percent against 7.4, 14.6 and 4 percent | small | C 2.6 |
| F8 | `SNOW_LINE_AU_AT_1_LSUN = 2.7` | It is the 170 K blackbody value (derivation in G425 2.12); real snow lines move by a factor of 2 to 3 (about 0.2 dex), so a fixed value per star is a simplification | design | C 2.1, D 2.4 |
| F9 | Thorngren citation | arXiv:1511.07442 in notes is not Thorngren et al. 2016; the paper is arXiv:1511.07854 | docs | C 2.6 |
| F10 | `validation.py` re-run | Recalculation after orbit validation redraws albedo, greenhouse and the mix (trace 2026-10-17); relevant here because composition must not be redrawn either | design | trace |
| F11 | `generate_asteroid_composition`, `ASTEROID_COMPONENTS` | Belt composition is four names drawn uniformly from an 84-name list mixing elements, major minerals, rare minerals and aromatic hydrocarbons, independent of the star, distance or snow line; belt density is a separate uniform draw (details in 14.5) | large | D 2.3 to 2.5, R1 |

## 7. Differences found between sources, and the call made

| Topic | Sources | Call |
|---|---|---|
| Closed-form radius for rock plus iron plus water | B (structure solver, fast fit max error 0.40 percent) and C (closed form (1.0688 − 0.2317 CMF (1−x)) M^(1/3.7) (1 + 0.4492x − 0.1597x²), rms 1.2 percent, max 3.6 percent) agree to 0.1 to 0.8 percent where they overlap | One solver (B) extended with C's ice layers is the reference; both closed forms stay as the fast path and as tests; documents carry each one's validity range (B 0.1 to 20 Earth masses, C 0.3 to 20) |
| Zeng 2016 Eq. 2 | calibrated on PREM (Earth's seismic model) | Not usable as a law without the "nothing copied from the Solar System" caveat; use as a test target only (B, G425) |
| Fortney 2007 pure water | 12 percent bigger than Zeng 2016 | Use the newer water EOS (G425 correction 3, C 2.3) |
| MgO bulk modulus | Seager 177 GPa against SLB 161.4 GPa | Choose and cite one (G425 correction 7); B uses SLB via burnman (2ND-HAND) |
| Thiabaud Fe/Si, Mg/Si, C/O box | A shows it is a compilation envelope wider than a homogeneous catalogue (88 percent of a solar-neighbourhood draw inside it; metal-poor high-alpha stars fall below its Fe/Si edge) | Use A's drawn distributions, keep Thiabaud only as a plausibility box |
| Slope of planet versus star iron | Plotnykov 2025: 0.94 (+1.02/−1.07); Brinkman 2024: 1.3 ± 1.0 (OLS) and 5.6 ± 1.6 (ODR) on the same data | 1:1 with scatter as default (Thiabaud 2015 agrees); keep the slope as a stored tunable |
| Where small-body mixing stops | G425 zero-pressure rule; D shows it gives 0.29 to 0.33 ice for Ganymede against about 0.45 from interior models | Simple mixing with porosity below about 1e20 kg; B/C structure above; D 2.10 |
| Chondrite and bulk heat | D (chondrites, 3.4 to 5.4 pW/kg now) and E (bulk silicate Earth) use the same isotope constants (Dye 2012, within 0.1 percent for U and Th) | One heat function (cluster 4, section 10) |

## 8. What still needs Boss, ranked, and what is newly reachable

Only items no reachable source could replace. Detail and "number it gives" are in each part's section 5.

| Rank | Reference | Gives | Parts |
|---|---|---|---|
| 1 | **Lodders 2003**, ApJ 591, 1220 (full text; READ on 2026-10-23 through the ADS link gateway, see 5.1; the IOP page is a flaky bot check) | Condensation temperatures and the water share beyond the snow line; K condensation; the exact solar set the earlier notes quote | C, A, E, D |
| 2 | **Stixrude and Lithgow-Bertelloni 2011 and 2024** originals (Oxford Academic refused) | Replaces the 2ND-HAND copy of mantle mineral EOS taken from burnman | B, E |
| 3 | Geochemistry compilations behind Elsevier and AGU: Palme and O'Neill 2014, Rudnick and Gao 2014, McDonough and Sun 1995, Gale 2013, Workman and Hart 2005 | Turns the oxide and K/Th/U tables from 2ND-HAND to READ | E |
| 4 | **Giant-planet interiors**: Helled 2022 review, Wahl 2017, Debras and Chabrier 2019; Thorngren and Fortney 2018 inflation table | Core mass and diluted-core placeholder; hot-Jupiter inflation efficiency | C |
| 5 | **Moon formation**: Canup 2012 (Science), Cuk and Stewart 2012, Mosqueira and Estrada 2003 | Composition of impact-formed moons; regular-satellite composition against distance | D |
| 6 | **Comet and ice references**: Mumma and Charnley 2011, Bardyn 2017, Fray and Schmitt 2009 (ice densities), Bierson and Nimmo 2019 full text | Comet volatile and refractory tables; solid CH4, NH3, H2S densities; compaction law for ices | D, C |
| 7 | **Stellar and planet pairing**: Brewer and Fischer 2016 tables via VizieR; Ca and Al per star; mantle FeO against star and water (Elkins-Tanton and Seager 2008) | Joint Fe/Mg/Si distribution; planet-minus-star scatter | B, A |
| 8 | Water EOS tables: French 2009, French and Redmer 2015, Haldemann 2020 (AQUA) | Superionic phase and a single water table; replaces a 300 K isotherm plus polytrope | C |
| 9 | Space weathering darkening, frost and dust reflectance, Hapke 2012 | Albedo of weathered regolith and frosts; lab-to-Bond conversion | E |
| 10 | Rogue-planet volatile stripping at ejection | None found anywhere | C |
| 11 | **Fray and Schmitt 2009** itself (Planet. Space Sci. 57, 2053; ScienceDirect refused) | Range of each vapour-pressure polynomial, the choice between the two sets, and an independent cold-end check of the working table in part F | F, C, D |
| 12 | **Sossi 2020 beyond Earth**: Gaillard 2014 and 2021, Hirschmann 2012 | How magma-ocean oxidation state varies across planets (needed for the drawn range in item F2) | F, B |
| 13 | **ICRP Publication 103** | Neutron weighting function and tissue table for the surface-dose conversion | F, E |

**Newly reachable and not yet extracted** (a short second pass can take them, so they do not need Boss): UNSCEAR 2008 Annex B (radon and soil Ra-226, E item 4), VizieR front pages (Brewer and Fischer tables, if the catalogue pages open), and a Zenodo data record for steam-world/migration populations (record 10527170; C item 2 wants planet-by-planet ice fractions). Probed at about 12:55 AM EDT, front pages only.

## 9. Verification state

Counts are table rows (or label mentions where a part said so) as each part reported them; they are not strictly comparable between parts but show the balance.

| Part | READ | 2ND-HAND | COMPUTED | RECALLED | NOT REACHED | SNIPPET |
|---|---|---|---|---|---|---|
| A stars | 92 (+8 targets) | 9 (+1) | 103 (+7) | 9 | 1 row, 2 items | 0 |
| B rocky | 33 | 46 | 31 | 2 | 1 row, 9 gap items | 0 |
| C volatile/giants | 78 | 8 | 33 | 7 | 5 | 7 |
| D small bodies/moons | 167 | 22 | 49 | 5 (+ prose items) | 10 references | 4 |
| E ground/radiogenic | 178 | 83 | 171 | 34 | 27 | 4 |

Every part ran at least one numeric check against a published table or a known body, and states the error: A fits to catalogues (88 percent inside the Thiabaud box, Sirius B mass-radius to 0.1 percent); B's solver 0.20 percent in radius against a scipy reference over 90 cases, the stoichiometry reproduces Unterborn 2023 (core fraction 0.2917 against 0.29); C's closed form rms 1.2 percent; D reproduces Carry's compaction transition near 1e20 kg; E recomputes isotope heat to 0.1 percent (U, Th) and 2.2 percent (40K). Calibration targets in the parts are test-only (Earth, Mars, Moon, Galilean moons, Sirius B): **no Solar System value is used as a default anywhere**.

## 10. Natural work clusters (input to the collision and consolidation sweep)

Each part proposed its own items (A 17, B 12, C 21, D about 24, E 16; the lists are at the end of each part). Where the same job appears in several parts, one worker should do it once:

| # | Cluster | Merges | Why one worker |
|---|---|---|---|
| 1 | **One schema and migration set** | A-S1, B.9, C item 18, D schema and back-fill, E13 | One revision after 0105; shared columns (`cmf`, `water_mass_fraction`, ice fractions, K/Th/U); the choice between columns and the child table `body_composition` (D) is made once: scalar fields as columns, species lists as the child table |
| 2 | **One materials table and loader** | A-T1, B.2, B.3, C items 5 and 6, D material, meteorite, taxon and radiogenic tables, E3 to E5 | GEN.427 already plans it; E's CSV is the base with a `labels` column kept |
| 3 | **Per-system snow line** | C item 8, D zone kernel and per-system scatter | Same draw (0.2 dex around the 170 K blackbody position), stored once per system |
| 4 | **One radiogenic module** | A-D5 (birth Eu/Th/U/K), E6, E8, E16, D heating | Same isotope constants, same decay law, Al-26 and long-lived terms together |
| 5 | **One structure solver and radius fits** | B.7, C items 12 to 14, D hand-off to B and C | Same ODE and EOS tables; one reference radius level |
| 6 | **Giant and brown-dwarf radius** | C items 7, 10, 11, 14, A-T3 BD limits | Same Fortney/Baraffe table interpolators and burning limits |
| 7 | **Star-side abundance draws** | A-D1, A-D2, B.1 (mass-to-molar conversion) | B reads A's ratios; the conversion is one function |
| 8 | **Hydrosphere upgrade** | C item 15, E7 (ground type), D ocean moons | One rewrite of `hydrosphere.py` against `water_mass_fraction` |
| 9 | **Ground, albedo and colour hand-off** | E5, E11, D taxon albedo, MAP.244 / GEN.418 | Albedo by band feeds the atmosphere-colour chain once |
| 10 | **Surface radiation** | E9, E10, `radiation.py` changes, E15 docs | Replaces `a_rad`, adds the column rule |
| 11 | **Tests** | A-V1/V2, B.11, C 19/20, D invariants, E14 | One sweep (TEST.210) with the invariants below |
| 12 | **Design docs** | A-X1, B.12, C 21, D docs, E15 | Four notes: stellar composition (A), bulk composition (B plus C, extending `bulk-composition.md`), small bodies and moons (D), ground, radiogenic and radiation (E) |
| 13 | **Research that blocks nothing** | A-R3, E1, D research group, C items 1 to 3 | Read-only; can start at once; the water-share decision (5.2) blocks only C's draw rule |

## 11. Invariants a test can assert (from the parts)

- Mass fractions of each body sum to 1 (rock, metal, ice species, envelope); X + Y + Z = 1 for stars; mantle oxides sum to 1.
- Round trip: composition to radius and back (route 1 inverse) returns the composition inside the stated uncertainty; radius is monotone in mass at fixed composition and in water fraction at fixed mass.
- Stored density equals mass over volume of stored radius; surface gravity from the structure equals the stored value.
- Core mass fraction of a star-matched planet equals the stoichiometric value from the star's Fe/Si and Mg/Si within the draw scatter; Fe/Si, Mg/Si, C/O of a drawn population keep the quantile targets of A section 4.
- Heavy-element mass of a giant never exceeds its mass (Z at most 1); brown-dwarf limits shift with Z and Y in the stated direction.
- Radiogenic heat of a body of age t equals the isotope-wise sum; chondrite and bulk-silicate values both reproduce their calibration target within tolerance.
- Hydrosphere: the sealed high-pressure-ice threshold follows from the water fraction (about 6e-3) and base pressure is stored.
- Small bodies: porosity is zero above the compaction pressure; the mixing rule hands over to the structure solver at the stated mass.

## 12. Next steps

1. Boss's decisions in section 5 are taken (1:24 AM EDT).
2. Hand this file and the five parts to the TODO thread: it builds the trees, back-scatters the building blocks to 8.1, runs the collision sweep with the clusters of section 10, then clusters the issues table and updates the documents.
3. A short second extraction pass (section 8, newly reachable) before the design is frozen.

## 13. Additions from Boss's research note (1:22 AM EDT)

Checked in [part F](body-composition-parts/F-boss-note-check.md). What changed:

- **Water-ice sublimation** (Feistel and Wagner 2007 form with the note's coefficients) is correct: COMPUTED against IAPWS R14-08 to 0.06 percent from 60 K to 273.16 K.
- **Vapour pressures of the cold ices** (N2, CO, CH4, CO2, NH3, H2S, SO2) now have a working table in bar (Fray and Schmitt 2009 form), transcribed from the image of the NASA/Icarus supplement, and checked against the NIST WebBook triple points (within 3 percent; the 4 percent case is the low-temperature CO2 polynomial). The cold end below 50 K and the range of each polynomial are not independently verified; the primary paper is still unreachable. The older two-parameter fits are far off in the cold (CO six times too high at 40 K).
- **Magma-ocean redox** (Sossi 2020) is confirmed by reading it (0.037 plus or minus 0.005, IW + 0.5, f_H2O/f_H2 = 1.05, f_CO2/f_CO = 0.21 at 2173 K). Two corrections: about 35 percent of the CO2 stays dissolved (the note says carbon goes to the gas), and everything is Earth's inventory, so the generator needs a drawn oxidation state per rocky planet (item F2) with Earth as a test target only.
- **ICRP 103 table**: weights are right but tissue weights are not tied to particle types, and the neutron weighting spans 2.5 to about 20.7, not 5 to 20 (RECALLED; ICRP 103 not read).
- **Comet ratios** are consistent with part D's READ values; **Thorngren scatter** stays 0.26 dex (the note's 0.3 to 0.5 is not in the paper as read); the **Mars dust row** attributes Lemmon 2004 and Tomasko 1999 to SO2 aerosols, which is wrong (atmosphere lane).
- New work items F1 to F5 (vapour-pressure function and test, magma-ocean redox state, lava-planet vaporisation, clathrates and brines, data licence column) attach to clusters 2, 7, 8 and 11 of section 10.
- Everything else in the note (Tiers 2 and 3 atmospheric optics, haze, photochemistry, trapped particles) belongs to the atmosphere lane.

## 14. Weighted possibility tables: the draw format (Boss, 1:40 AM EDT)

Requirement as typed by Boss: for every composition category there is a table of possibilities with weighted probabilities; each row lists the object classes it supports (the gate); the draw uses only the rows that pass the gate and the probabilities add up to 100 percent. The coordinator agreed with two refinements: weights may depend on the star and the formation distance (route 2), and every weight carries its source. The parts already contain most of the content; this section restates it in that format and says which tables are ready.

### 14.1 Table format (one file per category, loaded like the materials table)

| Column | Meaning |
|---|---|
| `table`, `row` | Category and possibility (for example `ground_type`, `basaltic`) |
| `gate_kinds` | Object kinds the row supports: `star_ms`, `star_evolved`, `white_dwarf`, `neutron_star`, `brown_dwarf`, `rocky_planet`, `ice_water_world`, `sub_neptune`, `ice_giant`, `gas_giant`, `rogue_planet`, `asteroid`, `comet_kbo`, `dwarf_planet`, `moon_rocky`, `moon_icy`, `moon_volcanic`, plus the existing class letters as an alias list while GEN.90 reworks the classes |
| `gate_conditions` | Optional physical predicates on stored quantities: mass range, equilibrium temperature, position relative to the snow line, water state, age, [Fe/H]. Predicates are preferred to class letters because the classes are being reworked |
| `weight` | A number, or the name of a weight function `w(star, formation distance, mass, temperature, age)` returning a value of at least 0 |
| `weight_source`, `label` | Citation and verification label (READ, 2ND-HAND, COMPUTED, RECALLED) or **PARAMETRISATION** (a flagged stand-in, not a published number) |
| `payload`, `payload_source` | The composition the row stands for (mass fractions, oxide table, ice species set, and so on) with its own source |

Draw rule: take the rows whose gate passes, evaluate each weight, normalise so the weights sum to 1, draw one row with the seeded wrapper (`draw.py`). Two table kinds exist and each category says which it is:

1. **Choice table** (pick one row): ground type, mantle branch, asteroid taxon, moon origin, white-dwarf type, abundance branch.
2. **Mixture table** (all passing rows contribute, proportions add to 100 percent): the share of each taxon in a belt, the ice species set of a body, crust and mantle end-members. Proportions come from the weights directly (deterministic) or from a Dirichlet draw centred on them with a stated concentration (a parametrisation, flagged).

Invariants a test asserts: every table has at least one row for every object kind that uses it (a gate that leaves no row is a load-time error, never a silent fallback); weights are at least 0 and normalise to 1 after gating; every row has `weight_source` or the PARAMETRISATION flag; payload fractions sum to 1; the table file version is part of the generator epoch (any changed weight changes the epoch, as the plan already requires for generator changes). The continuous quantities (core mass fraction, X, Y, Z, water share, heavy-element mass) stay distributions drawn from stored inputs, not tables; where Boss wants them tabulated they can be binned into rows with a weight function per bin.

### 14.2 Inventory: every category in the parts, in this format

Status codes: **G+S** gateable now, weights sourced; **G+P** gateable now, weights are a flagged parametrisation (a number exists, no published weight); **G+N** gateable, no weights yet (research or Boss needed); **D** deterministic from stored inputs (a rule, not a draw); **C** continuous distribution (not a table).

| # | Table (category) | Rows | Gate kinds | Weight or weight function | Source and label | Status | Part |
|---|---|---|---|---|---|---|---|
| T1 | Abundance branch of a star | thin, high-alpha, halo-low, CEMP | star_ms, star_evolved | P(high-alpha) = 1/(1 + exp(([Fe/H] + 0.350)/0.132)) times a population factor (young and intermediate: thin only; old, bulge: logistic; halo: high); CEMP 0.25 for -4 < [Fe/H] < -2 rising to 1 below -4 | Catalogue fits (Bensby 2014, Nissen 2014, Battistini and Bensby 2016): COMPUTED; CEMP fractions 20 to 30 percent below [Fe/H] of -2: READ, the stepwise rule: PARAMETRISATION | G+S (CEMP step G+P) | A 2.7, 3.2 |
| T2 | Stellar surface state, evolved stars | normal, carbon star, WN, WC, post-dredge-up | star_evolved | thresholds in mass, Z and phase (carbon stars need C/O above 1 and a minimum mass) | Boothroyd and Sackmann, Crowther, carbon-star minimum mass: READ and 2ND-HAND | D (thresholds sourced, no weights) | A 2.10 |
| T3 | White-dwarf type | DA, DB, DC, DQ, DZ | white_dwarf | DA about 0.80 of the population, pollution (DZ) 25 to 50 percent of cool WDs; dependence on Teff and age not sourced | A 2.11 (READ for the fractions as printed); age dependence NOT REACHED | G+S for the share, G+N for the temperature dependence | A 2.11 |
| T4 | White-dwarf core | He, C/O, O/Ne/Mg | white_dwarf | by progenitor mass through the initial-final mass relation and mass thresholds | IFMR READ (Kalirai, Cummings); core C/O split RECALLED | D (rule), C/O split G+N | A 2.11 |
| T5 | Rocky-planet iron branch | star-matched core, super-Mercury, coreless | rocky_planet | star-matched 1 - p_SM - p_coreless; p_SM 0.05 to 0.23; coreless a small fraction | Slope data Plotnykov 2025 and Brinkman 2024: READ; the branch probabilities: PARAMETRISATION | G+P | B 2.8, 3.1 |
| T6 | Mantle mineral assemblage | pyroxene-dominated, olivine and pyroxene, olivine plus oxide, carbide-bearing | rocky_planet, moon_rocky | fixed by Mg/Si (0.86 and 1.73 mass-ratio thresholds, molar 1 and 2) and C/O (0.6 mass, 0.8 molar) of the star | Thiabaud 2015: READ, thresholds re-derived: COMPUTED | D | B 2.5, 2.7 |
| T7 | Core light element | S, Si, O (each 0 to 10 percent of the core) | rocky_planet | weights not published; a range exists | Hakim 2018, Unterborn 2023: READ for the ranges; the draw: PARAMETRISATION | G+P | B 3.2 |
| T8 | Mantle oxidation state (FeO, Fe3+/total Fe, delta-IW) | reduced to oxidised | rocky_planet, moon_rocky | FeO 2 to 12 wt percent (median 7.7); delta-IW range IW - 4.5 to + 3.1 explored | Putirka and Rarick 2019 (READ); Sossi 2020 (READ); no distribution across planets exists | G+N (the distribution is Boss's research item 12) | B, F |
| T9 | Formation branch: born beyond the snow line or migrant | beyond, migrant, in situ dry | rocky_planet, ice_water_world, sub_neptune | beyond if a is at least r_snow; migrant probability 0.14 to 0.67 above a minimum mass 3 (M_star/Msun)^0.48 Earth masses, 0 below | Ranges from a population synthesis (READ, figure reads); the function: PARAMETRISATION | G+P | C 3.3 A |
| T10 | Ice species set of a body | H2O, CO2, CH3OH, NH3, H2S, CH4, CO, N2 | all bodies formed beyond the matching ice line | mass shares from the cumulative condensation sequence (mg per mol of H2: 10.72, 5.237, 2.858, 0.708, 0.406, 0.573, 3.333, 1.165), included when the formation temperature is below the species' condensation temperature | Marboeuf 2014 abundances and saturation curves: READ; condensation temperatures recomputed: COMPUTED | G+S (mixture table; the water share convention is a switch, decided) | C 2.1, 2.2 |
| T11 | Water delivered inside the snow line | dry, delivered (1.2e-4 to 5e-3), migrated icy core | rocky_planet | p_deliv 1 if a giant sits outside the snow line, else 0.5 | Raymond and Izidoro: READ for the range; the probabilities: PARAMETRISATION | G+P | C 3.3 B |
| T12 | Giant interior | diluted core, compact core | ice_giant, gas_giant | no weights; heavy-element mass is continuous (57.9 (M/MJ)^0.61 Earth masses, 0.26 dex) | Thorngren 2016: READ; core mass placeholder log-uniform 5 to 20: PARAMETRISATION | C (core structure G+N) | C 2.6, 3.3 F |
| T13 | Asteroid taxon shares in a belt | Bus-DeMeo classes (C, S, X, D, P, V and subclasses) | asteroid | mixture weights as a kernel in x = a/a_snow, per-system snow-line scatter 0.2 dex | DeMeo and Carry 2013 Table 6: READ; scatter: Mulders 2015. The kernel is the measured Solar System zone mix rescaled by a/a_snow, so it is a Solar System-derived shape. **Decided (Boss, 1:43 AM EDT): stand-in only**, flagged as such in the table; replacement item R1 | G+S as a stand-in | D 2.3, 2.4 |
| T14 | Meteorite analogue of a taxon (rock family) | CI, CM, CO/CV, H, L, LL, enstatite, iron, achondrite | asteroid, small moons | taxon-to-analogue map; within-taxon weights not sourced | Carry 2012, DeMeo and Carry: READ for the mapping | D for the map, G+N for within-taxon weights | D 2.3 |
| T15 | Comet formation flavour | with or without refractory organics | comet_kbo | factor 0.67 on the ice-to-rock ratio, p = 0.5 | Marboeuf 2014 (READ) for the two values; the probability: PARAMETRISATION | G+P | C 3.3 B |
| T16 | Comet volatile abundances (relative to H2O) | CO, CO2, CH4, NH3, H2S, ... | comet_kbo | medians and ranges, not weights (CO 3.3 percent, range 0.3 to 26; CO2 12, range 2 to 30) | Harrington Pinto 2022, Rubin 2019: READ | C (could be binned) | D 2.6 |
| T17 | Moon origin | disk, impact, captured, ring, binary, in situ | moon_rocky, moon_icy, moon_volcanic | rules by parent (giant: disk; terrestrial: impact or capture; ring only inside the Roche zone); probabilities not published | Canup and Ward scaling READ; the weights: PARAMETRISATION | G+P (rules D) | D 2.8, 2.9 |
| T18 | Moon ice fraction | ramp with r / r_ice | moon_icy | r_ice = 20 R_J (M_p/MJ)^(2/3) (tau/5 Myr)^(-1/3); ramp width not sourced | Sasaki 2010: COMPUTED from the paper's scaling; far from 1 MJ untested | D (G+N for the ramp) | D 2.8 |
| T19 | Dwarf-planet and KBO rock fraction | rock-rich to ice-rich | dwarf_planet, comet_kbo | about 70 percent rock, spread 0.58 to 0.82 | Bierson and Nimmo (abstract only); own COMPUTED values | G+P | D 2.7 |
| T20 | Ground type | ultramafic, anorthositic, basaltic, felsic, ice, sulfur, carbon, evaporite, clay, oxidised, metal | rocky_planet, moon_rocky, moon_icy, moon_volcanic, asteroid, dwarf_planet (each row lists its own subset) | decision rules on mobile-lid chance, ocean, ice inventory, heat flux, surface temperature; weights themselves not published | Crust genealogy Hu 2012: READ; the weights: PARAMETRISATION | G+P | E 3.2 |
| T21 | Grain state | slab, crushed, powder | any solid surface | not sourced; band albedo differs 2 to 4 times between slab and powder | USGS splib07a and Paragas 2025: COMPUTED band means; the weights: NOT REACHED | G+N | E 2.6 |
| T22 | Crust and mantle end-members (oxides) | the reservoirs in E 2.3 (8 Earth reservoirs) | rocky_planet, moon_rocky | mixture proportions from T20 and the differentiation rules | Oxide values 2ND-HAND (pyrolite package); primary geochemistry NOT REACHED. **Decided (Boss, 1:43 AM EDT): stand-in only** (measured Earth rock types used as a type library, flagged as such); replacement item R2 | G+S for the rows as a stand-in, G+N for the proportions | E 2.3 |
| T23 | Radiogenic budget | K, Th, U concentrations, Al-26 ratio | rocky_planet, moon_rocky, asteroid | Th and U follow Eu and the star (A); K depletion distribution not sourced; Al-26 log-uniform 5.25e-6 to 5.25e-4 | A 2.9 (READ and COMPUTED); Al-26 range: D (COMPUTED from Lodders 2021, Ruedas 2017) | C for Th and U, G+N for K, G+P for Al-26 | A 2.9, D 2.5, E 2.8 |

### 14.3 Reading the inventory

- **Ready as tables today (G+S): T1, T3 (share), T10, T13, and the row libraries of T22.** These have published weights or weight functions.
- **Gateable with flagged stand-ins (G+P): T5, T7, T9, T11, T15, T17, T19, T20, T23 (Al-26).** A number exists but is a parametrisation; each is marked in the file and listed for replacement when a source is found (the part files' "to verify" lists).
- **Gateable, no weights yet (G+N): T3 temperature dependence, T8, T12 (core structure), T14 (within taxon), T18 (ramp width), T21, T22 (proportions), T23 (K).** These need research or Boss before they can be tables; T8 and T21 are the largest holes.
- **Deterministic (D), no table needed:** T2, T4, T6; keep them as rules with their thresholds in the tuning data so the thresholds also carry sources.
- **Continuous (C):** T12 (mass), T16, T23 (Th, U) stay distributions; they can be binned if Boss wants every category in the same format.
- **Decided (Boss, 1:43 AM EDT, card tap relayed by the coordinator): mixture proportions are the table weights directly where a sourced kernel exists (T10, T13), and a Dirichlet draw where only a range is known.** The concentration of each Dirichlet is a stored, flagged parametrisation.
- **Weights that depend on the star and the formation distance** (the coordinator's refinement) are already the shape of T1 (on [Fe/H], age, population), T9 and T10 (on the snow line, formation temperature and mass) and T13 and T18 (on position relative to the snow line and the parent's mass).

### 14.4 Decided stand-ins and the research that replaces them

Boss's ruling (1:43 AM EDT): the Solar System measured shapes (asteroid taxon shares, Earth rock-type library) are **stand-ins only**: they stay in the tables flagged `STAND-IN` and each has a research item to replace it with a derived or general distribution. The same flag applies to the other table content that comes from measured Solar System objects (the meteorite-group library of T14, the Al-26 range of T23 and the Earth, Mars and Moon calibration rows); those are listed with the two named ones so the sweep does not miss them.

| Item | Stand-in today | Replacement research | Output | Depends on | Size |
|---|---|---|---|---|---|
| R1 Belt composition derived, not measured | DeMeo and Carry 2013 zone mix (T13) rescaled by a/a_snow | Planetesimal composition as a function of formation temperature (condensation sequence, C and D parts), then dynamical mixing by giant-planet migration and scattering (Walsh 2011 and Raymond and Izidoro 2017 type models); fit a general kernel `P(taxon | a/a_snow, giant-planet history)` and test that it returns the measured zone mix as one realisation, not as the input | A drawable, derived kernel with stated uncertainty | T9, T10 sources; per-system snow line (cluster 3) | large |
| R2 Rock types derived from bulk composition | Eight Earth reservoirs (T22, E 2.3) as the oxide library | Derive crust and mantle end-members from the body's bulk oxides (B's mantle oxides) by melting and differentiation: partial-melt and fractional-crystallisation models with an open petrology code (check what is reachable: pyMelt, burnman, Perple_X type tools on PyPI or raw GitHub) and the magma-ocean crust literature; the Earth reservoirs become test targets | Oxide table per body from its own bulk, with the rock-type names as labels of the result | B stored oxides; E crust types | large |
| R3 Meteorite groups as outputs, not inputs | CI, CM, CO/CV, H, L, LL, enstatite, iron, achondrite grain densities, porosities, chemistry (T14, D 2.2) | Express the groups as regions of composition space (oxidation state, metal fraction, water, refractory ratio) reachable from the drawn bulk, so a new system's asteroids are generated from its own chemistry; keep the measured groups as calibration | A composition-space map and a nearest-group label | R1, R2, F2 (redox state) | large |
| R4 Al-26 and Fe-60 | Log-uniform range 5.25e-6 to 5.25e-4 from the Solar System's own value and range (T23) | Derivation from stellar nucleosynthesis sources and the cluster the star formed in (supernova and Wolf-Rayet yields, birth-environment draw) | A drawn ratio tied to the star's formation environment | A (population, cluster), superbubble chain | medium |
| R5 Calibration rows kept but fenced | Earth, Mars, Moon, Galilean moons, Sirius B targets in the parts | Move to the test files with a `calibration_only` flag; assert no generator table reads them | Test that no table row is sourced only from a measured Solar System body | cluster 11 | small |

Until each lands, the table row carries `STAND-IN` in `label`, the loader reports the count of stand-in rows at start-up, and the TODO tree gets R1 to R5 as items under the 9.0 composition headline.

### 14.5 Existing code: the asteroid belt and field composition draw (Boss, 1:45 AM EDT)

Read from repo main (`7ae024f`, 2026-10-17): `src/planetgen/generation/belt.py` 29-60 (`generate_asteroid_composition`, called by `AsteroidBelt.__init__` at line 140), `tuning.py` 196-215 (`ASTEROID_COMPONENTS`) and 2900-2966 (`ASTEROID_FIELD_COMPOSITIONS`), `generation/phenomena/asteroid_field.py` 69-75 and 150-157, `db/models.py` 20-34 (`asteroid_belt_composition`).

**What it draws today.**

- **Belt composition is four components, not three.** `generate_asteroid_composition()` takes `draw.sample(ASTEROID_COMPONENTS, k=4)` (a uniform draw without replacement from a flat list of **84** names), shuffles them, and assigns the levels high, moderate, small, trace in that order. The result is stored as (component, concentration) rows in `asteroid_belt_composition` (`position`, `component`, `concentration`, concentration limited to those four words) and as the text `composition_summary` ("high concentrations of iron, moderate concentrations of nickel, ...").
- **The list mixes everything with equal weight**: 16 elements (including gold and the platinum-group metals), 13 rock-forming and metal minerals (olivine, pyroxene, kamacite, troilite, magnetite, ...), 8 rare-earth minerals, 6 platinum-group minerals, 7 uranium and thorium minerals, 7 titanium species, 14 oxides, salts, carbonates and gem or industrial minerals (borosilicate, sodium chloride, calcium fluoride, ...), 4 garnets, 2 micas, chrysoberyl, and six polycyclic aromatic hydrocarbons (naphthalene to coronene). A belt is as likely to have gold, monazite or coronene as its "high" component as olivine.
- **Nothing is physical.** The draw does not read the star, the belt's distance, the snow line, the system's [Fe/H] or the belt's mass. The only other belt quality, density (`dense`, `sparse`, `typical`), is a separate uniform draw (`belt.py` line 138). No quantity is stored, only the four words.
- **Standalone asteroid fields** (`phenomena/asteroid_field.py`) already have a weighted table: eight families drawn by `frequency` (carbonaceous 30, stony 20, icy 20, mixed 10, metallic 8, dust 5, basaltic 4, collisional 3; that is 30, 20, 20, 10, 8, 5, 4 and 3 percent), each with its own 5 to 9 component list sampled uniformly (the `mixed` and `collisional` families sample the same 84-name list as a belt). The class letters (A to U) come from family and density. The frequencies carry no source in the code.

**Fold-in, in the section 14 format.**

| # | Table (existing) | Rows | Gate kinds | Weight today | Status in the new plan |
|---|---|---|---|---|---|
| E1 | Belt components | 84 names | asteroid (belts bound to a star) | uniform (1/84 each, four drawn) | **Replace.** Gate by object kind; weights from T13 (taxon mixture by a/a_snow) and T14 (analogue chemistry) so the components are the ones the taxon mix implies. Flag the existing list `STAND-IN` until R1 lands; the list stays only as a vocabulary of names |
| E2 | Belt density | dense, sparse, typical | asteroid | uniform | G+N: no source; should follow belt mass (stored mass per area, from the disc) rather than a draw. Link to R1 and to the planetesimal-mass work in cluster 3 |
| E3 | Field families | the 8 above | asteroid fields (standalone) | frequencies 30, 20, 8, 20, 4, 10, 5, 3 | **G+P**: the table already is a weighted, gateable table; the weights are a parametrisation without a source. Flag as a stand-in; replacement is R1 for fields (a field's composition should come from the system whose belts or bodies it was stripped from, or from a derived population, not a fixed table) |
| E4 | Component lists per family | 5 to 9 names each, sampled uniformly | asteroid fields | uniform | **Replace** with the T14 analogue chemistry of the family's taxon; keep the family names as labels |

**Proposed replacement for the stored data** (compatible with part D 3.6, which keeps `composition_summary` as derived text): a belt stores the taxon mass shares (T13) and the derived mass fractions of its main components (iron, nickel, olivine, pyroxene, phyllosilicates, water, carbon, metal), and the four-level words (high, moderate, small, trace) are produced from those fractions by fixed thresholds (for example above 20 percent of the belt mass, 5 to 20, 1 to 5, below 1) rather than drawn. Exotic minerals (rare-earth, platinum-group, radioactive, gemstones, polycyclic aromatic hydrocarbons) become accessory species included only when the chemistry supports them (a ppm threshold from part D's chondrite chemistry), not as 1/84 draws. A belt keeps at most four named lines in the text so the existing description and database columns still fit; the thresholds and the cap are stored, flagged parametrisations.

Tie to the decided items: this is the code that R1 replaces. Until R1 lands the current draw stays as the stand-in, but it should be flagged `STAND-IN` in the table and in the loader's start-up count, and its non-physical equal weighting recorded as finding **F11** (below).

| # | Where | Finding | Size | Part |
|---|---|---|---|---|
| F11 | `generate_asteroid_composition`, `ASTEROID_COMPONENTS` | Belt composition is four names drawn uniformly from an 84-name list mixing elements, major minerals, rare minerals and aromatic hydrocarbons, independent of the star, distance or snow line; belt density is a separate uniform draw | large | D 2.3 to 2.5, R1 |

### 14.6 One component vocabulary for every body (Boss, 1:47 AM EDT)

Boss's point: anything an asteroid field or belt can be made of, a planet or moon can be made of too. The design consequence is that the **component list is shared**, not per object type.

- **One components table** (it is the GEN.427 materials table, extended): every element, mineral, ice, organic and alloy any body can contain has one row with its formula, density, equation-of-state constants and the rest of the material properties (E 2.10, B 2.1, C 2.3, D 2.1). The 84 names of `ASTEROID_COMPONENTS`, the family lists of `ASTEROID_FIELD_COMPOSITIONS`, the class composition strings of planets and moons ("silicon, iron, and magnesium", "frozen hydrocarbons and ice", ...) and the comet component list all become views of this one table. Names that exist only as text today (anthracite, sulfides, mineral dust, hydrocarbons) get rows or are mapped to rows.
- **Every body composition is a mixture of rows of that table**, stored as component mass fractions (section 14.1 mixture tables). Belts, fields, asteroids, comets, planets, moons and dwarf planets use the same storage; the difference between a belt and a planet is only which rows pass the gate and what the weights are.
- **The gate does the physics, not the object type.** A row does not say "asteroid" or "planet"; it carries `gate_kinds` that can list all of them plus formation and condition predicates: condensation temperature at the formation distance (ices, organics, hydrated silicates), mantle oxidation state (native iron, sulfides, oxides), water and carbon inventory (carbonates, clays, serpentine), star abundances (carbides when C/O is high, rare-earth and uranium minerals only as trace accessories scaled to the star's Eu and Th). So olivine or kamacite are possible for an asteroid, a rocky planet or a moon alike, while coronene needs the organic-rich cold branch and gold is only a trace accessory.
- **Accessory minerals** (rare-earth, platinum-group, uranium and thorium minerals, gemstones, polycyclic aromatic hydrocarbons) enter any body at trace level from the body's bulk composition and a ppm threshold, never as an equal-weight draw; the same rule serves a belt and a planet's crust.
- **Descriptions** ("high concentrations of iron, moderate concentrations of nickel...") are produced from the stored fractions by fixed thresholds for any body, so a planet's atmosphere-and-crust text and a belt's text use the same code (GEN.430 and GEN.415 read the stored mixture).
- **Effect on the existing classes:** the planet classes' fixed composition strings are replaced by the drawn mixture; the class label becomes a summary derived from it (as D 3.6 already proposed for belts and as GEN.426 states for planets). The class letters of asteroid fields remain labels of family and density.

Work items: one components-table loader and name mapping (cluster 2), one mixture storage for all bodies (cluster 1, replacing `asteroid_belt_composition` and the planet composition text), one gate evaluator with the condition predicates (cluster 7), one thresholded text writer (cluster 12); the stand-ins of 14.4 apply to the vocabulary's first weights.

## 15. Boss's second research note, body claims (2:11 AM EDT): checked

Four read-only checks, one per topic, each with a verdict table, the figures to keep and an unreachable list: [G1 comets](body-composition-parts/G1-comets.md), [G2 radiogenic ratios and Thorngren](body-composition-parts/G2-radiogenic-thorngren.md), [G3 magma ocean and lava worlds](body-composition-parts/G3-magma-ocean-lava.md), [G4 ices, brines and licences](body-composition-parts/G4-ices-brines-licences.md). Scripts and tables: `parts/G-files/`. The atmosphere lane has the remaining sections of the note (optics, haze, particles, photochemistry).

| Check | Rows | Confirmed | Partly | Contradicted | Unverifiable |
|---|---|---|---|---|---|
| G1 comets | 31 | 17 | 9 | 2 | 3 |
| G2 radiogenic, Thorngren | 17 | 3 | 12 | 1 | 1 |
| G3 magma ocean, lava | 26 | 5 | 16 | 3 | 2 |
| G4 ices, brines, licences | 26 | 10 | 10 | 3 | 3 |
| **Total** | **100** | **35** | **47** | **9** | **9** |

A row is CONFIRMED only when a source that was read gives the figure or one inside the claimed range. The note's `[cite: n]` tags point to nothing, so it is a pointer list, not a source.

### 15.1 Corrections to the note that change the design

| Topic | Note says | Found | Label | Part |
|---|---|---|---|---|
| NH3 photolysis lifetime at 1 AU | 2.5e4 s | quiet Sun about 6.5e3 to 6.7e3 s (Heays 2017 rate 1.5e-4 per s; Biver 2007); a quiet-and-active average of 1.8e-4 per s gives 5.6e3 s (atmosphere lane GEN.496, from Drozdovskaya 2023); literature range 5.4e3 to 7.7e3. Both lanes read the same data; the two figures differ by 1.2 because one is the quiet Sun and one is a quiet/active mean. Store the quiet-Sun rate plus an activity factor | READ + COMPUTED | G1 |
| NH3 sub-mm line | 1.15 mm | none; ground rotational line 572.498 GHz (0.524 mm), inversion lines about 1.3 cm | READ | G1 |
| Comet ranges | CO 0.5 to 30, CO2 2 to 20, HCN 0.05 to 0.5 percent | CO 0.3 to 26 (0.15 for 103P), CO2 2 to 30, HCN from 0.03 | READ | G1 |
| Lifetimes | one value per species | quiet-Sun values; active Sun shortens H2O by 1.7 and CO2 by 2.35; "scale by stellar UV" is species-specific, not one scalar (cool-star CO falls by 3.5e-4, CH4 by about 80, NH3 barely) | COMPUTED | G1 |
| Comet class | dynamical class replaced by chemistry | class is a weak predictor: only CO separates (28 percent of variance, p = 0.035); continuous log-normals plus a carbon-chain-depleted flag is the parametrisation | READ + PARAMETRISATION | G1 |
| Venus K/U | 10 to 13, "Earth-like, confirmed by Venera" | not supported: Venera 8 and Vega give K/U 5.9 to 15 (x10^3), Th/U 2.3 to 3.0, errors 40 to 130 percent | 2ND-HAND | G2 |
| Mercury Th/U | 3.5 to 4.0 | 2.4 +/- 0.9 (K 1150, Th 220 ppb, U 90 ppb; abstract only) | SNIPPET + COMPUTED | G2 |
| Earth and Moon ratios | MORB Th/U 3.7 to 4.0; Moon K/U 2.5 to 3.0 | MORB Th/U about 3.0 (part E keeps 3.04), K/U 12.7e3; Moon K/U about 2.3e3 | 2ND-HAND | G2 |
| Stellar Th, U | metal-poor stars have higher [Th/Fe], [U/Fe] | Th/Fe rises as [Fe/H] falls only inside +/- 0.1 dex of solar (slope -0.338 +/- 0.049); U is not measured in disc stars | READ | G2 |
| Thorngren 2016 | Z_p/Z_star exponent 0.45 to 0.55, "atmospheric" metallicity | exponent -0.45 +/- 0.09, prefactor 9.7; Z is the bulk heavy-element fraction (upper limit on atmosphere Z); sample starts near 19 Earth masses | READ | G2 |
| Mass-Z law | (part C) 57.9 (M/MJ)^0.61, 0.26 dex | confirmed twice; own re-fit 55.6 (M/MJ)^0.60 with 0.32 dex raw rms | READ + COMPUTED | G2 |
| QFM against IW | Delta-IW +2 equals QFM | QFM is IW + 3.1 to 3.6 (1 bar, 2200 to 1500 K) and +3.4 to 4.0 (1 GPa); CO2 > CO crossover is IW + 1.0 to 1.4, so the "+2" boundary is right but is not QFM | COMPUTED | G3 |
| Carbon in reduced melt | experiments Delta-IW -1.4 to -4.7, CH4 | the cited papers are models; -4.7 appears only in N2 experiments; reduced melts hold mainly CO, CH4 needs water; gas CH4 is minor above 2000 K | 2ND-HAND | G3 |
| Solubility laws | linear in NBO/T, exp(-Delta-IW) | exponential in NBO/O (Iacono-Marziano 2012), mass action gives 10^(-0.5 Delta-IW); the note's forms are schematic | 2ND-HAND | G3 |
| Lava-world vapour | Na and K dominate to 2200 K; mean mass 40 to 60 g/mol | Na dominates to 2700 to 2800 K; mean molar mass about 25 near 2000 K and 33 near 3000 K | COMPUTED | G3 |
| Sossi 2020 scale | Earth at IW + 0.5 | the same gas ratios are Delta-IW 0.0 on the Hirschmann 2021 scale; the IW reference must be fixed in the generator, +0.5 is a Sossi-scale calibration target only | COMPUTED | G3 |
| Cold-ice vapour pressure | Blakley 2024, Grundy 2024, factor 1.6, "9.2 kJ/mol" | Blakley 2024 (O2, NH3 only), Grundy 2024 (N2, CO, CH4); a 2026 Grundy preprint multiplies all five by 1.19 +/- 0.03; 9.2 percent is the enthalpy excess (about 10.1 kJ/mol), not 9.2 kJ/mol | READ | G4 |
| CO2 solid EoS | Jager and Span: 0 to 800 K, 12 GPa | that range belongs to Trusler 2011; use Span and Wagner eq. 3.12 for sublimation (842.65 Pa at 150 K, 103971 Pa at 195 K, test targets) | READ | G4 |
| H2 | non-condensable above 10 K | H2 freezes below 13.96 K; Tc 33.145 to 33.18 K; HD Tc 35.908 K | READ | G4 |
| CO2 hydrate Q1 | 271.60 K, 1.044 +/- 0.002 MPa | 271.60 K and 1.044 MPa is Semenov 2023 with +/- 0.021 MPa; Q2 283.13 K, 4.494 MPa | READ | G4 |
| Licences | NIST, JPL 19-5 free | NIST SRD is copyrighted (permission needed); JPL 19-5 is "protected by copyright, US government purpose use rights"; only the 2013 Mainz snapshot states CC BY 4.0; Leiden, Bremen, Karkoschka state none | READ | G4 |

### 15.2 What changes in this draft and in parts A to F

1. **Part D 3.4.2, comet volatile draw:** replace "log-normal around 67P, sigma 0.5 dex" by the comet-population table of G1 (CH3OH population geometric mean about 1.9 percent against 0.21 at 67P; CH4 0.68 against 0.34; real scatter 0.22 to 0.37 dex, CO 0.64 dex). 67P stays as a test target. Interstellar-object CO/H2O range widens (3I/ATLAS 1.65). Table T16 takes the G1 ranges.
2. **Part F item F1 (vapour pressures):** for N2, CO, CH4, O2 and NH3 take the cold end from the QCM-corrected tables of G4 (F&S is up to about 15 times too high for CO near 25 K); CO2 from Span and Wagner eq. 3.12 (F&S set 1 is 7 percent off at 125 K and 41 percent at 100 K). The part F transcription of N2-1, CO-1 and CH4-1 agrees to 0.5 percent with the F&S curves in the Grundy 2024 data package (relabel READ + COMPUTED). Constants C and D2 to D4 for N2, CO and CH4 are curve-fit reconstructions (COMPUTED, valid only inside the fitted range).
3. **Part F item F2 (magma-ocean redox):** the generator fixes one IW scale (the Hirschmann 2021 one) and stores Delta-IW on it; Sossi's +0.5 for Earth is on his own scale. There is no published per-planet Delta-IW distribution; G3 F8 gives a labelled PARAMETRISATION built from Doyle 2019, Deng 2020, Guimond 2023 and the Rubie to Siebert range, and the stored quantities needed.
4. **Part F item F3 (lava planets):** use VapoRock-type chemical-potential vapour equilibria; species and the 25 to 33 g/mol mean mass from G3 replace the note's range. LavAtmos is arXiv:2210.10463 (van Buchem, Miguel, Zilinskas, van Westrenen; 9 oxides, GPL-3.0); VapoRock is AGPL-3.0. Both are tools to read, not code to copy into the repo.
5. **Part F item F5 (licences):** the claim that NIST, PubChem and USGS data are public domain is wrong for NIST SRD and unproven for PubChem and USGS. (Boss, 2026-10-22: the licence is not a column of any table; it is recorded in `DATA-LICENSE.md`.) derived tables from NIST SRD, JPL 19-5, Leiden, Bremen and Karkoschka need permission or a re-derivation from the primary measurements. The full table is in G4.
6. **Part E 3.3 (radiogenic draw):** keep x_Th sigma 0.22 dex as an upper bound (solar-twin scatter 0.05, pooled Th/Mg 0.10, Eu 0.16); add a mean shift with [Fe/H] (Th/Fe slope -0.34 dex per dex, Eu/alpha -0.16 to -0.23). `d_K` log-uniform 0.03 to 1 stays (observed bodies lie at 0.03 to 0.28). Keep part E's N-MORB Th/U 3.04. K/U and Th/U by reservoir are calibration targets with the scatters of G2, not presets.
7. **Part C 2.6 (giants):** M_Z law unchanged; state that Z is bulk, not atmospheric; atmosphere Z is f Z_p with f about 0.2 for Jupiter and Saturn (Thorngren and Fortney 2019); add the raw 0.32 dex rms next to the 0.26 dex intrinsic figure; the core-mass placeholder (log-uniform 5 to 20 Earth masses) is supported by Chachan 2025 (14.7 Earth masses).

### 15.3 New or changed work items

| Title | Scope | Attaches to | Size |
|---|---|---|---|
| G-1 Comet population table | Species ranges and log-normal scatters (G1) replace the 67P-centred draw; carbon-chain-depleted flag | cluster for D | small |
| G-2 Photolysis lifetime table | Species lifetimes at 1 AU for the quiet Sun, active-Sun factors, species-specific stellar-spectrum ratios; constants of the species table, not stored per body | species table | small |
| G-3 Cold-end vapour pressure | QCM-corrected N2, CO, CH4, O2, NH3 and Span and Wagner CO2; extends F1 and its test (triple points plus the 1.19 correction as a switch) | F1, cluster 2 | small |
| G-4 One IW scale | Fix the reference (Hirschmann 2021), store Delta-IW, convert Sossi and Frost buffers; test: CO2/CO and H2O/H2 ratios reproduce Sossi Eq. 4 | F2 | small |
| G-5 Vapour mean mass | Lava-world vapour species and mean molar mass from T and composition | F3 | medium |
| G-6 Licence audit (removed with GEN.500, Boss, 2026-10-17, 2:35 PM EDT) | Licence column filled from G4; permission list for NIST, JPL, Leiden, Bremen, Karkoschka | F5 | small |

### 15.4 Still unreachable (Boss must research), ranked

1. Fray and Schmitt 2009 and Grundy 2024 main text with its constants C to D4 (G4); Jager and Span 2012 full text.
2. Hirschmann 2016 Table 2 and Eguchi and Dasgupta 2018 (carbon solubility constants), Armstrong 2015 and Ardia 2013 originals, Duncan 2017 (G3).
3. Venera 8, 9, 10 originals (Surkov), Peplowski 2011 full text, Hypatia download, Frank 2014, Delgado Mena 2019 (G2).
4. Dello Russo 2016, Mumma and Charnley 2011, Bockelee-Morvan and Biver 2017, Huebner and Mukherjee 2015 (G1; the lifetime values come second-hand).
5. Gough 2011 and Nuding 2014 (deliquescence relative humidity); the JPL copyright page (G4).

## 16. Where the TODO thread placed the findings and the clusters

Findings of section 6: F1 and F2 (radiogenic heat 37 to 50 percent high, the ground dose independent of the heat) are GEN.459 (8.4 phase 1); F3 (class D density) is a requirement and test target of GEN.428, GEN.429 and TEST.210; F4 (hydrosphere depth) is GEN.461 (8.4 phase 2); F5 (neutron-star radius floor, lifetime) are GEN.463 and GEN.467 (8.1 phase 1, with the stellar tables of GEN.407, Core 3 proposed); F6 is GEN.464 (solids law, 8.1 phase 3) and a note on GEN.405 (gradient); F7 (giant radius law) is GEN.460 (8.4 phase 1); F8 (snow line) is GEN.458 (8.4 phase 1); F9 (the Thorngren citation) is not in the repository; F10 is GEN.422.

Clusters of section 10: 1 schema is DB.69 and DB.70; 2 materials table is GEN.427; 3 snow line is GEN.458; 4 radiogenic module is GEN.459; 5 structure solver is GEN.429; 6 giant and brown-dwarf radius is GEN.460; 7 star-side draws are GEN.439, GEN.440 and GEN.464; 8 hydrosphere is GEN.461; 9 ground, albedo and colour is GEN.441, GEN.462 and GEN.418; 10 surface radiation is GEN.444 and GEN.453; 11 tests are TEST.210 and TEST.211; 12 design docs are GEN.426, GEN.435, GEN.436 and GEN.457 (small bodies and moons); 13 research is done (GEN.432 and GEN.433 retired), with the second extraction pass GEN.465 and the water-share research GEN.468. Moons are GEN.466. Part D's drawings of taxa, comets and the compaction law are inside GEN.442. The supporting scripts and data of part C stay in /mnt/project-files/design-drafts/parts/C-volatile-and-giants-files/ and are not committed. Owners are proposals: Core 4 for the body chain, Core 3 for GEN.463 and GEN.467.

Part F and section 13 (Boss's research note, checked by the lane at 1:28 AM EDT): F1 (the vapour-pressure function with its test) is GEN.470 (8.1 phase 2, with the tables) and GEN.452 reads it; F2 (the magma-ocean redox state per rocky planet) is folded into GEN.451 with the stored field in DB.69; F3 (lava-planet vaporisation) is GEN.456, which now follows GEN.451; F4 (clathrates and brines) is GEN.471 (8.4 phase 3); F5 (the licence column) is a requirement of GEN.427, GEN.410 and GEN.438. The three gap rows (Fray and Schmitt 2009 itself, Sossi 2020 beyond Earth with Gaillard and Hirschmann, ICRP 103) are on Boss's research list (section 8 rows 11 to 13). The supporting script and images of part F (`F-files/`) stay in /mnt/project-files/design-drafts/parts/ and are not committed.

Section 14 (the weighted possibility tables, Boss, 1:40 AM EDT) is GEN.472 (the one gated format and loader) and TEST.214 (the sum test over every registered table), recorded in [composition-model.md](composition-model.md) section 4b. Its two decisions (Boss, 1:43 AM EDT, "Weights, Dirichlet where unsure", and 1:43 AM EDT, "Stand-in only") are in the same note, and the five stand-in replacements R1 to R5 of section 14.4 are GEN.473 (R1, the belt composition derived), GEN.474 (R2, rock types derived from the bulk oxides), GEN.475 (R3, meteorite groups as outputs), GEN.476 (R4, Al-26 and Fe-60 from the formation environment) and TEST.215 (R5, the calibration rows fenced into tests). The T1 to T23 inventory of section 14.2 is the list of tables GEN.472's loader receives, with the gate and weight of each. The asteroid belt and field code of section 14.5 (finding F11, tables E1 to E4) is GEN.477, with E3's field families folded into GEN.473; its stored shares are DB.70. Section 14.6 (one component vocabulary for every body) is the shared table of GEN.427 with its loader and name mapping, the one mixture storage of DB.69 and DB.70, the gate evaluator with the condition predicates in GEN.472, and the one thresholded text writer of GEN.430 and GEN.446 (the notes on those items). Section 15 (the check of Boss's second note, 100 claims) is recorded on GEN.451, GEN.456, GEN.459, GEN.460, GEN.470 and the licence column of GEN.410, GEN.427 and GEN.438, with six new small items: G-1 is GEN.495 (the comet population table), G-2 is GEN.496 (the photolysis lifetime table), G-3 is GEN.497 (the cold-end vapour pressures), G-4 is GEN.498 (one IW scale), G-5 is GEN.499 (the lava-world vapour mean mass) and G-6 is GEN.500 (the licence audit, removed by Boss at 2:35 PM EDT).
