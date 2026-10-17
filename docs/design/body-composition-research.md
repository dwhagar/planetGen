# Body composition: research findings and the composition design (draft for the 9.0 headline)

**Status (2026-10-17):** research result of GEN.432 and GEN.433 (the Body composition research lane), merged here by the TODO thread so that every lane has it in the repository. It feeds the designs GEN.426 (bulk), GEN.435 (stars), GEN.436 (radiation) and GEN.457 (small bodies and moons), the tables GEN.427, GEN.438 and GEN.439, and the generation chain of [composition-model.md](composition-model.md). Times in it are as the lane's clock reported. Section 14 says where the TODO thread placed each finding and cluster. The five decisions of section 5 are taken (Boss's card taps of 1:21 to 1:24 AM EDT, as the coordinator relayed them).

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

`docs/design/atmosphere-colour.md` section 5 already holds GEN.425 (done), GEN.426 (this design), GEN.427 (materials table), GEN.428 (fractions drawn), GEN.429 (density, mass, radius derived), DB.69, GEN.430, TEST.210, DOC.65, plus GEN.411 (molar density from partial pressures; needs gravity from GEN.429). The research widens three of them:

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
| 2 | Water share of solids beyond the snow line | **Decided: keep the switch** (Marboeuf-derived 0.14 to 0.39 or the 0.5 convention). The default setting was not stated; the draft carries the Marboeuf-derived sequence as default until the Lodders 2003 full text is read |
| 3 | Radius reference level | **Decided: 1 bar** (scale-height correction for envelope tables) |
| 4 | Reference abundance scale | **Decided: Asplund 2021 proto-solar** log eps as the scale for ratios (not a preset); spread among A21, Lodders 2009 and Lodders 2019 recorded as the uncertainty on ratios |
| 5 | The two stellar figures | **Decided: replace both**: neutron-star radius from the 11.6 to 13.1 km band (no 10 km floor), main-sequence lifetime from the star's own mass and composition instead of a fixed 10 Gyr (A 2.12, 3.4) |

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
| 1 | **Lodders 2003**, ApJ 591, 1220 (full text; IOP page is a bot check) | Condensation temperatures and the water share beyond the snow line; K condensation; the exact solar set the earlier notes quote | C, A, E, D |
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
- **Vapour pressures of the cold ices** (N2, CO, CH4, CO2, NH3, H2S, SO2) now have a working table in bar (Fray and Schmitt 2009 form), transcribed from the image of the NASA/Icarus supplement, and checked against the NIST WebBook triple points (within 3 percent; the 4 percent case is the low-temperature CO2 polynomial). The cold end below 50 K and the range of each polynomial are not independently verified; the primary paper is still unreachable. The older Prialnik 2004 two-parameter fits are far off in the cold (CO six times too high at 40 K).
- **Magma-ocean redox** (Sossi 2020) is confirmed by reading it (0.037 plus or minus 0.005, IW + 0.5, f_H2O/f_H2 = 1.05, f_CO2/f_CO = 0.21 at 2173 K). Two corrections: about 35 percent of the CO2 stays dissolved (the note says carbon goes to the gas), and everything is Earth's inventory, so the generator needs a drawn oxidation state per rocky planet (item F2) with Earth as a test target only.
- **ICRP 103 table**: weights are right but tissue weights are not tied to particle types, and the neutron weighting spans 2.5 to about 20.7, not 5 to 20 (RECALLED; ICRP 103 not read).
- **Comet ratios** are consistent with part D's READ values; **Thorngren scatter** stays 0.26 dex (the note's 0.3 to 0.5 is not in the paper as read); the **Mars dust row** attributes Lemmon 2004 and Tomasko 1999 to SO2 aerosols, which is wrong (atmosphere lane).
- New work items F1 to F5 (vapour-pressure function and test, magma-ocean redox state, lava-planet vaporisation, clathrates and brines, data licence column) attach to clusters 2, 7, 8 and 11 of section 10.
- Everything else in the note (Tiers 2 and 3 atmospheric optics, haze, photochemistry, trapped particles) belongs to the atmosphere lane.

## 14. Where the TODO thread placed the findings and the clusters

Findings of section 6: F1 and F2 (radiogenic heat 37 to 50 percent high, the ground dose independent of the heat) are GEN.459 (8.4 phase 1); F3 (class D density) is a requirement and test target of GEN.428, GEN.429 and TEST.210; F4 (hydrosphere depth) is GEN.461 (8.4 phase 2); F5 (neutron-star radius floor, lifetime) are GEN.463 and GEN.467 (8.1 phase 1, with the stellar tables of GEN.407, Core 3 proposed); F6 is GEN.464 (solids law, 8.1 phase 3) and a note on GEN.405 (gradient); F7 (giant radius law) is GEN.460 (8.4 phase 1); F8 (snow line) is GEN.458 (8.4 phase 1); F9 (the Thorngren citation) is not in the repository; F10 is GEN.422.

Clusters of section 10: 1 schema is DB.69 and DB.70; 2 materials table is GEN.427; 3 snow line is GEN.458; 4 radiogenic module is GEN.459; 5 structure solver is GEN.429; 6 giant and brown-dwarf radius is GEN.460; 7 star-side draws are GEN.439, GEN.440 and GEN.464; 8 hydrosphere is GEN.461; 9 ground, albedo and colour is GEN.441, GEN.462 and GEN.418; 10 surface radiation is GEN.444 and GEN.453; 11 tests are TEST.210 and TEST.211; 12 design docs are GEN.426, GEN.435, GEN.436 and GEN.457 (small bodies and moons); 13 research is done (GEN.432 and GEN.433 retired), with the second extraction pass GEN.465 and the water-share research GEN.468. Moons are GEN.466. Part D's drawings of taxa, comets and the compaction law are inside GEN.442. The supporting scripts and data of part C stay in /mnt/project-files/design-drafts/parts/C-volatile-and-giants-files/ and are not committed. Owners are proposals: Core 4 for the body chain, Core 3 for GEN.463 and GEN.467.

Part F and section 13 (Boss's research note, checked by the lane at 1:28 AM EDT): F1 (the vapour-pressure function with its test) is GEN.470 (8.1 phase 2, with the tables) and GEN.452 reads it; F2 (the magma-ocean redox state per rocky planet) is folded into GEN.451 with the stored field in DB.69; F3 (lava-planet vaporisation) is GEN.456, which now follows GEN.451; F4 (clathrates and brines) is GEN.471 (8.4 phase 3); F5 (the licence column) is a requirement of GEN.427, GEN.410 and GEN.438. The three gap rows (Fray and Schmitt 2009 itself, Sossi 2020 beyond Earth with Gaillard and Hirschmann, ICRP 103) are on Boss's research list (section 8 rows 11 to 13). The supporting script and images of part F (`F-files/`) stay in /mnt/project-files/design-drafts/parts/ and are not committed.

