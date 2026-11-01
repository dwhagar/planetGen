# Small bodies and moons: the composition draws, compaction, moon origins and the ice fraction

**Status (2026-10-21):** design note for GEN.457, written from the Body composition research lane's part D ([D-small-bodies-moons.md](body-composition-parts/D-small-bodies-moons.md)), its comet check ([G1-comets.md](body-composition-parts/G1-comets.md)), the master note [composition-model.md](composition-model.md) and the bulk note [bulk-composition.md](bulk-composition.md) (GEN.426). It changes no stored value and no draw; the build items do (section 11). It settles what part D left to the design: which kinds of body draw their composition first and which derive their density, the compaction law and its hand-over to the structure solver (with a finding that changes the way Boss's 2:10 PM EDT ruling of 2026-10-17 can be met), the taxon, comet, icy-body and moon-origin tables in the gated format of the GEN.472 loader, what is stored and what is derived, the parametrisations with their replace-when-sourced rule, a worked example for five bodies and the test targets. Times are EDT as the lane that reported them gave them, on 2026-10-17 unless a date is shown. Verification states follow the rule of 3:02 AM EDT: **verified**, **partly verified** (with the limit named), **open** (no source; the seeded range used meanwhile is shown), **stand-in** (a measured Solar System shape or a placeholder, flagged, with the item that replaces it). Where two sources disagree the rules of composition-model.md section 4a item 12 decide.

The numbers marked COMPUTED were made for this note with the two small programs in `/mnt/project-files/design-drafts/small-bodies/` (`handoff.py`, `taper.py`: a uniform-composition sphere in hydrostatic balance, Bierson and Nimmo's two porosity laws mixed by volume, third-order Birch-Murnaghan grains: ice K0 9.2 GPa, K' 5.5, silicate K0 110 GPa, K' 4.0, which are my stand-ins for the solver's EOS). They are a check of the shape of the design, not a replacement for the reference solver of GEN.429.

**Implementation status (2026-10-23, GEN.442):** sections 3 to 5 and the compaction of section 4 are built in `physics/compaction.py` (with the table `physics/data/compaction_table.json`, made by `scripts/build_compaction_table.py`) and `physics/small_body_composition.py` (the taxon table, the comet, icy-body and ring draws). Belts, comets and interstellar comets read them (`AsteroidBelt` takes its star for the stored snow line; a comet stores `mass_kg` from its makeup; `tuning.COMET_COMPOSITION` and `tuning.COMET_NUCLEUS_DENSITY_KG_M3` are gone). Section 1 below describes main before that work. Stored composition rows are DB.70; the condition-based classes of lane D read the stored rows. Lane K's interface note: `/mnt/project-files/reports/lane-k-contracts-2026-10-23.md`.

**Implementation status (2026-10-24, GEN.466 wiring):** every generated moon draws its origin (GEN.714), and its metal, rock, water and ice shares, rock porosity and mass come from the origin's makeup in `bulk_draw` (`_draw_moon`): a giant's moons share one disk system (`draw_disk_system`, drawn once before the moons), the formation distance of a disk moon is its present distance (moons migrate, so a stand-in), the rock's FeO, Al2O3 and Mg/Si still come from the star, and a moon lighter than 1e21 kg takes its mass from its class radius by the tapered brittle law (`moon_composition.mass_from_radius_kg`) while a heavier one goes to the structure solver with no macroporosity. The moon's class letter still decides its radius, atmosphere and life until the class cut-over (GEN.646), so a class and a composition can disagree until then. The disk system's partition of the total moon mass among the large moons is drawn but the room-filling loop of `generate_moons` still sets the number and size of the moons. The origin and disk are held while the moon is generated and are not stored (DB.70).

## 1. What exists at main

- **Belts and fields.** `AsteroidBelt` (generation/belt.py) draws four names from the flat list `tuning.ASTEROID_COMPONENTS` and gives them the words high, moderate, small and trace; the density word (dense, sparse, typical) is a separate uniform draw; the mass is stated as zero (NAV.62). A standalone field draws one of eight families by frequencies with no source. The replacement of that draw is GEN.477; this note gives it the taxon table it reads.
- **Comets.** `Comet` (generation/comet.py) draws a nucleus diameter of 0.5 to 20 km, an orbit and a list of ten components; `bulk_mass.comet_nucleus_mass_kg` takes the mass from one tuning density. No volatile abundance, ratio or porosity is stored.
- **Dwarf planets, class D.** A class D body (GEN.1091: small, airless (surface pressure under 10 Pa) and at least 1 percent of its mass water or ice; size is not in the rule; the intended scale is a dwarf planet or a moon, radius 50 to 500 km) draws the rocky density of every terrestrial class, 3.93 to 5.37 g/cm3 (finding F3).
- **Moons.** `generate_moons` (physics/planets.py) fills the room around a planet with independent `Planet(is_moon=True)` bodies: class from the terrestrial classes except `MOON_BLACKLIST`, radius up to the parent's divided by 10^(1/3), mass at most a tenth of the parent's, distance log-uniform from a Roche-like bound. Nothing about the parent's composition, the circumplanetary disk or the age enters.
- **The loader and the snow line.** The gated table loader exists (GEN.472, `generation/possibility_table.py`); `formation.snow_line_au` is 2.34 AU times the square root of the luminosity (a 182 K blackbody value, [R2726]; 2.7 AU at the 170 K first used); the per-system snow line with its 0.2 dex scatter is GEN.458.

## 2. The decisions

| # | Decision | State | Used meanwhile, source, limit |
|---|---|---|---|
| S1 | New small bodies and moons draw their composition first and derive the bulk figures (route 2); a stored body (an old row, an edit, an imported density) takes route 1, the inverse of D 3.3 | Decided (Boss, 1:21 AM, card "Composition first"); part D's advice to keep route 1 for bodies whose density is already drawn is answered by this table | Section 3 says kind by kind what is the stored size and what is derived |
| S2 | The bulk density of a small body is the zero-pressure mixing rule with a porosity, `rho = rho_g (1 - phi)`, `rho_g = 1 / sum(w_i / rho_i)`, and the porosity is the Bierson and Nimmo brittle law tapered to zero between 1e19 and 1e21 kg (section 4) | Law verified (Bierson and Nimmo Table 1 and Eqs. 1 and 2, re-read); taper is a parametrisation | The paper's law alone cannot close the jump at the solver (section 4.2): a finding that needs Boss's card |
| S3 | Ice and silicate grain densities are 920 and 3500 kg/m3; other ices are close to water ice | Verified (Bierson and Nimmo Table 1; Fray and Schmitt stays on Boss's list, the NIST and PubChem bracket is the labelled stand-in, Boss, 3:07 PM on 2026-10-18) | N2 950, NH3 850, CO 1000, CO2 1000 to 1600, CH4 500 kg/m3 (bulk-composition.md section 6) |
| S4 | Belts and fields draw their taxon mass shares as a Dirichlet draw around a kernel in `x = a / a_snow,sys`; the kernel and the shares are flagged `stand-in` rows | Decided (Boss, 1:43 and 1:43 AM); stand-in until GEN.473 | Section 5.1; the concentration is open |
| S5 | Comets, Kuiper belt objects and dwarf planets draw a refractory-to-ice ratio and, for comets, volatile percentages of water from the population table of GEN.495 | Verified for the ranges ([R0776, R3835], Harrington [R3535, R3782] database); the log-normal widths from G1 | Section 5.2; the carbon-chain-depleted flag and the organic share are parametrisations |
| S6 | Icy bodies draw a rock mass fraction from a truncated normal with centre 0.70 and a collisional tail | Verified for the centre (Bierson and Nimmo: 70 percent matches 14 of 17 Kuiper belt objects within 2 sigma); the 0.08 width and the 0.45 to 0.95 limits are not in the paper | Decision 26 of the bulk note: 0.30 to 0.82 range, most permissive; section 5.3 |
| S7 | A moon takes an origin type from a gated table (disk, impact, captured, ring, binary) and its composition from the origin | Weights all open; the rules of each origin are part D 2.8, 2.9 and 3.4.3 to 3.4.5 | Section 5.4 |
| S8 | A disk moon's ice fraction rises with its formation distance past the disk's ice line: `f_ice = f_inf clamp(ln(r / (0.4 r_ice)) / ln 2.5, 0, 1)`, `f_inf` N(0.35, 0.10) on [0.1, 0.6]; the bracket over all moon systems is 0.5 to 0.85 with [R0804] at the centre | Decided (Boss, contradiction sweep, taps 6:10 to 6:12 PM on 2026-10-18, fork 15) for the bracket; the ramp is a parametrisation | Section 5.5. [R0804, R3922] disagree on the ice fraction: the newest wins |
| S9 | The hand-over mass to the structure solver is a continuity condition per composition, not the constant 1e20 kg | Decided (Boss, card tapped 8:13 PM EDT on 2026-10-21): taper porosity to zero between 1e19 and 1e21 kg | Section 4.2 |
| S10 | One storage for every body's composition (DB.69 and DB.70), one thresholded text writer (GEN.430, GEN.446); belts keep their four-line text from the stored shares | Decided (composition-model.md section 4b) | Section 6 |
| S11 | Everything here rides the one generator epoch bump of the 8.4 chain (GEN.428 to GEN.430 and DB.69) | Decided (Boss, taps 6:08 to 6:10 PM on 2026-10-18, fork 1) | The table versions pin in `Registry.versions` |

## 3. Which kinds draw what first

The stored size is the quantity the generator draws today and keeps. Composition is drawn before density; the density, and the mass where the size is a radius, follow.

| Kind | Composition first | Stored size (unchanged) | Derived |
|---|---|---|---|
| Asteroid belt or field | taxon mass shares (5.1), then per taxon the bulk density, albedo and analogue | belt distance, limits, mass (NAV.62: the mass stays an input of the disc, GEN.477 E2) | mean bulk density `1 / sum(w_i / rho_i)`, mean albedo, ice, metal and water mass fractions, radiogenic heat, the density word and the text |
| Comet, interstellar object | refractory-to-ice ratio, organic share, volatile percentages (5.2) | nucleus diameter | porosity from its draw, bulk density, mass `rho V` (replaces the tuning density of `bulk_mass.py`) |
| Class D body, Kuiper belt object, dwarf planet | rock fraction (5.3), metal fraction | radius | grain density, porosity from the compaction law (4), bulk density, mass, gravity; above the hand-over the solver |
| Ring | as the icy body or the moon that made it (5.4) | ring extent | ice fraction and albedo from the composition |
| Moon | origin (5.4), then its row of 5.5 or the impact, captured, ring or binary rule | radius, mass (the existing limits stay the cap), distance | density, ice fraction, mix of rock and metal |

One exception to "the size stays": a moon of a disk system has its mass partition drawn from the disk's total (alpha over f, GEN.466); the existing room-filling loop of `generate_moons` stays for impact, captured and binary moons and caps the disk moons' total at `M_T`.

## 4. The mixing rule, the compaction law and the hand-over

### 4.1 The rule

Per body: mass fractions `f_ice`, `f_rock`, `f_metal` adding to 1 (tested to 1e-9) and the grain densities of S3. Then `rho_g = 1 / (f_ice/920 + f_rock/rho_rock + f_metal/rho_metal)`, `rho = rho_g (1 - phi)`, with the porosity `phi` of the law below. A bulk density under the ice grain density at the porosity ceiling (0.8 for a pristine comet, 0.6 otherwise) is refused and the porosity lowered, with a flag.

The law (Bierson and Nimmo, Table 1 and Eqs. 1 and 2, `P` in MPa): ice `phi_i = max(0.60 exp(-0.1 P), 0.20)`, silicate `phi_s = min(0.6 P^-0.11, 0.6)`; the two layers' porosities mix by volume. The earlier 1e6 Pa guess is dropped (decision 25 of the bulk note). A comet nucleus does not use the law: its porosity is drawn (5.2). The viscous law (their Eq. 3) prints a positive sign and is not coded.

### 4.2 The hand-over: a finding

Boss ruled at 2:10 PM EDT that the paper's law stays and the hand-over at 1e20 kg is shifted until the porosity jump (0.86 to 1.55 g/cm3) disappears. The computation for this note says that **shifting alone cannot do it**. The law has a porosity floor of 0.20 for ice and the silicate law falls as a slow power of pressure, so the mean porosity of a uniform body stays near 0.26 to 0.30 even at the mass of Ganymede, and the zero-porosity density of the solver stays 1.32 to 1.35 times the law's from 1e22 to 3e23 kg (COMPUTED, `handoff.py`):

| Rock fraction | Mass (kg) | Law radius (km) | Law density (kg/m3) | Zero-porosity density (kg/m3) | Ratio |
|---|---|---|---|---|---|
| 0.70 | 1e19 | 143 | 815 | 1901 | 2.33 |
| 0.70 | 1e20 | 283 | 1049 | 1902 | 1.81 |
| 0.70 | 1e21 | 573 | 1272 | 1907 | 1.50 |
| 0.70 | 1e22 | 1198 | 1389 | 1931 | 1.39 |
| 0.70 | 1e23 | 2515 | 1500 | 2020 | 1.35 |
| 0.55 | 1e23 | 2681 | 1238 | 1640 | 1.32 |
| 0.85 | 1e23 | 2326 | 1897 | 2602 | 1.37 |

The ratio 1.81 at 1e20 kg is the jump Boss named (0.86 to 1.55 is 1.8). A real body does not hold that porosity: the paper itself says viscous compaction decides the density above about 200 to 500 km radius and its brittle model is a lower bound, and the observations the design is tested against say the same ([R0666]: the compact-to-rubble transition lies at 1e19 to 1e20 kg and dwarf planets have near-zero macroporosity; Orcus, 6.3e20 kg, 1.6 to 1.75 g/cm3; Quaoar 1.99 g/cm3).

**Proposed means (parametrisation, flagged):** keep the paper's law and taper its porosity to zero in log mass, `phi_eff = w(M) phi_law`, `w = clamp((log10 1e21 - log10 M) / 2, 0, 1)`, between `M_a` = 1e19 kg (Carry's lower edge of the transition) and `M_b` = 1e21 kg (above it dwarf planets are near zero). Above `M_b` the body has no macroporosity and the solver of GEN.429 takes over with the same grains, so the density is continuous there by construction (the solver adds only its own compression). The computed densities (rock fraction 0.70, COMPUTED, `taper.py`): 815 at 1e19 kg, 1526 at 1e20, 1739 at 3e20, 1847 at 6.3e20 (Orcus: radius 433 km, density 1.85 against 1.6 to 1.75 measured), 1907 at 1e21, with no step. The endpoints are the one place the design leans on a measured Solar System shape, so they are `stand-in` until viscous compaction (Eq. 3 with its sign checked against [R0372] or the published paper) replaces the taper; the Galilean-class solver test of decision 25 still comes before the 1e21 kg to 0.1 Earth mass range is called validated.

Decided (Boss, card tapped 8:13 PM EDT on 2026-10-21): the taper, with the endpoints flagged `stand-in`.

### 4.3 The fast path

The law costs a hydrostatic integration. A body is generated millions of times, so the compaction module returns the mean porosity from a table indexed by log radius and rock fraction, built once by the reference integration, interpolated at generation and pinned by a table version (`Registry.versions`, composition-model.md section 4b); a test says the interpolation is within 2 percent of the integration on a fixed grid. The table has the taper in it. A small-body radius is the stored size, so the table is indexed by radius and the hand-over mass follows from the radius and the density.

## 5. The tables, in the gated format

Each row has the columns of composition-model.md section 4b: possibility, supported classes (the gate, which is the object kind and never a class letter), weight or weight function, source and licence (`open` where none is found), state, and the notes. Weights that depend on the star or the formation distance read the stored `snow_line_au_sys` (GEN.458) and the stored [Fe/H] (GEN.440). Every table has a version; a change of row or weight is a version change inside the one epoch bump.

### 5.1 The taxon table (belts and fields)

Possibility: a taxon of Bus and DeMeo, with its analogue and its grain figures. Gate: object kind belt or field. Weight: the zone kernel at `x = a / (a_snow 10^s)`, where `s` is the system's stored snow-line scatter (N(0, 0.2) in log10, GEN.458). State: **stand-in** for every row; replaced by GEN.473.

The kernel is [R1065: Table 6] (READ), in percent of the zone's mass, with the zone limits converted to `x` by the Solar System's 2.7 AU (the only place a Solar System number sets a boundary, and the reason the rows are stand-ins):

| Zone | x from to | Classes (percent of zone mass) |
|---|---|---|
| Hungaria | 0.66 to 0.76 | E 50, C 21, S 9, A 7, L 7, K 5 |
| Inner belt | 0.76 to 0.93 | V 69, S 21, C 6 |
| Middle belt | 0.93 to 1.04 | C 70, B 15, S 8, P 4 |
| Outer belt | 1.04 to 1.21 | C 52, P 15, B 13, M 10, S 5 |
| Cybele | 1.21 to 1.37 | P 84, C 13 |
| Hilda | 1.37 to 1.56 | P 71, D 15, C 14 |
| Trojans | 1.87 to 2.00 | D 67, P 26 |

Between zone centres the weights interpolate linearly in log `x` (the width of 0.3 in log AU of the open-gaps report, BC14, a placeholder); outside the table the nearest zone holds (no extrapolation). The kernel gives the sequence of dry and igneous inside, hydrated C in the middle, ice-rich P and D beyond, and nothing in the sources tests it for other stars.

The shares are a Dirichlet draw around the kernel (Decided, 1:43 AM). **The concentration is open** (no source gives the spread among belts; the Solar System has one belt). Meanwhile the draw uses 10, written in the row as `open` and listed in the report; the widest range reading (Boss, 1:12 PM) says a lower value is the safer reading. For each taxon: bulk density from [R0666: Table 3] (mean and scatter, truncated normal), albedo from DeMeo and Carry Table 4 times a log-normal scatter (the sigma is not in anything read: open), analogue chemistry and grain density from the meteorite group (Britt and Consolmagno). The aggregate follows by section 3.

The eight field families of today become a gated table of their own (GEN.477 E3): carbonaceous = C, B (and P), stony = S, Q and L, metallic = M and X, icy = P and D, basaltic = V, mixed = a Dirichlet over all, dust = no fixed taxon (CI and CM chemistry at porosity 0), collisional = fragments of its dominant taxon. Their frequencies have no source: stand-in, replaced under GEN.473.

### 5.2 The comet and icy-body volatile tables

The population table of GEN.495 gives the comet rows; this note fixes the interface and the figures it holds meanwhile.

- **Refractory-to-ice ratio `delta_RI`:** log-uniform on [0.5, 8] ([R0776]: Rosetta RSI 3 to 7, CONSERT above 3, Fulle 7.5, COSIMA 0.2 to 3). `f_ice = 1 / (1 + delta_RI)`. Verified range, parametrised shape.
- **Refractory split:** silicate and organics; the organic share of the refractories is N(0.60, 0.15) clipped to [0, 1], a labelled stand-in at the widest bracket (Boss, 1:23 PM: no measurement exists; COSIMA's 0.45 lies inside it). Silicate dust 3500 kg/m3 (a stand-in: [R0430] say the density of cometary refractories is not known and gives 2600 to 5000 for minerals); organics from the shared materials table, 1500 kg/m3, the middle of the same paper's 1000 to 2000 (GEN.712: the earlier 2000 was attributed to [R0776] without the paper being read).
- **Volatiles as percent of water** ([R3835] for 67P; Harrington [R3535, R3782] for the population, G1): CO median 3.3, range 0.3 to 26; CO2 median 12, range 2 to 30 (CO + CO2 about 18); CH4 0.15 to 2.2, geometric mean 0.68, sigma 0.37 dex; C2H6 0.09 to 0.82; NH3 0.125 to 0.9 with Mumma and Charnley's bound near 1.5; HCN 0.03 to 0.36; H2S 0.09 to 1.5; CH3OH 0.21 to 4.5, geometric mean 1.9; H2CO 0.09 to 1.1; C2H2 0.03 to 0.28. These replace part D's single-comet centre with sigma 0.5 dex (67P is a test target only, G1: its CH3OH 0.21 and CH4 0.34 are low for the population). N2, O2, Ar, Kr and Xe have no population (67P only): open, stand-in. Class by dynamical type predicts only CO (28 percent of the variance), so it is not a gate.
- **Interstellar objects:** the same draw with the CO/H2O ratio widened to include 2I/Borisov and 3I/ATLAS (0.01 to 1.65 log-uniform), one object each, so a parametrisation.
- **Porosity:** a comet nucleus draws `phi` from U(0.5, 0.8) (Choukroun 65 to 79 percent); a Kuiper belt object under 100 km U(0.3, 0.6); larger bodies the compaction law of section 4.

### 5.3 The icy bodies

Possibility: the rock mass fraction of a class D body, Kuiper belt object or dwarf planet. Weight: a truncated normal with centre 0.70 and width 0.08 on [0.45, 0.95], with a collisional tail toward 0.9 for the dwarf planets ([R0269]: Haumea, Eris). The centre is verified (Bierson and Nimmo); the width and the limits are the repository's and not in the paper (decision 26 of the bulk note keeps the most permissive range, 0.30 to 0.82, for the small bodies, with dwarf planets 0.78 plus or minus 0.12). Computed rock fractions of the measured dwarf planets (Barr and Schwamb densities, rock 2700 or 3500): Pluto 0.69 to 0.77, Charon 0.63 to 0.70, Orcus 0.58 to 0.65, Quaoar 0.73 to 0.82, Eris 0.86 to 0.96, Gonggong 0.65 to 0.72; they are test targets, not rows. The metal fraction is zero for a class D body and a small draw for a differentiated one.

### 5.4 The moon origin table

| Possibility | Supported (the gate) | Weight | Composition rule | Source and state |
|---|---|---|---|---|
| disk | parent is a giant (an envelope on it) | open | regular moons from the circumplanetary disk, section 5.5 | [R0643, R0804, R3922]; weight open |
| impact | parent is terrestrial (no envelope) | open | no metal phase in the condensate (iron enters as FeO), metal U(0.01, 0.05), moderately volatile elements depleted, water 1 to 1000 ppm, ice 0; the parent's mantle oxides are the source composition | [R2724] (verified for the rule); the metal range is a parametrisation; weight open |
| captured | any parent (GEN.714: the origin is drawn first and the orbit follows from it, so there is no inclination cut) | open | the icy-body rule of 5.3, then ice lost by a factor U(0, 0.5) | [R0269] (qualitative); the factor is a parametrisation; weight open |
| ring | the parent has a massive ring and the moon forms near the Roche radius | open | ice-rich, near the ring's composition; masses rise with distance | [R0920] (abstract only, RECALLED for the composition); rules not sourced |
| binary | the primary is a dwarf planet | open | the icy-body rule of 5.3; mass ratio q in 5e-4 to 0.11 | Barr and Schwamb Table 2; weight open |

Every weight is `open`: the open rows of one gate share the probability equally (a giant's moon is a disk moon or a captured one with equal chance, a terrestrial planet's an impact moon or a captured one) and each is in the loader's report until a source or a derivation replaces it.

**The orbit follows the origin (GEN.714, `moon_composition.draw_inclination_deg`).** A regular moon (disk, impact, ring, binary) lies near the parent's equator: uniform on 0 to 25 degrees, the generator's long-standing bound, wider than the regular satellites show (the Galilean moons lie within 0.5 degrees) and so a flagged stand-in. A captured moon takes an isotropic orbit (cosine of the inclination uniform) outside the Kozai-unstable band of 55 to 130 degrees ([R3263]; the bounds as used in [R1747], Shen and Tremaine give 60 to 130, partly verified), so it is prograde below the band or retrograde above it: 43 percent prograde and 57 percent retrograde by the width of the cosines. The distance range of a moon is unchanged (retrograde orbits are stable out to about 0.7 Hill radii, prograde to 0.4; the generator keeps 0.49 for both). The origin is held by the moon while it is generated and is not stored until DB.70. A gate that no row passes is an error that names the parent and the table.

### 5.5 Disk moons: the system, the ice line and the ramp

1. The disk lifetime `tau` log-uniform 1 to 10 Myr (Sasaki used 5 as reference); `alpha/f` log-uniform on [1e-6, 5e-4]; the total satellite mass `M_T / M_p = 2.5e-4 (alpha/f / 3e-5)^(1/3)` (Canup and Ward Eq. 4.15), a range of 0.8e-4 to 6.4e-4 (COMPUTED); the number `N` of large moons 1 to 5 and the partition a Dirichlet draw (parametrisations; Sasaki Table 1 is the calibration for the shape only).
2. The ice line, in kilometres and not in planet radii: `r_ice = 20 R_J (M_p / M_J)^(2/3) (tau / 5 Myr)^(-1/3)` with `R_J` = 71 492 km: 1.43e6 km at 1 Jupiter mass and 5 Myr (COMPUTED), 6.4e5 km at 0.3, 3.0e6 km at 3. The formula was fitted for a Jupiter-like disk; it is untested far from 1 Jupiter mass, so outside 0.1 to 3 the row carries a flag (partly verified, limit named).
3. The ramp of S8, with `r` the formation distance and not the present one (moons migrate): Io 0, Europa 0.061, Ganymede 0.239, Callisto 0.35 for `f_inf` = 0.35 at `r_ice` = 20 R_J (COMPUTED); the ordering is a test and the values are under the interior-model figures, as the simple-mixing failure of D 2.11 says they should be.
4. The rocky block: anhydrous inside 0.4 `r_ice`, partly hydrated between 0.4 and 1.0 (serpentine, density 2585), C-like or CI-like beyond; the metal fraction by chondrite group (H about 17, L 8, LL 3 percent, RECALLED) or from the star's iron minus the oxidised share.
5. Density by the rule of section 4.1, and above the hand-over by the solver; a regular moon of 1e22 kg and more is differentiated, so the zero-pressure rule under-predicts the ice of Ganymede by 10 to 15 points and is not used there.

## 6. What is stored and what is derived

The rule of composition-model.md section 7.4: store what was drawn, cannot be derived or is expensive and exact; derive the rest and test that a cache equals the recomputation. A recalculation after an orbit validation or a zone change never redraws a composition (finding F10).

| Body | Stored (drawn) | Derived on demand |
|---|---|---|
| Belt or field | `snow_line_au_sys` and its scatter (once per system), taxon mass shares (child table `asteroid_belt_taxa`: belt, taxon, mass fraction), the main component mass fractions (DB.70) | mean bulk density, mean albedo, mean porosity, ice, metal and water fractions, radiogenic heat, density word and the four-line text (thresholds stored as parametrisations) |
| Comet | `refractory_ice_ratio`, organic share, `porosity`, the volatile percentages of water (child table `comet_volatiles`: comet, species, percent) | `f_ice`, grain density, bulk density, mass |
| Icy body or moon | rock, ice and metal fractions (the one mixture storage of DB.69 and DB.70: body, component, mass fraction), `origin`, `formation_distance_ratio` (`r / r_ice` or `a / a_snow`), the surface porosity draw `phi0`, `drawn` or `inverse`, the hand-over flag | grain density, porosity, bulk density, mass or radius, gravity |

Names, units and ranges are those of part D 3.6; the columns belong to DB.69 and DB.70, not to this note. A route-1 row (an edited or imported density) stores `inverse` and the range of fractions that fit, as the bulk note section 7 says.

## 7. The parametrisations and when each is replaced

Each is a stored, flagged constant; a better published distribution replaces it and bumps the table version.

| Quantity | Value now | Basis | Replaced by |
|---|---|---|---|
| Snow-line scatter | log10 N(0, 0.2) | [R3193, R3196] as the coordinator quoted it (RECALLED) | GEN.458 when read |
| Taxon kernel and shares | section 5.1 | [R1065: Table 6] (stand-in) | GEN.473 |
| Dirichlet concentration of the taxon shares | 10 | none (open) | a spread among belts, when one is found |
| Albedo scatter within a taxon | none found | open | a measured spread |
| Compaction taper | 1e19 to 1e21 kg | [R0666] (stand-in) | viscous compaction (Eq. 3 sign checked) |
| Surface porosity `phi0` | C-like U(0.25, 0.55), S-like U(0.10, 0.35), comets U(0.5, 0.8) | [R0544, R0666, R0776] | a measured distribution |
| Rock fraction of icy bodies | N(0.70, 0.08) on [0.45, 0.95] | Bierson and Nimmo (centre only) | a derived distribution |
| `delta_RI`, organic share | log-uniform [0.5, 8]; N(0.60, 0.15) | [R0776]; Boss 1:23 PM stand-in | a measurement of the organic share |
| Comet volatiles | G1 ranges and sigma | [R3782] database, [R3835] | GEN.495 and a larger database |
| Moon system | `alpha/f`, `tau`, `N`, partition, ice ramp, `f_inf` | [R0643, R0804, R3922] | [R3180] when read |
| Impact-moon metal, captured-moon ice loss | U(0.01, 0.05); U(0, 0.5) | [R2724], Barr and Schwamb (qualitative) | [R0638, R0958] when read |
| Al-26 ratio per system | log-uniform [0.1, 10] x 5.25e-5 | [R2678] | GEN.476 |

## 8. Worked example: five bodies, as test targets only

Nothing here is a preset (TEST.215); each is a body some draw inside the allowed distributions can reach.

1. **A belt in the middle zone (`x` about 0.97).** Kernel C 70, B 15, S 8, P 4, renormalised to 0.722, 0.155, 0.082, 0.041. With the [R0666] bulk densities 1570, 2190, 2660, 1450 kg/m3 the mean bulk density is `1 / sum(w/rho)` = 1696 kg/m3 and the mass-weighted albedo 0.094 (COMPUTED). A draw moves the shares by the Dirichlet scatter; the sequence (C-dominated here, more P and D at larger `x`) is the test.
2. **A comet.** `delta_RI` = 3 (inside Rosetta's 3 to 7): `f_ice` 0.25; organics 0.45 of the refractories, silicate 3500 and organics 1500 kg/m3 give a refractory grain density of 2188 and `rho_g` = about 1630 kg/m3; porosity 0.70 gives about 490 kg/m3 (COMPUTED), against 537.8 for 67P (test tolerance 100). With the organic share at the stand-in centre 0.60 the result is about 460; at porosity 0.65 it is about 570, at 0.79 it is about 345: the porosity range alone spans 67P's measurement.
3. **A Kuiper belt object near 6.3e20 kg and rock fraction 0.70.** Taper weight 0.10, radius 433 km, density 1847 kg/m3 (COMPUTED, section 4.2), against Orcus at 430 to 450 km and 1.6 to 1.75 g/cm3: inside the 0.3 g/cm3 tolerance; the same fraction at 1e19 kg (143 km) has density 815.
4. **A regular moon system at 1 Jupiter mass and 5 Myr.** `r_ice` = 1.43e6 km; formation distances of 5.9, 9.4, 15.0 and 26.3 Jupiter radii give ice fractions 0, 0.061, 0.239, 0.35 (COMPUTED): Io, Europa, Ganymede and Callisto in order, the ordering exact.
5. **A captured moon.** Rock fraction 0.70, ice lost by a factor 0.3 of the ice: rock 0.79, grain density 2203 kg/m3, a mass above 1e21 kg so no macroporosity: 2203 against Triton's 2061 kg/m3, inside 0.3 g/cm3 (COMPUTED).

## 9. Test targets and invariants

Invariants (TEST.210 style; part D 3.7): fractions sum to 1 within 1e-9; the stored bulk density is reproduced from the stored fractions, grain densities and porosity within 1e-6 relative; the grain density lies between the smallest and largest component density; the ice fraction of a disk moon is 0 inside 0.4 `r_ice` and does not fall with distance; porosity does not rise with mass for a fixed family; a class D body has density under 3000 kg/m3 when its rock fraction is under 0.9 (F3); the long-lived heating of a belt is within 10 percent of the taxon's chondrite value; Al-26 exceeds the long-lived heating by more than 1000 at formation; a body above the hand-over carries the flag for a solver-sourced density; the density is continuous in mass across 1e19 and 1e21 kg within 5 percent of a step of 10 percent in mass (section 4.2, the test that Boss's ruling asks for); the same system gives the same shares and moons at any worker count; a recalculation moves nothing stored (F10); the table sums of the taxon, icy-body and origin tables are 1 within 1e-12 for every gate (TEST.214).

Targets from part D section 4 (each met for some draw, never for one): ordinary chondrite bulk density 3.42, 3.36, 3.22 g/cm3 within 0.2; 67P 537.8 kg/m3 within 100; mean Kuiper belt object 0.77 and comet 0.47 g/cm3 within 0.3; rock fractions 0.58 to 0.82 for Pluto, Charon, Quaoar and Orcus inside the draw; compaction transition near 1e20 kg within a factor 3; Galilean ice order exact and densities 3.53, 3.01, 1.94, 1.83 g/cm3 within 0.3 (Ganymede needs the solver); total regular-satellite mass ratio 2.1e-4 (Jupiter) and 2.5e-4 (Saturn) within a factor 3; Io 3.53 g/cm3 within 0.2; dwarf-planet mass ratio q from 5e-4 to 0.11 covered; Near-Earth object mix S/Sq/Q about 50 and C about 20 percent within 15 points; comet NH3 up to 1.5, H2S 0.25 to 1.5 and SO2 0.2 percent of water inside the draw; Makemake and 2007 OR10 near 1800 and MU69 near 750 kg/m3 reachable.

## 10. Open weights and gaps

- Open weights (the loader lists them): the five weights of the moon origin table; the Dirichlet concentration; the albedo scatter; the share of captured moons (an equal share with the other rows); the Kozai band's edges; the comet N2, O2, Ar, Kr and Xe; the field family frequencies (stand-in).
- No validated solver from 1e21 kg to 0.1 Earth masses (bulk note section 12 item 1): Ganymede-class moons and dwarf planets sit in the gap and run flagged extrapolated until the Galilean test passes.
- Not read (part D section 5): [R0638, R0958] (the impact moon's iron), [R3180], the full text of Bierson and Nimmo (Eq. 3), [R0920, R3206] (the ring rules), and a metal weight percent per chondrite group.
- The 1e19 to 1e21 kg taper endpoints are two numbers read off [R0666]'s statements, not a fit.

## 11. The items that follow

GEN.442 (belts, comets, rings, dwarf planets) reads sections 3, 4 and 5.1 to 5.3; GEN.466 (moons) reads 5.4 and 5.5; GEN.477 reads 5.1; GEN.495 (comet population table) fills 5.2; GEN.473 replaces the stand-ins; DB.69 and DB.70 store section 6; the compaction module and its table (section 4.3) are part of GEN.442's build with the test of 4.2. All ride the one epoch bump of the 8.4 chain.

## 12. Decisions for Boss

1. **Decided (taper, 8:13 PM EDT on 2026-10-21). Close the porosity jump at the solver.** Boss ruled at 2:10 PM EDT on 2026-10-17 to shift the 1e20 kg hand-over until the jump disappears. It does not disappear: the paper's law leaves the density 1.32 to 1.35 times low up to 3e23 kg (section 4.2). Recommended: keep the paper's law and taper its porosity to zero between 1e19 and 1e21 kg (continuous, matches Orcus and Triton, endpoints flagged `stand-in`). Alternatives: keep the step and move it to about 1e23 kg with a 1.35 times jump left; or leave the step at 1e20 kg until viscous compaction is coded.
2. **The taxon and field-family rows stay `stand-in` with an open Dirichlet concentration of 10.** Recommended yes (nothing sourced exists; GEN.473 replaces it).
3. **The organic share of comet refractories N(0.60, 0.15).** Recommended yes (the widest bracket around COSIMA's 0.45, as ruled at 1:23 PM EDT).

## 13. The derived taxon kernel (GEN.473, lock-free step)

`physics/taxon_kernel.py` replaces the Solar System's zone table (`small_body_composition.ZONES`, [R1065: Table 6]) with a kernel built from what makes the zones, in four steps (registry TK1 to TK6, section 6.41):

1. **Four families by disk temperature.** Dry (S, V, A, M, K, L), hydrated carbonaceous (C, B), primitive (P) and ice-rich (D). The dry family ends at the water line of [R2726] (182 K, the temperature a body must be colder than to accrete water ice), the primitive one at the ammonia line (131 K); the line between hydrated and primitive has no condensation line and is fitted (166 K, tweak L3). The temperature at a distance is the star's blackbody one, so the zones scale with the star's own luminosity and nothing is in the Solar System's 2.7 AU.
2. **A prediction, not a fit.** Ammonia hydrate condenses at 1.93 snow-line distances, and the Trojans (5.2 AU, 1.9 snow lines) are 67 percent D in Table 6: the one zone boundary that no number was chosen for. The C onset at the snow line (x of 0.93 to 1.04 in the table) and the P onset at about x 1.2 are the water line and the fitted one.
3. **Mixing.** 8 percent of each family's mass is exchanged with its neighbours (scattering by giant planets), fitted to the C share of the inner belt and the S share of the outer one (tweak L2). [R3681, R4603] are not read; the wiring step can give a system with no giant none.
4. **Melting from the star's own aluminium-26.** The dry family splits into the undifferentiated taxa (S, K, L) and a melted share's core (M), mantle (A) and crust (V). A body melts when its aluminium-26 heating exceeds the heat to its solidus: the ratio `Q / (X_Al / M_Al N_A E_decay)` is 1.1e-5, and with accretion times spread over 3 Myr after the isotope formed the melted share is 0.55 at the Solar System's 5.2e-5 and zero below the melting ratio. So the star's birth environment (GEN.476) decides how igneous its belt is: round a typical star (a median 26Al ratio of a few percent of the mean, `birth_radioisotopes`) a belt has no V, A or M at all.

The measured zones come back as one realisation: the families summed from Table 6 agree with the model to 0.16 rms in the share by zone (middle belt worst: the table has C 88 percent at x 0.98 where the model has the water line at 50 percent). What is not derived: the E taxon (a reduced, igneous one; zero here), the hot refractory K and L taxa (no temperature of their own), and the shares within each family (TK6, a stand-in at tweak L4; the meteorite-group map of GEN.475 is the replacement). The wiring step replaces `small_body_composition.kernel(x)` with `taxon_kernel.kernel(temperature, ratio)` for belts and fields, gives the belt the star's aluminium-26 ratio from `birth_radioisotopes`, and moves the taxon shares of every belt (fill stage, no migration).

**Wired (GEN.473, step 2, train 7).** `small_body_composition.draw_belt(distance_au, luminosity_w, al26_ratio)` takes the blackbody disk temperature at the belt's distance (`condensate_shares.disk_temperature_k`) and the aluminium-26 its star was born with (`radiogenic.star_al26_ratio`, shared by every body of the star, GEN.476) and draws the Dirichlet shares around `taxon_kernel.kernel`; `kernel(x)` over the Solar System's zones is gone and `ZONES` stays as the measured table the tests compare with. A belt's `detail` has `temperature_k` and `al26_ratio` in place of `x`. The star's snow line no longer decides the taxa (the water line is Lodders's 182 K at the star's luminosity, 2.2 AU for the Sun, against the stored snow line of GEN.458 which is 2.7 AU for it); the taxon rows stay stand-ins for their mixes of materials (GEN.475). Effect on generated systems: belt taxa and so belt density, albedo and water shift (a dry belt is as igneous as its star's aluminium-26 allows; the middle belt is hydrated carbonaceous, the outer primitive and ice-rich), the belt's stored composition words change.

## 14. The mass left in a belt (GEN.477 E2, lock-free step)

The belt density (dense, sparse, typical) was a uniform draw and the belt's mass was stated as zero (NAV.62); the item's E2 asked for a mass per area from the disc with a depletion factor that had no source. `physics/belt_mass.py` gives both from one law (registry BM1 and BM2, section 6.43):

- **What formed** is the disc's solids over the belt's annulus (`planet_mass.solids_surface_density_gcm2` integrated): 6.9 Earth masses at 2.2 to 3.4 AU round a Sun-like star.
- **What is left** follows [R4785] ([R4785], READ): a planetesimal belt in collisional equilibrium falls as `M0 / (1 + t / tc)` with `tc` falling with the mass, so the late-time mass does not depend on the mass that formed (Eq. 14 and 15): `M_max = 3.8e-6 rho r^3.5 (dr/r) Dc / (M*^0.5 t) * ... / G(11/6, Xc)`, which reproduces their Eq. 17 (0.009 at their stated assumptions, tested). It depends on the belt's radius and width, the star's mass, the age, the largest body `Dc` and the strength of the largest bodies (`Q*D`: the generator's gravity-bound Leinhardt and Stewart threshold, 2.5e3 J/kg at 100 km and 2.5e5 at 1000 km, with Wyatt's 200 J/kg of km-size bodies as a floor) and the particles' eccentricity (the belt's width gives it: `e = (upper - lower) / (upper + lower)`).
- **The Solar System's belt is a test, never an input:** 2.2 to 3.4 AU, 4.5 Gyr, Ceres-size bodies, gives 0.5 thousandths of an Earth mass at 500 km and 2.6 at 1000 km, against the 0.6 inferred ([R2471], as Wyatt et al. quote it); Wyatt et al. themselves quote 0.4 for their belt at 3 AU. The primordial disc is over 1000 times the result: the depletion factor of E2 is not a constant but the collisional law's.

The mass goes as `Dc` to a power near 2.3, so the belt is only as good as its largest body: `Dc` comes from the size laws (`small-body-size-laws.md`) once GEN.566 is built, and a belt stirred by a giant has more eccentricity than its width gives (so less mass). The density word is the belt's mass per area over a reference by a factor of 3 each way (BM2, T9, tweak L4, burn-down line); the reference is the same belt with the median largest body, around a star of the Sun's metals. What E1, E3 and E4 asked (taxon shares, the four-level words by thresholds, the field families) were built with GEN.442.

**Wired (GEN.477, train 7).** `AsteroidBelt` takes its mass and density word from `belt_mass.belt_state` at the star's mass, snow line, metals and age (`bulk_draw.metal_ratio_to_solar`), and stores the mass in `asteroid_belts.mass_kg` (a stated zero before; the column existed, so there is no migration). No size law is built, so `Dc` is a stand-in (LB1, T9, tweak L3): log-uniform between 100 and 1000 km, drawn from the belt's own stream, so the run's draws do not move; GEN.1000 replaces it with a draw from the size law. At that range a tenth of a decade of `Dc` is a factor 1.3 in mass and the shares come out 38 percent typical, 34 sparse and 28 dense over 300 belts, about the thirds of the old uniform draw. A belt's page lists its mass. The sectors' summed mass (`sector_stats.total_mass_kg`) now includes the belts' real masses for new belts; belts stored before keep their stated zero. The composition thresholds still take the Dirichlet draw's shares (the mass is not an input there).

## 15. The meteorite groups as outputs (GEN.475, lock-free step)

`small_body_composition.TAXON_DATA` takes each taxon's mix from a meteorite analogue (CI and CM, ordinary, CV and so on): the groups are inputs. `physics/meteorite_groups.py` holds them as points of a space instead, so a derived body is given a group as a label (registry MG1 and MG2, section 6.44):

- **The axes** are four rows of [R4658: Table 1]: the olivine's mol percent fayalite (the oxidation state, `fayalite_mol_percent` from the body's FeO and MgO), and the volume percent of metal, of refractory inclusions (CAIs and amoeboid olivine aggregates) and of matrix (the fine grains that carry water and organics).
- **The groups** are the table's averages for CI, CM, CO, CV, CK, CR, CH, CB, H, L, LL, EH, EL, R and K. `nearest_group` returns the closest by a scaled distance over the axes both have (CO and CV have only an olivine range, so they are placed by the other three), or `none` past `NEAREST_LIMIT`: a metal-rich, ice-rich or differentiated body is not forced into a chondrite group.
- **What it cannot do:** EH and EL have the same coordinates here (they differ in the table by sulfide and silicon in the metal, which are not axes), and CV against CK differ only by oxidation. The oxygen isotopes, chondrule sizes and petrologic types are not in the space.

**Wired (GEN.475, train 7).** `small_body_composition.chondrite_group(components)` reads the four axes off the grains of a belt or field (`group_axes`: the olivine's fayalite against forsterite in mol percent, then the volume percent of the metal rows, the refractory rows (spinel, anorthite) and the matrix rows (serpentine, magnetite, organics, calcite)), asks `nearest_group`, and calls the mix `none` when it holds `ICE_MAX` (5 percent) of ice or water, since the groups are rocks (MG3, T9, tweak L1). The draw is unchanged and takes no random number: the label is written into the makeup's detail (`group`, `group_distance`), kept on `AsteroidBelt.chondrite_group` and `AsteroidField.chondrite_group`, stored in `asteroid_belts.chondrite_group` and `asteroid_fields.chondrite_group` (schema migration 0187, NULL on a row made before, which the sheet leaves out) and shown as "Nearest Meteorite Group". The map does not give the within-family shares of the taxa (TK6 stays a stand-in, now `no data` in the burn-down register): it names what a derived mix looks like, it does not yet build the mix.
