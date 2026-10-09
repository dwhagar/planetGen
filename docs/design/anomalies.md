# Anomalies for the starmap

Which anomalies from the two Star Trek documents Boss uploaded are worth adding to the starmap, each with its real-science basis, rate, placement and how it is drawn, plus the conversion from those documents' 20 light-year sectors to this program's 4 pc sectors and the errors found in them. The two sources are `Star Trek Anomalies and Science.md` and `Star Trek Anomaly Probabilities.md` in this folder (not edited); the asteroid field and belt render plan is in [nebula-and-asteroid-field-classes.md](nebula-and-asteroid-field-classes.md).

Informs: GEN.113, GEN.114, GEN.103, GEN.47, MAP.132, GEN.104, GEN.128, GEN.129, GEN.130

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [S] seen in a search result (URLs under Sources), [C] computed, [R] recalled and unconfirmed. The research environment could read search-result text only, not the papers; the [R] items are listed under Evidence notes.

## Decisions already taken

- Boss (2026-10-07 11:47Z): "Two additional files have been added, but these are part of adding more anomalies to the starmap so it'll be probably Phase 2 and 3 just to analyze those." GEN.113 asks for a list of anomalies worth adding, each with real-science basis, rate, placement and how it is drawn, reviewed by Boss, noting the documents' 20 ly sectors against our 4 pc (about 13 ly) sectors. GEN.114 adds the chosen ones, one kind per PR.
- Nothing has been chosen yet. The tiers below are a proposal.

## Summary

- The documents' rates are per 20 ly cube (8,000 ly^3). Our sector is a 4 pc cell: 13.05 ly edge, 64 pc^3 = 2,220 ly^3, which is 0.2776 of the documents' volume [C]. Point-like objects (black holes, magnetars, anything much smaller than a sector) take the documents' probability times 0.2776; extended objects (nebulae, bubbles, radius above about 10 ly) have a per-sector probability equal to the volume filling factor, which does not change with sector size. The documents mix the two. Details in "From 20 ly sectors to 4 pc sectors".
- Real science covers only a few of the classes, and fewer fit the data model cheaply. Everything in the "subspace" family (rifts, tears, interphase, sandbars, Omega, anti-time, implosion rings) has no physical basis and is a lore decision for Boss.
- Tier 1 (add now, real, cheap): magnetars as a subtype of the existing neutron star, about 30 to 300 active in the galaxy, placed with the young population; Einstein radius shown on black hole, neutron star and rogue planet pages; correcting the pulsar fraction of neutron stars, now 70% (`NEUTRON_STAR_PULSAR_CHANCE`) against about 1e-4 to 1e-3 real.
- Tier 2: Wolf-Rayet stars with their ring bubbles, a superbubble ("Local Bubble") cavity field beside `nebula_field`, drawn jets after GEN.104 gives spin axes, symbiotic binaries inside GEN.128 to GEN.130, a derived radiation-hazard score, and a Fermi-bubble overlay on the Galaxy Map.
- Tier 3 (skip, or a lore decision): all subspace anomalies, cosmic strings ("quantum filaments"), wormholes, the Dark Matter Anomaly as an object, dark matter nebulae, Thorne-Zytkow objects, quark stars, red novae, gamma-ray bursts, fast radio bursts, the Great Attractor and the Bootes void.
- Rate flag outside the documents: GEN.47's molecular cloud field fills 3% to 46% of volume against 0.5% to 1% observed (see [nebula-and-asteroid-field-classes.md](nebula-and-asteroid-field-classes.md), "The molecular cloud field").

## What to carry over from the two documents, and what is superseded

Still valid, and used below:

- The Poisson form `P = 1 - exp(-lambda)` with `lambda = V * rho * K` (Probabilities, "Sector Baseline"). The project already uses it per star (`phenomenon_rate_per_star`, [interstellar-object-rates.md](interstellar-object-rates.md)).
- The split between "pure scientific" and "lore-adjusted" probability, and the idea of a per-kind multiplier (kept as a per-kind `K_lore`, default 0, for fiction kinds only).
- Latitude scaling of gas-related anomalies with distance from the plane (Probabilities, "Sector Modifier Modulations"): kept in spirit, but the project's disk density model ([galaxy-disk-density.md](galaxy-disk-density.md)) and `nebula_field.gas_factor` already do it.
- Physics sections: the Alcubierre energy density is negative, so the weak energy condition is violated [R]; the straight cosmic string metric `ds^2 = -c^2 dt^2 + dz^2 + dr^2 + (1 - 4G mu/c^2)^2 r^2 dphi^2` with deficit angle `8 pi G mu / c^2`, no Newtonian pull and double images without magnification; the Kerr ergosphere radius `r = GM/c^2 + sqrt((GM/c^2)^2 - a^2 cos^2 theta)` and tidal force across a body `F ~ 2 G M m L / r^3`; Beer-Lambert sensor extinction `I = I0 exp(-integral sigma rho ds)`. These are correct standard results. The extinction law is what `NEBULA_CLASSES` A_V ranges already express; the string and Kerr formulas are background for the Tier 3 notes.
- Real-world nebula types in the catalog table: emission nebulae are lit by hot young stars and ionized, dust clouds absorb, which maps onto classes C to G and M to Q. "Mutara" and "Paulson" style dust clouds map to dark classes M to N with a high A_V; "Class 11 / sirillium" and "ionite" are lore composition tags, not placement logic (see Tier 3). Crab Nebula, Veil Nebula and Coalsack are listed in the document among emission and stellar-nursery nebulae, but the Crab is a supernova remnant with a pulsar wind nebula (class T or U), the Veil an old remnant (class V) and the Coalsack a dark cloud (class N) [R].

Superseded or not adopted:

- The 20 ly sector and its "6 to 10 star systems": replaced by the 4 pc sector (see errors).
- A global narrative multiplier `K_lore` on real kinds: not adopted (below).
- The path model `lambda_path = A_sensor L rho_sector`: not adopted; course avoidance is NAV.51 and keep-out, and the sensor-sweep area depends on a sensor model that does not exist. A 1 to 3 ly sensor radius sweeps only 1.8% to 17% of a 13 ly sector's volume on one traverse [C], so it could not use a per-sector rate in any case.
- The game-engine blueprint (voxel grid, vector field engine, ship physics layer, shear damage, FTL speed cap, sirillium ignition cascade, the translation matrix): the program is a generator and map with course planning, with no ship damage or sensor-lock simulation. Nothing in the repository reads those formulas, so none is adopted.
- The probability table's numbers: replaced by the table in "Phenomenon by phenomenon".

## Errors found in the two documents

| Claim | Check | Verdict |
|---|---|---|
| A sector is 8,000 ly^3 and holds 6 to 10 star systems | 8,000 ly^3 = 230.6 pc^3; at 0.10 to 0.14 stars/pc^3 that is 23 to 32 stars [C] | Wrong by 3 to 4 times. 6 to 10 matches our 64 pc^3 sector (6.4 to 9 stars). |
| Black holes: 1e8 in 1e12 ly^3, "0.08%" per sector | 1e8 * 8,000 / 1e12 = 0.8, not 0.0008 [C]. The real disk is about 7e11 pc^3 = 2.4e13 ly^3, so 3.3% per 20 ly sector at the galaxy average, 0.9% per 4 pc sector [C] | Arithmetic and volume both wrong. The "lore" 1 to 2% is not an inflation of the science; the generator's 0.9% per 4 pc sector (`black-hole` 1.4e-4 pc^-3 * 64) is already above the lore value once converted (1 to 2% * 0.2776 = 0.28 to 0.55%). |
| Gas nebulae fill 1 to 2% of the disk | Molecular clouds: about 0.5% (inner Galaxy), about 1% (disk, from 8,107 clouds of 30 pc radius) [S/C]. Hot ionized gas fills more but is not a nebula | Upper end optimistic for clouds; fine for all nebulae only if H II and planetary are added, which are far smaller. |
| Ion and plasma storms, f_V 0.2 to 0.3 | The warm ionized medium fills roughly 0.2 to 0.5 of the disk [R] | The volume is plausible but it is neither a storm nor a hazard; no ship feels the warm ionized medium. |
| Cosmic string mu about 1e18 kg/m | G mu / c^2 = 7.4e-10 [C]. PPTA limit 5.1e-10, NANOGrav-based 5.3e-11 [S]. GUT-scale strings (G mu 1e-6, mu 1.3e21 kg/m [C]) are excluded 1e4 times over | Already excluded by pulsar timing. |
| Quantum filament "shears hulls by tidal shear from the deficit angle" | A straight cosmic string is locally flat: the deficit angle gives lensing and a conical geometry, no tidal force [R, follows from the metric the document writes down] | Physically wrong; at most a lore hazard. |
| String "hundreds of metres long" | With that mu a 100 m loop radiates its energy in about 9 s (Gamma = 50 [R]) and would have a mass of 1e20 kg [C] | Cannot be a persistent object. |
| DMA ergosphere "up to 12 AU" shears planets and ships | An ergosphere of 12 AU needs M = 6.1e8 solar masses, 140 times Sgr A* (4.3e6 [R]); tidal acceleration across a 100 m ship at 12 AU is 2.8e-6 m/s^2 [C]. Tidal damage at the horizon grows as 1/M^2 | Backwards: a supermassive hole is gentle. Stellar-mass holes are the dangerous ones: 10 solar masses gives 2.7e5 m/s^2 across 100 m at 1,000 km [C]. |
| Alcubierre energy density written with `rho` for both density and radius | In the standard form the factor is `(y^2 + z^2)`, not `rho^2` [R] | Notation clash only. |
| Lore table applied per sector | Over 2.1e10 sectors in the default galaxy [S: CHANGELOG] the lore rates imply 3 to 5 billion nebulae, 2 to 3 billion plasma storms, 1e8 to 6e8 rifts and 2e7 to 1e8 topological defects in one galaxy [C] | The lore column is a per-encounter number, not a population. Lazy generation hides the totals. |

## From 20 ly sectors to 4 pc sectors

- Volume ratio: (4 pc = 13.046 ly)^3 / 20^3 = 2,220.5 / 8,000 = **0.2776** [C]. Inverse 3.60.
- Point-like or small objects (radius under about 1 ly: black holes, magnetars, neutron stars, Wolf-Rayet stars, symbiotic systems, jet sources): the Poisson mean scales with volume, `lambda_4pc = 0.2776 * lambda_20ly`, `P = 1 - exp(-lambda)`.
- Extended objects (radius at least a sector, 10 ly and up: nebulae, superbubbles, any "region" anomaly): the chance a random point is inside is the filling factor f_V whatever the sector size. A sector touches an object if its centre is within `r + half the cell diagonal` (`_db.sectors_reached_by`, the test `nebula_field` uses): `P_touch = n * (4/3) pi (r + 0.87 e)^3` with e = 4 pc. For r above 20 pc this is about `f_V (1 + 0.13 e/r)^3`, so take it as f_V.
- Per-star scaling is the project's own method: a rate quoted per 20 ly sector is first turned into a density per pc^3 (divide by 230.6 pc^3) and then into a per-star rate (divide by 0.14). Never multiply by 0.2776 and add the result to a per-star rate.
- Lore multiplier: if Boss wants one, make it a per-kind multiplier on Tier 3 kinds only, default 0, shown in the admin report like `PHENOMENON_RATE_SCALE`. Do not apply a global `K_lore` to real kinds: after conversion the real rates are at or above the documents' lore values (black holes 0.9% against 0.28 to 0.55%; neutron stars 4.5% per sector; terrestrial rogues about 50 per sector).

Converted rows of the documents' probability table:

| Documents' class | Documents' pure-science probability | Per 4 pc sector | Note |
|---|---|---|---|
| Gas nebulae (emission, dark, planetary) | 1 to 2% (f_V) | 1 to 2% (f_V, unchanged) | the project's dark-cloud field is at 13% on arms; observed about 1% |
| Volatile chemical nebulae | about 0.1% (f_V) | unchanged | a composition tag, no placement |
| Ion and plasma storms | 20 to 30% | not an object | the warm ionized medium; see radiation hazard |
| Singularities (black holes) | "0.08%" (computed wrongly) | 0.9% at the galaxy average, 0.9% at local density | generator already at this value |
| Subspace rifts and tears | 0% | 0 | lore; depends on warp traffic the project does not model |
| Topological defects | under 0.0001% | about 1e-17 pc of string per sector [C] | none |
| Catastrophic (Omega, anti-time) | 0% | 0 | lore |

## Phenomenon by phenomenon

Counts are whole-galaxy; "per sector" is the mean over the 2.1e10 qualifying sectors [S: CHANGELOG] using the 7e11 pc^3 disk of [interstellar-object-rates.md](interstellar-object-rates.md) [C]. Placement rules decide where a real kind actually goes.

| Phenomenon (docs' class) | Real basis | Count in the Milky Way | Per sector (mean) | Fits a 4 pc sector? |
|---|---|---|---|---|
| **Magnetar** (not in the documents; nearest to "particle fountain" and FRB lore) | Neutron stars with B about 1e14 to 1e15 G, spin 2 to 12 s, active about 1e4 yr [R]. Birth rate 0.1 to 6 per century [S] | 17 expected (Gill and Heyl) up to a few hundred; about 30 known [R] | 1e-9 to 1e-8 [C] | Yes, a point object |
| **Radio pulsar** | Spinning beamed neutron stars; density model n ~ R^2.35 exp(-R/1530 pc) exp(-abs(z)/300 pc) [S]; beaming about 21% for young ones [S] | tens of thousands beamed; young about 3,400 [S]; total active about 1e5 to 1.6e5 [R] | 2e-6 to 8e-6 [C] | Yes |
| **Wormhole** ("spatial displacement") | Morris-Thorne: needs exotic matter violating the null energy condition [S]; no detection | 0 observed | n/a | n/a |
| **Cosmic string** ("quantum filament") | Not observed; G mu below 1.3e-7 (CMB), 5.1e-10 (PPTA), 5.3e-11 (NANOGrav-based) [S] | expected string length in the Galaxy about 4e-7 pc [C] | about 1e-17 pc of string [C] | no map object |
| **Void** ("cosmic void boundaries") | Bootes void about 330 Mly across [S]; voids hold about 40% of volume [S] | n/a inside one galaxy | n/a | no; the galaxy-scale analogue is the Local Bubble, 100 to 200 pc, n_e about 0.005 cm^-3, hot [S]; simulated superbubble filling about 12% [S] |
| **Dark matter halo / "nebula"** | Local density 0.4 to 0.6 GeV/cm^3 = 0.011 to 0.016 Msun/pc^3 [S] | about 1e12 Msun [R] | 0.5 to 1.0 Msun per 4 pc sector [C] | smooth background only; collisionless, no drag |
| **Gravitational lens / Einstein ring** | `R_E = sqrt(4GM/c^2 * D_L D_LS / D_S)`: 10 Msun at 100 pc, source at 8 kpc: 2.8 AU [C]; 1 Msun at 2 pc: 0.13 AU [C]; bulge optical depth about 1.75e-6 [S] | n/a | n/a | an effect of existing objects, not an object |
| **Stellar merger / red nova** | Merger rate 0.5/yr (M_V < -3), 0.1/yr (brighter than V1309 Sco), 0.03/yr (V838 Mon) [S, Kochanek 2014]; ZTF finds the bright end 5 to 100 times lower [S] | transient | about 5e-12 per sector per year [C] | a transient |
| **Symbiotic star** | White dwarf with red giant; about 300 known [S]; estimates 1.2e3 to 1.5e4, 5.3e4, 3e5 [S] | 1e3 to 3e5 | 5e-8 to 1.4e-5 | Yes, a binary subtype |
| **Quasar / AGN jet** | The existing quasar kind; Fermi bubbles about 50 degrees, about 8 kpc [R] | 1 | n/a | exists at ring 0 |
| **Microquasar / particle fountain** | SS 433: precessing jets at 0.26 c, 162.3 d period [S]; several dozen known, perhaps 100 to 1,000 [S/R] | 1e2 to 1e3 | 5e-9 to 5e-8 [C] | jet 0.1 to 10s of pc spans sectors |
| **Protostellar jet / Herbig-Haro** | Up to 150,000 in the Galaxy (estimate) [S] | about 1e5 | 5e-6 [C] | Yes, inside class Q cores (already named in class Q) |
| **Wolf-Rayet star and bubble** | 642 known, 1,200 +/- 100 expected [S]; ring nebulae 1 to 70 pc across, NGC 6888 7.6 x 5.0 pc [S] | 1,200 | 5.7e-8 [C] | star yes; bubble 2 to 35 pc radius spans sectors |
| **Thorne-Zytkow object** | Neutron star inside a supergiant; predicted 20 to 200 [S]; HV 2112 disputed, none confirmed | 20 to 200 | 1e-9 to 1e-8 | as a star |
| **Quark / strange star** | Hypothetical, no candidate [R] | unknown | n/a | n/a |
| **Fast radio burst** | Galactic magnetar SGR 1935+2154 burst in 2020; 0.0036 to 0.8 per magnetar per year [S] | transient | n/a | a flag on a magnetar |
| **GRB** | Long GRB local rate 0.7 to 1.3 Gpc^-3 yr^-1; Milky Way 1e-6 to 3e-4 per year [S] | transient | about 1e-14 per sector per year [C] | no |
| **Supernova precursor / LBV** | Eta Carinae-like events about 0.094 (0.04 to 0.21) of the core-collapse rate [S]; Galactic core-collapse rate 1.9 +/- 1.1 per century [S] | transient | n/a | a star-type flag |
| **Great Attractor** | About 150 to 250 Mly away [R] | n/a | n/a | no; at most a line of flavour text |
| **Intergalactic stars** | Unbound stars exist only in clusters; the nearest analogues are the halo and hypervelocity stars (partly modelled) | halo about 1e9 Msun [R] | n/a | yes, halo above the plane |
| **Rogue objects** | Built: rogue planets, brown dwarfs, comets | see [interstellar-object-rates.md](interstellar-object-rates.md) | | |
| **Ion / plasma / proton storms** | Stellar flares and CMEs, the warm ionized medium, supernova shock fronts, cosmic-ray overdensities; no "storm" at interstellar scale | n/a | n/a | a derived hazard of existing objects |

Notes on the less obvious rows.

- **Magnetar numbers.** Active count equals birth rate times active lifetime: 0.1 to 6 per century [S] and about 1e4 yr [R] give 10 to 600, matching the 17-object figure at the low end. A magnetar older than about 1e5 yr has decayed to an ordinary high-field neutron star. A generator rule that needs no new table: a neutron star younger than 1e4 yr with a 10% birth probability (1 to 10% [S]; a 2026 preprint's "about half" is provisional [S]) is an active magnetar.
- **Pulsar fraction.** The project has about 5e8 isolated neutron stars (7e-4 pc^-3 * 7e11 pc^3) and `NEUTRON_STAR_PULSAR_CHANCE = 0.7`. With about 1e5 beamed radio pulsars [R/S] the active fraction is about 1e-4 to 1e-3 [C]. Radio activity ends after about 1e7 yr [R]; isolated neutron stars average a few Gyr. The 30% millisecond share (`PULSAR_MILLISECOND_CHANCE`) needs accretion from a companion, which an isolated star lacks (recycled pulsars that later lose the companion exist [R]). Suggested rule: pulsing probability `exp(-age / 1e7 yr)` with a floor for recycled pulsars in binaries (GEN.128 on).
- **Gravitational lensing** is not an object; the physical Einstein radius of every compact object is a number the pages can show (1 to 10 AU for a stellar-mass hole at 100 pc to 1 kpc [C]). A ring cannot be drawn at sector scale because the angular radius is below 0.5 arcsec for most sources.
- **Voids.** The existing zero-density sectors (GEN.102) are the galaxy's edge, not voids.
- **Dark matter nebulae.** Dark matter does not clump on 10 ly scales and exerts no drag [R]. A 100 ly "nebula" at local density holds about 1,200 Msun over 1.2e5 pc^3 [C], which is the background.

## Per-kind parameters, drawing and effects (Tier 1 and 2)

### Magnetar (Tier 1)

- Data: a subtype of `neutron-star` (`phenomenon_scatter` rows already carry `subtype`; `SCATTERED_KINDS` appends new kinds last so earlier draws stay unchanged). B log-uniform 1e14 to 1e15 G (existing young pulsars use 1e11 to 1e13, `PULSAR_MAGNETIC_FIELD_GAUSS_RANGE_YOUNG`), P 2 to 12 s, age 1e3 to 1e4 yr, a supernova remnant companion row when younger than about 3e4 yr, a cosmetic `fast_radio_burst_active` flag on a few percent.
- Rate: about 100 over the galaxy (range 30 to 300), weighted by the young-star density with arm contrast (`nebula_field.gas_factor` is the nearest existing weight), not by total density. Per star that is about 3e-10 [C], so `phenomenon_rate_per_star` needs a young-population weight; the generic per-star form would put them in the old bulge.
- Drawing: the Sector Map reuses the `neutronStar` recipe in `sectorscene.js` with a thin outer ring (`#9fb8ff`; `#b59cff` is taken by rogue planets) and glow `{power 1.0, strength 2.6, scale 1.5}`; the page view in `phenomenonrender.js` is the neutron star with a dipole field-line pair; the Galaxy Map gets a MAP.132 marker glyph.
- Effect: a giant flare of about 1e46 erg delivers 8.4e7 erg/cm^2 at 1 pc, 8.4e5 at 10 pc, 8.4e3 at 100 pc, 84 at 1 kpc and about 1 at 8.7 kpc [C], matching the 2004 SGR 1806-20 flare reaching Earth at about 1 erg/cm^2 [R]. The damage threshold for an atmosphere is unconfirmed, so keep it a flag ("within 10 pc of an active magnetar") feeding the habitability radiation term, not a ship hazard.

### Einstein radius numbers (Tier 1)

Derived from mass and distance, nothing stored. Show "Einstein radius against the galactic centre (AU)" on black hole, neutron star, rogue planet and quasar pages: `R_E = sqrt(4GM/c^2 * d * (1 - d/D_S))`.

### Wolf-Rayet bubble (Tier 2)

- Needs a W spectral class in the star model (`SPECTRAL_PROBABILITIES_LARGE_STAR` has O to M only). 1,200 galaxy-wide [S], 5.7e-8 per sector [C]; placement in massive clusters and arms, scale height about 50 to 70 pc [R].
- Bubble: nebula family `emission`, radius 2 to 35 pc (diameter 1 to 70 pc [S]), a thin shell with a hollow interior drawn with the existing `nebula_shape` mesh (shell profile from the appearance table in [nebula-and-asteroid-field-classes.md](nebula-and-asteroid-field-classes.md)); the host rule goes into `NEBULA_HOST_RULES`.

### Superbubble / Local Bubble cavity (Tier 2)

A seeded field like `nebula_field`, cell 50 pc, radius 50 to 500 pc; inside, the dark-cloud probability is zero and an "ionised, 1e6 K" tag applies. Fill is about 10% of the disk [S, simulation]; the Sun's real position inside such a bubble [S] is a calibration point. Drawn as a faint thin shell on the Galaxy Map only. It also gives the lever for the over-full cloud field (zero clouds inside).

### Jets (Tier 2)

Protostellar (class Q cores): length 0.05 to 1 pc. Microquasar: 1 to 100 pc, SS 433's 0.26 c and 162 d precession [S]. Pulsar jets. Parameters need a spin axis (GEN.104). Drawn as two opposed elongated glow sprites along the axis; no gameplay effect beyond a keep-out cone.

### Symbiotic binary (Tier 2)

A subtype of the wide-binary or compact-primary designs (GEN.128 to GEN.130): white dwarf with red giant at 1 to 10 AU, period 200 to 1,000 days, an outburst flag, a small circumbinary nebula. Rate per sector 5e-8 to 1.4e-5; use 1e-6 per system and flag the range as unsettled.

### Radiation hazard (Tier 2)

Not an object: a number per sector and system from flare-star class, distance to OB and Wolf-Rayet stars, supernova remnants and active magnetars. It replaces the documents' "ion storms" and ties to `Naturally Occuring Ionizing Radiation.md`. It is also the clean place for MAP.132's marker overlay (same budget and toggle).

## Tier 3 and the lore question

What remains is lore fiction: spatial implosions, distortion rings, subspace rifts, tears, fissures, interphase pockets, sandbars, Omega destabilization, anti-time eruptions, graviton ellipses, "particle fountains" in the documents' sense, the DMA, Class 11 sirillium nebulae, ionite and dichromic nebulae, dark matter nebulae. If Boss wants any:

- Do not call them physics. One new `lore-anomaly` kind with a `subtype`, per-kind `K_lore` default 0, off in the Generate page unless enabled, flagged "fiction" on the page.
- Rifts and tears depend on warp traffic and former war zones; the project has no traffic or war model, so any placement rule is arbitrary.
- The DMA as an object duplicates the quasar and black hole kinds; a `subtype` `artificial` on `black-hole` would do.
- Sirillium and similar chemistry could be a composition tag on existing molecular cloud classes with no placement logic.

## Placement by population (the GEN.103 rule set)

| Population | Follows | Kinds |
|---|---|---|
| Young thin disk and arms | young-star density with arm contrast (`nebula_field.gas_factor`) | O, B and Wolf-Rayet stars, magnetars, giant molecular clouds, supernova remnants, class Q cores, protostellar jets |
| Thin disk, kick-inflated scale height | thin-disk density, older objects thicker | pulsars, isolated neutron stars, open clusters (about a thousand known [R]) |
| Bulge and thick disk | bulge and thick-disk density | symbiotic stars, white dwarfs, low-mass X-ray binaries, old binaries |
| Halo | halo profile | globular clusters (about 150 [R]), the stellar halo (about 1e9 Msun [R]) |
| Between arms, where supernovae are dense | | Local Bubble style cavities |

The `phenomenon_scatter` draw is weighted by total stellar density per ring and angle bin (`kind_mean`). Young-population kinds need the young-star density; old-population kinds need the bulge and thick-disk density.

## How it fits the code

- New objects go through `tuning.PHENOMENON_DENSITY_PC3`, `phenomenon_scatter.SCATTERED_KINDS` (append only), `objectId.KIND_CODES` (codes 12 to 63 are free), the `queryDb` readers, the class pages (`classref`) and the galaxy-map marker set. A subtype avoids a new kind and a new ID code: magnetars, symbiotic stars and the artificial black hole can all be subtypes.
- MAP.132 already plans markers for black holes, nebulae and habitable worlds with a per-kind toggle (MAP.123) inside MAP.116's budget; the new kinds join that list instead of a second overlay system.

## Ranked list and PR order (one PR per kind)

| # | Tier | Item | Why | Depends on |
|---|---|---|---|---|
| 1 | 1 | Magnetar subtype of neutron star and the age-dependent pulsar fraction | real, about 100 objects, existing table and scatter | a young-population weight for the scatter |
| 2 | 1 | Einstein radius on compact-object pages | display only, no storage | none |
| 3 | 2 | Wolf-Rayet class and ring bubble | real, spans sectors through the nebula machinery | W class in the star model |
| 4 | 2 | Superbubble cavity field | calibrates GEN.47; zero clouds inside | the `nebula_field` pattern |
| 5 | 2 | Jets as drawn features | cheap draw | GEN.104 spin axis |
| 6 | 2 | Symbiotic binaries | real binary subtype | GEN.128 to GEN.130 |
| 7 | 2 | Radiation-hazard score and MAP.132 markers for the new kinds | one derived field, one marker set | items 1, 3 |
| 8 | 2 | Fermi-bubble overlay on the Galaxy Map | draw only | none |
| 9 | 3 | `lore-anomaly` kind, `K_lore` 0 by default | only if Boss wants lore | Boss decision |

## Evidence notes

[R] claims to check once arXiv and Wikipedia are reachable:

- About 30 known magnetars, the 1e4 yr active lifetime, B 1e14 to 1e15 G, P 2 to 12 s, three Galactic giant flares since 1979, 1e46 erg flare energy, one magnetar near Sgr A*.
- Pulsar death after about 1e7 yr, 1e5 to 1.6e5 total radio pulsars (ATNF count not confirmed by search), isolated millisecond pulsars.
- Sgr A* mass 4.3e6 Msun; Gamma = 50 for loop decay; the cosmic-string scaling network (xi = 10, arithmetic of the research draft).
- Warm ionized medium filling 0.2 to 0.5; stellar halo about 1e9 Msun; about 150 globular clusters; about 1,300 open clusters; Fermi bubbles about 50 degrees and 8 kpc; Great Attractor distance; dark matter halo total; the giant flare fluence at Earth.
- The Alcubierre `(y^2 + z^2)` form, and the Crab, Veil and Coalsack classifications.
- Pinned by computation only: the 0.2776 ratio, the black hole arithmetic, GMC filling in the project, Einstein radii, flare fluences.

## Sources

Search results used (WebFetch was blocked; [S] claims rest on search-result text):

- Magnetar birth rate and fraction: https://arxiv.org/pdf/astro-ph/0703346, https://arxiv.org/pdf/0807.2106, https://arxiv.org/pdf/1903.06718, https://arxiv.org/html/2601.16159v1
- Cloud filling and catalogues: https://arxiv.org/pdf/1301.3905, https://arxiv.org/pdf/astro-ph/0701877, https://arxiv.org/abs/1602.02791, https://dc.g-vo.org/rr/q/lp/custom/CDS.VizieR/J/ApJ/834/57
- Red novae: https://web3.arxiv.org/abs/1405.1042, https://arxiv.org/pdf/2211.05141, https://arxiv.org/pdf/1906.00812
- Symbiotic stars: https://arxiv.org/pdf/2301.08201, https://arxiv.org/pdf/astro-ph/0208085, https://arxiv.org/html/2504.02090v2
- Wolf-Rayet: https://arxiv.org/abs/1412.0699, https://arxiv.org/pdf/astro-ph/0003053, https://arxiv.org/pdf/0909.0621
- Thorne-Zytkow: https://en.wikipedia.org/wiki/Thorne%E2%80%93%C5%BBytkow_object, https://arxiv.org/pdf/1601.05455, https://arxiv.org/abs/2407.11680v1
- Cosmic strings: https://arxiv.org/pdf/2205.07194, https://arxiv.org/pdf/1907.04960, https://arxiv.org/pdf/1206.2924
- Wormholes: https://arxiv.org/pdf/0807.2774, https://arxiv.org/pdf/1009.6084, https://arxiv.org/pdf/2405.05476
- GRB rates: https://arxiv.org/pdf/1609.09355, https://arxiv.org/pdf/astro-ph/0610043, https://arxiv.org/pdf/astro-ph/0310667, https://arxiv.org/pdf/astro-ph/0701748
- Supernova rate: https://arxiv.org/pdf/astro-ph/0603669, https://arxiv.org/pdf/1110.4105
- FRBs and magnetars: https://arxiv.org/pdf/2005.05283, https://openaccess.inaf.it/handle/20.500.12386/35137
- Microlensing: https://arxiv.org/pdf/1906.02210, https://iac.es/es/ciencia-y-tecnologia/publicaciones/microlensing-event-rate-and-optical-depth-moa-ii-9-yr-survey-toward-galactic-bulge
- Local Bubble and superbubbles: https://arxiv.org/pdf/0812.0505, https://arxiv.org/pdf/1010.5654, https://arxiv.org/pdf/2403.12135
- Dark matter density: https://arxiv.org/pdf/2012.11477, https://arxiv.org/abs/1404.1938
- Microquasars and jets: https://arxiv.org/pdf/astro-ph/0506008, https://arxiv.org/html/2506.01106v1; Herbig-Haro estimate via the HH-30 Wikipedia page
- Pulsars: https://arxiv.org/pdf/astro-ph/0412641, https://arxiv.org/pdf/2510.16618, https://arxiv.org/pdf/2205.05200
- Eta Carinae and LBV rates: https://arxiv.org/pdf/1010.3718, https://arxiv.org/abs/0905.3338
- Voids: https://en.wikipedia.org/wiki/Bo%C3%B6tes_Void, https://ar5iv.labs.arxiv.org/html/astro-ph/0702257

Project files read: `Star Trek Anomalies and Science.md`, `Star Trek Anomaly Probabilities.md`, [nebula-and-asteroid-field-classes.md](nebula-and-asteroid-field-classes.md), [interstellar-object-rates.md](interstellar-object-rates.md), [galaxy-coordinate-system.md](galaxy-coordinate-system.md), [galaxy-disk-density.md](galaxy-disk-density.md), [object-ids.md](object-ids.md), [reproducible-galaxies.md](reproducible-galaxies.md); `tuning.py`, `generation/phenomenon_scatter.py`, `galaxy/nebula_field.py`, `generation/phenomena/compact_remnant.py`.
