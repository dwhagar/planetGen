# Phenomena catalogue and the fill rules

One master list of every kind of object the generator has, has designed, or should have: what is observed, what the generator does now, whether the kind may have a star system and with what probability, what it must carry inside it, where it occurs in the sphere profile, whether it comes from a scatter or from a sector fill, and the fill-time rule that keeps every kind represented. It folds in the earlier notes (listed under "Notes folded in") and the whole-sphere design in [intergalactic-and-halo-contents.md](../galaxy/intergalactic-and-halo-contents.md).

Informs: GEN.100, GEN.113, GEN.114, GEN.128 to GEN.130, GEN.153 to GEN.155, GEN.159 to GEN.164, GEN.233 to GEN.249, and the new items in section 8.

Status: design, 2026-10-15; nothing here is built. The rules in section 1 are Boss's of 2026-10-15; the defaults chosen where they leave a gap are marked "Default" and listed in section 9.

Evidence tags: **[S]** seen in this repository or in a search-result summary (the papers themselves were not opened; the environment cannot reach arXiv or journals); **[C]** computed with the repository's own code in this session; **[R]** recalled and unconfirmed; **[I]** inferred, with the step shown. A figure keeps the tag of its weakest input. Times shown to Boss are US Eastern (EDT) only.

## Summary

- **Scatters never carry planets or moons.** Bright and mass stars, neutron stars, black holes, remnants, planetary nebulae and hypervelocity stars are scatter rows: a point, a seed and a kind, nothing around them. Nebulae need stars and phenomena, never planets.
- **A fill makes everything.** When a sector is filled, every object already in it that can have a system rolls for one, and every object a fill creates rolls for its own, with the observed probability of that kind (section 5). The roll uses the object's own seed, so the answer is the same whenever and in whatever order sectors are filled.
- **Today only stars get this.** The fill already builds a full system (companion, planets, moons) around each pre-placed bright star (`run_sector.py:515`), but a scattered black hole or neutron star becomes a bare phenomenon (`run_sector.py:449`, `:454`), and a system around one exists only through the `--anchor-system` option (`run_phenomenon.py:70`). That is the gap the system pass closes. The cost is tiny: about 0.002 extra systems per local sector, 0.03% of the 6.3 a sector holds [C].
- **One registry.** The kind list is now spread over five places (`SECTOR_PHENOMENON_KINDS`, `SCATTERED_KINDS`, `PHENOMENON_DENSITY_PC3`, `KIND_CODES`, the class pages). Section 6 puts every column of this catalogue in one table that the code reads, with a completeness test, so a kind cannot exist without a rate, a gate, a builder, an ID code and a page.
- **"All types represented" is audited, not forced.** Forcing every kind into every sector would distort the observed rates. Instead a statistical audit fills sample sectors in each kind of region and compares each kind's count with its expected count (section 7).
- **The catalogue** (section 4) has 57 entries in six families: 28 built, 20 designed and already filed in TODO (GEN.130, GEN.153 to GEN.164, GEN.233 to GEN.249), 7 new here, and 2 not objects at all; a further 11 kinds are deliberately skipped (section 4.7). Nine items are proposed in section 8.

## 1. The rules (Boss, 2026-10-15)

1. Scatters and backscatters (bright stars, massive stars, neutron stars, black holes and the rest) get no systems, no planets and no moons. Planets and moons are not generated when a sector is not being filled with everything.
2. Nebulae need stars and phenomena inside them, not planets.
3. When a sector is filled, or whenever a star system would be generated, all types of objects are represented. If the sector already holds phenomena that could have a system of planets, generate one if it gets one, with probabilities from real observations, like everything else.
4. When a new object that could have a system is generated during a sector fill, generate its system.
5. Go through all the notes about phenomena and rebuild the list and its generation: this note.
6. Standing rules that still hold: redesign freely (no backward compatibility), no change to anything but docs in this research thread, and TODO.md is edited only by the TODO thread.

## 2. How a fill works today [S]

`generate_sector` (`run_sector.py:531`) runs, for a galaxy-placed sector, in this order:

| Step | What | Where | Planets and moons? |
|---|---|---|---|
| 1 | every pre-placed bright or massive star gets a full system, built from the row's own seed, at its stored point | `_add_preplaced_systems`, `run_sector.py:515` | yes: companion, planets, moons |
| 2 | every scattered phenomenon row becomes its object: black hole, neutron star, supernova remnant, planetary nebula (with its own new central white dwarf system), hypervelocity star (an ordinary young system), the galactic nucleus | `add_scattered_phenomena`, `run_sector.py:402`; `_build_scattered`, `:433` | only the planetary-nebula white dwarf and the hypervelocity star are systems; black holes and neutron stars are bare |
| 3 | the sector's other systems, `--density` many, ages from the cell's population mix, below the bright-star threshold | `generate_sector`, `run_sector.py:531` onward | yes |
| 4 | the per-star Poisson rolls for kinds that are not scattered: rogue planets, brown dwarfs, comets, asteroid fields (rate 0), and in a stand-alone sector all of them | `generate_sector_phenomena`, `run_sector.py:146`, over `SECTOR_PHENOMENON_KINDS` (`:50`) | a rogue planet keeps a moon with probability 0.2 (`SYSTEM_PASS_ODDS`) |
| 5 | molecular clouds come from the seeded `nebula_field`, not from a roll | `generate_sector_phenomena` | none |
| 6 | emission and reflection nebulae grow around the sector's own hot stars | `add_star_hosted_nebulae`, `run_sector.py:262` | none |
| 7 | neutron stars and black holes below the scatter's mass cut, rolled in the sector | `add_below_cut_remnants`, `:463` | none |
| 8 | runaway and hypervelocity flags | `flag_fast_stars`, `:302` | n/a |
| 9 | the galactic nucleus, in one sector only | `add_galactic_nucleus`, `:357` | none (a quiescent supermassive hole, or a quasar with 10% chance) |

What is missing against the rules of section 1:

- **Step 2 gives no system to a black hole or neutron star.** The rules for their company are written (`multistar-and-compact-systems.md` 7.3 and 7.6, refined in `exotic-environments-planets-and-compact-binaries.md` 2.1) but GEN.129 and GEN.130 are open, and `StarSystem.__init__` skips binary generation whenever it is handed a remnant (`system.py:255`).
- **Steps 4 to 8 draw on the sector's one stream in a fixed order**, so adding a kind (say Wolf-Rayet nebulae) in the middle shifts every draw after it. The scattered objects avoid this by being built from their own seeds (`_seeded_build`, `run_sector.py:394`).
- **No step checks representation**: nothing says a kind with a rate and a builder actually appears at its rate.
- **The star census has no way to hold** a Wolf-Rayet star, a protostar, a brown dwarf companion, a pulsar-bearing binary or a symbiotic star. The model's living-star shares are 94.2% main sequence, 5.2% white dwarf, 0.31% subgiant, 0.27% giant, about 0 supergiant [C, 60,000 draws of `sample_living_star`], and 0.63% of all stars born have collapsed to a neutron star or black hole [C], which matched the 0.5% and 0.1% that `PHENOMENON_DENSITY_PC3` held before the GEN.100 retune (the retuned shares are 0.4% and 0.05%, section 5).

## 3. The fill, as it should run

Each stage draws on its own stream (`seed(sector) + stage name`), so adding a stage or a kind never moves another stage's draws. The new stages are marked **new**.

| Stage | What | Planets and moons |
|---|---|---|
| A | pre-placed stars: ordinary systems (as step 1 above) | yes |
| B | scattered rows become objects, each from its own seed (as step 2) | no |
| C **new** | **system pass**: every object from A and B whose kind has a system model rolls once, on its own seed, for a system (companions, planets, belts, moons, disk), and gets it (section 5) | by kind |
| D | the cell's other systems, the whole star-state spectrum of section 4.1 | yes |
| E | per-star kinds not scattered (rogue planets, brown dwarfs, comets) and the **system pass** for the ones that can have a moon or companion | by kind |
| F | field nebulae, then hosted nebulae, each **making its required stars as real systems inside its bounds** (GEN.241); a star made for a nebula has belts at most, as the planet rules for stars under 10 Myr already say (`system.py:1024`, `PLANET_MIN_STAR_AGE_GY = 0.01`) | belts only |
| G | below-cut remnants (same as B and C) | by kind |
| H | flags, then the nucleus | none |
| J **new** | the neutron stars and black holes the galaxy ejected (GEN.238), a count per kind from the sector profile; runs after G | none (a kick strips the planets) |
| I **new** | **audit hooks**: counters per kind per stage, read by the audit of section 7 | n/a |

Rules inside the stages:

- **Never planets for a scatter row alone.** A scatter row stores a point and a seed. Nothing with a planet or moon is stored until its sector is filled. The pages of an unfilled sector show the row, not a system.
- **Order of the system pass inside C.** Stars first (they already are systems), then compact remnants, then the rest, so that a remnant's companion is placed with the Hill spheres of the pre-placed stars in view, as `add_scattered_phenomena` already does by placing scattered objects before the sector's own systems (`run_sector.py:402`).
- **No room.** A system that does not fit drops its outermost body first, then its companion; the object itself always stays (the same rule `generate_sector` uses for a full cube).
- **The same seed gives the same system** whatever order sectors are filled, because the roll uses `draw.bound(row.seed)` with a stage-tagged suffix, the way `_add_preplaced_systems` does.

### As built (GEN.251)

- **Streams.** `generation/fill_stages.py` names the stages (`pre-placed`, `scattered`, `system-pass`, `systems`, `per-star-kinds`, `nebulae`, `below-cut`, `ejecta`, `flags`, and `clouds`, GEN.357: what the cloud a placed system sits in does to it, one stream per placement) and draws one 64-bit base from the sector's stream when `generate_sector` begins; a stage binds `draw.Stream(f"{base}:{stage}")`, and a stage that holds several independent parts binds one stream per part (`per-star-kinds:<kind>`, `nebulae:<kind or hosted>`, `flags:nucleus`, `flags:fast-stars`, `systems:fallback`, `system-pass:below-cut`). A piece called by hand on a sector no fill began on keeps the stream it has. The same galaxy seed gives a sector different contents than it did on one stream; redrawn sectors under a directive differ from each other because each attempt draws its own base.
- **Stage I.** `sector.stage_counts` holds the objects each stage made by `(stage, kind)`: star systems in A and D, scattered rows by kind in B, kinds in E to G, the fast-star count in H, and in C the objects rolled (`key`), those that became systems (`key:system`) and those whose planet model put no planet in them (`key:empty`). The representation audit (section 7, item 4) reads them; `fill_stages.describe` prints them.
- **The system pass.** `generation/system_pass.py` rolls the objects of B and G (compact remnants first) on `draw.Stream(f"{seed}:system-pass")`; `seed` is the scatter row's seed, or for a below-cut remnant `f"{below_cut.seed}:{n}"`, `n` its place in the stage's draws, counted whether or not it found room. The odds are the rows of the host tables (`generation/phenomena/hosts.py`, GEN.552) by object key (`remnant-satellites`: `neutron-star:millisecond` 0.7% planets; `planet-satellites`: `rogue-planet` 20% moons, rolled when the planet is made). A hit builds `StarSystem(compact_remnant=...)` on the same stream and replaces the bare object with it at its point; the system's galactic Hill sphere is its primary's, the one the object was placed with, and the planet model keeps every body inside its orbit ceiling, so every system fits and the drop-the-outermost-body rule never fires. A hit whose planet model puts no planet in the system keeps the bare object. The pulsar type used is that of the object as built from the row's seed; GEN.264 reconciles it with the stored members. Companions, belts, disks and the other models are not rolled here: companions are drawn at scatter time (GEN.262) and built from the stored members at fill (GEN.264); a planetary nebula's circumbinary disk and a pulsar's planet count (the planet model gives a recycled pulsar as many as ten) belong to GEN.130.
- **Built (GEN.320): the stars a nebula needs.** `generation/nebula_contents.py`. A hosted nebula's roll (chance, class, age) runs on a stream seeded from the host star's own mass, age and temperature (`hosted_seed`), so the same star gets the same nebula whether the sector pre-places it or a fill makes it. A field cloud draws its members from a seed of its own (`cloud_seed`: the galaxy seed and the cloud's centre in millipc): a Bok globule (P) holds one protostar in 30% of clouds, a star-forming core (Q) a zero-truncated Poisson count of mean 2 capped at 4 (`tuning.NEBULA_CLOUD_PROTOSTARS`), each of IMF mass under 8 solar masses and under 0.5 Myr old, placed in a Gaussian clump a third of the cloud's radius wide. Every sector the cloud reaches draws the same members and builds only the ones inside its own cell (`add_cloud_protostars`, stage F), so a cloud across a border is built once. Each member is an ordinary young system (GEN.242: protostar, disk, Perseus companions) with no planets. A planetary nebula's central star takes a close companion in 20% (`tuning.PLANETARY_NEBULA_COMPANION_CHANCE`, [R: 15 to 20% of central stars are close binaries]). A massive member of an emission region (an O or early B younger than 0.3 Myr) is in a young group at the catalogue's 5%, its mates built as systems of their own (GEN.328, `multi-star-systems.md` 13.1). The members of the C and G nebulae's companies belong to GEN.241, and a pre-placed star's nebula shown from its scatter row is a remainder (TODO).
- **Built (GEN.530): a waiting star's hosted nebula is read from its row.** `nebula_contents.row_nebula` rolls `hosted_nebula_of` on the four numbers a row holds (`row_star`: type, age, mass, temperature), the same a built star carries, with the `"nebula"` factor of the row's cell (`db/hostednebula.py`), so the nebula is the one the sector builds. The uncharted star's page names it (`hosted_nebula`) and the uncharted sector's scene draws it at the star (`hosted_nebulae`, linking to the star's page); a star built into a system lists none, the system's own phenomenon holds it. The Galaxy Map's tiles are not changed: they list scattered phenomena by their own rows, and a hosted nebula has none until it is built.
- **Built (GEN.321): the ionizing members of emission regions D and E.** `nebula_contents.emission_members`. The region's `hosts` text counts its stars of 5.5 solar masses and up by spectral type (`physics/hii.cluster_hosts`); the seed star is one of them, so the members are the text less one star of the seed's type. Each member is a mass from the IMF within its type's range (`hii.mass_range_of_type`, the log-midpoints to the neighbouring anchor types; a draw whose star would be dead at the nebula's age is made again, the text having counted survivors), the nebula's age, and a point in a Plummer sphere about the seed star of half-mass radius 1 pc (D) or 2 pc (E) (`tuning.NEBULA_CLUSTER_HALF_MASS_PC`) from the table and scale of the satellites' profile (`physics/component_profiles`, GEN.237), redrawn while it lies outside the nebula's radius. The seed is the galaxy's and the seed star's point in millipc (`emission_seed`), so the members are a function of the galaxy, the nebula and that point alone. A sector builds those whose points fall in its cube (`run_sector.add_emission_members`, stage F, each an ordinary young system as in GEN.320), from two lists: the regions its own stars light (`add_star_hosted_nebulae`'s `hosts_out`) and the regions the pre-placed hot stars (O and B0 to B2) of the sectors around it light, read from their `bright_stars` rows within the Plummer cut radius plus the sector's reach (`run_galaxy._emission_hosts`, `FillContext.emission_hosts`): each row is rolled as its own sector rolls it (`nebula_contents.row_host`, with that sector's `"nebula"` factor), whether or not that sector was filled. The hosted roll's seed now takes mass, age and temperature as floats, so a star made in a fill (whole-kelvin integer temperature) and its stored row (a float) roll the same. Not covered: a hot star the scatter did not place (a B1 or B2 star under the luminosity floor of its tier) is no row, so its region's members fall in its own sector alone. A region's seed star must still be the host the sector's fill made; its stars below 5.5 solar masses stay the sector's background density (the default of the design note).
- **Not built here.** The members of the C and G nebulae's companies (GEN.241) and the companion draw for objects made at fill (stages D and G, GEN.264).

## 4. The catalogue

Columns. **Status**: built, filed (TODO id), new (proposed in section 8), skip. **Observed**: the real count or rate with its tag. **Now**: what the generator does. **System?**: probability that the object has a bound companion or planets, with the rule. **Carries**: objects the kind must hold. **Profile**: what drives its rate in the sphere profile (`mix` populations, `gas`, `kind_factors`; section 5.3 of the whole-sphere note). **Where**: scatter (S), fill (F), field (D, derived from a seed), flag, or component.

Counts per local sector use a 64 pc^3 cell with 8 stars [S: `interstellar-object-rates.md`].

### 4.1 Stars and stellar states (what a system's stars can be)

| Kind | Status | Observed | Now | System? | Carries | Profile | Where |
|---|---|---|---|---|---|---|---|
| Main-sequence star O to M | built | the bulk | 94.2% of living stars [C] | companion by class: O 0.90, early B 0.84, mid B 0.76, B 0.65, A 0.55, F 0.47, G 0.44, K 0.40, M 0.26 (`BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS`) then shared between close, inner and outer separations by Moe and Di Stefano's lines (GEN.448); planets by the planet model | none | `mix` | S for bright and mass, F otherwise |
| Subgiant, giant, supergiant | built | giants about 1% of field stars [R] | 0.58% giant or subgiant, about 0 supergiant [C] | as above; planets inside twice the giant's radius engulfed (`GIANT_ENGULFMENT_RADIUS_FACTOR`) | none | old and intermediate | F; bright ones S |
| White dwarf | built | about 6 to 10% of stars [R] | 5.2% [C] | companion about 30% (22% wide) [S: multi-star 7.3]; surviving planets outside 1.5 AU (`WD_PROGENITOR_ENGULFMENT_AU`) with orbits widened by the mass ratio; metal pollution of 25 to 50% of cool ones [S] unbuilt | none | old share, so halo and clusters have them | F |
| Hot central white dwarf | built | 1 per planetary nebula | O3VII to B0VII, mass 1.1 to 1.44 (GEN.149 fixes it) | close companion about 20% [R: 15 to 20% of central stars are close binaries]; circumbinary disk only if close binary, never planets (`exotic-environments...` 1.3) | itself | old share | made with its nebula |
| Brown dwarf as a star's companion or primary | **new** | free-floating 1 per 4.7 stars [S]; companions to Sun-like stars under 1% inside 5 AU, a few % wide (the brown dwarf desert) [R] | only as a free-floating row | its own companion probability 0.15 for the pair [R] | none | `mix` | F |
| Pre-main-sequence star, T Tauri, Herbig Ae/Be | filed GEN.242 | every star still contracting: 1.1e-4 of stars are under 10 Myr (2.8 stars/yr x 10 Myr = 2.8e7) [C], and since GEN.529 removed the 10 Myr cap about 3.6% of the disk's draw are red dwarfs younger than their contraction time (0.1 to 0.7 Gyr) | none; the model has no such state | companion 0.33 (Class II, six systems) rising to 0.57 (Class 0) [S: Perseus survey]; disk fraction falls from about 0.8 at 1 Myr to 0.1 at 5 to 6 Myr [R]; no planets | a disk descriptor | young share | F, inside nebulae |
| Protostar, Class 0/I | filed GEN.242 | about 1.4e6 in the Galaxy (2.8 stars/yr x 0.5 Myr) [C, R inputs] against the 2.6e5 the P and Q cores of GEN.241 give: 5x short | none | multiplicity 0.57 Class 0, 0.23 Class I, companion star fraction 1.2 and 0.23 [S] | a disk | `gas` and young share | F, inside P and Q cores and embedded in M clouds |
| Wolf-Rayet star | filed GEN.243 | 1,200 +- 100 expected, 642 known [S]; the model reproduces 1,180 [C] | none | WN 0.52, WC 0.9 (capped from 0.96) [S: Dsilva et al.]; planets none (under 10 Myr) | its ring nebula 40% [R, I] | young share only | S via the mass pass; F |
| Luminous blue variable, hypergiant | filed GEN.246 | tens in the Galaxy [R]; Eta Carinae-class events 0.09 of the core-collapse rate [S] | none | as massive stars (binary 0.9) | an ejecta nebula | young share | S, rare |
| Be star, carbon star, Mira, S-type | filed GEN.249 | Be about 20% of B stars [R] | none | ordinary | a decretion disk (Be) | intermediate and old | F, attribute |
| Blue straggler | filed GEN.162 | 10 to 30% of cluster cores [R] | none | binary origin | none | `cluster` | F |
| Runaway star | built | 1 to 2% of stars above 30 km/s [S] | flag on an ordinary system | ordinary | none | `mix` | flag |
| Hypervelocity star | built | 1e3 to 1e4 [S] | scatter row, an ordinary young system at fill | ordinary | none | r^-2 law | S |
| Thorne-Zytkow object | skip | 20 to 200 predicted, none confirmed [S] | | | | | |

### 4.2 Compact remnants

| Kind | Status | Observed | Now | System? | Carries | Profile | Where |
|---|---|---|---|---|---|---|---|
| Neutron star, isolated, young pulsar | built | about 1e9 neutron stars [S]; 1e5 beamed pulsars [R] | 5.6e-4 pc^-3; pulsar chance 0.02 (`NEUTRON_STAR_PULSAR_CHANCE`) | companion 1 to 3%; planets none (needs a disrupted white dwarf) [S: multi-star 7.3, exotic 2.1] | none | `mix`; the thicker layer of `compact-remnant-regions.md` is a proposal, not built | S |
| Neutron star, millisecond pulsar | built | about 600 known; field MSPs mostly with a white dwarf companion [R] | 10% of pulsars (`PULSAR_MILLISECOND_CHANCE`) | companion 75 to 80%; planets 0.7% (upper bound 1%): disk planets 25%, captured giant 25% in clusters only, ablated companion 50% [S: exotic 2.1] | the companion | old and bulge, clusters | S |
| Magnetar | filed GEN.153 | 17 expected, about 30 known [S, R] | not built | none | a young remnant under 3e4 yr | young share, thin layer | S |
| Stellar black hole | built | about 1e8 [S] | 7e-5 pc^-3 (superseded by the metallicity-dependent density of GEN.563, black-hole-density.md); the stellar class is 5e-4 per star, 40% of it above 14 solar masses | companion 1.25%, not 10 to 30% (Boss, 2026-10-15, 6:08 PM EDT, multi-star-systems.md; Gaia BH1, BH2, BH3 [S, low confidence]); planets of the companion only inside the Holman-Wiegert limit; none in the hole's own disk; 0.1% have an accretion disk (`BLACK_HOLE_ACCRETION_DISK_CHANCE`) | none | thick layer, rising toward the core | S |
| Intermediate-mass black hole | built | a handful of candidates, none confirmed in the field [R] | 3e-6 of black holes (`BLACK_HOLE_INTERMEDIATE_MASS_CHANCE`, GEN.259: Boss, decision card, 2026-10-16, "hundreds, as the observations suggest"), about 300 in the Galaxy of 1e8 black holes; it was 0.1%, about 1e5, four orders above the candidates | none | none | `cluster`, halo | S |
| Supermassive black hole, quasar | built | one per galaxy | one in the core cell, quasar with chance 0.1 | none: no stellar system, no blanets [S: exotic 2.4] | the S-star cluster as text (GEN.155) | the nucleus | F, one sector |
| Active X-ray binary, microquasar | filed with GEN.130 and GEN.247 | several dozen to 1,000 microquasars [R] | none | by definition a binary: a donor at hours to a day | an accretion disk, jets | young for high-mass, old for low-mass | F |
| Nuclear star cluster | filed GEN.155 | one, 2.5e7 solar masses [R] | none | its stars are ordinary | the supermassive hole | the nucleus | F |

### 4.3 Small bodies and debris

| Kind | Status | Observed | Now | System? | Carries | Profile | Where |
|---|---|---|---|---|---|---|---|
| Rogue planet, six mass bins | built | 5.6 terrestrial per star, giants 0.028 [S] | about 58 per local sector | a moon with probability 0.2 | none | per star | F |
| Rogue brown dwarf | built | 0.03 pc^-3 [S] | built | none | none | per star | F |
| Interstellar comet | built | 1e15 pc^-3 real, 0.007 pc^-3 kept as rows [S] | rows plus a computed count | none | none | per star | F |
| Asteroid field | built, rate 0 | disperses in 1e6 to 1e7 yr [S] | hand-made only | none | none | none | hand-made |
| Debris belt in a system | built | all stars with planets | in the planet draw | n/a | n/a | the star | F |
| Debris disk around A and F stars | filed GEN.245 | about 20% of A and Sun-like stars [R] | none | n/a (a descriptor) | dust radius and mass | intermediate | F |
| Protoplanetary disk descriptor | **new**, inside GEN.242 | most stars under 3 Myr [R] | none | n/a | disk mass and radius | young | F |

### 4.4 Nebulae and gas

| Kind | Status | Observed | Now | Needs inside | System? | Profile | Where |
|---|---|---|---|---|---|---|---|
| Warm neutral and ionized medium, A and B | built as background | | not generated | none | none | `gas` | D |
| Emission C, D, E, G | built; stars real in GEN.241 | about 8,000 catalogued [S]; 67,500 hosted in the model (270,000 O stars x 0.25) [I] | host-first; members text only | O or early-B star; D and E their cluster members (42 and 549 on average) | members are stars with belts at most | young share | F, host first |
| Reflection F | built | about 150 within a few kpc [S] | host-first | a B or A star | belts at most | young and intermediate | F |
| Planetary H to L | built; GEN.149 fixes the star | about 20,000 [S] | own white dwarf, always | one hot white dwarf | see 4.1 | old share | S, made with its star |
| Giant molecular cloud M | built | about 8,100 to 9,700 [S]; the generator keeps 102,000 | field | none; embedded protostars not built, shortfall accepted (GEN.258) | none | `gas` | D |
| Dark cloud N | built | | field | none | none | `gas` | D |
| Bok globule P | built | 3.2e5 small globules [R] | field | protostar in 30% [I] | protostar systems | `gas` | D |
| Star-forming core Q | built | | field | 1 to 4 protostars | protostar systems | `gas` | D |
| Supernova remnant R to V, W | built | 400 to 3,000 visible [S] | 1.5e-9 pc^-3 | core-collapse: a neutron star or black hole in 60%; Type Ia none | none | young share (core collapse), old share (Type Ia) | S |
| Wolf-Rayet ring X | filed GEN.243 | tens of rings; 40% of WR stars [R, I] | none | its W star | the W star's companion | young | F, host first |
| Pulsar wind nebula | built as morphology | Crab and others [S] | morphology plerion | its pulsar | none | young | S |
| Superbubble, Local Bubble | **new** | Local Bubble 100 to 200 pc; superbubbles fill about 12% of the disk [S] | none | none | none | `gas`, young | D |
| High-velocity clouds, Magellanic Stream | filed GEN.239 | [S] | none | none | none | halo | component |
| Herbig-Haro and protostellar jets | **new** | up to 150,000 [S] | none | its protostar | none | `gas` | D |
| Microquasar and pulsar jets | **new**, needs the spin axis (GEN.104) | [S] | none | its source | none | per source | F |
| Fermi bubbles | **new**, map overlay | two lobes, about 8 kpc [R] | none | the nucleus | none | nucleus | map only |

### 4.5 Groupings and components

| Kind | Status | Observed | System? | Notes |
|---|---|---|---|---|
| Globular cluster | filed GEN.159 to GEN.164, GEN.209 | 157 in the Milky Way [R] | members are ordinary systems with the hard-orbit planet cut and type-B pulsar planets | `cluster` population |
| Young open cluster, OB association | filed GEN.244 | about 3,000 catalogued, total in the tens of thousands [R] | members ordinary, age gates | young; feeds D and E members |
| Satellite and dwarf galaxy | filed GEN.237 | about 60 within 300 kpc [R] | members ordinary | a component |
| Tidal stream | filed GEN.248 | Sagittarius and others [S] | | an elongated density term |
| Stellar halo, thick disk | filed GEN.234 | 1.5e9 solar masses halo [S] | members ordinary, old | halo population |
| Neighbour galaxies | filed GEN.157 | about 25 real [R] | n/a | catalogue file |

### 4.6 Flags, derived values and effects (not objects)

| Item | Status | Rule |
|---|---|---|
| Runaway and hypervelocity flags | built | flag on an ordinary system |
| Einstein radius on compact-object pages | filed GEN.154 | derived from mass and distance |
| Radiation hazard score | **new**, part of GEN.114 | derived from distance to O/B, Wolf-Rayet, remnants and magnetars |
| Fast radio burst flag | skip as an object | a flag on a magnetar |
| Interstellar debris count | built | computed on the sector page |
| Dark matter and gas medium | filed in the profile (whole-sphere note) | profile fields |

### 4.7 Skipped, with the reason

| Kind | Reason |
|---|---|
| Subspace anomalies (rifts, tears, interphase, sandbars, Omega, anti-time, implosion rings) | no physical basis; a lore decision for Boss, default off (`anomalies.md`, Tier 3) |
| Wormholes, cosmic strings, quantum filaments | not observed; strings excluded by pulsar timing [S] |
| Dark matter nebulae and a Dark Matter Anomaly object | dark matter does not clump on 10 ly scales [R]; the profile carries the background |
| Quark stars | no candidate [R] |
| Thorne-Zytkow objects | none confirmed [S]; Boss's 6:08 PM EDT widening puts them back as a flagged theoretical type (100 in the Galaxy), see 4.8 |
| Stellar mergers and red novae | transients; 5e-12 per sector per year [C] |
| Gamma-ray and fast radio bursts | transients with no persistent object |
| Great Attractor, Bootes void | beyond the galaxy; at most flavour text |
| Blanets (planets around a supermassive hole) | no confirmed one; no stellar system has a supermassive primary (`exotic-environments...` 2.4) |
| Kilonovae, compact-object mergers | out of scope (needs a binary-evolution code) |
| Free-floating asteroid fields | disperse in 1e6 to 1e7 years; rate 0 stands |

### 4.8 Member kinds and system types added by the catalogue of multi-star types

Boss (2026-10-15, 6:08 PM EDT) widened the multi-star plan to every observed and theoretical type of multi-star or star-like-object system. [multi-star-systems.md](multi-star-systems.md) sections 13 to 16 hold that catalogue (about sixty types, each with its observed or theoretical rate, where it occurs, and whether it is drawn at scatter time or at fill). It adds these to the kind list of the registry (GEN.250), with their rates in `SYSTEM_TYPE_RATES` (GEN.268):

| Kind or node | Status | Where it comes from |
|---|---|---|
| `stripped-helium` (hot subdwarf, helium star, helium donor) | filed GEN.269 | post-common-envelope pairs, AM CVn donors, Be companions |
| `thorne-zytkow` | filed GEN.269, flagged `theoretical`, 100 in the Galaxy | a giant with a neutron-star core, single |
| `wolf-rayet` as a star kind | filed GEN.269 on GEN.243's W state | young massive stars and their companions |
| `group` node (young non-hierarchical group) | filed GEN.273 | trapezia under their decay age |
| `compact` node (a tight compact pair promoted to the scatter) | filed GEN.260 | GEN.263 |
| interaction states on a pair (detached, semi-detached, contact, accreting, post-common-envelope, merger-product) | filed GEN.269 | the Roche classifier |

## 5. The system rule for every kind

A kind may have a system if its row in sections 4.1 to 4.4 gives a probability. The roll and what it builds:

| Object | P(system) | What the system holds | Source |
|---|---|---|---|
| Ordinary star | by spectral class (0.26 to 0.90), then wide 0.5 | companion, planets, belts, moons | existing model |
| Star under 10 Myr (nebula member, protostar, Wolf-Rayet, T Tauri) | as its class, higher when younger | companion, belts, a disk descriptor; **no planets** | `system.py:1024`; Perseus survey |
| White dwarf | 0.30 | main-sequence or white dwarf companion; surviving planets | multi-star 7.3 |
| Planetary-nebula central star | 0.20 | close companion; circumbinary disk, no planets | exotic 1.3; [R] |
| Young neutron star | 0.02 | Be or main-sequence companion | multi-star 7.3 |
| Millisecond pulsar | 0.78 | white dwarf or ablated companion; planets 0.007 | multi-star 7.3; exotic 2.1 |
| Stellar black hole | 0.0125 | Sun-like or giant companion at 1 to 20 AU; wide tertiary only for direct-collapse holes | multi-star 7.3; exotic 2.3 |
| Magnetar, intermediate-mass hole, supermassive hole, quasar | 0 | none | |
| Rogue planet | moon 0.2 | one captured moon | `hosts.PLANET_SATELLITES` |
| Rogue brown dwarf | 0 | none | |
| Comet, asteroid field | 0 | none | |
| Nebula | 0 for planets | its required stars (section 4.4) | rule 2 |

The expected cost of the new system pass for a local sector (64 pc^3) is: neutron stars 0.036 x 0.0215 plus black holes 0.0045 x 0.0125 plus a planetary-nebula term below 1e-6: about 0.0008 extra systems against 6.3 [C]. White dwarfs are ordinary stars of stage D and add nothing. In the halo and clusters the same fractions apply to a far smaller count.

Where the observed probability is poorly known (the black-hole companion fraction rests on three Gaia systems), the table gives 1.25% for the black-hole companion fraction (Boss, 2026-10-15, 6:08 PM EDT, multi-star-systems.md, in place of the earlier note's 10 to 30%) and the build keeps it in one tuning dict so it can be changed without code.

## 6. One registry

Today a kind is defined in `SECTOR_PHENOMENON_KINDS` (`run_sector.py:50`), `SCATTERED_KINDS` (`phenomenon_scatter.py:55`), `PHENOMENON_DENSITY_PC3` (`tuning.py:2726`), `KIND_CODES` (`names/object_id.py:36`, six bits, codes 12 to 63 free), the class pages (`classref`), the marker list of MAP.132, and the pre-placed passes. Proposed: one table, `phenomena/registry.py`, one row per kind, with these columns (the same columns as section 4):

`kind`, `family`, `status`, `rate` (density, scaling law, and the profile fields that drive it), `gate` (which populations and `gas` allow it), `mode` (scatter, fill, field, flag, component), `builder`, `host` (the host class, GEN.552: what can accompany and orbit it; `generation/phenomena/hosts.py`), `carries` (required objects and their builder), `id_code`, `page_class`, `marker`, `evidence_note`.

The scatter lists, the fill stages, the rate tables, the ID codes and the class pages are generated from it. Completeness test: every row has every column; every `builder`, `system_model` and `carries` resolves; no ID code repeats; no code path creates an object of a kind that is not in the table. Adding a kind is then one row, a builder and a test.

## 7. Representation: how "all types" is kept true

Interpretation (Default): "all types represented" means every kind in the registry has a working generation path at its observed frequency in every region where it occurs, and the generator can prove it. It does not mean every sector holds every kind, which would put Wolf-Rayet stars in the halo and a supernova remnant in every cell. If Boss means the stronger reading, the change is a per-kind minimum in the registry and I will cost it.

The audit, run in the test suite and as an admin report:

1. **Completeness.** The registry test of section 6.
2. **Expected counts.** For each region class of the profile (core, bar, thin disk arm, thin disk interarm, thick disk, halo, Magellanic-type satellite, globular cluster core, intergalactic floor) the registry gives the expected count of each kind in a sector from its rate and gate.
3. **Observed counts.** Fill a fixed set of seeded sectors in each region class (the same seeds every run) and count each kind by stage.
4. **Compare.** Each kind's total must lie within Poisson bounds of its expectation (four sigma), and any kind with an expected total of 20 or more over the set must appear at least once. The report lists kinds that never appear, kinds that appear where their gate forbids (a Wolf-Rayet nebula in the halo), and kinds whose count is out of bounds.
5. **Systems.** For every kind with a system model, the fraction of its objects that got a system must match the table of section 5.
6. **Required contents.** Every hosted nebula holds its host; every planetary nebula its white dwarf; every core-collapse remnant a compact core at the stated share; every D or E region its members inside its radius; every Wolf-Rayet nebula a W star.

The set is sized to run in a few minutes; the full-galaxy rates are checked once per release by the full run that already exists.

## 8. Proposed TODO items

Existing items stay as filed (GEN.153 to GEN.164, GEN.209, GEN.233 to GEN.249, GEN.130, GEN.114). New or re-scoped items; the TODO thread numbers them and assigns lanes. Order is build order.

1. **Phenomena registry** (Core; first): section 6, with the completeness test; the existing lists read from it; no behaviour change.
2. **Fill stages with their own streams and the system pass** (Core; after 1): section 3 and 5. Rescopes GEN.130: GEN.130's slices (a) to (c) become the contents of the system pass for compact remnants, while this item builds the machinery and the white-dwarf, planetary-nebula-star and rogue-planet rows.
3. **Compact-remnant systems** (Core; after 2; this is GEN.130 reworded): companions and pulsar planets by the probabilities of section 5, the kick-and-survival check, Holman-Wiegert zones.
4. **Representation audit** (Quality; after 2): section 7, with the admin report.
5. **Brown dwarfs as stars** (Core; later): an L/T/Y state for companions and primaries, from `ROGUE_BROWN_DWARF_*` cooling, with the brown-dwarf-desert companion rule. Built as GEN.252 (`brown-dwarf-state.md`): the state and the free-floating row's type; the companion rows are GEN.275's.
6. **White-dwarf pollution flag and cool-white-dwarf presets** (Core; later): 25 to 50% of cool white dwarfs polluted [S]; part of the GEN.130 white-dwarf row.
7. **The anomalies of GEN.114, one item per kind** (Core and Interface; later): superbubble field, Herbig-Haro and microquasar and pulsar jets (spin axes exist), radiation hazard score, Fermi-bubble overlay. GEN.114 is the umbrella that never got split.
8. **Protostars in giant molecular clouds** (Core; with GEN.241): the P and Q cores give 2.6e5 protostars against 1.4e6 expected [C]; embedded Q-type groups in M clouds make up the rest. Default: accept the 5x shortfall for the first release. **Decided (GEN.258, closed by the TODO thread's answer of 12:17 PM EDT 2026-10-23, the open-question defaults being Decided since 2026-10-17):** the shortfall of about 1.14e6 protostars (1.4e6 expected against 2.6e5, 5.4 times short) is accepted for the first release and the embedded groups are not built; building them would add about 11 protostars to each of the generator's 102,000 M clouds and a rescatter. The yield is measured by the representation audit of a whole fill (TEST.160 does not count the cloud's protostars against 1.4e6); a later item reopens this by changing `tuning.NEBULA_CLOUD_PROTOSTARS` or adding the groups in `nebula_contents.py`.
9. **Docs** (Quality): the catalogue pointers in `interstellar-object-rates.md`, `anomalies.md`, `multistar-and-compact-systems.md` and `nebula-and-asteroid-field-classes.md` so they say where each kind now lives.

## 9. Decisions and defaults

Defaults the build proceeds on, each reversible:

1. **Representation means audited frequency**, not a guaranteed instance in every sector (section 7).
2. **System probabilities** are the earlier notes' numbers (section 5); the black-hole companion fraction, 1.25% (Boss, 2026-10-15, 6:08 PM EDT; the earlier 0.20 is replaced), is the weakest.
3. **Planetary-nebula central star companion 0.20** [R], unsourced beyond recall.
4. **Protostar shortfall accepted** for the first release (item 8).
5. **Skips** of section 4.7 stand; the subspace anomalies stay a lore decision, default off.
6. **Registry first** in the build order, because every later item adds rows to it.

## Notes folded in

`interstellar-object-rates.md`, `compact-remnant-regions.md`, `multistar-and-compact-systems.md`, `exotic-environments-planets-and-compact-binaries.md`, `anomalies.md` (GEN.113), `nebula-and-asteroid-field-classes.md`, `nebula-density-vs-reality.md`, `globular-clusters.md`, `phenomenon-scatter-mass-cut.md`, `intergalactic-and-halo-contents.md`, the research files `11-nebula-planets-and-required-stars.md` and `20-anomalies-and-asteroid-fields.md` (in /mnt/project-files/research), and the TODO items named in the tables.

## Sources

See the sources of each folded note. Added for this note (search-result summaries only): [R4379], [Galactic Wolf-Rayet multiplicity](https://cris.iucc.ac.il/en/publications/the-multiplicity-of-galactic-wolf-rayet-stars/).
