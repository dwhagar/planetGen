# Interstellar asteroid fields: rogue-planet collisions and the other sources

Status: plan (items filed 4:51 PM EDT 2026-10-30; Boss decided 4:54 PM EDT: physics only with no game rate, retire a field past 1 ly, the leftover planet of a collision is its own rogue planet; nothing built). Written 4:40 PM EDT 2026-10-30 as the lane reported. No generator code changed. Asked for by Boss (4:40 PM EDT 2026-10-30): make asteroid fields the result of two rogue planets colliding, evolve them when a sector is generated, give large bodies ID names, and treat the field like an in-system belt. Builds on `belt-evolution-and-body-ledger-2026-10-30.md` (the "ledger draft": mass cut, L1/L2/L3 tiers, status words, `<belt ID>.<rank hex>` designations) and `belt-origins-scatter-and-class-d-2026-10-30.md`. Repo today: `generation/phenomena/asteroid_field.py`, `tuning.ASTEROID_FIELD_*`, `docs/design/systems/interstellar-object-rates.md`, `nebula-and-asteroid-field-classes.md`.

Evidence tags: [S] a source exists and was found by search this session (title and arXiv number seen, paper body not opened, so the building lane reads it before using a number); [C] computed this session from the stated inputs (script in the scratchpad, formulas below, reproducible); [R] recalled, not opened; [REPO] read in the repo. Every [R] number the build uses goes into GEN.725 with tier and tweak level.

## 1. Short answers

1. **Rogue-rogue collisions are real physics but almost never happen.** About **one collision in the whole life of the Galaxy** (range 0.1 to 10) [C]. A field made that way would be unrecognisable within about 60 thousand years, so the standing number in the Galaxy is about 1e-5. A generator that follows the physics makes none. That is what the repo's rate 0 already says, and the reason it gives (fields disperse) is right.
2. **Even when it happens, it rarely leaves asteroids.** At normal rogue speeds (30 to 60 km/s) two Earth-size rocks hit with peak shock pressure of 1 to 4 TPa, which is mostly melt and vapour that recondenses as droplets and dust. Boulders and 100 km fragments need a slower hit, roughly 15 to 30 km/s, a few percent of collisions [C]. Slower than about 15 km/s the pieces mostly fall back together into one planet.
3. **A field cannot evolve by grinding.** With no star holding it together the cloud flies apart in a straight line. Collisions stop within hours for big pieces and within years for grit [C]. After that "evolution" is only expansion (a pure function of age), then galactic shear. So bodies are **derived from a seed and an age**, not simulated.
4. **The best physical sources are not collisions.** The leftover belts of stars that went supernova are unbound by the mass loss: about 100 to 200 fields standing in the Galaxy [C with [R] inputs], roughly 4 orders of magnitude more than collisions, and still about 1 sector in 2e8. They are the right ones to host, because a young supernova remnant already exists in the generator.
5. **To have any fields at all we need a labelled game rate.** Recommendation in section 8: physical sources hosted where they belong, plus one stand-in "wanderer" rate at tier "from needed output", tweak L4, default 1 field per 10,000 sectors, registered next to the physical figure so research can replace it.
6. **Fit with the belt design is good.** A field is a belt with no star: same ledger, same mass cut, same designation rule, same planet-status test. Differences are: it expands, it has a hard lifetime, and its members are drawn once and moved by `position = v x age`.

## 2. How often do two rogue planets collide?

### 2.1 Inputs [REPO unless marked]
* Rogue density: 6.5 rogues per star, so 0.91 per pc3 at the solar neighbourhood, of which terrestrial 0.78 per pc3 (`tuning.ROGUE_PLANET_RATE_PER_STAR`, `PHENOMENON_DENSITY_PC3`). The observational root is the microlensing work (Mroz et al. 2017, 2020; Johnson et al. 2020, cited in tuning) [REPO]; the terrestrial number is already an extrapolation.
* Relative speed: rogues share the stellar velocity dispersion, plus the ejection kick. Relative 3-D speeds of 30 to 60 km/s in the old thin disc, 100 km/s and up in the bulge, 2 to 4 km/s inside a young cluster [R]. I carry 20, 40 and 80 km/s.
* Gravitational focusing: cross-section `sigma = pi (R1+R2)^2 (1 + v_esc^2 / v^2)` with `v_esc` at touching of the pair (11.2 km/s for two Earths). At 40 km/s focusing is 1.08, so it hardly matters.

### 2.2 Rate per rogue [C]
`Gamma = n sigma v`. For two Earth-size rogues at 0.78 per pc3:

| Relative speed | Gamma per rogue per year | Per 10 Gyr |
|---|---|---|
| 10 km/s | 9.6e-24 | 1e-13 |
| 20 km/s | 1.1e-23 | 1.1e-13 |
| 40 km/s | 1.8e-23 | 1.8e-13 |
| 80 km/s | 3.5e-23 | 3.5e-13 |

With the real mass function (most rogues are 0.1 to 0.5 Earth masses, radius 0.5 to 0.7 Earth) the number is 2 to 4 times lower.

### 2.3 Galaxy total [C]
`events/yr = 0.5 x sigma v x (6.5)^2 x integral(n_star^2 dV)`. Exponential disc (scale length 2.6 kpc, height 0.3 kpc, 0.14 stars per pc3 at the Sun, 8.4e10 stars) gives 3.4e10 star^2 per pc3 pc3; a 2e10 star Plummer bulge of 700 pc gives 5.1e10.

| Case | Collisions per year, Galaxy | Over 10 Gyr |
|---|---|---|
| Earth-size pairs, disc 40 km/s, bulge 120 km/s | 8.9e-11 | 0.9 |
| Mass-function pairs (0.7 Earth radius) | 4.3e-11 | 0.4 |
| Young clusters (n* 100 to 2e4 per pc3, 1.5 to 4 km/s, 10 Myr, 1 rogue per star): per rogue 5e-14 to 4e-12, so for 1e11 clustered stars | | 0.01 to 1 |

**Answer: order one per Galaxy lifetime, about one per 10 Gyr to one per 1 Gyr at most.** Standing fields at any moment: rate x field lifetime (section 5) = 1e-10 x 6e4 yr = **6e-6**. In a 3e10-sector Galaxy, no sector holds one, ever, unless we put it there.

Cluster collisions are the only ones slow enough to leave real fragments, and they sit inside a cluster where the fragments are dispersed in a few Myr.

### 2.4 Related: tidal disruption of a rogue by an isolated neutron star or black hole [C]
Tidal radius for an Earth: 4.9e5 km (1.4 Msun neutron star), 9.5e5 km (10 Msun black hole). Focusing is large (escape speed 870 and 1670 km/s there). Per object 1.7e-18 and 7.8e-17 per year; with the repo's 1e9 neutron stars and 1e8 black holes, 30 percent of them in rogue-rich space, about **5e-10 and 2e-9 events per year**, 1 per 0.5 to 2 Gyr. Debris is thrown at about 100 km/s, so it leaves a 1 ly cloud within about 3 thousand years: standing number 1e-5. The bound half is a disc around the remnant (this is the pulsar-planet origin that already exists). Not a field source. Listed for completeness.

## 3. What a collision leaves behind

### 3.1 Energy and shock [C]
Specific impact energy of an equal-mass pair `Q_R = v^2 / 8`; gravitational binding energy of an Earth `3GM/5R = 3.7e7 J/kg`; complete vaporisation of silicate is of the order 1e7 J/kg for the heating share [R].

| v km/s | Q_R (J/kg) | Peak shock pressure, rock on rock [C, forsterite C 5.5 km/s, s 1.3] | Outcome |
|---|---|---|---|
| 5 | 3e6 | below 100 GPa | merger, small ejecta, nothing escapes |
| 10 | 1.3e7 | 200 GPa | partial disruption; mostly re-accretes (hit and run if oblique) |
| 15 | 2.8e7 | 400 GPa | near catastrophic threshold, large remnant plus fragments |
| 20 | 5e7 | 600 GPa | catastrophic: fragments from 10s to 1000s of km, some melt |
| 30 | 1.1e8 | 1.2 TPa | super-catastrophic: melt and vapour dominate, small remnant |
| 40 | 2e8 | 2.1 TPa | mostly vapour and droplets |
| 60 | 4.5e8 | 4.4 TPa | essentially all vapour and melt |

Outcome regimes and the largest-remnant law are in Leinhardt and Stewart 2012, arXiv:1106.6084 and 1109.4588 [S, to read]. I used only the threshold idea (disruption when `Q_R` exceeds `Q*_RD`, which for planet size is of the order of the binding energy); the building lane replaces my table with their law, `M_lr / M_tot` as a function of `Q_R / Q*_RD`, and registers any piece it cannot read as a stand-in.

**The window that makes an asteroid field is 15 to 30 km/s.** Fraction of collisions in it, for a Maxwellian relative speed with 1-D dispersion s [C]:

| s (km/s) | below 15 | 15 to 30 | 30 to 60 | above 60 |
|---|---|---|---|---|
| 25 | 3 percent | 17 percent | 60 percent | 20 percent |
| 50 | 0.2 percent | 1.6 percent | 15 percent | 83 percent |

So only about 1 to 17 percent of the (already rare) collisions leave rubble. The rest give a vapour and droplet cloud (still a "debris field" but made of mm to m material) or a merged planet.

### 3.2 What is in it
* **Mass:** the unbound part, not the sum. Unbound mass fraction runs from near 0 (merger) to near 1 (super-catastrophic); a remnant planet is a **separate rogue planet** with its own ID, left near the field's centre.
* **Composition:** from the two parents. A differentiated rogue has an iron core (about 0.3 of its mass) and a silicate mantle; water-rich rogues (GEN.26 `ROGUE_WATER_RICH_CHANCE`) add ice shards. The repo's `small_body_composition.draw_field("collisional")` picks one taxon; replace with shards: iron/M-type from the cores, V/A or S-type from the mantles, D/P-type from ice. The ledger draft's rule (differentiated bodies only if the parent was) is automatically met because the parents are differentiated.
* **Speed:** ejecta dispersion about 0.3 to 1 times the impact speed above `v_esc`: stand-in 5 to 15 km/s (L1, nudged inside the physics bracket; to be derived in item F1).
* **Sizes:** a steep tail in the catastrophic regime. Largest fragment at 15 to 30 km/s: tens to hundreds of km (stand-in; from LS12 once read). Below a few km the vapour/droplet share dominates. The L1 ledger tier (mass cut rule E, floor 1e19 kg) holds the large pieces.
* **Shape:** a sphere is acceptable for a collision. Supernova-unbound belts (section 4) are a thin expanding **annulus/pancake**, not a ball; store the shape flag.

## 4. Other plausible sources

Rates are stand-in orders of magnitude; the same lifetime rule gives the standing number.

| Code | Source | Mechanism | Numbers | Standing Galaxy-wide | Verdict |
|---|---|---|---|---|---|
| S1 | Rogue-rogue collision | section 2 | 0.4 to 1 per Galaxy lifetime [C] | 1e-5 | physics-true, essentially none |
| S2 | **Belt unbound by a supernova** | An A to B star keeps a debris belt (A-star 24 micron excess is common at 10 to 100 Myr; Su et al. 2006, arXiv:astro-ph/0608563 [S, to read]; stand-in 30 percent keep one). The core collapse removes more than half the mass, the belt (50 AU, orbital speed 13 km/s around 10 Msun, escape speed 7 km/s around the 1.4 Msun remnant) is unbound and leaves in a flat expanding ring. | core-collapse rate 0.02 per yr [R], belt fraction 0.3 [R, stand-in], leaves 1 ly in 25 to 37 kyr at 8 to 12 km/s [C] | **90 to 220** [C], 6e-10 to 9e-10 per star, 3e-9 to 8e-9 per 4 pc sector | the best physical source. Hosted by a young SNR (age under about 40 kyr). Shock heating of the belt by the blast: ice lost, rock kept [R]. |
| S3 | Young-cluster stripped swarms | Flybys truncate discs at about 0.3 of periastron (Holman and Wiegert style cut already in repo; Breslau et al. 2014 [R]); the cut-off planetesimals fly free. Cluster dispersion 1 to 3 km/s, so the swarm leaves in a few Myr, then joins the interstellar background. | cluster-limited | tens of clusters x 1e5-ish swarm members, no coherent field after 10 to 20 Myr | hook only: a young cluster (age under 20 Myr) may list a swarm. Low priority. |
| S4 | Rogue disrupted by NS/BH | section 2.4 | 5e-10 and 2e-9 per yr | 1e-5 | not a field source. Feeds the pulsar-planet disc. |
| S5 | Ejected belt of a disrupted system | A stellar flyby or instability ejecting a whole belt. A star keeps its belt through a flyby unless the encounter is closer than about 3 times the belt radius; the ejected share is part of the isolated-body background rather than one cloud. | [R] | none | goes to the interstellar-comet background (`INTERSTELLAR_DEBRIS_DENSITY_PC3` 1e15 per pc3, consistent with Do, Tucker and Tonry 2018 arXiv:1801.02821, about 0.1 to 0.2 per AU3 [S]), which is not a field. |
| S6 | Debris of a destroyed in-system planet | A giant impact inside a system: debris is **bound** to the star. | belongs to the belt `collisional-debris` origin (origins doc, section 3) | not interstellar | already a belt, not a field. |
| S7 | Disrupted giant | A gas giant torn up by a star or hole forms a stream of gas, little rock. | rare | none | skip. |
| S8 | Hand-made or facility | admin / future 8.6 edits | | | keep the type for it. |

Source of rogues themselves: planet-planet scattering ejects them with a kick of several km/s [R]; all of them already come from `RoguePlanet`, which matters for section 7.

## 5. How a field disperses and ages

### 5.1 No central star, so no orbit
Each body moves in a straight line relative to the field centre at its ejection velocity. Radius `R(t) = R0 + sigma_v t`. Galactic tidal shear (Oort A 15 km/s per kpc) adds a drift of `2 A x t` for a body at offset x: for x = 1 pc and t = 10 Myr that is 0.3 pc, small beside `sigma_v t` until about 100 Myr and 10 pc [C]. So shear is a second-order step for fields that matter.

### 5.2 Age at a given radius [C]

| sigma_v (km/s) | to 63 AU (the smallest field) | to 1 ly (the largest field) |
|---|---|---|
| 0.1 (tidal disruption of one rubble body) | 3 thousand yr | 3 Myr |
| 1 | 300 yr | 300 thousand yr |
| 5 | 60 yr | 60 thousand yr |
| 10 | 30 yr | 30 thousand yr |
| 15 | 20 yr | 20 thousand yr |

The repo's field radius range (63 AU to 1 ly, uniform draw) is thus a range of ages up to 20 kyr to 3 Myr. **Radius becomes a result of age**: births are uniform in time, so the age is uniform in `[0, t_life]` and radius follows. A field larger than 1 ly is no longer a field (spacing between 100 km bodies of an Earth mass cloud at 1 ly is thousands of AU; its column density is far below the interstellar-debris background of 100 m bodies) and is **retired**. This is the first object type whose end is physics rather than a hand-set removal; bring-to-epoch (GEN.1006) must handle it (item L5).

### 5.3 Grinding freezes out within hours to years [C]
Collisions need optical depth above 1: `tau = total cross-section / (pi R^2)`. For one Earth mass at 5 km/s ejection:

| Piece size | Radius where tau = 1 | Time to get there |
|---|---|---|
| 100 km | 6.5e-4 AU | 5 hours |
| 1 km | 6.5e-3 AU | 2 days |
| 1 m | 0.2 AU | 2 months |
| 1 mm | 6.5 AU | 6 years |

So the size distribution of the pieces above a metre is whatever the collision (and the first hours) made; **there is no Dohnanyi grinding in the field, and no Wyatt `1/t` mass loss.** The ledger's grinding and Yarkovsky steps are switched off for fields (no star, no luminosity), which is a correct simplification, not a hole. Grit and dust persist but contribute nothing to the ledger (L3 tail). Without a star there is no radiation pressure or Poynting-Robertson loss; interstellar gas drag and sputtering are slow [R].

### 5.4 Density per area
Column density `Sigma = M_field / (pi R(t)^2)`, falling as `t^-2`. A "dense/typical/sparse" word (the repo's three levels) is derived from it with the same thresholds the belt rule uses. The word changes as the field ages, never stored.

## 6. Storage, ID and bodies

* **Mass rule (GEN.1161, Boss approved 3:14 PM EDT 2026-10-31).** A field holds what the disc stretch it came from has left, by the belts' own rule, and its bodies add up to it: a supernova-unbound belt is the progenitor's belt ground for the progenitor's life (`belt_ledger.solve_belt`) and then kept; a collision field is the unbound share of the pair (`rogue_collision.unbound_fraction`); nothing is drawn from a table of typical masses. Built in `physics/field_mass.py`; `mass_total_kg` below stores its result.
* **One row per field** (the existing `asteroid_fields` row with its 12-byte `uid`) gains: `mass_total_kg`, `seed`, `birth_age_yr` (age at the sector's epoch), `sigma_v_kms`, `shape` (ball, annulus), `origin` (S1, S2, ...), `history` text, `parent_a` and `parent_b` (mass, class, speed; text), `remnant_id` (the rogue left behind, if any), `size_law` (slope, largest-fragment mass). Units go in `db/column_units.py`.
* **Position and size at any epoch are formulas of the row**, so the bring-to-epoch step is a pure function, with one exception (retirement, 5.2).
* **Bodies are derived.** L1 members (mass cut rule E, floor 1e19 kg, cap 300) come from the seed: rank k has mass from the order-statistic draw, velocity from `Normal(0, sigma_v)` per component, position `x_k = v_k x age`. Stable rank by mass drawn once at birth, as the belt draw. L2 (1e15 to 1e19 kg) is drawn on demand; L3 is the tail as total mass.
* **Designation:** `<field ID>.<rank hex>`, the same function as belt members (ledger draft section 9, `designation(parent_id, rank)`). The field's own ID is its existing uid. The older display name `AF C3-FE81000A2B-01` stays as the field's label. A remnant planet takes an ordinary rogue-planet ID (kind code 1); a pre-existing rogue is **not** reused as a parent, since the parents are drawn virtually (below).
* **Status words** per L1 body: the same `body_status.py` flags. Field members are on a free orbit, so "cleared" uses Soter mu against the field's own mass and is nearly always "minor body". The only rounded bodies would be a few big fragments above the roundness mass, which are "dwarf planet" by the repo's word and have no star to orbit (the IAU word needs a star). Store the numbers; the displayed word for free bodies is `free dwarf body` or `minor body` (Boss to confirm word in the ledger decision, section 11 of the ledger draft already stores both).

## 7. How the existing rogue-planet generation feeds a field

A collision cannot use existing rogue rows as the parents. At 0.78 per pc3 the nearest two rogues are about 0.5 pc apart, and one of them is not in the same sector as a rule; there is no event, only the draw. So the field draws **two virtual parents** from the same machinery the rogues use:
* mass: `tuning.ROGUE_PLANET_MASS_BINS` with the mass function exponent 0.65, terrestrial and sub-Neptune bins only (gas giants leave gas, not rubble);
* class and composition: `roguePlanetData.choose_rogue_planet_class` (classes C and D for rock and ice; water-rich chance `ROGUE_WATER_RICH_CHANCE`); core fraction from the existing small-body family composition;
* relative speed: Maxwellian from the sector's local velocity dispersion (`galactic_orbit`), restricted to the 15 to 30 km/s window for rubble;
* the history sentence: "two rogue planets of 0.8 and 0.3 Earth masses met at 22 km/s about 31 thousand years ago; the 0.5 Earth-mass remnant is 0.1 ly away".

Per-sector draw, not a global event list, as everywhere else in the generator.

## 8. Rate calibration

The generator's rate scale is "per star, times the sector's star count" (`phenomenon_rate_per_star`). A 4 pc sector holds about 9 stars; the Galaxy has about 3e10 sector-equivalents. Physical sources, per sector [C]:

| Source | Standing per sector | Sectors per field |
|---|---|---|
| S1 rogue collision | 2e-16 | 5e15 |
| S2 supernova-unbound belt | 3e-9 to 8e-9 | 1.2e8 to 3e8 |
| S4 tidal disruption | 3e-16 | far fewer |

Physics-only means a player never meets one. Three options, with what they mean for a player:

| Option | Rule | Standing in the Galaxy | Distance to the nearest field |
|---|---|---|---|
| A. Physics only (S2 hosted by young SNRs, S1 by hand) | no free draw | about 100 to 200 | thousands of pc |
| B. A plus a labelled game rate of **1e-4 per sector** (1.6e-6 per pc3, 1.1e-5 per star) | stand-in, tier "from needed output", tweak L4 | about 3e6 | about 85 pc |
| C. A plus 1e-3 per sector | same, 10 times higher | about 3e7 | about 40 pc |

Option B means births of about 50 fields per year Galaxy-wide, ten orders of magnitude above the collision rate. That mismatch is the honest price of a visible field, and it is why the stand-in is registered with the physical figures beside it. **Recommendation: B**, because "plausible rather than common" is a field every 10,000 sectors (a bubble about 85 pc across holds one) and because S2 stays the only physically rooted channel, which shows next to every young supernova remnant. A (nothing generated) is the strictly honest answer and stays the default if Boss wants no stand-in.

Calibration rule for the dial: `lambda_sector = lambda_hosted(S2, S3) + lambda_floor`, where `lambda_floor = FIELD_STORY_RATE_PER_SECTOR` (registered GEN.725 row, L4), scaled by sector star density like the other phenomena, and every drawn floor field carries origin `collision (stand-in rate)` so it is searchable and removable. A test checks the expected number per 10,000 generated sectors, as a distribution, and that none of them is older than `t_life`.

## 9. Stand-ins this design registers (GEN.725)
* `FIELD_STORY_RATE_PER_SECTOR` 1e-4 (tier "from needed output", L4; physical figure 2e-16 beside it).
* Rubble window 15 to 30 km/s (tier "from physics, educated guess", L1; to be replaced by the LS12 law).
* Ejecta dispersion 5 to 15 km/s (L1).
* Supernova belt fraction 0.3 (tier "Marginal Unverified", L1; replace after reading Su 2006 and the Herschel DEBRIS fractions in the origins doc).
* Largest fragment mass in the rubble regime (L1).
* `t_life`: not a constant. It follows from `sigma_v` and the 1 ly retirement radius.
All carry the ledger draft's Solar System rule: no Solar System input; test targets only.

## 10. Reading list
Leinhardt and Stewart 2012 (arXiv:1106.6084, 1109.4588); Su et al. 2006 (arXiv:astro-ph/0608563) and Rieke et al. 2005 [R] for A-star belts; Do, Tucker and Tonry 2018 (arXiv:1801.02821) for the interstellar background; Mroz et al. 2017 and 2020, Johnson et al. 2020 (already cited in tuning) for rogue density; Pierazzo and Melosh style shock-vaporisation work for the 400 GPa to 1 TPa vaporisation thresholds [R]; Breslau et al. 2014 for flyby truncation [R]. Boss cannot access some papers: use second-hand where closed.

## 11. TODO items

| Item | What | Owner (best guess, can move) |
|---|---|---|
| GEN.1098 | F5 research read | Aux 4 |
| GEN.1099 | F1 rogue_collision.py | Aux 4 |
| GEN.1100 | F2 field_dispersal.py | Aux 4 |
| GEN.1101 | F4 field_rates.py | Aux 4 |
| GEN.1102 | F3 composition from two parents | Aux 4 or Sources 1 |
| GEN.1103 | L1 columns, radius from age (lock) | Aux 4, train 19 |
| GEN.1104 | L2 supernova host rule (lock) | Aux 4, train 19 |
| GEN.1105 | L3 members, designations, pages (lock) | Core 2 with Aux 4, train 20 |
| GEN.1106 | L4 render (lock) | Core 2, train 20 |
| GEN.1107 | L5 bring-to-epoch (lock) | Bug 1, train 20 |

## 12. Decisions (Boss, 4:54 PM EDT 2026-10-30)
1. **Rate: option 1, physics only.** Fields come only from physical sources (supernova-unbound belts and the others); no labelled game rate, no floor dial. The same holds for rogue planets (GEN.1109 reviews them).
2. **A field that outgrows 1 ly: retire it.**
3. **The remnant planet of a collision: its own rogue planet.**

## 13. Boss decision (4:54 PM EDT 2026-10-30) and rogue-planet stand-in audit

**Decided by Boss: physics all the way. Option A of section 8 (no game rate, no `FIELD_STORY_RATE_PER_SECTOR`), and rogue planets are to be physics-driven as well.** My two small choices are taken as recommended (a field past 1 ly is retired; the collision remnant is its own rogue planet). Consequences for the items: drop the floor dial from F4 and the L4-tier row from section 9; fields are generated only by the hosted physical channels (S2 supernova-unbound belts, S3 young clusters) and by hand. Section 2's rate stays as the documented physical figure. The rogue-planet rework below is new work for TODO to number.

### 13.1 Where rogue-planet generation uses stand-ins or flat rates today [REPO, read in `tuning.py` lines 2915 to 3130]

| What | Today | Kind | Physics-driven replacement (to research) |
|---|---|---|---|
| Rogues per star | `ROGUE_PLANET_RATE_PER_STAR` 5.8, flat, same in every galaxy region | flat rate | Derive from the planet architecture: count planets the architecture scatters out (GEN.728 already moves and ejects giants) per system, by star mass, metallicity and age, plus cluster stripping; the 5.8 becomes a test target against microlensing (Mroz, Johnson et al.) |
| Mass function | one power law `dN/dlogM ~ M^-0.65`, 0.1 Earth to 13 Jupiter, four bins | fitted slope, bin edges by hand | the ejected-population mass function from the architecture's own ejections (light bodies are ejected more easily) |
| Terrestrial floor 0.1 Earth | hard lower edge | stand-in | ejection of the oligarch chain below `PLANET_MIN_MASS_EARTH`, and the debris class below it |
| Gas-giant threshold | `ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER` 0.031 (10 Earth masses, the bound chain's core limit, GEN.1133) | rule of thumb | derive from the envelope fraction (core mass and gas capture) |
| Age | `ROGUE_PLANET_AGE_RANGE_GY` 0.1 to 10, uniform | default | the age of the ejecting system plus time since ejection, from the galaxy's star formation history |
| Hydrogen envelope kept | `ROGUE_HYDROGEN_ENVELOPE_CHANCE` 0.1 / 0.5; pressure 10 to 1000 bar | defaults | the planet's own formation (nebula gas capture, GEN.8 architecture) and loss by photoevaporation before ejection |
| Water-rich | `ROGUE_WATER_RICH_CHANCE` 0.35 / 0.8, share 1e-4 to 0.1 log-uniform | defaults | formation zone relative to the stored formation snow line of the ejecting system |
| Frozen atmosphere | mass 0.3 Earth cut | default | outgassing and volatile inventory from the body's composition and impact history |
| Tidal heating | `ROGUE_TIDAL_FLUX_RANGE_W_M2` 1e-4 to 1e-2, log-uniform, for rocks that kept a moon | default | whether a moon survived ejection (Hill sphere against the ejection encounter) and its orbit (tides work in GEN.842) |
| Radiogenic inputs | thorium/uranium offset sigma 0.22 dex, potassium range, 26Al factor 0.1 to 10 | placeholders | the ejecting star's stored composition and formation environment (GEN.476) |
| Heat law anchors | Jupiter flux 5.4, Neptune 53 K, cooling `t^-1.3`, photosphere 0.3 bar, 0.1 bar, T cap 2800 K | anchored fits | keep (these are sourced), but the age feeding them must come from the real ejection history |
| Escape parameter | `ROGUE_MIN_ESCAPE_PARAMETER` 15 | rule of thumb | keep (Jeans criterion) |
| Speed | rogue planets take the generic galactic orbit (`set_galactic_motion`); interstellar comets have a flat `INTERSTELLAR_OBJECT_SPEED_KMS_RANGE` | no ejection kick | add the ejection speed from the architecture's encounter (escape speed at the scattering orbit) plus the host's velocity, because the collision rate and the dispersal both depend on it |
| Origin | rogue rows carry no parent | missing | record the ejecting system (or "formed free"), so a rogue has a history sentence and a composition that match the architecture |

The common thread: a rogue is today drawn from a rate table, not made by a system. The physics-driven version is **ejection as an output of system generation**: every ejection the architecture performs becomes a rogue (mass, composition, moons kept, speed, age, parent star), counted per sector with the sector's stars. The flat rate then becomes a check, not an input. That is a lock item (epoch, `--affects`) touching the architecture slot loop, so it comes after the planet architecture and clearance ledger items already on the chain.

### 13.2 Items from this decision
* R1 = GEN.1109 (research, non-lock): ejection yield per star from the architecture (mass, count, speed, moons kept), checked against microlensing rates; output report and stand-in rows.
* R2 = GEN.1111 (lock): rogue rows created from architecture ejections, with parent system, ejection speed, age and kept moons; the rate table becomes a test target. Epoch, migration.
* R3 = GEN.1112 (non-lock): replace the defaults above one by one with the derivations in the right column, each as its own small physics module.
* R4: the field parents (section 7) are drawn from R2's ejection population, and the collision speed window uses R2's speeds.
