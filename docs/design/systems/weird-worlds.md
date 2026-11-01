# Weird worlds: physically plausible oddities for the generator

Status: plan (research and design; Boss 4:15 PM EDT 2026-10-30 said yes to the named Odd line, tilt-driven climate and Venus-type spin first, and the weirdness dial; build items filed in section 6 as the weird worlds goal #4117). Written 4:09 PM EDT 2026-10-30 as the lane reported. No generator code changed. Asked for by Boss (4:09 PM EDT 2026-10-30): "account for weirdness, the very odd things that could happen, strange combinations of factors and surface conditions, thrown in based on the physics so that they aren't fantastical and are plausible."

Related and not repeated here: GEN.1034 (ocean tides, exotic liquids), `docs/design/habitability/life-classes.md` (nine life classes), `docs/design/composition/dry-carbon-model.md`, `docs/design/composition/planet-architecture.md`.

Evidence tags: [S] the paper or object is well known and the number is as published as far as recalled; [R] recalled, no paper opened in this research (web search returned only link lists, no text, so every rate below is [R] unless it comes from a repo file); [C] computed from repo values; [REPO] read in the repo this session. Source tiers use GEN.725 (T1 to T9, L0 to L4). Every number marked [R] must be confirmed by the lane that builds it, and any value that stays unconfirmed is registered as a stand-in (T2 or T3, L0 to L1).

## 1. The idea in one paragraph

Do not paint weirdness on the finished planet. Make each oddity a **cause** that enters at an input of the existing chain (orbit, spin, bulk composition, history, moon set), then let the solvers already built (interior heat, tides, air, hydrosphere, climate, class, life) produce the odd surface. A planet is then weird because of what happened to it, and the class letter, the air and the life scores follow from the physics the same way as for any body. The one exception is a small set of **bistable** states (snowball or warm, runaway or temperate, collapsed or kept air) where the same inputs allow two outcomes, and a history draw picks one.

## 2. What the generator already makes (so we only add what is missing) [REPO]

| Oddity | Already in the generator | Where |
|---|---|---|
| Rocky planet spin axis isotropic, so Uranus-like tilts are common | yes | `spin.py` (obliquity `sin(e)/2`), `impact_spin.py` |
| Tidal locking, 3:2 type resonances | yes | `tidal_spin.py` (GEN.731) |
| Tidal heating from stored eccentricity, melt, heat pipe, icy flux | yes | `tidal_heating.py` (GEN.737) |
| Captured and retrograde moons, impact and ring moons | yes | `moon_composition.py`, `moon_orbits.py` |
| Hot Jupiters, inward migration, scattering | yes | `planet_architecture.py` |
| Radius valley, envelope loss, stripped cores | yes | `radius_valley.py`, `envelope_escape.py` |
| Lava worlds and rock-vapour air | yes | `mixture_lava.py` |
| Carbon-rich dry worlds | yes | `dry_carbon.py` |
| Rogue planets, pulsar and white dwarf planets | yes | `rogue_surface.py`, `pulsar_planets.py`, `white_dwarf_planets` |
| Circumbinary zones, close binary types | descriptor and zone only | `multistar.py` |
| Stagnant vs mobile lid | yes | `interior_heat.py` (GEN.721) |
| Ammonia-water and silicon-rich bodies | yes | `strange_solvents.py` |
| Hycean-like ocean worlds | classes O and N | `planet_classes.json` |
| Ice-albedo climate feedback | not yet (GEN.972) | open item |
| Eccentricity drawn per body | yes, but no high-e tail from scattering or Kozai on purpose | GEN.727 |

Gaps that matter for weirdness: (a) **climate has no obliquity term** (the swing of the starlight follows from eccentricity only and is the stored `flux_variation`; class Q was retired in GEN.1135 [REPO, condition-based-planet-classes.md]), so a 98 degree tilt is stored but changes nothing a reader would see; (b) no Venus-type retrograde end state for thick-air slow rotators; (c) no collisional mantle stripping that gives iron-rich Mercury-type worlds; (d) no resonant-chain heating coupling; (e) no disintegrating or Roche-overflow bodies; (f) no super-puff tail; (g) no abiotic-oxygen path checked end to end; (h) no bistability.

## 3. The catalogue

Columns: what it is and the real examples; the trigger from inputs we already store; the rate (conditional on the trigger unless said otherwise); what it changes. "Cov" says whether the generator already covers it (Y), partly (P) or not (N).

### 3.1 Orbit and spin

| ID | Oddity and real cases | Trigger (stored inputs) | Rate [R] | Effect | Cov |
|---|---|---|---|---|---|
| W1 | **Extreme obliquity seasons.** Uranus 98 deg, Mars chaotic 0 to 60 deg over Myr (Laskar and Robutel 1993), Earth without the Moon (Lissauer, Barnes and Chambers 2012 say milder than first thought) | any rocky planet without a large moon, or in a spin-orbit secular resonance; tilt already drawn | tilts over 60 deg are about 13 percent of isotropic draws [C: `1 - cos 60`] already | seasonal temperature swing from the tilt, polar summer hotter than the equator, wet-bulb and ice limits by season | P (tilt exists, climate ignores it) |
| W2 | **Venus-type retrograde slow rotator.** Venus 243 d retrograde, thermal tide against core-mantle friction gives four spin equilibria for a thick-air planet (Correia and Laskar 2001, 2003, Nature 411 and Icarus) | rocky, air above about 10 bar, orbit where the torques balance, not locked | the equilibrium is reached by some fraction of such planets; the share is unknown, say 25 percent as a stand-in (T4, L3) | slow rotation, no tide lock, very long day on a thick air planet | N |
| W3 | **Highly eccentric orbit.** HD 80606 b e = 0.93, HD 20782 b e = 0.95, Mercury 0.2, Pluto 0.25 | giants after scattering, planets in packed chains, Kozai-Lidov with an outer companion | about 10 percent of RV giants beyond 1 AU have e above 0.5 [R, Winn and Fabrycky 2015 review] | years of 200 K swings, flash-heated atmosphere, the stored `flux_variation` already carries the flux swing (class Q retired, GEN.1135) | P |
| W4 | **Polar and retrograde orbits (stellar obliquity).** hot Jupiters round hot stars (above 6,250 K) are often misaligned, polar peak near 90 deg (Albrecht 2012, 2021); Kepler-56 ; HAT-P-7b retrograde | hot Jupiter, host above 6,250 K, or an outer stellar or giant companion | misaligned fraction of hot Jupiters about a third [R], most for hot hosts | orbital inclination to host equator; moons' stability (Hill); mostly flavour plus a stored angle | N |
| W5 | **Resonant chain with forced eccentricity.** TRAPPIST-1 seven planets in a Laplace chain, Kepler-223, TOI-178, Io-Europa-Ganymede | 3 or more planets with period ratios inside 5 percent of 3:2, 4:3, 5:3, 2:1 | about 1 to 2 percent of Kepler multis are in chains with an excess just wide of resonance [R, Fabrycky 2014; Goldberg and Batygin 2022] | forced e of 0.001 to 0.01 gives Io-like tidal heating for rock planets close to the star (feeds `tidal_heating` as it does for moons) | P (no resonance detection) |
| W6 | **Circumbinary world with two suns.** Kepler-16b, 34b, 35b, 38b, 47b c d, TOI-1338b | close binary host with period under about 40 d, planet beyond the Holman-Wiegert limit | at least 10 percent of close eclipsing binaries host a planet above 6 Earth radii inside 300 d if coplanar [REPO, multi-star-systems.md] | flux varies with the binary phase; double sunset text; zone from the pair | P |
| W7 | **Co-orbital (Trojan) companions and horseshoes.** Jupiter Trojans, Janus and Epimetheus, the asteroid 2010 TK7 | giant planet with a planet or moon at L4 or L5, mass ratio under 0.038 (Gascheau-Routh) | no confirmed exoplanet case; Solar System only | a stored companion at 60 deg; low priority | N |
| W8 | **Double planets and binary dwarfs.** Pluto-Charon (barycentre outside Pluto), Earth-Moon 1:81 | impact or capture origin moon with mass ratio above about 0.1 | rare | mutual tidal lock, both faces fixed | Y (origin "binary") |
| W9 | **Moon spiralling in, then a ring.** Triton retrograde decay, Phobos about 30 to 50 Myr from Roche, ring-moon cycling (Hesselbrock and Minton 2017) | `moon_orbits.py` result with decay inside the planet's life | follows from the existing evolution | a ring from a disrupted moon, a planet briefly with a ring | P |

### 3.2 Bulk composition and interior

| ID | Oddity and real cases | Trigger | Rate [R] | Effect | Cov |
|---|---|---|---|---|---|
| W10 | **Iron-rich Mercury type.** Mercury core 0.7 by mass fraction, K2-229 b, Kepler-107 b and c density contrast (Bonomo 2019), giant impact mantle stripping (Benz 2007, Marcus 2010) | rocky inner chain, high collision count in `planet_architecture`, mass under about 5 Earth | a few percent of super-Earths are iron-enhanced outliers [R, Plotnykov and Valencia 2020; Adibekyan 2021] | bulk density, thin crust, big magnetic field, little air | N |
| W11 | **Coreless or iron-poor world.** composition tracks the host's Fe/Mg (Schulze 2021) with scatter | host Fe/Mg tail | tail of the stellar distribution already drawn | low density, weak dynamo | P |
| W12 | **Carbon planet.** C/O above 0.8 in a few percent of stars (Brewer and Fischer 2016, already in repo), 55 Cnc e carbon claim now doubtful | host C/O tail | about 3 percent of hosts [R] | graphite, SiC, diamond crust, no water | Y |
| W13 | **Lava ocean planet with rock-vapour air.** 55 Cnc e, K2-141 b (supersonic Na, SiO winds, Nguyen 2020), CoRoT-7 b, TOI-561 b dense volatile air over magma (Hu 2024) | locked, equilibrium T above about 1,500 K | covered | night side freezes the vapour; day side evaporates | Y |
| W14 | **Disintegrating planet.** Kepler-1520 b (KIC 12557548 b), K2-22 b, KOI-2700 b: dust tails, mass loss up to about 1 Earth mass per Gyr (Rappaport 2012; Sanchis-Ojeda 2015) | rocky, under about 0.1 Earth mass, period under 1 d, magma ocean | three in a Kepler and K2 sample of order 1e5 stars, so about 1e-5 per star, 1 in 1,000 ultra-short-period planets [R] | tail of dust, deep transit, mass declining, planet dies in Gyr | N |
| W15 | **Super-puff.** Kepler-51 b c d (below 0.1 g/cc, 5 to 9 Earth radii, 2 to 8 Earth masses), Kepler-79, Kepler-87 c, WASP-107 b (Libby-Roberts 2020, Lee and Chiang 2016, Millholland 2019 tides on inflated young planets) | gas dwarf or ice giant, mass under 10 Earth, young age under 1 Gyr, orbit 10 to 100 d, hazes | under 1 percent of sub-Neptunes [R] | huge radius, density 0.03 to 0.1; fades as the star ages | N |
| W16 | **Helium-rich envelope.** after photoevaporation H escapes first, leaving He-dominated air (Hu 2015, Malsky and Rogers 2020) | sub-Neptune, strong XUV history, envelope loss above 90 percent | follows from the existing escape model if it fractionates | He-dominated air, no hydrogen | P (check) |
| W17 | **Very low or very high water.** Kepler-138 d, TOI-1452 b water-rich candidates; desiccated cores | disc ice share tail | already in mass fractions | ocean world vs dust world | Y |

### 3.3 Air, climate, surface

| ID | Oddity and real cases | Trigger | Rate [R] | Effect | Cov |
|---|---|---|---|---|---|
| W18 | **Eyeball and nightside ice-trap worlds.** Pierrehumbert 2011, Joshi 1997, Leconte 2013 (3D, thick air heat transport), Kopparapu 2016 | locked, cool star, ocean | all locked rocky planets in the zone | one open sea under the sub-stellar point, rest ice; or air freezes out on the night side | P (`day-and-night.md`) |
| W19 | **Atmospheric collapse.** N2 or CO2 freeze on the dark side of a locked planet when the air is thin (Joshi 1997, Wordsworth 2011) | locked, pressure under about 0.1 to 0.3 bar | all thin-air locked planets | air lost to nightside caps; day side bare | N |
| W20 | **Snowball and warm bistability.** Budyko 1969, Sellers 1969, Hoffman 1998, Pierrehumbert 2011; runaway greenhouse boundary (Kasting 1993) | ice-albedo feedback, flux within about 5 percent of the inner edge or outer edge of the ecosphere | width of the bistable window about 10 to 20 percent in flux for an Earth twin [R] | same inputs, ice-covered or ocean-covered by a history draw; exit needs a CO2 build-up | N (waits for GEN.972) |
| W21 | **Abiotic oxygen planet.** Luger and Barnes 2015 (water loss in the pre-main-sequence phase of an M dwarf leaves hundreds of bar O2), Wordsworth and Pierrehumbert 2014 (O2-Ar air), Schaefer 2016 (magma ocean oxidation) | M dwarf host, planet in the ecosphere now, pre-main-sequence luminous phase, low non-condensable gas | large share of M dwarf HZ planets lose water; fate depends on sinks [R] | O2 air with no life; oxygen is not a biosignature; life scores must not read O2 alone | P (`biogenic_oxygen`, `mixture_escape`; check path) |
| W22 | **Ultra-hot giant weather.** WASP-76 b iron rain (Ehrenreich 2020), WASP-121 b stratosphere and vanadium oxide, HD 189733 b silicate and glass clouds, HAT-P-7 b corundum clouds, KELT-9 b dayside 4,600 K with H2 dissociated and ionised iron, WASP-12 b egg shape and orbital decay (about -29 ms per yr) | giant at equilibrium T above 1,500 K, locked | all such giants | condensate lists, transit asymmetry, dark side cloud rain; WASP-12 b type Roche overflow at the rarest end | P (cloud table; check rows) |
| W23 | **Dark planet.** TrES-2 b albedo about 0.01 | hot giant with alkali and TiO absorption, no cloud | rare | very dark albedo, hot day side | P (`aerosols`) |
| W24 | **Volcanic, tidal heated worlds.** Io (about 2 W/m2 average, mean surface heat flow 2.5), Europa, Enceladus plumes, Triton nitrogen geysers, TRAPPIST-1 b c; Peale 1979 | forced e from resonance, a close giant, or W5 | follows from W5 and moons | sulfur or sulfate air, molten patches, plumes, sealed oceans | Y for moons, N for planets |
| W25 | **Titan and Pluto type air.** methane cycle on a cold moon, haze-covered | cold, mass and escape allow N2 | derived outputs already | covered | Y |
| W26 | **Hycean worlds.** K2-18 b, TOI-270 d (Madhusudhan 2021; the DMS claim is contested, Taylor, Welbanks and others) | sub-Neptune, 1.5 to 2.6 Earth radii, ocean under hydrogen air | about a third of observed sub-Neptunes in the ecosphere may qualify [R] | covered by O and N plus envelope | P |

### 3.4 Not weirdness, left out on purpose

Fiction rules: warp, subspace, exotic matter; detections still disputed with no physics (exomoon Kepler-1625 b-i is a candidate only, so no rate); Dyson swarms (KIC 8462852 is dust); phosphine on Venus (design flag, not an oddity).

## 4. The mechanism

### 4.1 One table of oddity rows, same shape as the existing tables

A data file `physics/data/oddities.json` read by `physics/oddities.py`, rows like `planet_classes.json` and `moon_composition.py`'s `Possibility` rows. Each row has:

- `id` (W-number kept), `name`, `kind` (planet, moon, giant, rogue, small body, star system);
- `requires`: conditions on stored values only (mass, orbit, age, host type, eccentricity, composition), the same style as `planet_classes.json`;
- `rate`: the chance given the trigger, a number with source, tier and tweak level, registered in GEN.725;
- `apply`: a **named effect** on one input of the chain (for example "set eccentricity to a draw from a Rayleigh tail", "strip mantle fraction f", "set spin to the retrograde equilibrium"), never an output;
- `limits`: the invariants that must hold afterwards (section 5);
- `source` and `text`: one reading sentence for the page ("This world turns backwards because ...").

### 4.2 Where it runs

At the end of the architecture and bulk steps, before the air solver: oddity rows modify inputs, then interior heat, tides, air, hydrosphere, climate, class and life are computed as they are for every body. Bistable rows (W19, W20) run after the climate solver: it hands back both branches and `oddities.pick_branch` takes one by a history draw. The class letter stays derived (`classify` is unchanged, no oddity names a class).

### 4.3 Rate rule: tails, not coin flips

An oddity is a **rate on a gated population**, never a global percentage. Example: W10 applies only to rocky planets under 5 Earth masses in a chain where `planet_architecture` recorded at least one merger, at the literature share of such planets. If the gate is not met the rate is zero. Where the trigger already makes the outcome (W5 resonance gives tidal heating), no extra draw is made.

### 4.4 Randomness and determinism

All draws go through `util.draw` with one **sub-stream per oddity id** keyed on the body, so adding or removing a row never changes any other draw (`test_order_independence` and the golden seed keep holding). A world-level **weirdness dial** (default 1.0, per-oddity override like `PHENOMENON_RATE_SCALE`) scales every rate, 0 turns all off. The generator epoch rises when rows are wired.

### 4.5 What is stored and shown

A text column `oddities` (comma list of ids) on planets and moons, new numeric columns in `db/column_units.py`, and a line in the object page: "Odd: retrograde rotation (Venus type)", with the cause sentence. Admin and API can force an oddity for a body (like `FORCE_CLASS`) for tests and for Boss's own checks.

### 4.6 Sources and stand-ins

Secondary sources accepted (Boss's rule of 8:16 PM EDT 2026-10-29). Every rate and threshold not read from a paper is a GEN.725 row with tier and tweak level. The weakest here: W2 share (T4, L3), W15 share, W20 window width (T6, L3), W4 fraction. Each carries the falsifier "replace when a survey statistic is read".

## 5. Tests: every weird outcome stays physical

A **plausibility harness** (`src/tests/oddity_plausibility_cli.py`, modelled on `physical_plausibility_cli.py` and `phenomena_plausibility.py`) forces each oddity on a sample of bodies (n 200 per row, fixed seed) and asserts the row's invariants plus the general ones:

General, for every body of every oddity:

1. Energy: surface temperature within the band the flux, albedo, greenhouse and internal heat allow; tidal power never above the orbital energy available over the body's age (`tidal_heating`).
2. Structure: radius inside the composition curves for its mass (mass-radius sample) except where the row says how far it may leave them (super-puff, W15, bounded by Lopez and Fortney 2014 envelope scaling).
3. Stability: Hill and mutual Hill separation hold; no orbit below the Roche limit unless the row is W9 or W14 and the body is flagged as dying.
4. Air: pressure under the retention cap (`surface_pressure`), mixture sums to 1, no species above its limit (`species_limits`).
5. Class is derived and total; the life scores do not read an oddity name.
6. No NaN, no negative mass, unit registry satisfied.

Row-specific: W2 spin period over 20 d and sign negative; W5 forced e under the value that keeps the chain stable; W10 core fraction under 0.85 (iron plus thin mantle; no body denser than pure iron at its mass); W14 mass loss rate times age under the body's mass; W16 H mass fraction under 1e-3; W20 ice albedo branch is energy-consistent, and a snowball needs CO2 to exit.

Rate tests (census, in the prevalence suite): the share of each oddity among its gated population within the literature range; and an **isolation test**: turning the dial to 0 changes no draw other than the oddity's own outputs.

**Built (TEST.490, TEST.494).** `src/tests/oddity_plausibility_cli.py` (the harness) and `src/tests/test_oddity_plausibility.py` (its tests). `force(row, number)` draws a gated body with the generator, calls the row's effect from `EFFECTS` with the row's own stream, and `run_row` checks the six general invariants (`general_failures`), the row's own `limits` (`oddities.limit_failures`) and the ring list of 9.3 (`ring_failures`: ring mass under 1e-4 of the planet, outer edge inside the Roche limit to the ring law's margin of 1.1 unless the body carries R4, optical depth within what the formula gives for the ring's mass, edges and kind, no ring inside 0.1 AU, ring plane in the equator, no ring field on a body with no ring). A module that implements a row registers its effect in `EFFECTS`; GEN.1047's wiring then replaces `force` with the generator's own forcing. The stub rows S1 to S3 pass; the broken rows B1 to B9 and R91 to R95 each break one invariant and the tests require that invariant to report. The isolation tests show a roll reads nothing from the run's stream at any dial, the dial at 0 turns the row off without moving its own draw, and one row's draw does not move when another is rolled. Two findings while building it: the generator stores air thinner than 1 Pa (0.015 Pa on non-G2V hosts), so the air invariant is "not negative" and the weight cap, not the 1 Pa floor of `test_planets`; and a gas-partial dict holds pascals that sum to the pressure, not fractions that sum to 1.

## 6. Proposed TODO items (for the TODO thread to file)

Lock means changes generated output and waits for a train (an epoch). Non-lock work merges alone.

### Non-lock (any waiting lane, one at a time on its own branch)

1. **W-core: oddity table and loader** (`oddities.py`, `oddities.json`, row schema, tier and tweak fields, registry rows in GEN.725), no wiring. Tests: schema, every row has a source, gate and limits.
2. **W-harness: oddity plausibility harness** (section 5), run on forced draws with a stub effect; the invariants for the general list first.
3. **W-rates: confirm the [R] rates** for W3, W4, W5, W10, W14, W15, W20 against papers (secondary accepted), write the literature table with sources, register the stand-ins.
4. **W2 module**: Venus-type spin equilibrium (Correia and Laskar 2003 torque balance), pure function, tests on Venus, Earth, Mercury.
5. **W10 module**: collisional mantle stripping (Marcus 2010, Benz 2007 scaling) to a core fraction, pure function, tests on Mercury and Kepler-107 b c.
6. **W14 module**: mass-loss and lifetime of a disintegrating rocky planet (Perez-Becker and Chiang 2013), pure function, tests on Kepler-1520 b.
7. **W15 module**: super-puff radius and fade with age (Lee and Chiang 2016, Millholland 2019).
8. **W5 module**: detect resonant chains from the stored periods and compute forced e (Goldreich and Soter, Murray and Dermott), tests on TRAPPIST-1 and the Laplace moons.
9. **W19 module**: nightside collapse threshold for a locked planet (Wordsworth 2011, Leconte 2013).
10. **W21 audit**: end-to-end check that M dwarf water loss can leave O2 air and that life scores do not read O2 alone (test only, a bug filed if it fails).
11. **W16 audit**: does `envelope_escape` fractionate He; add the fractionation if not (pure module).
12. **Climate obliquity term** (W1): a seasonal swing from the stored tilt in `regional_temperature`; module and tests lock-free, wiring below.
13. **W20 design**: the two-branch climate output (waits for GEN.972 ice-albedo).

### Lock (one train each, after the non-lock parts)

14. Wire the table and the dial; epoch; stored `oddities` column; page line; API field.
15. Wire W2, W10, W15 (bulk and spin inputs).
16. Wire W5 resonance heating and W3 and W4 tails in `planet_architecture`.
17. Wire W14 with its dying-planet flag, plus W9 rings from decayed moons.
18. Wire W1 climate term, W19 and W20 bistable branches.
19. Census re-baseline for every oddity share (TEST item) and the prevalence suite rows.

Overlap checks for TODO: GEN.1034 (ocean tides, exotic liquids) is separate; GEN.972 is a prerequisite for W20; GEN.1006 (overlap-free systems) is a prerequisite for W5 and W7; the planet-architecture weaknesses list (more warm super-Earths wanted) interacts with W10.

## 7. Open choices for Boss

1. Should rare oddities be visible on the page as a named "Odd" line, or only show in the property sheet? Recommended: a named line with the cause sentence, because it teaches the reader the physics.
2. Which first: W1 (climate from tilt, common) or W10 and W2 (rare, striking)? Recommended: W1 and W2 first, because they change every visit to a map.
3. Is a global weirdness dial wanted in settings? Recommended yes, default 1.

## 8. Added 4:17 PM EDT 2026-10-30 (Boss): the core levers, and giants

### 8.1 The core levers (what an oddity may change)

An oddity row may only change one of these inputs. Everything else follows from the solvers.

| Lever | What it is | Oddities that use it |
|---|---|---|
| L1 Spin axis and rate | tilt, retrograde sign, period | W1, W2, G1, G2, G6 |
| L2 Orbit | eccentricity, inclination to host equator, resonance, binary host | W3, W4, W5, W6, W7 |
| L3 Bulk composition | core fraction, water, carbon, ice mix, envelope share | W10, W11, W12, W15, W16, G3, G4 |
| L4 Interior heat and structure | internal flux, tidal power, dynamo, lid regime | W5, W24, G2, G5, G7 |
| L5 History | collisions, migration, mass loss, age at which an event hit | W10, W14, W16, W21, G1, G2 |
| L6 Moon and ring set | moon origin, retrograde moons, rings | W8, W9, G1, G8 |
| L7 Air state (bistable only) | which branch of snowball or warm, collapsed or kept | W19, W20 |

The dial (section 4.4) scales how often a qualifying body gets an oddity. The gates are: host type and mass, body mass and kind, age, orbit, and the history record. So "the core things" are the seven levers plus those gates.

### 8.2 Giant planet oddities (Neptune, Uranus, Jupiter, Saturn type)

| ID | Oddity and real cases | Trigger | Lever | Rate [R] |
|---|---|---|---|---|
| G1 | **Giant-impact tilt.** Uranus 98 deg, with its tilted ring and moon system aligned to the tilt (Slattery 1992; Morbidelli 2012 says a moon-disrupting or two-impact path; Kegerreis 2018) | ice giant or gas dwarf, mass 5 to 30 Earth, outer chain with a recorded merger | L1, L5, L6 | tilt above 60 deg is a minority for giants (Rayleigh 15 deg today); a tail of about 1 in 5 ice giants in the Solar System, unknown elsewhere, stand-in T4 L3 |
| G2 | **Cold or hot ice-giant interior.** Uranus heat flux about 0.04 W/m2 against Neptune 0.43 (Uranus emits almost no internal heat, Neptune about 2.6 times what it absorbs); explanations: layered interior and suppressed convection (Kurosaki and Ikoma 2017; GEN.1033 is already filing the heat spread) | ice giant | L4, L5 | spread already owned by GEN.1033, do not duplicate; oddity row only sets the extreme tail |
| G3 | **Non-dipolar, off-centre field.** Uranus and Neptune fields tilted 59 and 47 deg and offset by 0.3 and 0.55 radii, from a dynamo in a thin ionic-water shell | ice giant | L4 | follows from shell structure; one row for the field shape |
| G4 | **Superionic and diamond-rain interiors.** carbon rain and superionic water (Millot 2018) | ice giant, deep layers | L3, L4 | derived, a text line only |
| G5 | **Fast winds and giant storms.** Neptune 2,100 km/h, Great Dark Spot, Jupiter's Red Spot, Saturn's hexagon and Great White Storm | any giant | L4 | derived flavour from heat and spin, not a draw |
| G6 | **Seasonal giants.** Uranus season 21 years with 98 deg tilt; Saturn tilt 27 deg with ring shadow | tilted giant | L1 | follows from W1 climate term |
| G7 | **Hot Neptune desert and chthonian cores.** few Neptunes inside 3 d (Mazeh 2016), stripped giant cores | giant inside the desert | L3, L5 | covered by `envelope_escape`; oddity is the survivor outlier (LTT 9779 b survives in the desert) |
| G8 | **Odd ring and moon systems.** Saturn's rings (young or old, debated), Chariklo rings, J1407b giant ring system, Neptune's Triton captured retrograde moon, Uranus' dark rings | any giant | L6 | rings partly covered (`rings.py`), captured moons covered |
| G9 | **Super-Jupiters and brown-dwarf edge giants.** HAT-P-2 b, Kepler-39 b, 10 to 13 Jupiter masses, direct-imaged HR 8799 and Beta Pic b in wide orbits | high-metallicity or massive disc | L3 | within existing mass tail; check ceiling |
| G10 | **Helium rain, Saturn extra heat.** Saturn's excess heat from helium sedimentation (Stevenson and Salpeter 1977) | gas giant 0.2 to 1 Jupiter mass, age over 2 Gyr | L4 | derived from interior model |
| G11 | **Dark and reflective giants.** TrES-2 b dark, Kepler-7 b bright clouds | hot giant | L3 | W23 |

Most G rows are derived outputs, not draws. Only G1 and the survivor outliers need a rate row.

## 9. Rings (Boss, 4:17 PM EDT 2026-10-30)

### 9.1 What is built [REPO: `physics/rings.py`, `planetary-rings.md`, GEN.732]

One summed ring system per body (kind, origin, inner and outer radius, mass, optical depth, ice share). Roche limit and snow line are physics; **whether a body has a ring is a draw**: a giant always has one (kinds dusty : narrow : broad = 1 : 2 : 1, the Solar System's four), a rocky or dwarf planet 1 percent, rogues too. All numbers are stand-ins BC20 to BC23. The painter (MAP.275) draws from the stored fields.

### 9.2 What is missing, and the fix: rings follow from history, as oddities

The ring draw is independent of the planet's history, so every giant has a ring and rings never explain anything. Make a ring an **outcome of the L6 lever (moon and ring set)** and let oddities use it:

| ID | Ring oddity and real cases | Trigger | Rate [R] |
|---|---|---|---|
| R1 | **Ring from a disrupted moon.** Saturn (Canup 2010; young rings 100 to 400 Myr per Cassini mass, contested), Phobos about 30 to 50 Myr from disruption, ring-moon cycling (Hesselbrock and Minton 2017) | moon inspiral below the Roche limit (W9) | follows from `moon_orbits`; replaces the always-a-ring rule |
| R2 | **Ring from a captured comet.** Dones 1991 | comet passing inside the Roche limit | rare |
| R3 | **Rings aligned with a tilt.** Uranus rings and moons lie in the 98 deg equator (G1) | giant-impact tilt | by construction: ring plane = equator |
| R4 | **Rings beyond the Roche limit.** Quaoar ring at 7.4 radii (Morgado 2023), Haumea, Chariklo, Chiron | small body with a spin-orbit resonance to a satellite or fast spin | rare, no rate; stand-in T4 |
| R5 | **Shepherded and eccentric rings, arcs.** Uranus epsilon ring, Neptune's Adams arcs held by Galatea, Saturn's F ring | narrow ring plus an inner moon | derived flavour |
| R6 | **Dark or bright rings.** Saturn 90 to 99 percent ice, Uranus and Neptune dark | snow line, ring age, pollution | derived from ice share and age |
| R7 | **Giant ring systems on super-Saturns.** J1407b about 0.6 AU (Mamajek 2012); rings as a cause of super-puff radii (Piro and Vissapragada 2020, now disfavoured for Kepler-51) | low-mass young giant with many moons | a few named cases, no rate |
| R8 | **Ring shadow climate.** ring shades the winter hemisphere and the moons | broad ring plus tilt | derived from W1 plus the ring |
| R9 | **Ring lifetime near a star.** rings of close planets are short lived (Schlichting and Chang 2011) | planet inside about 0.1 AU | gate: no ring |
| R10 | **Propeller moons and gaps.** Saturn's propellers, Encke gap | broad ring with an embedded moon | painter flavour |

### 9.3 Changes this implies

- The ring draw becomes conditional on moon history (R1, R2), so a giant with no disrupted moon or capture has no ring, and Saturn-like broad rings are the tail. "Four of four in the Solar System" is a single-system anchor, so BC20 is registered T4 L3 and replaced by Cassini ring age and disruption rates.
- A ring stores its **age**; ice share and darkness follow from pollution over time (R6).
- Ring plane is the equator, so oddities that move the spin (G1, W1, W2) move the ring.
- Invariants for the harness: ring mass under 1e-4 of the planet, ring inside the Roche limit unless R4 names its resonance, optical depth from the formula (an output), no ring inside 0.1 AU.

**Built (GEN.1057, wiring).** R1, R2, R4 and R9 are wired in `rings.apply` (a ring from `moon_history`'s fallen moons and `ring_origin`'s comet passages; a dwarf planet's resonance ring; no ring inside 0.1 AU); a ring stores its age (`ring_age_yr`, schema v187). BC20 and BC22 are replaced (kind from mass, `rings.KIND_MASS_LIMITS_KG`). R6 (GEN.1055) is wired in `rings.apply` (`rings.darken`): the stored `ring_ice_fraction` is the share the ring has now, after the dark dust that fell on it over its age (`physics/ring_darkening.py`, stand-ins RGD1 to RGD4, section 6.148 of the gap-filling plan), so a thin old ring is dark and Saturn's rings keep over 90 percent ice; the reflectance follows from the share (`ring_darkening.albedo`). R8 (GEN.1055), `ring_darkening.shadow_transmission` (the share of sunlight that reaches a latitude through the rings), is read by `latitude_bands.ring_shade`: each band of a ringed, tilted rocky planet keeps the insolation-weighted mean of that share over its latitudes and the year, and its temperature follows. The world's own mean is not lowered (it is solved before the ring is drawn), so the shade moves heat between the bands, cooler behind the ring and warmer at the poles; a moon's ring does not shade its bands, because a moon's bands are made before its ring is drawn.

### 9.4 Added TODO items

Non-lock: (20) ring origin from moon decay and comet passage, pure function, tests on Saturn and Phobos; (21) ring age and pollution darkening; (22) R4 resonance rings for small bodies, Quaoar and Haumea tests; (23) ring invariants in the plausibility harness. Lock: (24) wire conditional rings (epoch, replaces the BC20 draw); (25) wire ring age, R4, R9; (26) census re-baseline for ring shares.
