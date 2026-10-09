# Population and politics: species, civilizations and territories

Design for POP.1 to POP.4 (`docs/TODO.md`, "Population and Politics"),
extended on 2026-10-09 with the design for POP.7 to POP.10: tech levels
(the tech-level draw, its parameters, band names and the age definition),
the limiting-domain and lopsided-species handling, and facility types with
ratings.

Informs: POP.7, POP.8, POP.9, POP.10 (new sections); POP.1 to POP.6 are built.

**Status (2026-10-01, checked against 7.58.2):**

| Part | Built in |
|---|---|
| The pass, the four tables (schema v44) and the read API | 7.49.0, PR #169 |
| 1 in 1,000 civilization chance | 7.49.0, PR #169 (commit b62d0aa) |
| Galaxy Map Territories overlay | 7.54.0, PR #179 |
| Overlay hidden until a polity exists | 7.58.1, PR #181 |
| Optional and off by default; `GET /api/population` | 7.58.2, PR #180 |
| Web pages (Species list and page, polity page, "Dominant species", "Territory of ...") | not built yet |
| Tech levels (POP.7 to POP.9) and facility types (POP.10) | designed here (2026-10-09), not built |

Boss's decisions are recorded in "Decisions", and the reasons for the main
choices in "Why it works this way".

## What existed before

- Every habitable-zone planet (never a moon) gets an evolutionary
  timeline (`evolution.get_evolutionary_timeline`), stored as paragraphs
  in `planet_evolutionary_paragraphs`. The furthest milestone it reached
  is one of `abiogenesis`, `photosynthesis`, `complex_cells`,
  `multicellularity`, `technological_civilization`, recoverable with
  `evolution.life_stage_from_paragraphs`. Planet classes cap the
  milestone (`PLANET_CLASS_MAX_LIFE_STAGE`).
- Each timeline has a pace (`fast`, `normal`, `slow`); the milestone ages
  for each pace live in `program_constants.EVOLUTIONARY_TIMELINES`, and
  the host star's age is stored on `stars.age_gy`.
- A system's galactic position is its sector's `center_*_pc` plus the
  system's own `position_*_mpc` / 1000.

Before v44 nothing named a species, dated a civilization or claimed a
system.

## The model

### Life worlds and their dominant species (52)

A **life world** is a planet whose milestone is `multicellularity` or
`technological_civilization` (`tuning.LIFE_WORLD_STAGES`): complex
life. As first built (7.49.0) every life world got a named dominant
species. Since GEN.80 (PR #448) only a world with a technological
civilization gets a species row (see "Civilization age"); a world with
complex life and no civilization, and the simpler biospheres (microbes,
mats, single cells), keep their life chemistry and get none. The pass reads the milestone and both ages (the
system's age and the milestone's) back out of the stored paragraph text
(`population.parse_timeline`); a planet whose text can't be read is
skipped.

Each world that gets a species gets exactly one row in `species`: a generated name, the
homeworld planet and system, its biochemistry (`life_chemical`), the
milestone it reached, and a few traits derived from the homeworld so a
species reads like it belongs to its planet:

| Trait     | From                                   | Values                              |
|-----------|----------------------------------------|-------------------------------------|
| `build`   | surface gravity                        | gracile (< 0.7 g), medium, robust (> 1.4 g) |
| `climate` | surface temperature                    | cold (< 250 K), temperate, hot (> 320 K) |
| `size`    | a random draw weighted by `build`      | small, medium, large (a robust build leans small, a gracile one large) |

A missing gravity counts as 1 g and a missing temperature as 288 K.

Names come from the star-name generator (`generate_phoneme_salad_name`
over `STAR_NAMES`, `STAR_PREFIXES` and `STAR_SUFFIXES`, at most 10
letters, checked by `is_name_valid`, which applies `offensive_words.txt`).
Species names are unique across the galaxy, like stars and systems: after
50 colliding draws a number is appended ("Voranthis 2").

### Civilization age (54)

The timeline's `technological_civilization` milestone only says the
window for a civilization has opened, and it opens on about one generated
system in seven (measured on 150 generated sectors: 235 such worlds among
1,639 systems). Taken literally that would fill the galaxy with empires, so a world past
the milestone has a civilization now only with `CIVILIZATION_CHANCE`
(default 1 in 1,000, moved by the generating run's intelligent-life
prevalence, `system_configs.prevalence_intelligent_life`), or always when
its system was generated with intelligent life forced on
(`system_configs.intelligent_life`). That gives about one civilization per
7,000 systems. Every other world past the milestone gets no species: since
GEN.80 (PR #448) species, polities and population are made only for
worlds with a technological civilization, and a population pass removes
any stored species without one.

**Spacing.** At the disk's local density (about 0.1 systems per cubic
parsec) one civilization per 7,000 systems is 1.4e-5 civilizations per
cubic parsec [C]. That is a mean spacing of 134 ly (n^(-1/3)), a mean
nearest-neighbor distance of 74 ly, an 82% chance that the nearest other
civilization is within 100 ly, and 1.7 other civilizations inside a 100 ly
sphere on average. This document (and the docstring of
`tuning.CIVILIZATION_CHANCE`, which still says so) used to say neighbors
sit "roughly 150 to 200 ly apart" and that territories at the 100 ly cap
"just about meet". That was optimistic: 150 to 200 ly is above the mean spacing, and territories at
the cap overlap on average instead of touching. The overlap is resolved by
the weighted Voronoi split below (`reach / distance`). Only the spacefaring
civilizations (about 83%) hold territory, and a polity reaches the cap
only at about 800,000 years (25 ly at 50,000 years, 50 ly at 200,000), so
younger polities are smaller; about 46% of civilizations (the Ancient and
Elder eras) sit at the cap [C, from the 10,000-species run under "Tech
levels"].

A civilization is younger than its window. Its age is a log-uniform draw
between 100 years (`CIVILIZATION_MIN_AGE_YEARS`) and the window (`system
age - milestone age`, both read back from the stored paragraph), so ages
spread evenly across orders of magnitude and none is older than its star
allows. A window shorter than 100 years (the stored ages are rounded to
10 million years, so this is rounding) gives 100 years.

Age is counted in years since the industrial transition, the point where
the Energy index crosses about 1.5 to 2 ("Tech levels", "Age: the
definition"). The era thresholds and the tech-level curve both depend on
that reading.

Age maps to an **era**, which is what the pages show and what drives
everything else:

| Era            | Age (years)          | Spacefaring | Notes |
|----------------|----------------------|-------------|-------|
| Industrial     | < 300                | no          | one world, spaceflight just beginning (Earth is 266 years past its industrial transition) |
| Interplanetary | 300 to 2,000         | no          | settles its own system |
| Interstellar   | 2,000 to 50,000      | yes         | first colonies on nearby stars |
| Established    | 50,000 to 1 million  | yes         | a stable territory |
| Ancient        | 1 million to 100 million | yes     | at the reach cap |
| Elder          | > 100 million        | yes         | at the reach cap; reads as legendary |

The era thresholds live in `program_constants.CIVILIZATION_ERAS` so they
can be tuned without a migration: every pass recomputes `species.era` and
`species.spacefaring` from the stored age.

### The species database (53)

**Spacefaring** means era Interstellar or later. The species database is
the `species` table itself, filtered by `spacefaring = 1`; non-spacefaring
species and plain life worlds stay in the same table so a later pass can
promote them if eras are retuned.

### Governments and territories (51)

Each spacefaring species founds one **polity** (`polities` table): its
name is the species name plus a form of government picked stably from
the species id out of 14 (`program_constants.GOVERNMENT_FORMS`:
Hegemony, Concord, Republic, Union, Directorate, Dominion, Assembly,
Commonwealth, Collective, Compact, Federation, Sovereignty, Ascendancy,
League), its capital is the homeworld system, and its map color is a
golden-ratio hue from the species id, so neighboring ids differ.

A polity's **reach** grows with the square root of its age, like a
diffusing frontier:

    reach_ly = min(TERRITORY_REACH_CAP_LY, 5 * sqrt(age_years / 2000))

(5 is `TERRITORY_BASE_REACH_LY`, 2,000 years the start of the
Interstellar era), so a new interstellar civilization holds its neighbors
within 5 ly, one 50,000 years old 25 ly, one 200,000 years old 50 ly, and
anything older than about 800,000 years sits at the cap (default 100 ly).
Reach is 0 before a species is spacefaring.

Every generated system within some polity's reach is **owned** by the
polity with the strongest claim, `reach / distance` (a weighted Voronoi
split: at the same distance the polity with the longer reach wins). An
exact tie goes to the lower polity id. A capital always owns its own
system, and a polity whose capital was never placed in the galaxy (a
standalone system, or a sector made with no galaxy) holds only its
capital. Ownership is stored per system in `system_owners`
(`star_system_id` primary key, `polity_id`, `distance_ly`), which is what
the 3D territory overlay draws: each owned system is a point colored by
its polity, and the polity's reach is a sphere around its capital.

### When it runs

`planetgen population` is a separate pass over what is already stored
(`population.run_pass`), so it works on a galaxy generated before v44
with no regenerate:

1. `scan_life_worlds`: scan planets above a watermark on `planets.id`
   (kept in `population_state`), 5,000 ids per query, recover each one's
   milestone, and add `species` rows for new life worlds. Traits, the
   civilization draw and the age come from a generator seeded by the
   planet id; names come from the shared name generator.
2. `refresh_civilizations`: recompute every species' era from its stored
   age (cheap; picks up retuned thresholds), dissolve the polity of any
   species no longer spacefaring, found polities for newly spacefaring
   species, and update every polity's reach.
3. `refresh_territories`: recompute ownership from scratch. Each polity
   looks only at the sectors its reach touches
   (`_db.sectors_reached_by`), then the systems in them, so the cost is
   the sum of the territories, not polities times all systems.

Flags: `--rescan` (forget every species, and with them every polity and
territory, and scan every planet again) and `--territories-only` (only
step 3). The two can't be combined.

Population is optional and off by default (Boss, 2026-10-01): nothing
runs it unless asked. `planetgen galaxy` and `planetgen sector` run
the whole pass after they save only with `--population` (this replaced
7.49.0's `--no-population`); the admin Generate page's jobs don't pass
it. `install.sh` (and `install.ps1`) offer to
run it after the database step (`offer_population_pass`,
`Invoke-OptionalPopulation`), y/N with a 30-second timeout defaulting to
No and skipped with no terminal; `POPULATION=1` (`-Population` on
Windows) runs it without asking. Nothing in system generation itself
changes.

## Tech levels (POP.7, POP.8, POP.9)

Not built. This section is the POP.8 design: how the six domain indices
are drawn from a species' age, era and world. POP.9 builds it. Boss's
source is [Technological Assessment.md](Technological%20Assessment.md)
(uploaded, not edited); "TA" below cites its sections. Evidence tags: [S]
seen in search results, [C] computed, [R] recalled and unconfirmed (see
"Evidence notes").

### What the source says and what the research changed

Carried over unchanged from TA:

- Six domains, each scored as an index from 0 to 7 by named technology
  bands (TA, "Domain Index Specifications"): Energy (EI), Materials (MI),
  Information (II), Medical (MeI), Propulsion (PI), Defense (DI).
- The weighted formula (TA, "Domain Weighting System"), weights summing
  to 1: `TL = 0.25 EI + 0.20 MI + 0.20 II + 0.15 MeI + 0.10 PI + 0.10 DI`.
- The reason for the weights (TA, "Weighting Rationale"): Energy sets the
  limit of all work; Materials is the structural bottleneck; Information
  governs analysis and automation; Medical is mastery of organisms;
  Propulsion and Defense are applications downstream of Energy and
  Materials. TA's worked example (E4 M3 I6 Me5 P4 D3 gives TL 4.25) was
  recomputed and is right [C].

Changed or added by the research (2026-10-09):

| Topic | TA | Now |
|---|---|---|
| Weights | Stated as a priority order | Kept as the headline formula. They barely matter in practice: equal weights move TL by 0.03 on average (0.27 at most) and keep the same ranking (Spearman 0.9986) on 10,000 simulated species [C]. No calibration work is needed. |
| Index values | A band per index, no rule for fractions | Stored as integers 0 to 7. A fractional value in the draw means: index k is routine technology in the leading society (TRL 9), and k + f means a fraction f of band k + 1's headline items are demonstrated (TRL 6 or higher). Admins can use the same rule by hand. |
| Calibration | None | Earth in 2026 scores E 3.0, M 3.0, I 3.4, Me 3.2, P 3.0, D 3.2, so TL about 3.1 [judgment against the bands, R]. TA's example (TL 4.25) is a 22nd to 23rd century Earth, not an exotic alien. Index 7 is a rarely reached ceiling. |
| Age | Not mentioned | `civilization_age_years` is defined as years since the industrial transition (below). |
| What the sum hides | Says weights guard against "jump gates with no metallurgy" | The sum charges that species only 0.6 TL points. Handled by a limiting-domain value and a lopsided flag (below). |
| Text problems | Pasted LaTeX debris in the title and first paragraph ("TL}TL\text{TL}TL"), "an species", Defense 5 (Gauss guns exist as prototypes) is ahead of readiness of Defense 4 (lasers being fielded) | Source left alone; the ladders below are cleaned copies. Earth sits at D about 3.2 and Defense is not strictly ordered by readiness. |

### The six scales

Condensed from TA, "Domain Index Specifications" sections 1 to 6. The
source file keeps the full lists.

| Index | Energy (TA 1) | Materials (TA 2) | Information (TA 3) | Medical (TA 4) | Propulsion (TA 5) | Defense (TA 6) |
|---|---|---|---|---|---|---|
| 0 | Muscle, open fire | Wood, stone, bone, hides, raw clay | Oral tradition, tally marks, cave art | Poultices, bone setting, cautery | Walking, pack animals, rafts | Wood and stone weapons, hide padding, wicker shields |
| 1 | Draft animals, waterwheels, windmills | Bronze, forged iron, high-carbon steel, fired brick | Alphabetic writing, arithmetic, scrolls, codices | Surgery without sterility, early pharmacology | Wheeled carts, multi-masted ships | Metal blades, mail and plate, siege engines |
| 2 | Coal steam, DC generation, petroleum engines | Structural steel, cast iron, vulcanized rubber, early plastics | Movable type, telegraph, analog radio | Anesthesia, germ theory, antibiotics, blood typing | Locomotives, motor vehicles, submersibles, atmospheric flight, chemical rockets | Black powder, rifled artillery, ironclads, early armored vehicles |
| 3 | Fission, AC grids, gas turbines, chemical cells | Precision plastics, titanium alloys, industrial composites | Digital computing, packet networks, cellular | Transplants, CT/PET/MRI, antivirals, diagnostic beds | Jets, nuclear vessels, multistage orbital rockets | Automatic firearms, guided missiles, Kevlar |
| 4 | Controlled fusion, micro fuel cells, room-temperature superconductors | Carbon nanotube structures, grown supermaterials | Early true AI, VR networks, neural interfaces | Gene selection, tissue engineering, basic cybernetics | Transit grids, space elevators, gravitic lift, sublight interplanetary ships | Electrolasers, directed-energy cannons, powered battlesuits, combat robots |
| 5 | Gravity-induction and mass reactors, He-3 fusion networks | Molecularly woven composite armor, nanotech matrices | Hyper-immersion networks, system-wide neural grids, biodroids | Cross-species splicing, artificial wombs, full bionic replacement, automed pods | Gravitic hovercraft, non-tactical warp, static jump gates | Gauss guns, particle beams, personal force screens, nanotech armor |
| 6 | Antimatter annihilation, micro-fusion cells, local energy taps | Liquid-state alloys, crystal carbon armor, isotope-encoded matter | Atomic-scale memory, self-directing AI, instant subspace relays | Cellular rejuvenation, mind uploading, automed regeneration, nanite symbionts | Tactical warp, quantum slipstream, short-range transporters | Compact phasers, molecular disruptors, ship-wide deflectors |
| 7 | Zero-point extraction, cosmic energy taps, mass-to-energy conversion | Living metal, programmable matter, forcefield hybrids | Subatomic computation, quantum god-minds, temporal-gate processing | Regeneration rays, metamorphosis, age reversal, immortality | Space-folding drives, chronal engines, planetary wormhole networks | Disintegrators, nuclear dampers, black-hole bombs, living-metal plating |

Cross-checks against real frameworks [S unless marked]:

- Energy against Kardashev/Sagan: K = (log10 P - 6) / 10 with P in watts;
  Earth, 620 EJ in 2023, is 19.6 TW and K = 0.729 [C]. Using Cook's
  per-capita ladder (about 4,000, 12,000, 26,000, 77,000 and 230,000
  kcal/day from hunter-gatherer to 1970s USA [R, secondary summaries]),
  a factor of about 3.9 per index fits Energy 0 to 3; extrapolated to 7
  that is 2.7 MW per person, 2e16 W for 8e9 people, K = 1.03. So the
  top Energy index is only a Type I civilization. Old ("Elder")
  civilizations are physically far beyond that but saturate at 7 here. If
  they should feel different, add a derived power display (Sagan K), not
  more indices.
- The Energy ladder describes the leading society, not the world average
  (Earth's mean is 2.4 kW per person, below the Energy 2 value).
- Barrow's scale (smallest scale a civilization controls) agrees with
  Materials: 4 is atomic-precision manufacturing, 6 to 7 subatomic.
- The UNDP Human Development Index is the precedent for combining
  normalized dimensions with a geometric mean, so a very low dimension
  drags the total (see the limiting domain below).

### Age: the definition

`civilization_age_years` is a log-uniform draw from 100 years up to the
window since the milestone, and the era thresholds only make sense if age
zero is the industrial transition: Earth is 266 years past 1760, at the
Industrial/Interplanetary boundary with spaceflight just starting. So:

> **Civilization age is the number of years since the industrial
> transition, the point where the Energy index crosses about 1.5 to 2.**

On the fitted curve below, Energy reaches 1.5 at about 21 years and 2.0 at
about 51 years [C]. Pre-industrial species never occur (minimum age 100
years), which matches the shipped behavior. The curve still works down to
about 10 years, if admin-made pre-industrial species (TL 0 to 2) are ever
wanted.

### Drawing the indices

Each index is a logistic in log10 of effective age, because time spent
per technological tier falls by orders of magnitude (about 3.3 million
years of Stone Age, about 2,100 of Bronze Age, 600 to 700 of Iron Age,
about 80 for the Industrial Revolution [S]; about 5,800 years at Energy 1
against about 180 at Energy 2 [R]; world population doubling times of
about 123, 47 and 48 years between 1804, 1927, 1974 and 2022 [C from S
dates]). A logistic in linear time has the wrong shape.

For each species, from its stored age `A` in years, using
`draw.Stream(f"tech:{planet_id}")` (a stream of its own, see "Fit with the
code"):

1. **Effective age.** Collapse events are a Poisson process, 0.10 per
   decade of age beyond 100 years. At each event the effective age is
   multiplied by a draw between 0.01 and 0.30 and then keeps growing with
   real time. With 12% probability progress freezes at a plateau age drawn
   between 10^2.7 and 10^6 years.
2. **Shared pace.** `z ~ N(0, 0.30)` decades, the whole species being
   quicker or slower, plus trait shifts (below).
3. **Raw index** for domain `d`:
   `7 / (1 + exp(-k * (log10(A_eff) - m - shift_d)))`, where
   `shift_d = offset_d + z + t(4) * 0.09`.
4. **Soft ceiling.** A species ceiling `c ~ N(6.2, 0.6)` clipped to 4 to
   7, a domain ceiling `c + N(0, 0.4)`, and a soft minimum of the raw
   index and the ceiling (sharpness 5). Physical limits are universal; what
   varies is how much of the technology tree a species finds.
5. **Prerequisite clamps**, in this order: Information at most
   Materials + 1.5; Medical at most min(Information, Materials) + 1.5;
   Propulsion at most min(Energy, Materials, Information) + 1.0;
   Defense at most min(Energy, Materials) + 1.0.
6. **Lopsided species.** With 3% probability one domain is raised by 1.5
   to 3.0 after the clamps (capped at 7). This is the only way large
   spreads appear.
7. Round to integers 0 to 7; `tech_level` is the weighted sum.

Domain offsets (decades; negative = earlier) are solved so the median
vector at 266 years equals Earth's: Energy +0.066, Materials +0.066,
Information -0.171, Medical -0.053, Propulsion +0.066, Defense -0.053.
Information and Medicine lead the calendar and Energy, Materials and
Propulsion lag, as Earth's history shows; the prerequisite clamps give
Energy and Materials the structural lead.

| Parameter | Default | Meaning |
|---|---|---|
| `k` | 0.97 | Index gained per decade of log age at the midpoint |
| `m` | 2.655 | Centre, log10 years (about 450 years); the curve is 7 x logistic, so 3.5 there |
| `sigma_shared` | 0.30 decades | Species-wide pace spread (a factor of 2 in time per sigma) |
| `sigma_domain` | 0.09 decades, t(4) | Per-domain noise, heavy-tailed, sd about 0.13 decades |
| `ceil_mean`, `ceil_sd`, `ceil_dom_sd`, `soft` | 6.2, 0.6, 0.4, 5.0 | Per-species ceiling |
| `p_stagnate`, plateau range | 0.12, 10^2.7 to 10^6 years | Frozen progress |
| `collapse_rate`, `collapse_keep` | 0.10 per decade, 0.01 to 0.30 | Collapses, share of effective age kept |
| `p_spike`, spike size | 0.03, +1.5 to +3.0 | Lopsided species |
| Prerequisite margins | I: M + 1.5; Me: min(I, M) + 1.5; P: min(E, M, I) + 1.0; D: min(E, M) + 1.0 | Clamps |

Trait hooks are all 0 until the trait exists (today `species` stores
only `build`, `climate` and `size`): homeworld harshness from the
habitability tier (GEN.89, GEN.92; -0.15 decades x harshness 0 to 1,
default 0.5 with no effect until it exists), resource richness (-0.10 on
Energy and Materials), tool-use dexterity (-0.08 on Materials and
Propulsion), militarism (-0.12 on Defense, where "defense lags in peaceful
species" lives) and lifespan (-0.06 on Medical). These coefficients are
design levers with no empirical basis, small on purpose (at most a 40%
change in effective time).

**Fit.** Seven anchors (Earth 2026 at age 266 plus the era boundaries)
against the plain logistic, before noise and ceilings, RMSE 0.17 TL [C]:

| Age (years) | 100 | 266 | 300 | 2,000 | 50,000 | 1 million | 100 million |
|---|---|---|---|---|---|---|---|
| Anchor TL | 2.2 | 3.1 | 3.3 | 4.8 | 6.0 | 6.5 | 6.9 |
| Curve | 2.4 | 3.1 | 3.2 | 4.6 | 6.2 | 6.7 | 7.0 |

The soft ceiling then lowers the old-age end to about 6.1; the anchors at
6.5 and 6.9 are deliberately sacrificed so old species are not all 7.

**Results** on 10,000 species using the real `civilization_age` (windows
drawn uniformly from 0 to 7.5 Gy, rounded to 10 My like the stored text;
seeds 1 to 3 agree within 1 percentage point per band) [C]:

| Era | Share of species | TL p5 | TL median | TL p95 |
|---|---|---|---|---|
| Industrial (< 300 y) | 6.3% | 2.03 | 2.92 | 3.92 |
| Interplanetary | 10.7% | 2.92 | 4.02 | 4.98 |
| Interstellar | 18.8% | 4.29 | 5.41 | 6.11 |
| Established | 17.8% | 4.98 | 6.04 | 6.54 |
| Ancient | 27.8% | 4.96 | 6.10 | 6.73 |
| Elder | 18.6% | 4.98 | 6.12 | 6.77 |

Median index by age (3,000 draws each) [C]:

| Age (y) | E | M | I | Me | P | D | TL |
|---|---|---|---|---|---|---|---|
| 100 | 2.46 | 2.45 | 2.84 | 2.63 | 2.44 | 2.63 | 2.59 |
| 266 (Earth) | 3.10 | 3.10 | 3.50 | 3.30 | 3.09 | 3.29 | 3.23 |
| 1,000 | 4.08 | 4.08 | 4.45 | 4.26 | 4.08 | 4.27 | 4.20 |
| 2,000 | 4.54 | 4.53 | 4.88 | 4.70 | 4.53 | 4.71 | 4.64 |
| 5,000 | 5.06 | 5.06 | 5.32 | 5.20 | 5.05 | 5.20 | 5.15 |
| 50,000 | 5.87 | 5.87 | 5.98 | 5.93 | 5.86 | 5.92 | 5.88 |
| 1,000,000 and older | 6.1 | 6.1 | 6.1 | 6.1 | 6.1 | 6.1 | 6.1 |

The median at 266 years sits above Earth's 3.13 because the curve is
convex there and the noise is symmetric in time. By integer-index band the
shares are: TL 2 to 3 3.0%, 3 to 4 7.5%, 4 to 5 11.6%, 5 to 6 34.0%, 6 to 7
40.6%, nothing below 2. The set is top-heavy because the log-uniform age
puts 64% of civilizations above 50,000 years (next section), not because
of the tech curve. Indices correlate 0.92 to 0.94 across the whole set and
0.76 to 0.79 at a fixed age. 30% of species have had a collapse (50% of
Elder), 8% are stagnated, and 5% of Elder end below TL 5. Adding 25%
fast-pace windows (milestone at 0.08 Gy, window under about 20 My, often
rounding to age 100 years) raises Industrial to 12% and cuts the
spacefaring share from 83% to 76%; that tail is an artifact of the 10 My
rounding.

Propulsion reaches 4 (sublight interplanetary ships) at a median age of
1,040 years and 5 (warp, jump gates) at 4,600 years. The shipped
Interstellar era (from 2,000 years) therefore corresponds to Propulsion
about 4.5, slow sublight colonization, consistent with a 5 ly reach that
then grows with the square root of age. No era retune is needed, but warp
drive is not what defines "Interstellar" here.

### Is the age prior reasonable?

Share of civilizations per era under physically motivated priors, on the
same windows [C]:

| Age prior | Industrial | Interplanetary | Interstellar | Established | Ancient | Elder |
|---|---|---|---|---|---|---|
| Log-uniform (shipped) | 6.6 | 11.1 | 18.8 | 17.7 | 26.2 | 19.7 |
| Uniform births, no extinction | 0.1 | 0.0 | 0.0 | 0.1 | 7.0 | 92.9 |
| Uniform births, mean life 100 ky | 0.4 | 1.7 | 37.8 | 60.1 | 0 | 0 |
| Uniform births, mean life 1 My | 0.1 | 0.2 | 4.6 | 57.9 | 37.2 | 0 |

With a constant birth rate and long or no lifetimes, an average
civilization is very old and Earth-like ones are rare coincidences. The
shipped prior guarantees variety and is kept. Consequences: 83% of
civilizations are spacefaring; the top-heavy TL distribution follows from
the prior; and a designer who wants fewer godlike species should lower
`ceil_mean` (not simulated at other values), not change the age prior,
which also drives territory sizes.

### The weighted sum hides a spike: limiting domain and lopsided species

TA motivates the weights by saying a species with jump gates and no
metallurgy is "fundamentally constrained by its resource base", but a flat
weighted sum does not enforce that. Weighted, geometric and harmonic
means (the last two on index + 1, minus 1) and the minimum [C]:

| Species | Weighted | Geometric | Harmonic | Min |
|---|---|---|---|---|
| TA example (E4 M3 I6 Me5 P4 D3) | 4.25 | 4.14 | 4.04 | 3 |
| Earth 2026 | 3.13 | 3.13 | 3.12 | 3 |
| Jump gates, no metallurgy (E3 M0 I3 Me3 P5 D3) | 2.60 | 2.16 | 1.55 | 0 |
| One 7 among 2s (Information) | 3.00 | 2.65 | 2.43 | 2 |
| All 3 | 3.00 | 3.00 | 3.00 | 3 |
| E, M at 2, the rest at 6 | 4.20 | 3.78 | 3.38 | 2 |
| E, M at 6, the rest at 2 | 3.80 | 3.39 | 3.04 | 2 |

The last two rows are the problem: advanced applications on a primitive
base outrank the opposite. In the generated set it rarely arises: the
weighted-to-geometric gap has mean 0.006 and maximum 0.13 TL, and only
0.9% of species spread 2 or more between their highest and lowest index,
all from the 3% lopsided ones. The prerequisite clamps do the structural
job the weights cannot, so this matters mostly for admin edits and the
lopsided species.

Handling (recommended):

- `tech_level` is TA's weighted sum, unchanged (it is Boss's formula).
- **Limiting domain**: the lowest of Energy, Materials and Information
  (the foundation domains; ties resolve in that order). Derived from the
  stored integers when displayed or served, not stored.
- **Lopsided**: true when the highest index minus the lowest is at least
  3. Also derived.
- Optional `tech_balanced`: the weighted geometric mean of (index + 1),
  minus 1, shown beside TL if Boss wants one alternative number. The
  geometric mean follows the HDI precedent and the Xenobiology doc's use
  of a geometric mean for its habitability scores.
- Do not replace the sum with the minimum: on the generated set that loses
  0.38 TL on average and lets whichever domain was randomly lowest decide
  the ranking.

### Tech bands

Named bands, centered on an integer so "TL 3 to 4" is a band (names are
admin-editable constants):

| TL | Name | Typical content | Earth analog |
|---|---|---|---|
| under 0.5 | Primal | Fire, stone, hides | to about 10,000 BCE |
| 0.5 to 1.5 | Agrarian | Draft animals, bronze and iron, writing, sailing | to about 1760 |
| 1.5 to 2.5 | Industrial | Steam, steel, telegraph, antibiotics | 1760 to about 1940 |
| 2.5 to 3.5 | Atomic | Fission, digital networks, jets, orbital rockets | about 1940 to now (Earth 3.1) |
| 3.5 to 4.5 | Fusion | Fusion, true AI, nanotube materials, in-system ships | TA's example |
| 4.5 to 5.5 | Gravitic | Gravity induction, immersion networks, warp or gates | |
| 5.5 to 6.5 | Antimatter | Annihilation power, self-directing AI, tactical warp | |
| 6.5 to 7 | Transcendent | Zero-point power, programmable matter | |

The age era says how long a civilization has existed and the tech band
says what it can do, so the species page shows both, with the tech band as
the primary label for "what are they like".

### Fit with the code (POP.9)

- A pure function beside `species_traits` and `civilization_age` in
  `planetgen/population/model.py`: `draw_indices(age_years, rng, traits)`
  returning the six integers. It uses its own seeded stream
  (`draw.Stream(f"tech:{planet_id}")`). It must not draw from the
  `draw.Stream(planet_id)` that `scan_life_worlds` uses for traits, the
  civilization draw and the age, or every stored species' age changes on
  the next pass. `draw.Stream` already has `random`, `uniform` and `gauss`;
  Poisson counts (Knuth's method) and t(4) variates (a gauss divided by
  `sqrt(-0.5 * ln(u1 * u2))`) are three lines each, so no numpy or scipy
  is needed. The simulation used numpy for those two draws only.
- `refresh_civilizations` already recomputes era from stored age on every
  pass; it also fills the indices of any species whose indices are NULL,
  which delivers POP.9's "existing species get theirs on the next
  population pass".
- When an admin edits an index by hand, `tech_source` is `'admin'` and no
  pass redraws it.
- GEN.92 and GEN.89 change which worlds reach `technological_civilization`
  and the window, not the tech curve, as long as `civilization_age_years`
  keeps its meaning. The draw needs only a 0 to 1 harshness number from
  the habitability score.
- Storage is under "Storage", below.

### Notes for a later head-count

`species` has no population yet, so nothing conflicts today. For when one
is added [S unless marked]: Earth had about 8 billion people on 15
November 2022; growth has fallen almost continuously since the 1960s, so
use a logistic to a capacity with a falling rate, not an exponential.
Published Earth carrying-capacity estimates run from about 2 billion to 40
billion, with outliers of hundreds of billions from photosynthesis limits,
and there is no consensus, so derive capacity from habitable area times a
density parameter by tier, not from energy. Energy per capita does not cap
population before Energy 5 or 6: 8e9 people at Energy 3 draw about 2e13 W,
0.01% of the 1.7e17 W of sunlight on an Earth-size planet, and about 3% at
Energy 6 [C]. Leave it out.

## Facility types and ratings (POP.10)

Not built. The research design for "facility types that are programmable
and can be selected from a dropdown", with affiliation and
Green/Yellow/Red ratings for crime, housing, resources, maintenance and
health (Boss, 2026-10-03).

### What exists

`facilities` (schema v42, `population/facilities.py`,
`tuning.FACILITY_KINDS`, `tuning.FACILITY_RULES`) has a fixed `kind`
(colony, outpost, mining-colony, station, starbase) with a CHECK, a
placement (terrestrial, orbital, asteroid, standalone) and host columns.
`FACILITY_RULES` decides which kinds go on which host. The placement form
(`web/system_facilities.py`, `static/facilityform.js`) has no type list,
affiliation or ratings.

### Data model

| Option | Fits "values in columns" | Search | Admin-defined fields | Cost |
|---|---|---|---|---|
| Fixed enum (today) | yes | yes | no | none |
| Type table + JSON column of values | no | generated columns, different syntax on MySQL 8.4 and MariaDB 11.4 [R] | yes | low |
| Type table + JSON Schema driven forms | no | as above | yes | a validator and a form generator (@jsfe/shoelace 0.4.0 was last released about two years ago [S], so do not depend on it) |
| **Type table + typed field definitions + typed value rows** | **yes** | **yes, one join** | **yes** | moderate, a few hundred lines |
| One table per type | yes | yes | no, DDL per type | not for admin-defined types |

Chosen: the fourth row. A **type refines a kind** ("Agricultural colony"
is a colony), so `kind` keeps driving `FACILITY_RULES` and the placement
checks, and those are untouched.

    facility_types        id, name (unique), kind, icon, description, sort_order, active
    facility_type_fields  id, type_id FK CASCADE, key (regex ^[a-z][a-z0-9_]{0,31}$), label,
                          datatype CHECK IN ('int','decimal','text','bool','choice'), unit,
                          required, min_value, max_value, choices (newline-separated text),
                          default_text, sort_order, UNIQUE(type_id, key)
    facility_field_values facility_id FK CASCADE, field_id FK RESTRICT (retire fields, do not delete),
                          value_int BIGINT, value_dec DOUBLE, value_text VARCHAR(255), value_bool TINYINT,
                          PRIMARY KEY (facility_id, field_id)
    facilities (+)        type_id FK RESTRICT NULL, affiliation_polity_id FK SET NULL NULL,
                          affiliation_name VARCHAR(255) NULL,
                          rating_crime, rating_housing, rating_resources, rating_maintenance,
                          rating_health  TINYINT UNSIGNED NULL CHECK (NULL or 1..3),
                          ratings_source VARCHAR(16), ratings_note TEXT, ratings_set_at TIMESTAMP

Values are validated with a Pydantic model built from the field
definitions by `pydantic.create_model` (Pydantic 2.13.5 is in
`requirements.lock`). Field definitions export to JSON Schema mechanically
if that is wanted later; `jsonschema` 4.25.1 is the last release found
that lists Python 3.9 (4.26.0 and Pydantic 2.14.0 need 3.10 or newer [C,
PyPI JSON API]), so check the lock's Python markers before adding it. The
rating columns are not custom fields: they sit on every facility, are
searched ("show red health everywhere") and need indexes.

### How programmable

| Level | Meaning | Decision |
|---|---|---|
| 0 | Fixed list in code (today) | replaced |
| 1 | Admin-defined types: name, icon, parent kind, defaults | do |
| 2 | Typed custom fields, validated server side | do |
| 3 | Declarative rules (which fields show when, default rating thresholds, from a closed list of operations) | maybe later |
| 4 | Expression engine for computed fields | only if asked; `simpleeval` 1.0.8, MIT, Python 3.9 or newer [C, PyPI] |
| 5 | Admin scripts | no: arbitrary code from a web form is a code-execution hole and nothing here needs it |

Safety for levels 1 and 2: admin scope only (the API already has scopes),
keys restricted by the regex, labels rendered through Jinja autoescape, at
most 24 fields per type, bounded choice lists, icon names checked against
a vendored allowlist.

### Icons

Shoelace's default icon library is Bootstrap Icons, MIT licensed
(`@shoelace-style/shoelace` 2.20.1 ships the licence; `bootstrap-icons`
1.13.1 has 2,078 SVGs, 1.2 MB [C]). The project vendors only `gear.svg`
plus five system icons under `static/vendor/shoelace/assets/`, registers
resolvers in `static/shoelaceicons.js`, and its CSP is `default-src
'self'`, so icons must be local. Vendor a curated subset (60 to 100 icons,
about 60 KB) with the licence file beside it. The admin picks from the
allowlist in a dropdown with a preview; the stored value is the icon's
name. Names confirmed present in 1.13.1 [C]: `check-circle-fill`,
`exclamation-triangle-fill`, `x-octagon-fill`, `dash-circle-fill`,
`shield-check`, `shield-exclamation`, `people-fill`, `house-door`,
`box-seam`, `wrench-adjustable`, `heart-pulse`, `building`,
`rocket-takeoff`, `cash-stack`, `fuel-pump`, `moon-stars`. Not checked
against the vendored Shoelace version beyond its licence file.

### Ratings

Boss's wording is Green/Yellow/Red; project management usually says
Red/Amber/Green [S]. Keep Boss's wording in the UI and database, and reuse
the chip component of the habitability tiers (Blue, Green, Yellow, Red;
[habitability-index.md](habitability-index.md)). Green is the good end for
all five dimensions.

| Dimension | Green | Yellow | Red | Icon |
|---|---|---|---|---|
| Crime | Low | Elevated | High | `shield-check` / `shield-exclamation` |
| Housing | Ample | Tight | Overcrowded | `house-door` |
| Resources | Sufficient | Strained | Critical | `box-seam` |
| Maintenance | Sound | Deferred | Failing | `wrench-adjustable` |
| Health | Healthy | At risk | Crisis | `heart-pulse` |

Storage: five nullable TINYINT columns, 1 = Green, 2 = Yellow, 3 = Red,
NULL = not assessed; the API speaks "green", "yellow", "red". An ordinal
is what Boss asked for, sorts and aggregates (`MAX(rating)` is
worst-of), and a numeric score would invite false precision. An "overall"
is the worst of the five, derived and never stored.

**Accessibility.** WCAG 1.4.1 (Level A) forbids colour as the only means
of conveying information [S]. A chip is one of three distinct icon shapes
(circle-check, triangle-exclamation, octagon-x), plus the word ("Yellow"),
plus an `aria-label` with the dimension ("Crime: Yellow, elevated").
Palette (Okabe-Ito [S]): bluish green `#009E73`, orange `#E69F00`
(labelled "Yellow"), vermillion `#D55E00`. Measured [C]:

| | Green | Yellow | Red |
|---|---|---|---|
| Contrast on white / on #121212 | 3.4 / 5.5 | 2.2 / 8.3 | 3.9 / 4.8 |
| Black text on the fill | 6.1 | 9.3 | 5.4 |

Black text on the fills passes 4.5:1 for all three and the fill does not
depend on the theme, so the colour goes in the chip's fill, with a 1 px
dark ring for 3:1 non-text contrast on pale surfaces. Under simulated
colour blindness the CIELAB distance for classic `#2e7d32`/`#c62828` is 19
(protan) and 16 (deutan), against 37 and 53 for `#009E73`/`#D55E00` [C,
recalled simulation matrices, indicative]. Yellow against red under
deutan is only 18, the weakest pair, so the icon shape and word carry the
distinction. Check the final design in greyscale.

### Setting ratings

Admin-set by default, because nothing in the data model gives a true value
(no population count, crime figure or housing stock). A "suggest" action
fills unset ratings and marks `ratings_source = 'generated'`, so admin
edits are never overwritten. A first rule set (design suggestion, not
run): each score 0 to 1, Green at or above 0.66, Yellow at or above 0.33,
Red below, with seeded noise of +/- 0.1 from the facility id so reruns
are stable.

| Rating | Inputs | Score |
|---|---|---|
| Health | species Medical index, host habitability (breathable, radiation, gravity against build), orbital or surface | 0.5 Medical/7 + 0.5 habitability |
| Housing | host habitability, Materials index, kind | 0.5 habitability + 0.3 Materials/7 + 0.2 kind size |
| Resources | host richness (belt against bare orbit), Energy index, distance to the owner's capital | 0.4 host richness + 0.3 Energy/7 + 0.3 (1 - distance/reach) |
| Maintenance | Materials and Energy index, distance to capital, harshness | similar |
| Crime | Information index, polity government form, distance from capital, unaffiliated penalty | |

Inputs that exist today: the host (type, body type, gravity, temperature),
distance (`system_owners.distance_ly`), and after POP.9 the six indices;
habitability arrives with GEN.89. Reds on harsh hosts are legitimate and
give the habitability tier an immediate use.

### Affiliation

A nullable FK to `polities` with `ON DELETE SET NULL`, plus
`affiliation_name`, a text snapshot written whenever the FK is set.
`planetgen population --rescan` deletes every species and with them every
polity and ids are reassigned; polities exist only after the optional pass
has run; and a facility may belong to nobody ("independent"). A bare FK
would lose every affiliation on a rescan. The system's owner
(`GET /api/systems/<id>/owner`) is the default when placing a facility,
editable by the admin. A MySQL CHECK cannot use a column that is in an
`ON DELETE SET NULL` FK [R], so this is validated in code like
`_db.add_facility` does for hosts.

### API and UI sketch

`/api/facility-types` (admin CRUD, public read). `POST /api/facilities`
accepts `type_id`, `affiliation_polity_id` or `affiliation_name`, the five
ratings as strings, and `fields` as an object keyed by field key.
`GET /api/facilities` filters by type, affiliation and any rating. The
system page's facility form gets a type `sl-select` (UX.49 is moving form
fields to Shoelace) and builds the custom-field inputs from the type
definition; the Facilities panel adds the five chips as a compact strip
and sorts worst rating first.

## Storage (schema v44)

- `species`: id, name (unique), homeworld_planet_id (unique),
  star_system_id, life_chemical, life_stage (`multicellularity` or
  `technological_civilization`), build, climate, size,
  civilization_age_years, era, spacefaring (every stored species has a
  civilization since GEN.80).
- `polities`: id, name (unique), species_id (unique), capital_system_id,
  government, color (`#rrggbb`), reach_ly.
- `system_owners`: star_system_id (primary key), polity_id, distance_ly.
- `population_state`: one row (`id = 1`), the planet-id watermark.

`_migrate_v43_to_v44` creates the four tables empty. Deleting a planet or
system cascades to its species, the species to its polity, and the polity
to its ownership rows. Since GEN.80 `life_stage` is always
`technological_civilization`.

### Planned columns (POP.9, POP.10; one new schema version each)

- `species` (POP.9): `energy_index`, `materials_index`,
  `information_index`, `medical_index`, `propulsion_index`,
  `defense_index` (TINYINT UNSIGNED, CHECK 0..7, NULL until drawn),
  `tech_level` DECIMAL(3,2) computed from the integers and indexed for
  range search, and `tech_source` (`'generated'` or `'admin'`, so a pass
  never redraws a hand-edited species). The CHECKs are fine on MySQL
  because these are plain columns, not columns of a cascading FK [R].
  Rounding left 105 distinct `tech_level` values among 10,000 simulated
  species [C]. Limiting domain and the lopsided flag are derived, not
  stored.
- `facility_types`, `facility_type_fields`, `facility_field_values` and
  new columns on `facilities` (type, affiliation, five ratings, rating
  source): the layout under "Data model" in the facility section.

## Presentation

- **API** (`src/planetgen/api/population.py`, its own blueprint):
  `GET /api/population` (whether a pass has run and whether any species,
  polity or owned system exists, so pages can hide themselves),
  `GET /api/species` (paged, `?spacefaring=`), `GET /api/species/<id>`,
  `GET /api/planets/<id>/species`, `GET /api/polities` (paged),
  `GET /api/polities/<id>` (with a page of its systems),
  `GET /api/systems/<id>/owner` (404 for an unknown system), and
  `GET /api/territories` (at most 20,000 owned systems, nearest their
  capitals first, with each polity's capital and reach). Documented in
  `docs/api.md`.
- **Galaxy Map:** a Territories button draws each polity's reach as a soft
  ball of its color around its capital and its systems as dots, with a
  legend naming each polity, its government and its system count (page
  endpoint `/galaxy/territories`, which merges `/api/territories` with
  `/api/polities`). It draws the same at every zoom. The button and
  legend appear only when `/api/polities` counts at least one polity.
- **Web pages (not built yet):** a Species list (paged, spacefaring
  filter), a species page, a polity page listing its systems, "Dominant
  species" on a life world's planet card and "Territory of ..." on an
  owned system's page. With POP.9 the species page shows the six indices,
  the tech level with its band name, the limiting domain and the age era
  side by side, and the species list can be searched by `tech_level` range
  and by each index.

## Decisions

Boss accepted every recommended default on 2026-10-01:

1. Named dominant species on multicellular worlds and up. Superseded by
   GEN.80 (PR #448, a bug report from Boss's 2026-10-07 list): only worlds
   with a technological civilization get a species.
2. Species names unique across the galaxy.
3. One polity per spacefaring species.
4. Territory reach cap 100 ly.
5. Territories recompute automatically after each fill. Superseded the
   same day: population is optional and off by default everywhere,
   including fills and install/update (see "When it runs").
6. Civilizations on 1 in 1,000 worlds past the milestone.

Decisions already taken for POP.7 to POP.10:

7. Boss (2026-10-03 05:38Z): generate a tech level for every technological
   species using "Technological Assessment.md" (POP.7). The six domains,
   the 0 to 7 scales and the weights 0.25 / 0.20 / 0.20 / 0.15 / 0.10 /
   0.10 are Boss's own document and stay as written.
8. Boss (2026-10-03 05:38Z): facility types that are programmable and
   chosen from a dropdown when adding a facility; facility data includes
   an affiliation and a place for Green/Yellow/Red ratings of crime,
   housing, resources, maintenance and health (POP.10).
9. Boss (GEN.80, 2026-10-07): species, polities and population are made
   only for worlds whose life reaches a technological civilization.

Everything else in the tech-level and facility sections (age definition,
logistic curve and its parameters, limiting domain, band names, the
declarative type model, palette, the suggest-ratings rules) is a
recommendation until Boss decides.

## Why it works this way

- **A pass over stored rows, not part of generation.** Everything it
  needs (the milestone and ages in the evolution text, each system's
  position) is already in the database, so an existing galaxy can be
  populated without a regenerate and system generation stays untouched.
  Territories depend on neighboring systems, which a single `StarSystem`
  cannot see while it is being built.
- **A watermark plus derived tables.** Only new planets need scanning;
  eras, polities and ownership are derived from stored ages and
  positions, so they are cheap to rebuild from scratch on every run.
  Retuning `CIVILIZATION_ERAS` or the reach constants needs no
  migration.
- **Seeded by planet id.** Traits, the civilization draw and the age come
  out the same on every run and on every server, so a rerun does not
  reshuffle the galaxy. Names are the exception (they use the shared
  generator), and a polity's government and color follow its species id,
  so `--rescan` changes those too.
- **1 in 1,000.** The milestone alone made nearly every seventh system a
  polity on a 150-sector test fill (commit b62d0aa). The chance was
  picked to keep civilizations scarce; the spacing it gives is 134 ly on
  average, so cap-sized territories overlap (see "Civilization age").
- **Log-uniform ages, square-root reach.** Ages spread across orders of
  magnitude, so all six eras occur; reach grows like a diffusing frontier
  and stops at a cap so the oldest polities don't swallow the galaxy.
- **`reach / distance`.** A weighted Voronoi split gives every system one
  owner and lets an older polity hold more ground than a young neighbor,
  with no overlap to store.
- **A logistic in log age, integers 0 to 7, a stream of its own.** Tier
  durations fall by orders of magnitude, so a logistic in linear time has
  the wrong shape. TA defines each index as a band of named technologies,
  so indices are small integers. A separate seeded stream means adding the
  tech draw cannot change any stored age.
- **Types refine kinds.** Keeping `kind` leaves `FACILITY_RULES` and the
  placement checks alone, and typed value rows keep the project's rule of
  values in columns.
- **Optional and off by default.** Boss's choice on 2026-10-01; no
  further reason is written down. It keeps scheduled `update.sh` runs and
  plain fills from spending time on the pass, and the pages and overlay
  hide themselves when it has not run.

## Open points found in this review

- `planetgen population --rescan`'s help says it gives "new names,
  ages and borders". Names are new, and so are governments and colors
  (they follow the new species ids), but ages, traits and the
  civilization draw are seeded by the planet id and come out the same,
  so borders stay the same unless the constants were retuned.
- The module docstring of `planetgen/population/model.py` still says the
  pass runs "after every `galaxy`/`sector` run"; since 7.58.2 it does
  only with `--population`.
- The docstring of `tuning.CIVILIZATION_CHANCE` repeats the "150-200 ly /
  just about meet" wording corrected under "Civilization age".
- "Technological Assessment.md" has pasted LaTeX debris in its title and
  first paragraph and "an species"; it is Boss's upload and was left alone.

## Evidence notes

The research environment could read only search-result text, not the
papers or pages themselves. Items to verify when paper access is allowed:

- [R] Cook's per-capita energy figures and Morris's index components came
  from secondary summaries (P2P wiki, an arXiv table, a Grokipedia
  summary); the factor of 3.9 per Energy index is a fit by eye. Maddison
  levels other than the 1 CE and growth-multiple figures are not quoted for
  that reason.
- [R] Earth's domain vector for 2026 and for 1760 to 1920 is a judgment of
  TA's band descriptions; the dates (fission 1942 to 1954, steam 1712 to
  1769) and the tier durations (about 5,800 and 180 years) are from memory.
  Boss should review the Earth scoring, because it anchors the whole curve.
- [R] Traveller and GURPS levels (TL0 stone age, TL10 to 11 interstellar)
  were not verified; no mapping to TA is offered. Star Trek and
  Civilization tech trees were not reviewed.
- [R] The colour-blindness distances used Machado et al. (2009) matrices
  typed from memory, applied in linear RGB; they are indicative. Contrast
  ratios are exact WCAG computations [C].
- [R] That a MySQL CHECK cannot use a column in an `ON DELETE SET NULL` or
  `CASCADE` FK is from memory of the manual (it matches the note in
  `docs/database-schema.md`); verify against the MySQL 8.4 and MariaDB 11.4
  manuals, as well as the JSON and generated-column syntax differences.
- [R] The Newman and Sagan diffusion models were not checked against the
  5 ly reach at 2,000 years rising with the square root of age (about
  0.0025 c of frontier speed at the start).
- Simulation: exact percentages move by about a point per numpy seed
  (seeds 1 to 3 tested). The real distribution of `window_years` in a
  generated galaxy was approximated (uniform windows, with an optional
  fast-pace mixture), not read from a database. `ceil_mean` was not
  simulated at other values. The suggest-ratings rules were not run.
- [S] Every URL below was seen in search results only; statements
  attributed to them are the search tool's summaries. The Sagan formula
  and K = 0.73 were confirmed by computation. Packages marked "read
  directly" were fetched from their registries.

## Sources

- Sagan interpolation, Earth K 0.73: scienceabc.com/nature/universe/what-is-kardashev-scale, kardashev1.com/scale, en.wikipedia.org/wiki/Kardashev_scale
- Morris energy capture: pzacad.pitzer.edu/~lyamane/ianmorris.pdf, wiki.p2pfoundation.net/Energy-Capture_Per-Capita_Index, eh.net review of "The Measure of Civilization"
- Cook energy stages: wiki.p2pfoundation.net/Evolution_of_Energy_Capture_Throughout_Human_History, arxiv.org/pdf/1406.6328
- HDI geometric mean: data.un.org/_Docs/FAQs_2011_HDI.pdf, roiw.org/2019/n4/roiw12370.pdf
- Maddison: oecd.org "The World Economy" (2001), h-net.org/reviews/showrev.php?id=30336
- Barrow scale: centauri-dreams.org/2023/12/12/seti-musings-on-the-barrow-scale
- TRL: esto.nasa.gov/files/trl_definitions.pdf, nasa.gov technology-readiness-levels
- Traveller and GURPS: wiki.travellerrpg.com/Technology_Level/meta, gurps.fandom.com/wiki/Tech_Level
- Era durations: en.wikipedia.org (Stone Age, Bronze Age), britannica.com/event/Industrial-Revolution
- Population: en.wikipedia.org/wiki/World_population_milestones; Energy Institute Statistical Review 2023 (620 EJ)
- Colour: w3.org/WAI/WCAG22/Understanding/use-of-color.html, developer.mozilla.org Use_of_color, sci-draw.com Okabe-Ito reference, rebelsguidetopm.com RAG
- Forms and JSON: github.com/json-schema-form-element/jsfe, npmjs.com/package/@jsfe/shoelace, pypi.org/project/jsonschema
- Read directly: pypi.org/pypi/jsonschema/json, pypi.org/pypi/pydantic/json, pypi.org/pypi/simpleeval/json, registry.npmjs.org/bootstrap-icons (1.13.1), registry.npmjs.org/@shoelace-style/shoelace (2.20.1), raw.githubusercontent.com/twbs/icons/main/LICENSE
- Repo: `Technological Assessment.md`, this document, [habitability-index.md](habitability-index.md), `docs/database-schema.md`, `docs/html-interface.md`, `src/planetgen/population/model.py` and `facilities.py`, `src/planetgen/tuning.py`, `src/planetgen/util/draw.py`.
- Simulation scripts (`fit.py`, `techsim.py`, `analysis2.py`, `priors.py`, `rag.py`, `wsens.py`) were run in the research session's scratchpad and are not in the repo; the results are the tables above.
