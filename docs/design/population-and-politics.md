# Population and politics: species, civilizations and territories

Design for POP.1 to POP.4 (`docs/TODO.md`, "Population and Politics").

**Status (2026-10-01, checked against 7.58.2):**

| Part | Built in |
|---|---|
| The pass, the four tables (schema v44) and the read API | 7.49.0, PR #169 |
| 1 in 1,000 civilization chance | 7.49.0, PR #169 (commit b62d0aa) |
| Galaxy Map Territories overlay | 7.54.0, PR #179 |
| Overlay hidden until a polity exists | 7.58.1, PR #181 |
| Optional and off by default; `GET /api/population` | 7.58.2, PR #180 |
| Web pages (Species list and page, polity page, "Dominant species", "Territory of ...") | not built yet |

Boss's decisions are recorded at the end, and the reasons for the main
choices in the last section.

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
`technological_civilization` (`program_constants.LIFE_WORLD_STAGES`):
complex life, so there is a dominant species worth naming. Simpler
biospheres (microbes, mats, single cells) keep their life chemistry but
get no species row. The pass reads the milestone and both ages (the
system's age and the milestone's) back out of the stored paragraph text
(`population.parse_timeline`); a planet whose text can't be read is
skipped.

Each life world gets exactly one row in `species`: a generated name, the
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
1,639 systems). Taken literally that would fill the galaxy with empires,
so a world past the milestone has a civilization now only with
`CIVILIZATION_CHANCE` (default 1 in 1,000), or always when its system was
generated with intelligent life forced on (`system_configs.
intelligent_life`). That gives about one civilization per 7,000 systems;
at the disk's local density (about 0.1 systems per cubic parsec)
spacefaring neighbors sit roughly 150 to 200 ly apart, so territories at
the 100 ly cap just about meet. Every other world past the milestone
still gets its named dominant species, with no civilization.

A civilization is younger than its window. Its age is a log-uniform draw
between 100 years (`CIVILIZATION_MIN_AGE_YEARS`) and the window (`system
age - milestone age`, both read back from the stored paragraph), so ages
spread evenly across orders of magnitude and none is older than its star
allows. A window shorter than 100 years (the stored ages are rounded to
10 million years, so this is rounding) gives 100 years.

Age maps to an **era**, which is what the pages show and what drives
everything else:

| Era            | Age (years)          | Spacefaring | Notes |
|----------------|----------------------|-------------|-------|
| Industrial     | < 300                | no          | one world, no spaceflight beyond orbit |
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

`generate.py population` is a separate pass over what is already stored
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
runs it unless asked. `generate.py galaxy` and `generate.py sector` run
the whole pass after they save only with `--population` (this replaced
7.49.0's `--no-population`); the admin Generate page's jobs don't pass
it. `install.sh`/`update.sh` (and `install.ps1`/`update.ps1`) offer to
run it after the database step (`offer_population_pass`,
`Invoke-OptionalPopulation`), y/N with a 30-second timeout defaulting to
No and skipped with no terminal; `POPULATION=1` (`-Population` on
Windows) runs it without asking. Nothing in system generation itself
changes.

### Planned: tech levels and facility types (2026-10-03, 2026-10-07)

- Species, polities and population are made only for worlds whose life
  reaches a technological civilization.
- **Tech levels** ("Technological Assessment.md"): six domain indices,
  each 0 to 7 (Energy, Materials, Information, Medical, Propulsion,
  Defense), combined as TL = 0.25 EI + 0.20 MI + 0.20 II + 0.15 MeI +
  0.10 PI + 0.10 DI. Each index is drawn from the species' civilization
  age, era and world, and stored in columns.
- **Facility types**: admins define facility types and pick one from a
  dropdown when placing a facility; each facility stores an affiliation
  and Green/Yellow/Red ratings for crime, housing, resources, maintenance
  and health.

## Storage (schema v44)

- `species`: id, name (unique), homeworld_planet_id (unique),
  star_system_id, life_chemical, life_stage (`multicellularity` or
  `technological_civilization`), build, climate, size,
  civilization_age_years (NULL without a civilization), era (NULL
  likewise), spacefaring.
- `polities`: id, name (unique), species_id (unique), capital_system_id,
  government, color (`#rrggbb`), reach_ly.
- `system_owners`: star_system_id (primary key), polity_id, distance_ly.
- `population_state`: one row (`id = 1`), the planet-id watermark.

`_migrate_v43_to_v44` creates the four tables empty. Deleting a planet or
system cascades to its species, the species to its polity, and the polity
to its ownership rows.

## Presentation

- **API** (`src/html/api/population.py`, its own blueprint):
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
  owned system's page.

## Decisions

Boss accepted every recommended default on 2026-10-01:

1. Named dominant species on multicellular worlds and up.
2. Species names unique across the galaxy.
3. One polity per spacefaring species.
4. Territory reach cap 100 ly.
5. Territories recompute automatically after each fill. Superseded the
   same day: population is optional and off by default everywhere,
   including fills and install/update (see "When it runs").
6. Civilizations on 1 in 1,000 worlds past the milestone.

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
  picked so that neighbors at the 100 ly cap just about meet.
- **Log-uniform ages, square-root reach.** Ages spread across orders of
  magnitude, so all six eras occur; reach grows like a diffusing frontier
  and stops at a cap so the oldest polities don't swallow the galaxy.
- **`reach / distance`.** A weighted Voronoi split gives every system one
  owner and lets an older polity hold more ground than a young neighbor,
  with no overlap to store.
- **Optional and off by default.** Boss's choice on 2026-10-01; no
  further reason is written down. It keeps scheduled `update.sh` runs and
  plain fills from spending time on the pass, and the pages and overlay
  hide themselves when it has not run.

## Open points found in this review

- `generate.py population --rescan`'s help says it gives "new names,
  ages and borders". Names are new, and so are governments and colors
  (they follow the new species ids), but ages, traits and the
  civilization draw are seeded by the planet id and come out the same,
  so borders stay the same unless the constants were retuned.
- The module docstring of `stellarObjects/population.py` still says the
  pass runs "after every `galaxy`/`sector` run"; since 7.58.2 it does
  only with `--population`.
- `docs/TODO.md` still lists the territory overlay as open; it shipped in
  7.54.0.
