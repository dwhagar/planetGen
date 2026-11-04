# What one master "mass object" class would save, with a 10 m floor

**Status:** plan (decided). Code analysis written after 11:39 PM EDT on 2026-10-17 from `main` at bc97f75; it adds to [unified-star-like-objects.md](unified-star-like-objects.md), which designs the behaviour (the recursive host chain), and this note measures the data and class side. Decided by Boss on cards E1 to E4 (11:46 to 11:49 PM EDT, 2026-10-17, as the coordinator relayed): **E1** the hybrid of section 10 (three family classes StellarObject, PlanetaryObject and SmallBody under MassObject, state and parent as data, one core table plus three profile tables; Boss typed his approval at 11:48 PM EDT); **E2** the 10 m floor is a representability floor (any body of 10 m or more can be a row, generated rows start at each population's own cut, belt members are drawn on demand from the belt's seed); **E3** built inside the host chain (with items U4, U7, U8 and U11 of the unified note), not as a separate refactor; **E4** the core is planned to take ships and stations in 9.3 and nothing is built for them now. Object IDs stay stable (card, 11:49 PM EDT): see DB.75. The items are GEN.565 (U4a), DB.74 (U11a), TEST.240 (U13a), GEN.566 (U15a) and DB.75 (stable IDs); section 7 below names the first four.

**Boss, 11:39 PM EDT:** "Code analysis, how much potential code savings could we get by making a master 'objects with mass' class, then make everything down to the lowest planet limit (class D I think is 500meters, might be as low as 50m) but the smallest object becomes the cutoff becomes say 10 meter diameter (5 meter radius), and that's the smallest object we'll track in this way."

Tags: [M] measured here by counting the repository; [E] my estimate (judgment, shown with its arithmetic); [R] recalled, not verified.


**Standing design rule (Boss, 5:20 PM EDT on 2026-10-18, typed in the project chat, relayed by the coordinator):** "For all objects are objects with mass in the system. I want to make sure that galactic motion is not the only thing that they handle because the trajectory of all motion, including planetary and lunar, should both be built into the object mass functions (since they're all the same rules)". So the mass object core and its family classes carry the trajectory of every kind of motion with one set of rules: galactic motion, a planet's and a moon's orbit about its parent, a comet's Kepler state and a star's orbit about a barycentre are all the same state vector relative to the parent (see the kinematic record of GEN.584), advanced by the same functions. GEN.565a (#2274) moved only the galactic-motion fields and the save and load into the core; GEN.565b and GEN.584 must bring the planetary and lunar trajectory into the same functions, not leave it in the family classes or the orbital update.

## 1. Short answer

- **Net lines saved: about 800 to 1,800, which is 10 to 20% of the roughly 9,000 lines the change touches** [E]. Most of the textual duplication is not in the classes themselves (under 10% of their code lines repeat across files [M]); it is in the serializers, the insert functions, the nine near-identical table definitions and the per-kind lists and branches scattered through 24 files.
- **The bigger payoff is not lines.** Eight of the 49 migrations (1,454 lines, 16%) changed 8 to 13 kind tables each [M]. With one shared core, a cross-cutting change (a spin column, the object ID, the point columns) is one change and not nine. And the master class is what removes the divergent copies of the same rule (17 copies of five rules, section 3) that make a rogue planet age differently from a star.
- **Alone, it is not worth the rewrite.** The rewrite touches every file that names a kind table (24 files for `rogue_planets` alone, 36 test files mention the rogue classes [M]). Done as part of the host chain of the earlier note (items U4, U7, U8, U11), most of the cost is already paid by work that is wanted anyway. **Recommendation: build it as a core plus profiles inside those items, not as a separate refactor.**
- **The 10 m floor is right as a floor of what the class can represent, and wrong as a floor of what is generated as rows** (section 4): the generator has no individual body below 50 m diameter in the code, and a belt holds an estimated 1e10 to 1e12 bodies above 10 m.

## 2. What exists, measured

### The classes (in memory)

| Family | Classes | Lines (file) | Notes |
|---|---|---|---|
| Star-like | `Star` 1,524; `BinaryStarProxy` 478; `WideBinaryPair` 371; `CompactRemnant`, `BlackHole`, `NeutronStar` 677 | 3,050 | `CompactRemnant` already subclasses `Star` and calls it with `_skip_property_init=True`, then fills the Star attribute surface by hand. A partial master class exists already. |
| Free planet-like | `RoguePlanet` (also the brown dwarf), `InterstellarComet` | 585 | Own surface model `rogue_surface` (284). |
| Planet | `Planet` (planets and moons) | 665 | Plus `physics/planets.py` 1,263; 109 serialized fields. |
| Small | `Comet` 419; `AsteroidBelt` 225; `AsteroidField` 251 | 895 | Belts and fields are statistics, not bodies. |
| Total | | 5,195, of which 2,079 are docstrings | About 3,100 lines of code. |

Serialized fields: 313 field slots across the nine classes, 179 distinct names; 53 names occur in two or more classes and 38 in three or more, among them name, mass, radius, spin axis (3), axial tilt, rotation period, the four galactic-orbit fields and the three reflex offsets [M].

Serializer boilerplate: `to_dict`, `from_dict` and `__str__` take 440 lines in the nine families [M].

### The tables

Ten kind tables (`stars`, `planets`, `moons`, `asteroid_belts`, `comets`, `black_holes`, `neutron_stars`, `rogue_planets`, `interstellar_comets`, `asteroid_fields`) take 910 lines of `schema.sql` and have 518 column declarations. Of those, 454 are in columns that two or more tables share and 274 are in 49 columns that three or more share (`id`, `uid`, `velocity_x/y/z_kms`, `epoch_unix`, `name`, `next_update_due`, `mass_kg`, spin axis, tilt, galactic orbit period and phase, rotation period, `radius_km`, `sector_id`, `center_x/y/z_pc`, `quadrant`, `inside_nebula_id`, `inside_remnant_id`, timestamps) [M]. `models.py` mirrors them.

### The code that knows the kinds

- `store.py`: ten insert functions of 558 lines for the mass-object kinds (star 60, planet 37, moon 34, belt 58, comet 64, black hole 85, neutron star 65, rogue planet 41, interstellar comet 49, asteroid field 65), 725 with nebula, remnant and quasar [M].
- The same list of kind tables is written out 22 times in 9 files for `rogue_planets` alone (`store.py` seven times, `query.py`, `edits.py`, `fingerprint.py`, `point_columns.py`, `sector_paths.py`, `api/ids.py`, `cli/orbits.py`, `web/lib/slowpages.py`) [M].
- Tests: 512 test files; 36 mention the rogue classes, and 60, 46, 45 lines of tests mention black holes, neutron stars and rogue planets by table name [M].

## 3. What a core class would remove

Design assumed: **`MassObject` core plus profiles**, not one fat class. The core holds what every body has (identity and ID, name, mass, radius, state vector relative to its parent with `parent_type` and `parent_id`, age, spin, galactic motion, sector and place, timestamps, serialization from a field table). Profiles hold what only some kinds have: a stellar profile (state, luminosity, temperature, habitable zone, activity), a planetary profile (class, composition, atmosphere, 109 fields), a small-body profile (nucleus, activity, composition). A profile is a component of the object, not a subclass, so a rogue planet is a `Planet` profile with no parent and a remnant is a stellar profile in the "remnant" state.

| Layer | In the code [M] | Removed [E] | Basis of the estimate |
|---|---|---|---|
| Serializers (`to_dict`, `from_dict`, `__str__`) | 440 | 300 to 360 | One table-driven pair, about 80 lines, replaces nine pairs. |
| Constructor boilerplate: galactic orbit, spin, activity, Hill sphere, name set-up, description scaffold | about 30 to 50 lines in each of 9 classes | 170 to 350 | The `_finish_init` of `CompactRemnant` reproduces the tail of `Star.__init__` by hand. |
| Copies of the same rule | 17 copies of 5 rules (below) | 100 to 200 | Lines; the behaviour gain is larger. |
| `store.py` inserts | 558 | about 350 | One insert driven by a column list, about 200 lines. |
| Per-kind lists and branches (`query`, `edits`, `fingerprint`, `point_columns`, `sector_paths`, `api/ids`, `cli/orbits`, `slowpages`, loaders) | about 22 sites per kind | 300 to 600 | Replaced by one registry lookup. Low confidence: not each site was read. |
| `schema.sql` plus `models.py` | 274 declarations of shared columns in each | about 450 in both | 274 declarations become about 49 in a shared table, in each file. |
| **Gross** | | **about 1,700 to 2,300** | |
| New code: core, profile registry, generic insert and load, migration | | +500 to +900 | |
| **Net** | | **about 800 to 1,800** | |

The 17 copies of 5 rules (from the earlier read) [M]: the white dwarf initial-final relation in four places (plus two side models for hot and cool white dwarfs); the neutron star and black hole mass draw in five places; three age rules (stars by population, remnants uniform 1e-5 to 10 Gy, rogue planets uniform 0.1 to 10 Gy); the Schwarzschild radius in two; the Burrows brown-dwarf luminosity in two (one with a hard-coded 0.05). Add the comet composition sample in two files and the rogue's galactic orbit called without a distance.

What it does not remove: the physics modules (`planets.py` 1,263, `rogue_surface` 284, `stellar_evolution` 621, `activity` 355), the descriptions, the web pages' per-kind wording, and the tests' per-kind expectations. The earlier note's host chain does not shrink them either; it changes what calls them.

## 4. The 10 m floor

### What the code has at the bottom

| Thing | Size | Source |
|---|---|---|
| Class D (smallest planet class) | a conditions rule (GEN.1091): surface pressure under 10 Pa and at least 1 percent water or ice by mass (`physics/data/planet_classes.json`); size is not in the rule. Its intended scale is a dwarf planet or a moon, radius 50 to 500 km, so 100 to 1,000 km across (not metres); it covers belt members and dwarf planets as well as moons | `planet_classes.json`, `tuning.PLANET_CLASSES["D"]` [M] |
| Star-bound comet nucleus | diameter 0.5 to 20 km | `BOUND_COMET_NUCLEUS_DIAMETER_RANGE_KM` [M] |
| Interstellar comet nucleus | diameter 0.05 to 5 km, so **50 m is the smallest individual body generated** | `INTERSTELLAR_COMET_NUCLEUS_DIAMETER_RANGE_KM` [M] |
| Asteroid belt, asteroid field | one statistical row each: a density and up to four named components | `belt.py`, `asteroid_field.py` [M] |
| Moons | smallest are class D, as above | `moon_class_options` [M] |

So the figure of 50 m in the ask matches the interstellar comet nucleus, and 500 m is close to the star-bound comet's lower end; class D is a small airless icy body whose intended scale is 50 to 500 km.

### What a 10 m (5 m radius) floor means

- A 10 m diameter rock at 2,000 kg/m^3 is 1.05e6 kg [C]. The heaviest object in the generator is a central or quasar black hole of up to 1e10 solar masses, 2.0e40 kg [S: `QUASAR_BLACK_HOLE_MASS_RANGE_SOLAR`]. The class would span about **34 orders of magnitude in mass**, which doubles hold. The class D size is 1e12 times heavier than the floor.
- **Positions must be relative to the parent, never galaxy-absolute.** 10 m at a galactocentric distance of 262 kpc (8e21 m) is a relative precision of 1e-21, far below the 1e-16 of a double. The scatter tables store positions as integer millipc (about 3e13 m a step) [S]. Planets and moons already store state relative to a parent (`parent_type`, `parent_id`) [S]; every small body must do the same, and a wanderer's position is the position of its cell plus an offset.
- **The floor cannot be a floor of generated rows.** An illustration with assumed inputs: if a belt has 1e6 bodies above 1 km and the cumulative size distribution falls as D^-2 to D^-2.7 [R], it has 1e10 to 1e12 bodies above 10 m [C]. Interstellar space is worse: the real density of interstellar debris in the code is 1e12 per pc^3 [S]. So 10 m is a **representability floor**: any body of 10 m or more can be a row (a named asteroid, a captured interstellar object, a mining claim, a ship or a station), and nothing smaller exists as a row.
- **What generates rows is a cut per population**, as with the neutron star and black hole mass cut that Boss approved on 2026-10-09: a planet or moon always, a comet above its size cut, belt members on demand. On-demand belt members fit the seed scheme: the k-th largest member of a belt is a pure function of the belt's seed and k, so it can be drawn when somebody asks and give the same body every time, without storing 1e10 rows.
- **The floor also meets version 9.3** (ships and stations): the `facilities` table holds stations, starbases and outposts on a host but has no mass or size column [M: schema], and a ship or station of 10 to 1,000 m would be a mass object with a parent. If the core is built with the floor, 9.3 joins it without a new object model. That is a reason for the floor, not a requirement.

## 5. Risks

- **A fat class.** The planet profile has 109 serialized fields and the star profile 37. One class that holds all of them is a bag of optional fields; the core plus profiles shape avoids it and keeps each profile's validation local.
- **The orbital update.** Planets and moons are stepped as state vectors (`kepler.universal_step`, GEN.289); stars, remnants and wanderers keep a galactic-orbit phase; comets keep their own Kepler state [S]. A real saving exists if all of them move as one state vector, but it was not measured here and the ordering of a nightly batch is the easiest place to break reproducibility.
- **Cost in effort, not lines.** The rewrite touches every kind table, so it needs a full database wipe and rescatter (allowed on the test server: no compatibility code). It is the same wipe as the host chain's.
- **Net savings are small beside the risk if done alone.** Hence section 6.

## 6. Choices for Boss (one card each)

**E1. Shape of the master class.**
1. One class with every field optional. Largest immediate line saving, worst for validation and the planet profile's 109 fields.
2. Core plus profiles (components). **Recommended:** 800 to 1,800 net lines, no bag of optional fields, and the same shape as the shared parent link of D1 in the earlier note.
3. Keep separate classes and share only through mixins. Smallest saving (the serializers and constructor boilerplate), lowest risk.

**E2. What the 10 m floor means.**
1. A representability floor: any body of 10 m or more can be a row; generated rows start at each population's own cut; belt members drawn on demand from the belt's seed. **Recommended.**
2. Generate every body down to 10 m. Not possible at 1e10 to 1e12 per belt.
3. A floor of 1 km, with smaller bodies as statistics only. Safer, but then ships, stations and small captured objects have no home.

**E3. Do it alone or inside the host chain.**
1. Inside the host chain, as the core and profiles of items U4, U7, U8 and U11 of the earlier note. **Recommended:** most of the rewrite is already paid for by work that is wanted anyway.
2. As a separate refactor first (behaviour-preserving, then the host chain). Cleaner history, but two wipes and a second pass over every kind.
3. Not at all. Keep the nine classes and fix only the 17 divergent copies. Cheapest; the cross-cutting cost stays.

**E4. Ships, stations and facilities.**
1. Plan the core to take them in 9.3, do not build them in this change. **Recommended.**
2. Move facilities into the core in the same change.
3. Leave them outside.

## 7. Proposed additions to the earlier item list

For the TODO thread to number, if E1 to E3 go the recommended way:

- **U4a:** the `MassObject` core (identity, mass, radius, parent-relative state vector, age, spin, galactic motion, field-table serialization) with the stellar, planetary and small-body profiles. Replaces the serializers and constructor boilerplate of the nine classes. Part of U4 and U8; Core 4.
- **U11a:** one generic insert and load driven by the column list, one registry of kind tables replacing the 22-site lists; the shared core table (D1 option 2). Part of U11; Core 1.
- **U13a:** a test that every kind's stored row round-trips through the core, and a test that fails if a kind table is named in a list that the registry should own.
- **U15a (research):** the size-frequency law of small bodies, for the on-demand belt members and the comet size cut (shared with the research of U10).

Order is as in the earlier note; none of this touches 8.1 beyond GEN.315.

## 8. Evidence and limits

All [M] figures were counted from the repository at bc97f75 with line counts, an AST walk (serializer lengths, docstring share), a column scan of `schema.sql` and a grep of table names; none was run against a database. All [E] figures are my arithmetic from those counts and are the weakest part: the per-site branch removal (300 to 600 lines) was not read site by site. The belt population estimate in section 4 is an illustration from an assumed count and an assumed slope, not a sourced result; a research lane must source the size-frequency law before any number from it is used.

## 10. Classes and subclasses, or a core plus profiles: what each looks like (the answer to E1)

Written 11:45 PM EDT on 2026-10-17; Boss chose the hybrid below (E1).

### The two shapes

**A. Classes and subclasses (inheritance).** A kind is a class; a more specific kind is a subclass.

```python
class MassObject:            # mass, radius, position and velocity relative to a parent, age, spin, name, ID
class StellarObject(MassObject):    # luminosity, temperature, habitable zone, activity
class Star(StellarObject): ...
class BrownDwarf(StellarObject): ...
class Remnant(StellarObject): ...
class WhiteDwarf(Remnant): ...
class NeutronStar(Remnant): ...
class BlackHole(Remnant): ...
class PlanetaryObject(MassObject):  # class letter, composition, atmosphere
class Planet(PlanetaryObject): ...
class Moon(PlanetaryObject): ...
class RoguePlanet(PlanetaryObject): ...
class SmallBody(MassObject): ...
class Comet(SmallBody): ...
```

Good: reads naturally, `isinstance` works, each subclass overrides only what differs, and it is close to the code (`CompactRemnant` already subclasses `Star`). Bad: an object that changes kind has to change class. A star becomes a white dwarf, a neutron star or a black hole as it ages; a planet becomes a rogue planet when it is unbound (GEN.291); a brown dwarf and a planet overlap in mass. Python can change an object's class, but a database row cannot, so every change of kind needs a table move or a hack. The tree also multiplies: kind times state.

**B. Core plus profiles (composition).** One object class holds what every body has; the kind-specific data hangs off it as profiles.

```python
class MassObject:
    id, name, mass, radius, parent, state_vector, age, spin, galactic_motion
    stellar:   StellarProfile | None    # state (protostar, main sequence, giant, white dwarf, neutron star, black hole, brown dwarf), luminosity, temperature
    planetary: PlanetaryProfile | None  # class letter, composition, atmosphere
    small:     SmallBodyProfile | None  # nucleus, activity, composition
```

Good: a change of kind is a change of data (a star's `state` goes from `giant` to `white dwarf`; a planet's `parent` goes to none), the database maps one to one (a core table and three profile tables), and one table-driven serializer covers all. Bad: `obj.stellar.luminosity` is one hop longer, and nothing stops an object from holding two profiles unless a check does.

### What I recommend: the two together, split by one rule

Use a small class hierarchy for the **families**, and use **data (state, parent) for anything that can change or is a matter of degree**.

- Classes and subclasses for the three families that never swap into each other: `StellarObject`, `PlanetaryObject`, `SmallBody`, all subclasses of `MassObject`. Behaviour lives here (evolve a star to an age, build a planet's atmosphere, step a comet's activity).
- A **state field** for what an object can become: a `StellarObject` has a state of protostar, pre-main-sequence, main sequence, giant, white dwarf, neutron star, black hole or brown dwarf. These are not subclasses. The same for a pulsar type or a black hole mass class.
- A **parent field** for where an object sits: a `PlanetaryObject` whose parent is a planet is a moon, one whose parent is a star is a planet, one with no parent is a rogue planet. Moon and rogue planet are not subclasses.
- Groups (the nodes of a multi-star tree) are not bodies, so they stay an `OrbitNode`.

The rule: subclass only when behaviour differs and the object never changes kind; use state or parent for the rest.

Against the code: `Star` becomes `StellarObject` (the 37 fields it already serializes), `CompactRemnant`, `BlackHole`, `NeutronStar` fold into it as states (the `_skip_property_init` workaround goes), `Planet` becomes `PlanetaryObject` (it already uses `is_moon`, which is the parent rule in disguise), `RoguePlanet` and the brown dwarf row fold into the two, and `Comet`, `InterstellarComet` and the member of a belt become `SmallBody`. In the database it is one core table plus three profile tables, not one table per kind (choice D1, option 2, of the earlier note).

### Cost and savings

The saving is the same as in the savings note, about 800 to 1,800 net lines, because the serializer, insert and column duplication is removed either way. The hybrid costs a little more than pure profiles (three family classes with real behaviour, so about 100 to 200 lines that pure profiles would not have) and a little less than pure inheritance in migration work (no table moves when an object changes kind). The risk is the same as before: it touches every file that names a kind table, so it belongs inside the host chain items and not as a separate refactor.

### Choice for the card

**E1 (replaces the earlier wording).**
1. Pure inheritance: a subclass per kind. Closest to the code, but a change of kind is a change of class.
2. Pure core plus profiles: one object class, no subclasses. Simplest storage, behaviour spread over profiles.
3. Hybrid: three family classes, state and parent as data. **Recommended.**


## 11. Motion lives in the core (GEN.565b)

Boss, 5:20 PM EDT on 2026-10-18: every object with mass follows the same rules, so the trajectory of all motion, planetary and lunar as well as galactic, is built into the one set of mass object functions. `MassObject` (`generation/mass_object.py`) now holds it:

- **State relative to the parent.** `spatial` (a `SpatialPosition3D`, "system" frame = the offset from the primary), `set_position_au`, `set_velocity_kms`, `orbit_mu_au3_per_year2` (each family says what it orbits), `orbit_from_vector` and `projected_orbit_au`. These were the `HoldsOrbitPosition` mixin in `physics/position.py`, which is gone.
- **One mover.** `advance(seconds)` steps a body that has a state vector with `kepler.universal_step` and the pull of its primary (a planet about its star, a moon about its planet, a comet about its host), and advances the phase of a body on a circular galactic orbit (a star, a remnant, a wanderer) by its period. The batch updates in `store.py` do the same arithmetic over rows; `tests/test_state_kernel.py` holds the batch step to the object's.
- **Boss's 6:09 PM EDT card (the full S16 record for every mass object)** is the stored form of the same thing: mass with mu, position and velocity relative to the parent, spin axis and rate, epoch and next-due. Its columns replace the per-class ones in the DB.73 revision (Core 1); the per-class fields that become derived (speed, axial tilt, rotation period, galactic period and phase) are dropped there. Until then they remain the stored fields; nothing here aliases them.

## 12. The kinematic record (GEN.565c)

The record Boss's S16 card gives every mass object is declared in code as `KinematicRecord` (`generation/mass_object.py`) and read by `MassObject.kinematic_record(scale)`:

| Field | Read from | Unit |
|---|---|---|
| `position_au`, `velocity_kms` | the point's "system" frame (relative to the primary) | AU, km/s |
| `spin_axis`, `spin_rate_rad_s` | `spin_axis_x/y/z`; `2 pi / rotation_period_hours` (a non-rotating black hole is 0) | unit vector, rad/s |
| `mass_kg`, `mu_m3_s2` | the point | kg, m^3/s^2 |
| `epoch_unix`, `next_due_unix` | the point's epoch; epoch + `position.update_interval_s` at the speed | Unix s |

What S18 makes derived is readable from it today, and each equals the stored value to the precision of the computation (`tests/test_mass_object.py`): `speed_kms()` (the velocity's length, or the circular galactic speed), `axial_tilt_from_state()` (the angle between the spin axis and `r x v`, or the galaxy's pole when the body has no state vector; `spin.orbit_normal_from_state`) and `spin_rate_rad_s`. Core 1's DB.73 and DB.74 store this record and drop the per-class columns; this note is their target. Not here: the spin moving onto the point itself (GEN.584 after the point class GEN.578), the tidally locked rate, and a black hole's a* from angular momentum.

## 13. Built (DB.74, first slice): the registry of kind tables

`db/kind_tables.py` lists the tables that hold one kind of object once, with their groups (`uid`, `id_block`, `timestamped`, `sector_filed`, `system_filed`, `star_filed`, `phenomenon`, `containable`, `path`). The fifteen lists in `store.py`, `fingerprint.py`, `sector_paths.py`, `edits.py`, `api/ids.py`, `galaxy/objectref.py` and the slow-page list read it, with the same members and order as before. Left as data because each carries more than a list: `query._PHENOMENON_TABLES` (label, descriptor and radius expression per table) and `objectref.PHENOMENON_KINDS` (an order that sorts search results); `tests/test_kind_tables.py` checks both against the registry. The generic insert and load (next slices) and the core and profile tables (the 8.4 revision) read the same registry, so the move to three profile tables changes the rows of this one file, not its twelve readers.

## 14. Built (DB.74, second slice): one insert for every kind

`db/object_rows.py` writes a row from a list of `Field`s (`attr`, `optional`, `context`, `computed`); `store.py` holds one list per kind table (`STAR_FIELDS`, `PLANET_FIELDS`, `MOON_FIELDS`, `BELT_FIELDS`, `COMET_FIELDS`, `BLACK_HOLE_FIELDS`, `NEUTRON_STAR_FIELDS`, `ROGUE_PLANET_FIELDS`, `INTERSTELLAR_COMET_FIELDS`, `ASTEROID_FIELD_FIELDS`). The spin, the galactic orbit, the placement and the movement vector are shared fields, so the 49 columns three or more tables share are each spelled once in the code as well. The 8.4 revision changes the lists, not the machinery: a column that moves to the core table is marked core, one that moves to a profile is marked with its profile, and the insert writes the core row, then the profile row. Nebulae, supernova remnants and quasars (not mass objects) keep their own inserts.

## 15. Built (DB.74, third slice): one load for every kind

A `Field` also carries `load(row)`, the `from_dict` keys its stored columns fill, so `object_rows.load_dict(fields, row)` reads a row back through the list that wrote it (`scaled` divides kilometres back to AU and names what a NULL means, `each` and `attr` read the column as it is, `_flag` makes a `bool`). `_star_row_to_dict`, the black hole and neutron star loaders (which read only the satellite row's own columns, the rest being on the `stars` row), `_planet_or_moon_row_to_dict` (from `BODY_FIELDS` and, for a planet, `PLANET_ONLY_FIELDS`), `_belt_row_to_dict` and `_comet_row_to_dict` are built this way; the child rows (reflection spectrum, composition) stay beside them. The star loader now returns the stored `metallicity_fe_h`, which the old loader left out. Free objects (rogue planets, interstellar comets, asteroid fields) and the per-sector queries keep their own readers until the 8.4 revision moves the tables.
## 16. The StellarObject family class (GEN.565d)

`generation/stellar_object.py` holds `StellarObject(LazySaladName, MassObject)`, the family class the hybrid of section 10 names. What a main-sequence star, a giant, a white dwarf, a binary's proxy, a neutron star and a black hole share lives there: the persisted field list (`SERIALIZABLE_FIELDS`) and the class-level defaults a reloaded object needs, the save and load of `habitable_zone` and `lifespan`, `planet_age_requirement_gy`, `age_ceiling_gy`, `adjust_age_for_planets`, `calculate_system_perimeter` and `calculate_galactic_orbit`, and the two helpers `get_star_spectral_class` and `get_star_evolutionary_profile`. `Star` is the drawn star on top of it (its draw, its population model and `from_params`); `CompactRemnant` (and `BlackHole`, `NeutronStar`) and `BinaryStarProxy` subclass `StellarObject` directly, so the `_skip_property_init` workaround, which told `Star.__init__` not to draw, is gone. A remnant is no longer an instance of `Star`. No field, stored value or draw changes. Still to come: folding the remnant classes into states of the one class, once the stored shape of DB.74 is in.

`SmallBody` (`generation/small_body.py`, GEN.565e) is the second family class: `Comet`, `InterstellarComet` and `AsteroidField` subclass it and share the constructor preamble (the configuration, a name drawn on first read, whether the name was given) and the save and load of a composition kept as a list of strings (the field keeps its pairs). Nothing about a comet's orbit, a field's radius or the draws changes. `PlanetaryObject` (`generation/planetary_object.py`, GEN.565f) is the third: `Planet` (a planet or, with `is_moon`, a moon) and `RoguePlanet` subclass it and share `is_giant` and `parent_kind` (star, planet or none). They share almost no constructor code, so nothing else moved; stored names (`body_type` and `planet_type`, `mass` and `mass_kg`) stay until the DB.74 table revision.

## 17. Built (DB.75, first part): a body's ID follows its slot

The object ID (`galaxy/object_uid.py`) is the system's ID with a body number, and the number was given in insertion order, so any redraw of a system's content renumbered every body. A body now has a slot (`store._body_slots`): a star's is its role; a planet's or belt's is the role of the star it orbits (`system` for a circumbinary one) and its place in that star's list of orbits; a moon's is its planet's slot and its place among the planet's moons; a comet's is its star and its rank among the system's comets. `replace_star_system_content` records the IDs by slot before it deletes the old bodies, gives the new body of a slot the ID of the old one, and numbers the other new bodies from the system's counter, which is first brought past every old body so a number a body once held is never given to another. A slot is a kind and a place, so a planet redrawn into a belt's slot is a new body; and the identity does not claim that the redrawn body is the same world, only that the public address of "the third orbit of this star" keeps working.

What is not here: a body that changes parent. Today nothing moves a body between tables; the ejection of a planet (GEN.291) and the shared parent link of the 8.4 revision (DB.73) do, and in the revision the ID is a column of the one core row and the parent another, so a move keeps the ID by construction. The slot-path form of GEN.556 (`populate(host)` seeds each orbiter by `host_seed:orbit:slot`) is the same slot with the host's path in front; when it lands, `_body_slots` reads the path instead of the orbital index.
