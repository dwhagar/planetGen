# Nebula, remnant and asteroid field classes

The design for `docs/TODO.md` items GEN.10 to GEN.14, recorded 2026-09-30 from Boss's
notes and the reference document Boss shared ("Astrophysical
Architectures and Speculative Mechanics of Nebulae, Stellar Remnants, and
Interstellar Collisions"). The physical ranges come from that document and
standard references.

**Status: built.**

| Item | What | Version | Code |
|---|---|---|---|
| 28, 31 | Letter classes for nebulae, remnants and asteroid fields | 7.19.0 (schema v38, PR #138) | `program_constants.NEBULA_CLASSES` (A-W), `ASTEROID_FIELD_COMPOSITIONS` |
| 29 | What sits inside a nebula | 7.24.0 (schema v39) | `inside_nebula_id` / `inside_remnant_id` columns |
| 27 | Nebulae generated with the stars they need | 7.30.0 (PR #144) | `generate.add_star_hosted_nebulae`, `_add_planetary_nebula`, `NEBULA_HOST_RULES`, `SUPERNOVA_KICK_SPEED_RANGE_KMS` |
| 30 | Standard names and designations | 7.31.0 (schema v40, PR #145) | `_db.reserve_system_name`, `cometData.comet_designation`, `roguePlanetData.interstellar_comet_designation`, `asteroidFieldData.asteroid_field_designation` |
| 29 | Heliopause squeezed inside a cloud | 7.39.0 (PR #161) | `starData.compressed_heliosphere_radius`, `queryDb.system_detail` |
| | Class shown on pages, search by class | 7.28.0, 7.29.0 | |
| | Nebulae and remnants drawn on the Galaxy Map | 7.36.0 | `queryDb.galaxy_clouds_in_box` |
| | Inside badges on system and phenomenon pages, "Inside <name>" in sector Contents, `inside` in `GET /api/systems/<id>` | 7.41.0 (PR #148) | `queryDb.system_detail` |
| | See-through cloud volumes on the Sector Map | 7.41.1 (PR #148) | `static/sectormap.js` |
| | Class reference pages for nebula, remnant and asteroid field classes; class labels link to them | 7.46.0 (PR #167) | `html/lib/classref.py`, `web/class_pages.py` |
| GEN.47 | Molecular clouds as a galaxy-scale field that spans sectors | | `nebulaField`, `_db._insert_field_nebulae` |

Habitability does not use the squeezed heliopause yet.

Boss's asks, in Boss's words:

- "Nebulae and Remnants should be generated and placed on the map.
  Research if we need stars at the center of these or what kind of star,
  etc, so we can make them."
- "We also want Nebulae and Stellar remnants to be classed like planets
  ... and we'll need to come up with What's IN the Nebulae ... develop a
  letter-class system similar to planets (A to Z) based on contents of
  the nebulae."
- "Add a DB field for if any stellar object (including systems) exist
  within a nebulae or similar (not asteroid fields, that wouldn't work)
  or stellar remnants if necessary."
- "Asteroid fields and comets should be named using a method that tells
  something about them by their name in letters and numbers in a
  standardized way. Nebulae should get names the same as star systems
  do, as do neutron stars, quasars, black holes, etc."
- "Asteroid fields should also have classes (A to Z) based on
  composition and density and size."

## What's in a nebula (by family)

From the reference document's table. nH is particle density; the gas to
dust mass ratio is about 100:1 throughout. Even a molecular cloud core
(10^6 cm^-3) is a harder vacuum than a laboratory vacuum (~10^7 cm^-3),
so a ship inside a nebula feels no drag and, inside an emission or
planetary nebula, sees an ordinary dark starry sky.

| Family | Dominant species | nH (cm^-3) | Temperature (K) | Optical extinction |
|---|---|---|---|---|
| Diffuse neutral (H I) | H, He, trace ions | 0.1-10 | 6,000-10,000 | transparent (A_V << 0.1) |
| Emission (H II) | H+, e-, O2+, N+, S+ | 10-10^4 | 8,000-12,000 | transparent; line emission ([O III] 500.7 nm, H-alpha 656.3 nm) |
| Reflection | silicate grains, PAHs, carbon soot, ices; neutral gas | ~10^2-10^3 | cold | partly; blue scattered light (∝ λ^-4) |
| Planetary | ionized H, enriched C, N, O, Ne | 10^2-10^5 | 10,000-20,000 | transparent; high surface brightness |
| Molecular (GMC, dark cloud, Bok globule) | H2, He, CO, PAHs, silicates, organics | 10^2-10^6 | 10-30 | opaque (A_V ~10-100 mag) |
| Supernova remnant (blast zone) | ionized ejecta, Fe, Si, S, ambient gas | 0.1-10^2 | 10^5-10^7 | transparent; synchrotron and soft X-ray |

## What sits at the center (research for GEN.10)

| Family | Central object the generator must place |
|---|---|
| Diffuse H I | Nothing. |
| Emission (H II) | One or more O or early-B stars (roughly 15 Msun and up), young (under ~10 Myr): only they emit enough ultraviolet below 91.2 nm to ionize hydrogen. Large regions hold a young cluster. |
| Reflection | A B or A star near or inside the cloud, too cool to ionize it. |
| Planetary | Exactly one hot central star: the exposed core of a 0.8-8 Msun star after the asymptotic giant branch, 0.55-0.9 Msun, 30,000-200,000 K, becoming a white dwarf (`starData` Yerkes class VII/D). Often with a close companion, which shapes bipolar nebulae. Visible for ~10,000-25,000 years. |
| Molecular | No illuminating star. Dense cores and Bok globules may hold protostars or T Tauri (pre-main-sequence) stars; giant clouds may hold embedded young clusters. |
| Supernova remnant, core collapse | A neutron star (most) or black hole, as `SupernovaRemnant.compact_remnant` already does, offset from the center by its birth kick (a few hundred km/s times the remnant's age). |
| Supernova remnant, thermonuclear (Type Ia) | No compact object; sometimes a runaway surviving companion star. |

## Nebula and remnant classes (GEN.11, built in v38, version 7.19.0)

One letter per class, like `PLANET_CLASSES`. I and O are left unused so
they aren't mistaken for 1 and 0; X-Z are reserved. Nebulae use A-Q
(`nebulae.nebula_class`) and remnants R-W (`supernova_remnants.remnant_class`).
`nebulae.nebula_type` is now the class's family (`diffuse`, `emission`,
`reflection`, `planetary`, `dark`; `diffuse` was added in v38), and
`SUPERNOVA_REMNANT_MORPHOLOGIES` (shell, plerion, composite) still sets a
remnant's morphology. The old `NEBULA_TYPES` constant is gone.

| Class | Name | Contents | nH (cm^-3) | Center |
|---|---|---|---|---|
| A | Diffuse neutral cloud | neutral H, He | 0.1-10 | none |
| B | Diffuse ionized gas | H+, e- | 0.1-1 | none nearby |
| C | Compact H II region | ionized gas still inside dusty birth cloud | 10^3-10^4 | one young O/B star |
| D | Classical H II region | H+, O2+, N+, S+ | 10-10^3 | small O/B cluster |
| E | Giant H II complex | as D, many shells, tens to hundreds of ly | 10-10^3 | young massive cluster |
| F | Reflection nebula | dust (silicates, PAHs, ices), neutral gas | 10^2-10^3 | B or A star |
| G | Emission-reflection mix | ionized core, dusty rim | 10^2-10^3 | early-B star |
| H | Young planetary nebula | ionized H, C, N, O, Ne; compact, bright | 10^4-10^5 | hot central star |
| J | Evolved planetary nebula | as H, expanded and faint | 10^2-10^3 | white dwarf |
| K | Carbon-rich planetary nebula | C-rich dust, PAHs | 10^2-10^4 | central star, 1.5-3 Msun progenitor |
| L | Nitrogen-rich bipolar planetary nebula | N and He enriched, bipolar lobes | 10^2-10^4 | central star, 3-8 Msun progenitor, often binary |
| M | Giant molecular cloud | H2, He, CO, PAHs, silicates | 10^2-10^6 | none; embedded clusters possible |
| N | Dark cloud / filament | H2, CO, cold dust | 10^3-10^5 | none |
| P | Bok globule | H2, CO, organics, dust; small, dense | 10^4-10^5 | none or one protostar |
| Q | Star-forming core | H2, outflows, Herbig-Haro jets, masers | 10^5-10^6 | embedded protostars |
| R | Young ejecta-dominated remnant | Fe, Si, S, O ejecta, 10^7 K | 0.1-10^2 | neutron star or black hole |
| S | Shell remnant | shocked gas and ejecta, 10^6-10^7 K | 0.1-10 | neutron star, black hole or none |
| T | Pulsar wind nebula (plerion) | relativistic electrons, magnetic field | low | pulsar (required) |
| U | Composite remnant | shell plus pulsar wind nebula | 0.1-10 | pulsar |
| V | Old radiative remnant | cooling shell merging with the ISM | 1-10^2 | neutron star, far off-center, or none |
| W | Thermonuclear remnant | Fe-rich, no hydrogen | 0.1-10 | none |
| X, Y, Z | Reserved | (hypernova remnant, tidal-disruption debris, fiction) | | |

Each class carries, like a planet class: description, composition,
radius range, density range, temperature range, extinction, central
object rule, and a relative frequency.

Boss (2026-09-30): "stellar remnants" means supernova remnants only; the
compact objects keep their own tables and need no letters.

## Objects inside a nebula (GEN.12)

- As built (schema v39), star systems, rogue planets, interstellar
  comets, black holes, neutron stars, asteroid fields and nebulae each have
  two nullable columns, `inside_nebula_id` and `inside_remnant_id`, naming
  the innermost nebula and the innermost supernova remnant containing them.
  The plan's single `nebula_id`, and its entries for stars in wide pairs and
  stand-alone facilities, were not built: a wide pair shares its system's
  row, and facilities have no such column. Asteroid fields can be inside a
  nebula but never contain anything.
- Nebulae are light-years across (emission up to 200 ly radius), far
  larger than a 13 ly sector, so containment is a 3D distance test
  against every nebula whose sphere reaches the object's sector, not a
  same-sector check. Computed at generation, when a nebula or remnant is
  placed (so a later sector inside an existing cloud sees it), by the v39
  migration for existing rows, and on each correlative update
  (`updateOrbits.py`, since 7.37.0).
- What being inside means, from the reference document: the stellar
  wind's heliopause shrinks as `R_HP ∝ (ρ_ISM v² + P_ISM)^-1/2`; in a
  dense cold cloud (nH ~3,000 cm^-3) the Sun's would sit at ~0.22 AU,
  exposing its planets to cosmic rays, mesospheric ozone loss, cooling
  and comet showers. Inside an emission nebula the sky stays dark and
  starry; inside a molecular cloud background stars vanish. This feeds
  system text, habitability, and the navigation hand-off radius. Built
  for the system text and navigation: `starData.compressed_heliosphere_radius`
  scales the stored open-space radius by `(P_ISM / P_cloud)^1/2`, with
  `P_cloud` the cloud's ram pressure at 26 km/s plus `n k T`;
  `queryDb.system_detail` returns it as `heliopause_au`. At nH 3,000
  cm^-3 that puts this program's Sun (~85 AU in open space) at ~0.6 AU.
  Habitability doesn't use it yet.

## The molecular cloud field (GEN.47)

Boss (2026-10-01): "No nebulae are being created at all." Molecular
clouds (classes M-Q) were rolled per sector at 5e-6 per pc^3 and placed
inside it, about 3e-4 per 13 ly sector, so a cloud tens of light-years
across showed up only if the one sector holding its center was
generated. Now they belong to the galaxy (`planetgen/galaxy/nebula_field.py`):

- The galaxy is cut into 50 pc cells (`NEBULA_FIELD_CELL_PC`). Each cell
  draws its clouds from its own seed (`galaxySeed.seeded`, kind
  `nebula-cell`, address `i/j/k`): a Poisson count with mean
  `PHENOMENON_DENSITY_PC3["molecular-cloud"]` (times its
  `PHENOMENON_RATE_SCALE`) times the cell's volume times the gas factor,
  each cloud's center uniform in the cell and its class, size and
  contents drawn as before. A galaxy planned without a seed uses a
  fixed one, so its sectors still agree.
- The gas factor is the young stars' density (they share the gas's thin
  layer and radial profile) with the gas's own arm contrast
  (`NEBULA_FIELD_ARM_AMPLITUDE` 0.6, about 4:1), relative to the average
  around the solar circle (2.82 disk scale lengths), to the power 1.4
  (`GMC_GAS_DENSITY_EXPONENT`), capped at 3
  (`NEBULA_FIELD_MAX_GAS_FACTOR`). At the solar circle that is about
  2.3 on an arm's crest, 0.9 midway and 0.2 between arms, and almost
  nothing 200 pc off the plane.
- A galaxy-placed sector takes every cloud whose sphere reaches it
  (center within its radius plus the sector's half diagonal, the same
  test as `_db.sectors_reached_by`), instead of its own
  `"molecular-cloud"` roll. `_db.insert_sector` stores each one under
  the neighbor lock unless a nebula is already stored at its center: the
  first sector saved that it reaches is its home sector, and containment
  is refreshed in every stored sector it reaches.
- Result at the solar circle: about 13% of arm sectors and 4% of
  between-arm sectors sit inside a dark cloud. Emission, reflection and
  planetary nebulae still come with their stars.

## Naming (GEN.13, built in v40, version 7.31.0)

- Nebulae, supernova remnants, neutron stars, black holes, quasars and
  rogue planets are named the way star systems are: through the
  system-name registry (`_db.reserve_system_name`/confirm, uniqueness,
  diminutive decoration, offensive-word filter) instead of today's
  unregistered `generate_phoneme_salad_name`.
- Comets and asteroid fields get standardized designations, modeled on
  IAU prefixes:
  - Star-bound comet: `P/<system>-<n>` for periodic (period under 200
    years), `C/<system>-<n>` otherwise, numbered in order around each
    star. A wide binary's comets carry their own star's name. Example:
    `P/Veranthi-2`. They follow a rename of their system or star.
  - Interstellar comet: `I/<sector designation>-<n>`. Example:
    `I/4F2A1-3`. A sector outside the grid is named instead.
  - Asteroid field: `AF <class><size>-<sector designation>-<n>`, where
    class is the letter from GEN.14 and size is floor(log10(radius in
    AU)) (a 0.001-1 ly field is 1-4). Example: `AF E3-4F2A1-02`.

## Asteroid field classes (GEN.14, built in v38, version 7.19.0)

Letter from composition and density (today's `sparse`/`typical`/`dense`);
size goes in the designation digit above. Composition families follow
asteroid taxonomy (C carbonaceous, S stony, M metallic, D/P icy
primitive, V basaltic).

| Class | Composition | Density |
|---|---|---|
| A / B / C | carbonaceous | sparse / typical / dense |
| D / E / F | stony (silicate) | sparse / typical / dense |
| G / H / J | metallic (iron-nickel) | sparse / typical / dense |
| K / L / M | icy, volatile-rich | sparse / typical / dense |
| N / P | basaltic (differentiated crust) | sparse or typical / dense |
| Q / R / S | mixed | sparse / typical / dense |
| T | dust-dominated, few large bodies | any |
| U | collisional family (fragments of one parent body) | any |
| V-Z | reserved | |

Boss (2026-09-30): size is "a digit and part of the class", so a full
class reads `C3` (`asteroid_fields.field_class`). The icy family's
components come from the cometary ices list.

## Why it works this way

- **Letter classes like planets.** Boss asked for nebulae and remnants to be
  "classed like planets ... A to Z based on contents", and asteroid fields
  "based on composition and density and size" (quoted above). The planet
  scheme was the model, including skipping letters that read as digits.
- **Classes carry contents, not just names.** Each class stores species,
  density, temperature, extinction and a center rule, because Boss asked
  what is *in* each nebula and because generation needs the center rule to
  place the right star (section "What sits at the center").
- **Supernova remnants only.** Boss said on 2026-09-30 that "stellar
  remnants" means supernova remnants; black holes and neutron stars keep
  their own tables.
- **Size as a digit in the asteroid field class.** Boss said size is "a digit
  and part of the class", so a class reads `C3`.
- **Stars placed by physics.** Only O and early-B stars ionize hydrogen, so
  H II regions are grown around them, and every planetary nebula needs a hot
  white dwarf at its center (`NEBULA_HOST_RULES`, CHANGELOG 7.30.0).
- **Containment by 3D distance.** Clouds are far larger than a 4 pc sector,
  so a same-sector test would miss most members.
- **Names through the system-name registry.** Boss asked that nebulae, neutron
  stars, quasars and black holes be named "the same as star systems do";
  sharing the registry also keeps every name in the galaxy unique.
- **Heliopause squeeze.** The reference document gives the pressure balance
  `R_HP ∝ (ρ v² + P)^-1/2`; the code applies it so navigation uses the real
  edge of a system inside a cloud.

### Alternatives not taken

- A single `nebula_id` column (built as two columns, one for nebulae and one
  for remnants; the reason is not recorded beyond the column comments).
- Generating diffuse gas (classes A-B): treated as background.
- Letters for black holes and neutron stars: not wanted (Boss, 2026-09-30).

## Planned changes (2026-10-03 and 2026-10-07)

- **Shapes**: a nebula gets an irregular, bulbous shape: 4 to 8 centres
  scattered in an anisotropic ellipsoid, evaluated as polynomial
  metaballs, with coordinates warped by domain-warped 3D simplex noise;
  marching cubes at an isovalue gives a mesh (a low-poly level for the
  Galaxy Map), and the same field gives an inside test. Boss also listed
  fBm, deformed spherical harmonics, diffusion-limited aggregation and
  level sets as alternatives.
- **Placed galaxy-wide first**: nebulae, black holes, neutron stars and
  quasars are placed across the galaxy when it is created, just before
  the bright-star scatter, and a sector fill keeps them and only adds.
- **Volume backfill**: a nebula that needs certain stars triggers a
  backfill of its volume down to 750 L_sun with the draw skewed to those
  types; existing backfilled stars inside are redrawn in place when a
  compatible type and brightness exist, and left alone when not.
- **Names**: every nebula gets a unique name from its ID (object-ids.md).
- **On the maps**: drawn from the mesh on every map, shaded over unfilled
  sectors too, with a show and hide toggle; the nebula page shows the
  whole shape with the dimmed galaxy around it.
- **Planets inside nebulae**: a feasibility study per class (disk
  photoevaporation near O and B stars, 26Al and 60Fe heating, extinction,
  cosmic rays) turns into generation rules; generation must know a
  system's surrounding cloud while it runs, not only on load
  (habitability-index.md, section 6).
- **Asteroid fields and belts as object systems**: still one object for
  positions and orbital speed, but rendered as many bodies seeded from
  the field's ID (a plan first).
