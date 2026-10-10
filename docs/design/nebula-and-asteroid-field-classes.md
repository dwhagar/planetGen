# Nebula, remnant and asteroid field classes

The design for `docs/TODO.md` items GEN.10 to GEN.14, recorded 2026-09-30 from Boss's
notes and the reference document Boss shared ("Astrophysical
Architectures and Speculative Mechanics of Nebulae, Stellar Remnants, and
Interstellar Collisions"). The physical ranges come from that document and
standard references.

**Status: built (classes, naming, containment, molecular cloud field, shapes). Planned and
researched, not built: the planet rules, the required-star backfill and the asteroid field
render (sections "Corrected facts" to "Asteroid fields and belts as rendered object
systems", research of 2026-10-09; evidence tags [S] seen in a search result, [C] computed,
[R] recalled and unconfirmed).**

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
| | Class reference pages for nebula, remnant and asteroid field classes; class labels link to them | 7.46.0 (PR #167) | `planetgen/web/lib/classref.py`, `web/class_pages.py` |
| GEN.47 | Molecular clouds as a galaxy-scale field that spans sectors | | `nebulaField`, `_db._insert_field_nebulae` |

Habitability does not use the squeezed heliopause yet. Anomalies other than nebulae (magnetars,
Wolf-Rayet stars, jets and the rest) are in [anomalies.md](anomalies.md).

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

The radius and density of classes C to G are drawn independently of the stars
placed inside them. That is physically off (see "Corrected facts": an H II
region's radius follows from the ionizing photon rate of its stars) and is
the first thing GEN.99 should change.

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
  (`planetgen.cli.orbits`, since 7.37.0).
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
  around the solar circle (3.15 disk scale lengths), to the power 1.4
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
- Rate check (research 2026-10-09, recommendation): `PHENOMENON_DENSITY_PC3["molecular-cloud"]`
  is 5e-6 pc^-3, the upper end of Kennicutt and Evans 2012. With the four dark
  classes at equal weight and uniform radii the mean cloud volume is 3.0e4 pc^3
  (class M alone 1.2e5) [C], so gas factors of 0.2, 0.9, 2.3 and 3.0 put 3%, 13%, 35%
  and 46% of the volume inside a cloud [C]. Observed filling is about 0.5% in the
  inner Galaxy and about 1% for the thin disk [S/C], and the catalogues hold about
  8,100 to 9,700 clouds (Miville-Deschenes) or 1,064 massive ones (Rice), which is
  1.2e-8 to 1.4e-8 per pc^3 of disk, 400 times fewer than 5e-6 [C]. Part of that gap
  is units (the catalogues count only clouds above a few Msun/pc^2). The field is
  roughly 10 to 40 times too full. Recommended: about 1e-7 pc^-3 for class M, with
  small dark clouds and Bok globules (N, P, Q) as their own higher-density,
  small-radius rows, or a lower `NEBULA_FIELD_MAX_GAS_FACTOR` and class M radius
  range. The 5e-6 is Boss's chosen research value (2026-09-30), so this stays a
  recommendation until he decides.
  Correction (research 2026-10-10, [nebula-density-vs-reality.md](nebula-density-vs-reality.md)):
  the four classes are drawn 1 : 3 : 2 : 1, not at equal weight, so the volume inside a
  cloud is 2%, 8%, 19% and 24% at those gas factors, not 3%, 13%, 35% and 46%; the
  metaball shape fills only 12% of its bounding sphere; and the recommended rates are in
  that note.

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
  compatible type and brightness exist, and left alone when not. The
  stars each class needs and the algorithm are in "Required stars per class
  (GEN.99)" and "Nebula volume backfill (GEN.99)".
- **Names**: every nebula gets a unique name from its ID (object-ids.md).
- **On the maps**: drawn from the mesh on every map, shaded over unfilled
  sectors too, with a show and hide toggle; the nebula page shows the
  whole shape with the dimmed galaxy around it.
- **Planets inside nebulae**: a feasibility study per class (disk
  photoevaporation near O and B stars, 26Al and 60Fe heating, extinction,
  cosmic rays) turns into generation rules; generation must know a
  system's surrounding cloud while it runs, not only on load
  (habitability-index.md, section 6). Studied in "Planets in and around
  nebulae"; the rule table is "Rule table (GEN.94)" and the generator hook is
  "Applying the rules at generation (GEN.95)".
- **Asteroid fields and belts as object systems**: still one object for
  positions and orbital speed, but rendered as many bodies seeded from
  the field's ID (a plan first). The plan is "Asteroid fields and belts as
  rendered object systems (GEN.112)".

## A nebula's shape (GEN.75)

`planetgen.galaxy.nebula_shape` draws each nebula a shape of its own: 4 to
8 metaball centres in an anisotropic, turned ellipsoid, a domain warp from
three channels of gradient noise, and an isovalue; a point is inside where
the warped field reaches it. `NebulaShape.contains` is the containment test,
`mesh("low" | "full")` a marching-cubes triangle mesh. Lengths are in units
of the nebula's radius, scaled so the surface's farthest point is at 1.0, so
`radius_ly` stays a bounding sphere. A shape is plain JSON (`to_dict`), drawn
from a seeded `random.Random`, so equal seeds agree on every worker. Storage,
the mesh API and the containment switch are the next two parts of GEN.75.

## Decisions already taken

All from Boss (2026-10-03 05:38Z) unless dated otherwise; everything in the sections below this one is a recommendation.

- GEN.93: "Account for nebula temperature and conditions when generating planets and calculating surface conditions", with a feasibility study of whether a planet could form in each nebula class and what it would do to the planet.
- GEN.99: nebulae that need certain stars trigger a backfill of the nebula's volume "down to 750 Solar Luminosities" with probability skewed to the right types; existing backfilled stars are regenerated in place (same position, similar luminosity) when the required class is compatible, and left alone when not. The volume is the nebula's shape.
- GEN.112: asteroid fields and belts stay one object for position and orbital speed; only rendering generates a system of objects.
- GEN.110: two colliding rogue gas giants merge and a star that forms "should get no planets and instead should be a planetary nebula" (wording discussed in "Merger of two gas giants").
- NAV.51 (2026-10-09 01:02Z): "Courses should avoid asteroid fields."
- 2026-09-30: the supernova remnant is the only "stellar remnant" with letters; black holes and neutron stars keep their own tables (recorded above).

## Corrected facts

The research of 2026-10-09 checked the classes, the host rules and the star draws against standard astrophysics. Corrections, with what the code does now:

| # | Current statement or code | Problem | Evidence |
|---|---|---|---|
| 1 | `NEBULA_HOST_RULES`: a B0 to B2 star has a 50% chance of a nebula of class C, D or G; an O star 100% of C, D or E | One B0V star ionizes a Stromgren sphere of 3.5 ly radius at n = 100 cm^-3 (16 ly at n = 10); one B2V star 0.8 ly. Class D (10 to 100 ly radius, n 10 to 1,000) needs a total ionizing photon rate Q from 1e47 to 1e54 s^-1, median 3e50, about 50 O7 stars; class E (50 to 200 ly) needs 1e49 to 1e55, median 1e52. One O star (Q at most 4e49) cannot make E and rarely D. | [C] |
| 2 | Class radii drawn independently of the stars | An H II region's radius follows from Q and n (next section). | [C], textbook |
| 3 | `PLANETARY_NEBULA_CENTRAL_STAR_TYPES` = O3VII to B0VII built by the generic hot-white-dwarf branch of `generation/star.py` | That branch draws mass uniformly from `HOT_WHITE_DWARF_MIN_MASS_SOL` (1.1) to the Chandrasekhar limit (1.44) and luminosity 0.1 to 100 Lsun, with class O temperature capped at 60,000 K (`TEMP_RANGES`). Real central stars have a median mass of 0.59 Msun, most below 0.64 [S, arXiv 1606.05489, 2408.06448], luminosity 1e2 to 1e4 Lsun and temperature 30,000 to 200,000 K [R]. The class table above already says 0.55-0.9 Msun; the code does not follow it. | [S] [R] |
| 4 | `SEDOV_TAYLOR_RADIUS_COEFFICIENT_LY = 0.35` (ly per yr^0.4) | Equals a Sedov blast in an ambient density of about 220 cm^-3 [C]; the remnant classes list n from 0.1 to 100. The physical form is `R = 1.03 ly (E51/n)^0.2 (t/yr)^0.4`, which gives 5.0 pc at 1,000 yr for n = 1. | [C] |
| 5 | "Visible for about 10,000 to 25,000 years" (planetary nebulae) | Fine: a population simulation uses 0 to 30,000 yr and a minimum central-star mass near 0.53 Msun [S, arXiv 1910.08748]; at 25 km/s, 10,000 yr is a 0.83 ly radius [C], consistent with class H (0.1 to 0.5 ly) and J (0.5 to 3 ly). | [S] [C] |
| 6 | `surrounding_cloud` set only at load (`db/store.py`, `StarSystem.__init__` sets `None`) | Confirmed. The dict also lacks extinction, radius, age and the nebula's hot stars. | code |
| 7 | Backfill floor 750 Lsun against the pre-placed bright stars at 1,000 (`BRIGHT_STAR_MIN_LUMINOSITY_SOL`) | Not an error, but a nebula backfill at 750 adds stars the galaxy-wide tiers would not. Recorded as an open question. | code |
| 8 | `sample_bright_stars(..., population=None)` | A draw from the whole-disk population is 76% K/M giants at or above 750 Lsun; the "young" population is 84% B dwarfs and 3.7% O dwarfs (table in the required-star section). | [C] |
| 9 | "Cloud-shielded planets run cooler" and "H II gas heats planets" (the intuition behind GEN.93) | Not supported: see "Planets in and around nebulae". | [C] |

### H II radius from the ionizing photon rate (replaces independent radii for C to G)

Stromgren radius, alpha_B = 2.6e-13 cm^3 s^-1 at 1e4 K:

`R_S = (3 Q / (4 pi n^2 alpha_B))^(1/3) = 3.1 pc (Q / 1e49 s^-1)^(1/3) (n / 100 cm^-3)^(-2/3)` [C].

Q(H) per star (main sequence). Martins, Schaerer and Hillier 2005 (A&A 436, 1049) lowered the older Vacca et al. 1996 Lyman fluxes by 0.2 to 0.8 dex [S, arXiv astro-ph/0503346]; Sternberg, Hoffmann and Pauldrach 2003 cover Teff 25,000 to 55,000 K [S]. Neither table could be read, so the log Q values are recalled (+/- 0.15 dex for O, +/- 0.3 for B):

| Type | log Q [R] | Older calibration seen [S] | Blackbody upper bound [C] | ZAMS mass (Msun) [R] | L (Lsun) [R] | R_S at n = 100 (ly) [C] |
|---|---|---|---|---|---|---|
| O3V | 49.64 | 49.85 (Schaerer and de Koter 1997) | 49.62 | 60 | 8e5 | 17 |
| O5V | 49.22 | 49.48 | 49.29 | 40 | 3.5e5 | 12 |
| O7V | 48.8 | 49.06 | 48.90 | 25 | 1.2e5 | 8.8 |
| O9V | 48.2 | 48.46 | 48.44 | 18 | 5e4 | 5.6 |
| B0V | 47.6 | | 48.08 | 15 | 2e4 | 3.5 |
| B1V | 46.8 | | 47.29 | 11 | 8e3 | 1.9 |
| B2V | 45.7 | | 46.48 | 8 to 9 | 4e3 | 0.8 |
| B3V | about 45 | | 45.92 | 7 to 8 | 1.5e3 | |
| B5V | about 43.5 | | 44.81 | 5.5 | 7e2 | |

The blackbody column is an upper bound (real atmospheres absorb at the H and He edges) and overshoots by about 0.5 dex for B0 to B2, so the recalled B values are the ones to use. The generator only needs this: Q falls by about a factor of 10 per subclass from B0 to B2, and B2 and later cannot ionize more than a 1 ly region. (The research draft listed 1.5 ly for B1V; recomputing from log Q = 46.8 gives 1.9.)

Expansion. The ionization front keeps growing as a D-type front, `R(t) = R_S (1 + 7 c_s t / (4 R_S))^(4/7)` with c_s about 10 km/s [R, Spitzer]: R is about 16 pc after 3 Myr for R_S = 3.1 pc (Q = 1e49, n = 100) [C]. Density falls as `n_0 (R_S / R)^(3/2)` [R]. Rule: `radius = min(class cap, R_S * expansion(age))`, with radius and density drawn jointly from Q.

Q needed by class, at the geometric middle of the class radius and density ranges [C]:

| Class | log Q needed, low / median / high | Meaning |
|---|---|---|
| C | 46.4 / 48.9 / 51.4 | median: one O7 to O9 star; the top of the range is a cluster, not class C |
| D | 47.0 / 50.5 / 54.0 | median about 50 O7 equivalents (the "small O/B cluster") |
| E | 49.1 / 52.0 / 54.9 | 1,000 to 1e5 O7 equivalents; a cluster of 1e4 to 1e5 Msun |
| G | 46.9 / 49.4 / 51.9 | one O9 to O5 or an early-B pair |
| H | 47.0 / 49.0 / 51.1 | a central star with Q of 1e47 to 1e48 (T about 1e5 K, L 1e3 to 1e4 Lsun) ionizes it |
| J | 45.1 / 47.2 / 49.4 | easily |

Cluster mass to supply Q, Kroupa IMF, about 5.4e46 s^-1 per Msun of stars [C]:

| Cluster stellar mass (Msun) | Stars (down to 0.08 Msun) | Stars at or above 750 Lsun (6.6 Msun) | O stars (18 Msun and up) | log Q total |
|---|---|---|---|---|
| 100 | 173 | 1.4 | 0.4 | 48.7 |
| 300 | 518 | 4.2 | 1.1 | 49.2 |
| 1,000 | 1,726 | 14 | 3.6 | 49.7 |
| 3,000 | 5,177 | 42 | 11 | 50.2 |
| 10,000 | 17,257 | 141 | 36 | 50.7 |
| 100,000 | 172,572 | 1,411 | 359 | 51.7 |

Check: the Orion Nebula Cluster has about 3,500 stars within 2.5 pc [S, arXiv 1409.2503, astro-ph/0603138], near the 1,000 to 2,000 Msun row. Its ionization is dominated by one O7V star (Q about 1e49) where a fully sampled IMF would give 1e49.7, ordinary sampling scatter. The Stromgren radius of theta1 Ori C is 0.21 pc [S, arXiv 2603.04521], matching class C's lower end.

### Supernova remnant phases (classes R to W)

The Cioffi et al. 1988 transition formulas are as quoted in search snippets of arXiv 1607.04654 and 1701.05942 [S]; the Sedov numbers are [C]:

| Phase | Time | Radius | Notes |
|---|---|---|---|
| Free expansion | 0 to about 200-2,000 yr (until swept mass equals the 1 to 10 Msun ejecta) | v t at 5,000 to 10,000 km/s, 0.005 to 0.01 pc per yr | Cas A at 350 yr: about 2.5 pc [R] |
| Sedov-Taylor | to `t_PDS = 1.33e4 yr E51^(3/14) n^(-4/7)` | `R = 1.15 (E t^2 / rho)^(1/5)`; n = 1: 5.0 pc at 1,000 yr, 12.5 pc at 1e4 yr; n = 100: 2.0 pc at 1,000 yr | shell 1e6 to 1e7 K |
| Pressure-driven snowplow | from `t_PDS`, at `R_PDS = 14.0 E51^(2/7) n^(-3/7) pc` (14 pc at n = 1, 5.2 at n = 10, 1.9 at n = 100) | | shell cools to 1e4 to 1e6 K |
| Merging with the ISM | 1e5 to 1e6 yr | 20 to 100 pc | no longer a remnant; the 1e5 yr cap is right |

The two formulas agree: the Sedov radius at `t_PDS` for n = 1 is 14.0 pc [C]. Rule for the code: use `1.03 ly (E51/n)^0.2 (t/yr)^0.4`, capped where the snowplow phase begins. A core-collapse supernova leaves a neutron star (most) or a black hole (about 15% in the project); a Type Ia leaves nothing, sometimes a runaway companion. The compact object's offset is kick speed times age: 400 km/s for 1e5 yr is 41 pc [C], so it can leave a 20 to 40 pc remnant, as class V already says.

### Lifetimes by class

How long an object looks like its class, which sets the host age limit in the next section [R unless marked]: C 0.1 to 1 Myr; D 3 to 10 Myr (bounded by the O stars' lifetimes); E 5 to 20 Myr; F 1 to 30 Myr (up to 50 for a late host); G 1 to 5 Myr; H up to about 5,000 yr at 25 km/s [C]; J, K, L to about 30,000 yr [S, population simulation]; M 10 to 30 Myr, destroyed by its own stars; N a few to 10 Myr; P 1 to 3 Myr; Q 0.1 to 1 Myr; R to W up to 1e5 yr.

## Required stars per class (GEN.99)

750 Lsun is a B5/B6 main-sequence star of 6.6 Msun and about 90 Myr lifetime [C, L = M^3.5, t = 1e10 yr M^-2.5]. "Compatible" means the required star is at or above 750 Lsun, so the backfill can draw it.

| Class | Required star(s) | Minimum luminosity of a required star | 750 Lsun compatible? | Host age limit |
|---|---|---|---|---|
| A, B | none | | nothing to backfill | |
| C | one O or B0-B1 star with Q above about 1e47 | about 8e3 Lsun (B1V) to 5e5 | yes | under 1 to 3 Myr |
| D | Q total 1e48 to 1e51: at least one O star (30,000 Lsun and up) plus the IMF's B stars, about 10 to 30 O equivalents at median size | O star 3e4+; B filler 750+ | yes: the B filler is exactly the 750 to 3e4 band | under 3 to 8 Myr (lifetime of the lowest-mass O star in the set) |
| E | cluster of 1e4 to 1e5 Msun: 36 to 360 O stars, a few O3 to O5 | O 3e4+ | yes | under 10 to 20 Myr |
| F | one B2-or-later or A star, Teff under about 20,000 K (Q below about 1e47) | about 20 Lsun (A2V) to 7e3 | only partly: window 750 to about 7,000 Lsun (B5 to B2); A and late-B hosts fall under the floor | up to 50 Myr |
| G | one early-B star (B0-B2) or late O | about 3e3 to 3e4 | yes | under 5 Myr |
| H | central star 0.5 to 0.7 Msun, 2e3 to 1e4 Lsun in its first few thousand years | 1e3 | yes while young | under 5,000 yr |
| J | white dwarf core, 1e2 to 1e3 Lsun and fading | 100 | no for most of J's life (L falls under 750 within roughly 10,000 to 15,000 yr [R, post-AGB tracks]); the star comes from `PLANETARY_NEBULA_CENTRAL_STAR_TYPES`, not the backfill, so no conflict | 5,000 to 30,000 yr |
| K, L | as H or J (L: progenitor 3 to 8 Msun, brighter and faster cores) | 3e3 to 1e4 | as H, J | |
| M, N | none (M: optional embedded cluster) | | | |
| P | none or one protostar or T Tauri star | 0.1 to 50 Lsun | no, below the floor | under 3 Myr |
| Q | at least one protostar; massive ones 1e3 to 1e5 Lsun | any | no for the typical one | under 1 Myr |
| R to W | compact object from the remnant's own tables (S, V: or none; W: none; T needs a pulsar, up to 20,000 yr) | | not a star-backfill matter | |

For the remnant classes the "required star" is a compact object that `run_phenomenon.py` already builds, so GEN.99 does not apply; a pulsar's optical luminosity is 1e-3 Lsun or less, so the floor is irrelevant.

Check of the draw, measured with the project's `star_population` (3,000 stars each) [C]:

| Population | K/M giants and supergiants | B V | O V | Other |
|---|---|---|---|---|
| young (age 0 to 100 Myr) | about 5% | 84% | 3.7% | B giants, A/F subgiants |
| intermediate | about 92% | under 1% | 0 | |
| old | 100% | 0 | 0 | |
| whole disk | 76% | 10% | 0.5% | |

Fractions of living stars at or above each floor [C]: whole disk 3.7e-4 at 750 and 1.0e-4 at 1,000; young population 4.1e-3 at 750 and 3.2e-3 at 1,000. The curve is not monotonic across populations because old giants brighten near the tip of the giant branch.

### How many stars, and how dense

A 4 pc sector is 64 pc^3 and holds about 6 to 9 stars; stars at or above 750 Lsun are about 5e-5 per pc^3 in the field, 0.003 per sector (one per 330 sectors) [C].

Expected counts at or above 750 Lsun [C]: a class D nebula of 15 pc radius (14,000 pc^3, 220 sectors) holds 0.7 from the field against about 79 from its cluster at median Q (20 O stars, a cluster of about 5,600 Msun); class E of 40 pc radius (4,200 sectors) 13 against 36 to 360 O stars; classes F and G (2 and 8 sectors) 0.005 and 0.03 against 1 and 1 to 2.

Spread uniformly, a class D nebula adds 0.36 bright stars per sector (79 over 220 sectors), 120 times the field rate but placeable. The Hill-sphere separation rule makes a real cluster impossible: Hill radius is 0.5 pc (1.6 ly) for 1 Msun, 1.07 pc (3.5 ly) for 10 Msun and 1.5 pc (5 ly) for 30 Msun [C, `MILKY_WAY_MASS = 1.15e12 Msun`, `calculate_hill_sphere`], the minimum separation is the sum, and random sequential packing of hard spheres tops out near 38% volume fraction, about 0.1 B stars/pc^3 and 0.03 O stars/pc^3, close to the field density (0.14). The Orion Nebula Cluster has about 54 stars/pc^3 on average [S]. Concentrating D's stars into the central 3 pc would need 0.7 per pc^3. Options: (a) spread the O/B stars through the whole nebula volume (recommended for D and E); (b) pass a flat `min_separation_ly` (about 0.5 ly) to `SpaceSector.add_system` for nebula members, giving up the "Hill spheres never overlap" guarantee inside nebulae (real clusters have truncated tidal radii of 0.1 to 1 pc anyway [R]); (c) stop at the packing cap and log the shortfall. A real H II region also holds about 100 low-mass stars per O star [R]; the 750 floor leaves them out, so a nebula reads as a thin association. That is a realism gap, not a bug.

## Planets in and around nebulae (GEN.94)

### Result

Almost no nebula class makes planets. Planets are made in discs that take about 1 to 10 Myr and survive only in low ultraviolet. Every class with stars inside it (C, D, E, G, Q and the embedded parts of M) is younger than about 10 Myr, and the project already gives stars younger than `PLANET_MIN_STAR_AGE_GY = 0.01` only belts. So GEN.95 needs three things, not twelve: (a) a per-system far-ultraviolet flux G0 from the nebula's hot stars, which decides whether outer discs, gas giants and volatiles survive; (b) a composition and radiation modifier per remnant class (26Al and 60Fe enrichment, cosmic-ray tier); (c) flags for the odd cases (white-dwarf survivors in planetary nebulae, rare pulsar planets).

A cloud does not cool or dim a planet that has a star. At 1 AU inside a cloud of n = 1e4 cm^-3 the dust column to the star is about 1.5e17 cm^-2, A_V about 8e-5 [C]. H II gas at 1e4 K is no heater either: its energy flux on a surface is of order 1e-7 W/m^2 [C, order of magnitude]. What matters is the sky (background stars, glow), the squeezed heliopause (more low-energy cosmic rays, `compressed_heliosphere_radius`) and, for starless rogue planets, the cloud dust temperature as the temperature floor.

### External photoevaporation

The relevant field is far-ultraviolet flux in Habing units, 1 G0 = 1.6e-3 erg cm^-2 s^-1 [S, arXiv 1909.04093]. Winter and Haworth 2022 (arXiv 2206.11910) call external photoevaporation the dominant environmental influence on discs [S]; Winter et al. 2018 (MNRAS 478, 2700) find tidal encounters matter only above a stellar density of about 1e4 pc^-3 over 3 Myr and that photoevaporation always dominates [S]. Thresholds from the literature, not confirmed here: significant mass loss above about 1e3 G0, discs of a few hundred AU destroyed in about 1 Myr above about 1e4 G0 [R, Johnstone et al. 1998; Haworth et al. 2018, whose FRIED grid spans 10 to 1e4 G0 [S]]; fluxes of 100 G0 and more already shift planet mass and migration in one growth model [S, snippet]. Orion anchors [S, astro-ph/0007044]: proplyd HST 182-413 sits 0.12 pc from theta1 Ori C, loses 4.1e-7 Msun/yr from a disc of about 100 AU and 0.04 Msun, a lifetime near 1e5 yr; measured rates are 1e-7 to 1e-6 Msun/yr.

FUV flux from a hot star: `G0 = f_FUV L / (4 pi d^2) / 1.6e-3`, with f_FUV (the 6 to 13.6 eV share) 0.45 at 40,000 K, 0.54 at 30,000 K, 0.46 at 20,000 K, 0.16 at 12,000 K [C]. Unattenuated:

| Star (L, Teff) | G0 at 0.1 pc | 1 pc | 10 pc | d for 1e2 G0 | d for 1e3 G0 | d for 1e4 G0 |
|---|---|---|---|---|---|---|
| B5V (800, 16 kK) | 530 | 5.3 | 0.05 | 0.23 pc | 0.07 | 0.02 |
| B2V (5e3, 21 kK) | 4.8e3 | 48 | 0.5 | 0.7 | 0.22 | 0.07 |
| B0V (2e4, 30 kK) | 2.2e4 | 220 | 2.2 | 1.5 | 0.46 | 0.15 |
| O9V (1e5, 35 kK) | 1e5 | 1e3 | 10 | 3.2 | 1.0 | 0.32 |
| O7V-like (2e5, 40 kK; theta1 Ori C scale) | 1.8e5 | 1.8e3 | 18 | 4.3 | 1.35 | 0.43 |
| O3V (8e5, 45 kK) | 6.4e5 | 6.4e3 | 64 | 8.0 | 2.5 | 0.8 |

Gas extinction in the embedded phase lowers these by roughly e^(-3 A_V) in the FUV [R]; ignored in the first version, so the table is an upper bound. For a stored system, sum G0 over the nebula's host stars. In D and E the O stars are few and far apart, so the exposed volume is small: one O7 gives 1e3 G0 within 1.35 pc, 10 pc^3 of a 14,000 pc^3 class D nebula (0.1%, plus clustering). In an Orion-like core nearly everything is exposed. Adams 2010 (Annu. Rev. Astron. Astrophys. 48) describes the Sun's birth cluster the same way: N = 4,300 +/- 2,800 stars to supply the radioisotopes, but a larger cluster would photoevaporate the disc before giants form [S, snippet of arXiv 1001.5444].

### Supernova blast and enrichment

- Blast stripping is weak. Ouellette, Desch and Hester 2007 find discs survive a supernova even at 0.1 pc with under 1% mass loss, shocked gas cushioning the disc; gas injection is too weak for the solar radionuclides, dust-borne injection is near 100% [S, arXiv 0704.1652]. Close and Pittard 2017 find at 0.3 pc mass-loss rates of 1e-7 to 1e-6 Msun/yr for about 200 yr, and only an edge-on low-mass disc keeps 50% [S, arXiv 1704.06308]; Chevalier 2000 is cited for survival to 0.2 pc [S]. A supernova within 0.3 pc of a disc is also very improbable. Rule: no blast-kill; the blast is a composition event.
- Enrichment is the real effect. Solar 26Al/27Al at CAI formation was about 5.25e-5, half-life 0.73 Myr [S, PMC4946628]. For a chondritic composition (1 wt% Al, 3.12 MeV per decay [R inputs]), 5.25e-5 gives 5.9e6 J/kg, an adiabatic rise of 5,900 K (cp = 1,000 J/kg/K), so melting; 1e-5 gives 1,100 K (borderline); 1e-6 gives 110 K (nothing) [C]. 60Fe/56Fe = 1e-8 adds only 9 K, 1e-7 gives 93 K [C]: 26Al dominates early heating by about 100 times and 60Fe is a long-timescale minor term (half-life 2.6 Myr [R]).
- The galactic mean is lower. With 2 Msun of 26Al in the Galaxy [R], 5e9 Msun of gas [R] and Al/H = 3e-6 [R], the mean 26Al/27Al is about 6e-6, a tenth of the solar value [C]. The typical system is "26Al-poor"; the Solar System's rich value needs a nearby massive star or supernova, tied by Adams to birth clusters of 1,000 to 10,000 stars [S]. More 26Al means drier planets; bulk water and radius anticorrelate with 26Al [S, arXiv 1902.04026, A&A 2021].
- Present-day heating from 26Al and 60Fe is zero (26Al gone after 7 Myr, 60Fe after about 26 Myr). The effect is composition (dryness, iron-core size), not surface temperature.

### After the star's death

- White dwarfs: planets survive the transformation [S, arXiv 2010.09747]. WD 1856+534 b is a Jupiter-sized candidate (at most 14 Jupiter masses, 1.4 day period) around a white dwarf of about half a solar mass, 80 ly away; common-envelope migration from 1.69 to 2.35 AU is one explanation [S]. Polluted white-dwarf atmospheres (a quarter to a half [R]) show rocky survivors. Orbits widen as M_initial/M_final (about 1.8 for 1 to 0.55 Msun; the law is verified for adiabatic mass loss, Jeans 1924 and Veras et al. 2011, see [exotic-environments-planets-and-compact-binaries.md](exotic-environments-planets-and-compact-binaries.md) 1.2). `WD_PROGENITOR_ENGULFMENT_AU = 1.5` suits 1 to 2 Msun progenitors; for 3 to 8 Msun progenitors (class L) the AGB radius is 3 to 6 AU [R], so use `1.5 + 0.8 (M_prog - 1)` AU as a first formula (design choice).
- Inside a planetary nebula a survivor is briefly baked: `T = 278 K (1 - A)^(1/4) L^(1/4) / sqrt(d)` gives 840 K at 5 AU for L = 3,000 Lsun (A = 0.3) and 1,140 K for 10,000 Lsun [C], plus extreme ultraviolet, for about 1e4 yr. The nebula gas is irrelevant to a planet.
- In-situ planets: second-generation formation from AGB ejecta or common-envelope debris is theoretical [R, Perets 2010]. PSR B1257+12 has three planets (two Earth-mass, one tiny) at 0.2 to 0.5 AU, probably from a fallback or merger disc, but fallback discs are short-lived [S, arXiv 1609.06409, 0709.2922]. PSR B1620-26 b in M4 (about 2.6 Jupiter masses, circumbinary, about 12.7 Gyr) formed through a stellar encounter [S]. Pulsar planets are rare: about 0.1% of known pulsars and about 0.7% of millisecond pulsars ([exotic-environments-planets-and-compact-binaries.md](exotic-environments-planets-and-compact-binaries.md) 2.1).
- Free-floating planets in H II regions: JWST found Jupiter-mass binary objects in the Trapezium (about 40, 25 to 400 AU apart, 0.6 to 13 Jupiter masses; Pearson and McCaughrean 2023) but the status is disputed (Luhman's reanalysis) [S]. Rogue planets inside nebulae stay ordinary rogue planets.

### Surface conditions

- **Extinction.** `A_V = N_H / 1.9e21 cm^-2` [R]: n = 100 over 10 pc gives 1.6; n = 1e3 over 5 pc gives 8; n = 1e4 over 1 pc gives 16; n = 1e5 over 0.3 pc gives 49 [C]. The `NEBULA_CLASSES` A_V ranges are consistent. It dims the night sky by 10^(-0.4 A_V) and nothing else for a world with a star.
- **Cosmic rays.** The molecular-cloud ionization rate is zeta_H2 about 1e-17 s^-1 (core values 5e-18 to 3e-17), falling with column density and a further 3 to 4 times by magnetic mirroring in cores [S, Padovani et al. 2009; Padovani and Galli 2011]; the diffuse ISM is about 1e-16 [R]. Inside a dense core there are fewer cosmic rays than outside. A planet is exposed when its orbit lies beyond `compressed_heliosphere_radius` (0.22 AU at n = 3,000, 0.6 AU for this program's Sun). A thick atmosphere stops cosmic rays below about 1 GeV (about 1,000 g/cm^2); the cloud column (N_H = 1e22 cm^-2 is 0.023 g/cm^2) is negligible shielding [C]. Design: radiation tier +1 when the orbit is beyond the compressed heliopause and the cloud density exceeds 1e3 cm^-3, +2 if the atmosphere is also thin or absent.
- **Young remnants.** Under 1e4 yr they are bright in X-rays and accelerate cosmic rays: factor about 10 in classes R, T, U under 5,000 yr [R, design], 3 for S, 1 for V and W. The hot gas (1e6 to 1e7 K, n 0.1 to 10) has 100 to 1e4 times the ISM pressure, so the heliopause shrinks by 10 to 100 (the existing formula covers the pressure term).
- **Temperature.** A starless world in a dark cloud sits at the dust temperature, 10 to 30 K, not 2.7 K (floor `max(2.7, T_dust)`). A rogue planet near an O star reaches `T = 278 K (L/Lsun)^(1/4) / sqrt(d/AU)`: 35 K at 0.1 pc and 110 K at 0.01 pc for L = 1e5 Lsun [C]. Young-star X-ray and EUV activity depends on the star's age, not the nebula, and belongs to the habitability index.
- **Oceans.** Oceans survive any cloud for a world with a star. They do not survive the planetary-nebula phase at 1 to 5 AU (300 to 1,100 K [C]); volatile loss in the UV phase gives "no surface liquid water in classes H to L; volatile inventory reduced".

## Rule table (GEN.94)

Boss's research of 2026-10-09 18:23Z confirms this table and adds three refinements (host-mass scaling of the photoevaporation cut, a circumbinary exception for H to L, a swallowed-giant flag); see [exotic-environments-planets-and-compact-binaries.md](exotic-environments-planets-and-compact-binaries.md) section 1.

Columns: *In situ* can a planet-forming disc exist and finish now; *Pre-existing* can a planet formed elsewhere or earlier be here; *Temp* change to a star-hosted planet's equilibrium temperature; *Radiation tier* added to the habitability index's dose tier (0 field values; +1 up to 3 times cosmic rays or UV; +2 up to 10 times; +3 above 10 times or ionizing X-rays; design numbers, [R]-based); *G0 rule* the tier from the photoevaporation section using the nebula's host stars.

| Class | In situ | Pre-existing | Temp | Radiation tier | Atmosphere / composition |
|---|---|---|---|---|---|
| A, B | n/a (normal) | yes | none | 0 | none (background only) |
| C | no for stars inside (age under 1 to 3 Myr: belts only); G0 above 1e4 within 0.4 pc of an O star destroys discs | rare (a passing older star) | none | +1 (G0 tier) | none; the O star itself: no planets |
| D | rare; stars inside are under 8 Myr (belts only); giants suppressed within 1e3 G0 | yes for older stars that wandered in (about 5% of the field) | none | +0 to +2 by G0 | gas giants removed or stripped above 1e4 G0; volatile fraction reduced above 1e3 G0 (outer disc evaporated) |
| E | as D | as D | none | as D | as D; enrichment chance higher (supernovae begin after 3 Myr) |
| F | yes: a star older than 10 Myr passing through a dust cloud keeps its planets; B/A host discs vanish early (under 3 Myr [R]) so the host has a debris system at most | yes | none | +0 (heliopause tier if n above 1e3) | none |
| G | no for the host (early B); other stars as D | rare | none | +1 to +2 by G0 | as D |
| H | n/a | survivors beyond the AGB engulfment radius | `278 L^0.25 / sqrt(d)`, up to 840 K at 5 AU | +3 UV/EUV | volatiles lost; oceans gone; scorched for about 1e4 yr |
| J | n/a | survivors | same, L falling | +2 | as H, partly recovered |
| K | as J | as J | as J | +2 | carbon-rich dust: slightly more C in the late veneer [R, design] |
| L | as J; a binary common envelope can make a rare close planet (WD 1856 b type) | as J | as J | +2 | as J; engulfment radius 3 to 6 AU |
| M | rare: embedded clusters make belts at most | no | dust floor 10 K (starless), none (with star) | 0 to +1 (heliopause tier) | cold volatile veneer: comets and ices rich |
| N | no | no | as M | +1 if n above 1e3 and orbit beyond heliopause | none |
| P, Q | disc phase only (age under 3 Myr): give a protoplanetary disc, not planets | no | as M | 0 | none |
| R | no (rare fallback or merger planets: 1% flag) | survivors only beyond the blast zone; the blast does not strip them | none | +3 (young) | surface irradiated; fresh 26Al/60Fe only in systems under 10 Myr |
| S | no (1% flag) | yes | none | +2 | as R, weaker |
| T | rare: pulsar planets about 1% of millisecond pulsars [R] | no | pulsar luminosity heating | +3 | |
| U | as T | as T | | +2 | |
| V | no | yes | none | +1 | |
| W | no | yes | none | +1 | no neutron star; iron-rich dust |

Enrichment modifier, applied to any system whose birth environment included a massive cluster or supernova, or that sits inside D, E or R to W (design, [R] basis): draw `k = (26Al/27Al) / 5e-5` log-normal, field median 0.1, sigma 0.5 dex. For D and E hosts and young remnants, `k` has median 1 with probability 0.15 (D), 0.3 (E), 0.5 (R, S, U), else the field value. Water factor for planets formed inside the snowline `w = 1 / (1 + k^2 / 0.25)` (k = 0.1: 0.96; k = 1: 0.2; k = 3: 0.03) and iron-core fraction up 10 to 30% when k is above 1. The link between more 26Al and drier planets is [S]; the constants are design values.

Playable consequences: planets inside nebulae are nearly all absent or young. The honest remaining options are discs and protoplanets in C, P and Q ("planets forming" as a text flag), wanderers that passed through, and survivors in planetary nebulae and remnants. White dwarfs with outer survivors should be common; pulsar planets are a curiosity at about 1%. The levers by which a cloud can honestly matter for play are the night sky, the cosmic-ray tier and the ice supply (comets); navigation inside a cloud has no drag.

## Applying the rules at generation (GEN.95)

Today `StarSystem.__init__` sets `surrounding_cloud = None` and `_db.load_star_system` sets it at load from `inside_nebula_id` / `inside_remnant_id`. Generation builds the system first, then `sector.add_system` picks a position; pre-placed bright stars have a position before the system exists (`add_preplaced_system`), field systems do not.

1. **Extend** the dict from `_db.surrounding_cloud` and the nebula record: `extinction_av`, `radius_ly`, `age_myr`, `host_stars` (a list of position, luminosity, Teff for the stars that matter) and `center`. Add `age_myr` and `q_total` to `nebulae`, and `age_years` to the remnants table (the class already has `age_range_years`).
2. **Find containers early.** In `run_sector`, for the sector's box call `store._placed_containers(conn, low, high, reach_pc)` once (it returns center and radius), build the shapes (`NebulaShape.from_dict`), and test each new system's position with `shape.contains` for the innermost nebula.
3. **Order of operations.** `_random_position` needs the star's Hill radius, so either split `StarSystem` into star plus placement and then planets, or (recommended, cheaper) generate as now, choose the position, then call a new `StarSystem.apply_cloud(cloud)` that discards and regenerates planets only if the cloud's rules change anything (G0 above 1e2, age rule, remnant class) and stamps the radiation tier, enrichment factor and composition flags. The same function serves GEN.99's regenerate-in-place.

```python
def apply_cloud(system, cloud, hosts):
    # cloud: dict from surrounding_cloud(); hosts: [(distance_pc, L_sun, teff_k)]
    rules = NEBULA_PLANET_RULES[cloud["class"]]          # the rule table, as data
    g0 = sum(frac_fuv(T) * L * L_SUN / (4*pi*(d*PC)**2) / 1.6e-3 for d, L, T in hosts)
    tier = rules["radiation_tier"] + g0_tier(g0)          # 0 below 1e2, +1 to 1e3, +2 to 1e4, +3 above
    cr = cosmic_ray_tier(cloud, system)                  # heliopause squeeze, orbit beyond it
    k = draw_26al_ratio(cloud, rules)                     # log-normal, as above
    if rules["in_situ"] == "no" or system.star.age_gy < tuning.PLANET_MIN_STAR_AGE_GY:
        keep = belts_only(system.planets)
    else:
        keep = system.planets
    if g0 >= 1e4:   keep = remove_gas_giants_and_cut_outer(keep, outer_au=10)
    elif g0 >= 1e3: keep = cut_outer(keep, outer_au=50); scale_volatiles(keep, 0.5)
    elif g0 >= 1e2: keep = cut_outer(keep, outer_au=200)
    system.planets = apply_water_factor(keep, k)
    system.cloud_effects = dict(g0=g0, radiation_tier=max(tier, cr), k_26al=k,
                                cloud_id=cloud["id"], note=text_for_class(cloud))
```

Edge cases: a system in nested nebulae uses the innermost for `class` and sums hosts of all; a system in a remnant uses the remnant row. Hooks: the `StarSystem(...)` calls in `run_sector`, `SpaceSector.add_system` and `add_preplaced_system`, and `habitability` (radiation tier). Existing sectors keep working: tiers are derived at load from `surrounding_cloud` when absent.

A wider question sits behind this: 70 to 90% of stars form in embedded clusters and only a few percent of clusters stay bound [R, Lada and Lada 2003], so the same photoevaporation and enrichment lottery applies at birth to nearly every star, wherever it is now [S, Winter and Haworth]. Applying it to all stars would change all planets, not only those inside nebulae; that is a decision for Boss (handoff). A cheap middle path is one per-system `birth_g0` roll with no effect on planets unless enabled.

## Nebula volume backfill (GEN.99)

1. **Cluster mass per nebula.** From class: D `log10 Q` uniform 48 to 51 (the lower half is more faithful to "small O/B cluster"), E 51.5 to 53.5, G 48.5 to 50, C and F one star only. `M_cl = Q / 5.4e46` Msun.
2. **Counts.** Expected stars at or above 750 Lsun `N750 = M_cl * 0.0141`, at or above 18 Msun `N_O = M_cl * 0.0036`. Draw Poisson, then force `N_O >= N_O_min` (C 1, D 3, E 20, G 1) by upgrading the most massive drawn stars.
3. **Draw from the young population with an age cap.** Required stars must be younger than the nebula (`age_myr`, host age limit column above): draw from the IMF above 6.6 Msun with age uniform in [0, min(nebula age, lifetime)], rejecting any star whose lifetime is under its age. This is the "skewed" probability, produced by physics instead of a weight, and it reproduces an O:B ratio near 1:4 at the floor.
4. **Place** uniformly in the nebula shape (`NebulaShape.contains`, rejection sampling in the bounding sphere), respecting Hill separation by default and recording a shortfall when a star cannot be placed in `SECTOR_MAX_PLACEMENT_ATTEMPTS` (1,000) tries.
5. **Existing stars inside the shape.** For each placed star with luminosity `L0` and type `T0`, pick a target type from the required set. It is compatible if `L0` is within a factor of 2 (0.3 dex) of a luminosity the target can have at its age, i.e. mass in `[(L0/2)^(1/3.5), (2 L0)^(1/3.5)]` for a dwarf, and that mass's lifetime exceeds the nebula age. If compatible, regenerate the star (same ID and position, `Star.from_params`) and rebuild its planets with `apply_cloud`: a 5 Gyr K giant that becomes a 3 Myr B star cannot keep its planets. If not compatible (an 800 Lsun K giant against a B V requirement: giants only map to supergiants or giants, which an H II region does not accept), leave it alone, as Boss said. Log the counts.
6. **Pre-placed bright stars (1,000+).** If nebulae are scattered first (the plan), the bright-star scatter should draw a nebula's volume from the nebula's population in the first place; for galaxies already scattered apply step 5.
7. **The 750 to 1,000 band.** The galaxy scatter stops at 1,000; the nebula backfill goes to 750 inside the nebula only (handoff question 1).
8. **Order.** After the nebulae and bright-star scatter exist, before the sector fill, so a sector's fill sees the new stars as pre-placed (`fill_context`).

Reuse `sample_bright_stars` (it takes `population`, `min_luminosity_sol`, `max_luminosity_sol`): add a `"nebula"` population with an age window `(0, age)` and a mass-slope argument rather than a second sampler. `_bright_table` is cached by `(threshold, population)`; make the age window part of the key or the cache is wrong.

## Appearance by class (GEN.47, GEN.75, MAP.113, MAP.142)

GEN.47 and GEN.75 are done and nothing in the research argues against them. Optional class-aware parameters for `nebula_shape` and the edge look:

| Family | Real shape | `nebula_shape` parameters | Edge look (MAP.142) |
|---|---|---|---|
| H II (C, D, E) | roughly spherical with one flattened blister side toward the cloud | 4 to 6 centres, low anisotropy, cut one side | soft edge; glow brighter at the rim (H-alpha red, [O III] teal) |
| Reflection (F) | irregular, wispy, around the star | 6 to 8 centres, high warp | very soft; blue |
| Planetary (H to L) | ring, or bipolar lobes for L | two lobes for L, ring for the others | thin bright rim, sharp outer edge, then fade |
| Dark clouds (M, N) | filamentary, elongated 3:1 to 10:1 | high anisotropy, long axis | occlusion (dark), soft fading edge |
| Bok, cores (P, Q) | compact, nearly round | 3 to 4 centres | sharper edge |
| Remnants (R to W) | shell, limb-brightened | shell of thickness 0.2 R | bright thin rim fading inward and outward |

MAP.142 hint: alpha from the field value, `alpha = smoothstep(iso, iso + w, f)` with `w` by family (H II 0.4, reflection 0.6, dark 0.3, planetary shell 0.15, remnant shell 0.1 of the field span). Shell classes should use a radial profile, not a filled ball.

## Merger of two gas giants (GEN.110)

Boss's wording makes the merged star "a planetary nebula". Physically a planetary nebula is the ejected envelope of a dying 0.8 to 8 Msun star around a hot core, so a merger product is not one: it is an M dwarf or brown dwarf at birth, with at most a short-lived debris shell like a luminous red nova [R]. The project's gas giants reach 13 Jupiter masses, so two make at most 26, a brown dwarf; hydrogen fusion needs about 0.075 Msun (78.6 Jupiter masses) [R], so only two rogue brown dwarfs (13 to 80 Jupiter masses) can reach it. Such collisions essentially never happen: a physical-radius cross-section (2 R_J) at the project's 0.9 rogue planets per pc^3 gives about 2e-21 per object per year, 2e-11 per object in 10 Gyr, and about 1e-13 for giants [C]. The handler is for extremely rare events (dense nebulae or hand-made scenarios). Recommended: do not call the result a planetary nebula (the class rule needs a white dwarf and a bright hot core); use the reserved class X as "Merger remnant nebula" (faint shell, 0.01 to 0.3 ly, 1e3 to 1e5 yr) or a "merger glow" text on the star page; state the fusion test as a 0.075 Msun threshold, brown dwarf below, new M-dwarf star at or above.

## Asteroid fields and belts as rendered object systems (GEN.112)

### Real populations

| Population | Number and mass | Evidence |
|---|---|---|
| Main belt, over 1 km | 7e5 (SDSS) to 1.2e6 +/- 0.5e6 (ISO) | [S]; 1.1 to 1.9 million quoted elsewhere is not confirmed |
| Main belt mass | 13.8 +/- 2.0 e-10 solar masses = 2.7e21 kg; 2.4e20 kg without the 300 largest | [S] Vinogradova 2015 |
| Main belt geometry | 2.1 to 3.3 AU, about 1 AU thick, about 20 AU^3 | [C] with i about 10 degrees, e about 0.15 [R] |
| Mean spacing of bodies over 1 km | 0.0257 AU = 3.8e6 km | [C]; matches "several million km" [S] |
| Size distribution | Dohnanyi steady state: differential q = 3.5, cumulative -2.5 (real values 3 to 4) | [S] astro-ph/0308467, 1407.3307, 1111.0667 |
| Mass concentration | Ceres about 1/3 of the belt, top four about 60% | [R] |
| Kirkwood gaps | `a = a_J (q/p)^(2/3)`: 3:1 at 2.50 AU, 5:2 at 2.82, 7:3 at 2.96, 2:1 at 3.28 for a giant at 5.2 AU | [C] reproduces the real gaps |
| Kuiper belt | 1e5 objects over 100 km at 30 to 50 AU, 1e9 over 5 km; 0.02 to 0.2 Earth masses | [S] |

Three facts drive the design. With q = 3.5, mass per logarithmic size bin grows as D^0.5 and cross-section falls as D^-0.5, so the largest bodies hold the mass and the smallest the optical depth: among 1e4 to 1e5 sampled bodies of 1 to 100 km the largest carries 12 to 15% of the mass and the top ten 37 to 45% [C], which justifies hero bodies, a mid layer and a dust haze. A single power law down to 1 m is not real, and mass in 1 to 100 km from one q = 2.5 slope is only 2% of the belt [C], so the largest bodies are explicit heroes. Free fields are far emptier than the main belt: at a belt-calibrated 4e-16 bodies per kg, one Earth mass inside a 0.001 ly (63 AU) radius gives 4% of the belt's number density (spacing 1.1e7 km); at 0.1 ly the spacing is 7.6 AU and at 1 ly 76 AU [C]. The project's radius range of 0.001 to 1 ly runs from belt-like to empty.

### What the code has now

- `generation/belt.py` `AsteroidBelt`: `distance`, `lower_limit`, `upper_limit` (AU), `density` (sparse, typical, dense), composition list. Stored as `asteroid_belts` with km limits and a `uid`; no mass, seed or body count.
- `generation/phenomena/asteroid_field.py` `AsteroidField`: class letter plus size digit, `composition_family`, `density`, `radius_ly` (0.001 to 1 ly), composition, galactic orbit. `asteroid_fields` has `uid BINARY(12)`, centre in parsecs, nullable sector; no mass or seed. Isolated fields have rate 0 (`PHENOMENON_DENSITY_PC3["asteroid-field"]`) and exist only by hand or after GEN.110.
- Drawing: `web/maps/systemmap.py` `_belt_band` draws a ring from the belt's real limits; `systemview3d.js` draws each belt as 700 `THREE.Points` (`BELT_PARTICLES`), lift +/-4% of radius; `sectorscene.js` `makeAsteroidTexture` is a tan sprite; `web/maps/phenomenonrender.py` has `_NO_VIEW = {"asteroid_field"}`, so a field page has no 3D view.
- Bug to fix in passing: the 3D belt is seeded from the database `id` (`belt.id * 2654435761`), which changes when the database is rebuilt. `reproducible-galaxies.md` says identity comes from the seed and the object ID, so the render seed must come from the `uid`.
- `galaxy/seed.py` has `unit_seed(galaxy_seed, kind, address)` (SHA-256), `short_seed` and `seeded`; `nebula_field` uses a fixed `FALLBACK_SEED` when a galaxy has no seed, and fields and belts should do the same. three.js is revision 186, vendored.

### Plan

**1. Storage.** Nothing is required to render. Inputs that exist: identity (`uid`), extent (`lower_limit_km`/`upper_limit_km` or `radius_ly`), class (`density`, composition, `field_class`, `composition_family`), host (for a belt, the star and any giant outside it, for gaps). Recommended addition: a nullable `mass_kg` on `asteroid_fields` (later belts), written by GEN.110 when a planet is destroyed; derive from class when null. Placeholder defaults, not research: fields `M = 0.1 Earth mass * f_density` with f = 0.1, 1, 10 for sparse, typical, dense; belts `M = 5e-4 Earth mass * (width / 1.2 AU) * (mid-radius / 2.7 AU)^2 * f_density`.

**2. Seed.** `seed = unit_seed(galaxy_seed, "asteroid-render", f"{render_version}/{uid_hex}")`, 128 bits split into four 32-bit words for a JS `sfc32` (or the existing `mulberry32`). The server computes the words and puts them in the view's `data-view` JSON (as `phenomenonrender.js` reads it), so there is one implementation, not two. `render_version` is frozen: changing the algorithm changes every field and is a deliberate version bump. Without a galaxy seed use the fixed fallback. Math functions may differ in the last bit between engines; at Float32 precision that is invisible, and byte-exact identity needs integer-only noise.

**3. How many bodies.** The rendered count is a budget; the real count is 1e6 to 1e10.

| Layer | Count | Draw | Memory | Shown when |
|---|---|---|---|---|
| Impostor | 1 sprite | existing sprite and glow | negligible | far (default sector and galaxy views) |
| Dust | 5,000 to 20,000 | one `THREE.Points`, size attenuation, 2 to 4 px | 0.06 to 0.24 MB positions plus colour | camera within about 3 field radii |
| Rocks | 1,000 to 5,000 (1e5 ceiling) | `InstancedMesh`, 4 to 6 variants, per chunk | 76 B per instance: 0.08 to 0.38 MB, 7.6 MB at 1e5 [C] | within about 1 field radius, or 3 times a belt's width |
| Heroes | 12 to 48 | individual meshes, detail 3 (320 triangles), name labels | KBs | within 0.2 radius |

Icosahedron detail d has 20 (d+1)^2 triangles in r186 (80, 180, 320 for detail 1, 2, 3) [C]; 1e5 instances at 20 triangles are 2.0 M triangles [C]. GPU frame times were not measured: treat 5,000 rocks at 20 to 80 triangles as the safe default and 1e5 as a high-quality option to test on the integrated GPUs the site targets. Use per-chunk `InstancedMesh` objects (a 4x4x4 grid) so three.js culls whole chunks. Measured CPU cost of the generation pipeline (inverse-CDF size, Plummer position, Euler rotation, matrix compose; Node 22, one core) [C]:

| Bodies | Sampling | Matrices | Instance matrix array |
|---|---|---|---|
| 1e4 | 17 ms | 53 ms | 0.64 MB |
| 1e5 | 103 ms | 285 ms | 6.4 MB |
| 1e6 | 1.2 s | 1.2 s | 64 MB |

Even the ceiling builds in under half a second; do it once per view in a worker or idle time and cache by `uid` + `render_version`.

**4. Size law.** Truncated power law by inverse CDF: `D = D_min (1 - u (1 - (D_min/D_max)^b))^(-1/b)`, `b = q - 1`, default q = 3.5. Defaults by class (design, not research): dust (T) q = 4.0; collisional family (U) q = 3.0 to 3.5 with the largest body a parent fragment at 0.2 of D_max; icy and carbonaceous 3.5; metallic 3.8. Heroes are the top ranks above D_hero, explicit; rocks run from D_hero down to the size under about 1.5 px on screen; dust below. Mass per body `(pi/6) rho D^3` with rho by family (C 1,300 to 2,000, S 2,700, M 5,000 to 7,000, icy 900 to 1,500 kg/m^3 [R]); the real count `N(>1 km) = 4e-16 * mass_kg` [C] is shown on the page ("about 2e9 bodies larger than 1 km"), never rendered.

**5. Spacing and distribution.** No enforced spacing (independent draws give the Poisson clumps real fields have), except one cheap relaxation for visible rocks: reject a rock whose centre is inside 0.6 times the mean spacing at that radius of an earlier rock (grid hash, O(N)).
- *Belt:* semi-major axis density rising linearly with a between the limits (equal bodies per unit area), cut by Kirkwood-style gaps from the system's real outer giant: for p:q in 3:1, 5:2, 7:3, 2:1 remove a width of 0.02 a around `a_J (q/p)^(2/3)`. Eccentricity Rayleigh sigma 0.1, inclination Rayleigh sigma 6 degrees, longitudes and mean anomaly uniform. This replaces the flat disc with +/-4% lift.
- *Field:* a Plummer sphere, density `(1 + r^2/a^2)^(-5/2)`, a = R/3, cut at R, flattened by a seeded ellipsoid (axes 1 : 0.5-1 : 0.2-1) and a seeded rotation; class U adds a dense core at a random offset; dust is the same law with larger a.
- Rings around planets (belt law with the Roche limit as inner edge) are out of scope.

**6. Shape.** Triaxial ellipsoids (b/a 0.6 to 1, c/a 0.4 to 0.8 [R]) scaled per instance on a shared noise-displaced icosahedron; 4 to 6 variants at detail 1, 2, 3 (0.01 to 0.02 MB each, 3 to 16 ms to build [C]), two octaves of value noise at amplitudes 0.35 and 0.15 of radius. Albedo by family [R]: carbonaceous 0.05-0.08, stony 0.15-0.3, metallic 0.1-0.2 with sheen, icy 0.5-0.7 blue-white, basaltic 0.3-0.4, via `instanceColor`. The system's star is the one directional light; a starless field gets dim ambient plus a faint rim light.

**7. Level of detail.** With camera distance d to the centre and radius R: d > 8R sprite only; 3R to 8R sprite fades out, dust fades in; d < 3R dust and rocks, rocks split by chunk distance into three meshes (80, 180, 320 triangles); heroes under 0.5R; within 0.05R switch dust to perspective-sized sprites. A body is a mesh above about 3 px, a point below 1 px, dropped below 0.1 px. True sizes are far below that (1 km body: 1.1 px at 1,000 km, 0.01 px at 1e5 km; 100 km body: 1.1 px at 1e5 km [C], 45 degree field of view, 900 px), so sizes are exaggerated: a factor that puts the largest hero near 6 px at the default distance (typically 1e3 to 1e6), applied to all bodies and stated on the page ("bodies drawn x N larger than life"), the convention `phenomenonrender.js` already uses ("illustrative, not to scale").

**8. Motion: still one object.** Updates (`planetgen.cli.orbits`) touch only the belt's or field's row. A belt moves with its star, a field with its galactic orbit. Belt bodies: the vertex shader solves Kepler from per-instance elements (a, e, i, Omega, omega, M0: 6 floats) and a `uTime` uniform (the view has `update(years)`), mean motion `sqrt(mu/a^3)`; nothing is stored per body. Field bodies are static relative to the centre or spin rigidly and slowly from the seed (real debris shears out over 1e6 to 1e7 yr [S, Raymond et al. 2020]). Picking selects the whole belt or field; do not raycast instances.

**9. Collision risk and keep-out (NAV.24, NAV.51).** Expected hits on a path L through bodies with a ship of radius r_s: `lambda = integral n(D) pi (r_s + D/2)^2 L dD`, `P = 1 - exp(-lambda)`. Main belt (n = 5.9e4 per AU^3 over 1 km, spacing 3.8e6 km), 100 m ship [C]:

| Smallest body counted | Expected hits per 1 AU of path |
|---|---|
| 1 km | 1.1e-11 |
| 100 m | 1.0e-10 |
| 10 m | 9.8e-9 |
| 1 m | 2.7e-6 |

A 1-Earth-mass field in 63 AU is 4% of the belt's density. The real hazard is dust at speed (a 1 g grain at 0.1 c carries 4.5e11 J [C]), so any risk model should be flux times speed squared, not a body count. NAV.51 ("Courses should avoid asteroid fields", Boss) is therefore a game rule, not physics. Code facts: `galaxy/keepout.py` and `queryDb.keep_out_radius` return `KeepOut(None, "none", PASS_THROUGH_NOTE)` for `asteroid_field` and `NO_KEEP_OUT` for `belt` (NAV.24 shipped "asteroid field: no keep-out"). NAV.51 must return `KeepOut(radius_ly * LY_KM, "radius", None)` for fields (margin 1.1 times radius is a cheap default), update the `keepout.py` docstring and tests, keep belts at none (a belt is inside a system and the system's own keep-out applies), and make `nav_page` and the planner treat the basis like `"radius"`. A 1 ly sphere is a small detour at galaxy scale. GEN.110: a new field from a destroyed terrestrial gets `mass_kg` = the planet's mass, `radius_ly` from collision speed times age (1 km/s for 1e5 yr is 0.1 pc = 0.3 ly [C]) and a class U composition (q 3.0 to 3.5, one hero fragment at 0.2 of D_max); a second encounter adds mass and widens it.

**10. Tests.** Same seed gives the same bodies (hash of the first 1,000 sizes and positions), a different `uid` different bodies; the size sample matches the truncated Pareto (KS test) and `N(>1 km)` matches `4e-16 * mass`; no rock overlaps, Kirkwood gaps empty, chunk culling leaves the picture unchanged; a browser frame-budget test for 5,000 instances.

**11. PR order**, each independent of the next: (1) `uid`-based render seed and a shared `asteroidfield.js` generating dust, rocks and heroes, replacing `BELT_PARTICLES = 700`; (2) field page 3D view (remove `asteroid_field` from `_NO_VIEW`); (3) instanced rocks with chunk LOD; (4) GPU Kepler motion for belts; (5) `mass_kg` column and the GEN.110 hook; (6) NAV.51's keep-out change.

## Evidence notes

The research environment could read search-result text only, not the papers, so these [R] items need checking once paper access is allowed:

- log Q(H) for O3V to B2V (Martins 2005 and Sternberg 2003 tables not retrieved; B1 to B2 need Lanz and Hubeny or a B-star calibration since Sternberg stops at 25,000 K); luminosities and ZAMS masses of O/B types (Martins 2005, Pecaut and Mamajek), +/- 30%.
- Planetary-nebula central-star luminosity (1e2 to 1e4 Lsun), temperature to 2e5 K, and the time a 0.6 Msun core takes to fall under 750 Lsun (10,000 to 15,000 yr; Miller Bertolami 2016 tracks).
- FUV thresholds (1e3 G0 significant, 1e4 G0 destroys outer discs in about 1 Myr; check Haworth et al. 2018 and Winter and Haworth 2022 sections 3 to 4), disc lifetimes around B stars (under about 3 Myr), Herbig Be discs.
- Galactic 26Al mass (2 Msun), gas mass, Al/H, the mean 26Al/27Al near 6e-6, 60Fe/56Fe, decay energies (3.12 and 3.0 MeV). Seen: the 5.25e-5 ratio and 0.73 Myr half-life.
- Lada and Lada 2003 cluster fractions (70 to 90% born in clusters, about 10% surviving); cosmic-ray enhancement of young remnants (factor 10); the fraction of core-collapse supernovae in associations.
- Spitzer D-type expansion and c_s = 10 km/s; typical lifetimes of compact and classical H II regions, Bok globules and GMCs.
- Pulsar planet fraction (about 1%), PSR B1257+12 details, white-dwarf pollution fraction, second-generation planets (Perets 2010), the 0.075 Msun limit (78.6 Jupiter masses), the luminous-red-nova remark, `N_H / 1.9e21`, Kroupa slopes and the 100 low-mass stars per O star.
- Cioffi et al. 1988 constants were used as quoted in snippets (`t_PDS = 1.33e4 yr`, `R_PDS = 14.0 pc`) and agree with a Sedov calculation; a 3.61e4 yr figure in the snippet is a different quantity (shell formation time) and was not used.
- Asteroids: main-belt e and i, the mass shares of Ceres and the top four, the size-frequency bump near 100 km, albedos, axis ratios, bulk densities, Saturn ring thickness; the 1.1 to 1.9 million main-belt count and the "one in a billion" probe figure were not confirmed (found 0.7 to 1.2 million). Pinned by computation only: belt spacing, hit probabilities, triangle counts, matrix memory, JS timings. GPU frame times were not measured.

## Sources

Search results used (contents were snippets, not full papers):

- O-star ionizing fluxes: https://arxiv.org/abs/astro-ph/0503346, https://arxiv.org/pdf/astro-ph/0312232, https://arxiv.org/pdf/astro-ph/9611068
- Disc photoevaporation: https://ar5iv.arxiv.org/html/2206.11910, https://arxiv.org/pdf/astro-ph/0007044, https://arxiv.org/pdf/1307.4368, https://arxiv.org/pdf/1702.00251, https://ar5iv.labs.arxiv.org/html/1909.04093; Winter et al. 2018, arXiv 1804.00013; Adams 2010, https://ar5iv.arxiv.org/html/1001.5444
- Supernova effects and 26Al: https://ar5iv.arxiv.org/html/0704.1652, https://arxiv.org/pdf/1704.06308, https://arxiv.org/pdf/0808.2506, https://pmc.ncbi.nlm.nih.gov/articles/PMC4946628, https://pmc.ncbi.nlm.nih.gov/articles/PMC4950964, https://ar5iv.labs.arxiv.org/html/1902.04026, https://www.aanda.org/10.1051/0004-6361/202039334
- White dwarf and pulsar planets, free-floating planets: https://ar5iv.arxiv.org/html/2010.09747, https://arxiv.org/pdf/1609.06409, https://arxiv.org/pdf/0709.2922, https://en.wikipedia.org/wiki/PSR_B1620%E2%88%9226_b, https://en.wikipedia.org/wiki/Jupiter-mass_binary_object, https://astrobites.org/2024/03/21/a-jumbo-surprise-challenges-theories-of-planet-formation/
- Cosmic rays in cores: https://arxiv.org/pdf/1305.5393, https://ar5iv.arxiv.org/html/1104.5445, https://www.aanda.org/10.1051/0004-6361/200911794/pdf, https://arxiv.org/pdf/2109.08169, https://arxiv.org/pdf/2003.05416
- Planetary nebulae and clusters: https://arxiv.org/pdf/1910.08748, https://arxiv.org/pdf/1606.05489, https://arxiv.org/pdf/2408.06448, https://arxiv.org/pdf/1409.2503, https://arxiv.org/pdf/astro-ph/0603138, https://arxiv.org/pdf/2603.04521; remnants: https://arxiv.org/pdf/1607.04654, https://arxiv.org/pdf/1701.05942
- Cloud filling and catalogues: https://arxiv.org/pdf/1301.3905, https://arxiv.org/pdf/astro-ph/0701877, https://arxiv.org/abs/1602.02791, https://dc.g-vo.org/rr/q/lp/custom/CDS.VizieR/J/ApJ/834/57
- Main belt and small bodies: https://arxiv.org/pdf/0903.4479, https://iaaras.ru/en/library/paper/1046/, https://arxiv.org/pdf/astro-ph/0308467, https://arxiv.org/pdf/1407.3307, https://www.emergentmind.com/papers/1111.0667, https://arxiv.org/pdf/2502.01744, https://ipac.caltech.edu/publication/2011ApJ...741...68M, https://www.scientificamerican.com/article/in-science-fiction-movies, https://www.astronomy.com/science/how-do-spacecraft-avoid-collisions-in-the-asteroid-belt/, https://faculty.epss.ucla.edu/~jewitt/papers/ESO/ESO.pdf, https://arxiv.org/pdf/1810.09771, https://en.wikipedia.org/wiki/Jupiter_trojan, https://ar5iv.arxiv.org/html/1803.00072

Project files read: this document, [interstellar-object-rates.md](interstellar-object-rates.md), [habitability-index.md](habitability-index.md) (section 6), `tuning.py`, `generation/run_sector.py`, `phenomena/nebula.py`, `phenomena/asteroid_field.py`, `generation/star.py`, `star_population.py`, `bright_stars.py`, `belt.py`, `galaxy/sector.py`, `keepout.py`, `seed.py`, `nebula_field.py`, `system.py`, `db/store.py`, `db/query.py`, `physics/constants.py`, `systemview3d.js`, `sectorscene.js`, `phenomenonrender.js`.
