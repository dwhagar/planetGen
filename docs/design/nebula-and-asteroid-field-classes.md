# Nebula, remnant and asteroid field classes (planned)

A draft for `docs/TODO.md` items 27-31, recorded 2026-09-30 from Boss's
notes and the reference document Boss shared ("Astrophysical
Architectures and Speculative Mechanics of Nebulae, Stellar Remnants, and
Interstellar Collisions"). The classes (items 28 and 31) are built in
schema v38: `program_constants.NEBULA_CLASSES` (A-W, with contents and
ranges) and `ASTEROID_FIELD_COMPOSITIONS` (letters plus a size digit).
Placement and central objects (27), containment (29) and naming (30) are
not built yet. The physical ranges come from that document and standard
references.

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

## What sits at the center (research for item 27)

| Family | Central object the generator must place |
|---|---|
| Diffuse H I | Nothing. |
| Emission (H II) | One or more O or early-B stars (roughly 15 Msun and up), young (under ~10 Myr): only they emit enough ultraviolet below 91.2 nm to ionize hydrogen. Large regions hold a young cluster. |
| Reflection | A B or A star near or inside the cloud, too cool to ionize it. |
| Planetary | Exactly one hot central star: the exposed core of a 0.8-8 Msun star after the asymptotic giant branch, 0.55-0.9 Msun, 30,000-200,000 K, becoming a white dwarf (`starData` Yerkes class VII/D). Often with a close companion, which shapes bipolar nebulae. Visible for ~10,000-25,000 years. |
| Molecular | No illuminating star. Dense cores and Bok globules may hold protostars or T Tauri (pre-main-sequence) stars; giant clouds may hold embedded young clusters. |
| Supernova remnant, core collapse | A neutron star (most) or black hole, as `SupernovaRemnant.compact_remnant` already does, offset from the center by its birth kick (a few hundred km/s times the remnant's age). |
| Supernova remnant, thermonuclear (Type Ia) | No compact object; sometimes a runaway surviving companion star. |

## Nebula and remnant classes (item 28, built in v38)

One letter per class, like `PLANET_CLASSES`. I and O are left unused so
they aren't mistaken for 1 and 0; X-Z are reserved. Today's
`NEBULA_TYPES` (emission, reflection, planetary, dark) and
`SUPERNOVA_REMNANT_MORPHOLOGIES` (shell, plerion, composite) map onto
these.

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

## Objects inside a nebula (item 29)

- Each object table that can sit inside one (star systems, stars in
  wide pairs, rogue planets, interstellar comets, black holes, neutron
  stars, asteroid fields, stand-alone facilities, other nebulae) gets a
  nullable `nebula_id` naming the innermost nebula or supernova remnant
  containing it. Asteroid fields can be inside a nebula but never
  contain anything.
- Nebulae are light-years across (emission up to 200 ly radius), far
  larger than a 13 ly sector, so containment is a 3D distance test
  against every nebula whose sphere reaches the object's sector, not a
  same-sector check. Computed at generation, when a later sector is
  generated inside an existing nebula, and on each correlative update.
- What being inside means, from the reference document: the stellar
  wind's heliopause shrinks as `R_HP ∝ (ρ_ISM v² + P_ISM)^-1/2`; in a
  dense cold cloud (nH ~3,000 cm^-3) the Sun's would sit at ~0.22 AU,
  exposing its planets to cosmic rays, mesospheric ozone loss, cooling
  and comet showers. Inside an emission nebula the sky stays dark and
  starry; inside a molecular cloud background stars vanish. This feeds
  system text, habitability, and the navigation hand-off radius (item
  33 uses a fixed ~120 AU heliopause today).

## Naming (item 30)

- Nebulae, supernova remnants, neutron stars, black holes, quasars and
  rogue planets are named the way star systems are: through the
  system-name registry (`_db.reserve_system_name`/confirm, uniqueness,
  diminutive decoration, offensive-word filter) instead of today's
  unregistered `generate_phoneme_salad_name`.
- Comets and asteroid fields get standardized designations. Draft,
  modeled on IAU prefixes:
  - Star-bound comet: `P/<system>-<n>` for periodic (period under 200
    years), `C/<system>-<n>` for long-period, numbered in order within
    the system. Example: `P/Veranthi-2`.
  - Interstellar comet: `I/<sector designation>-<n>`. Example:
    `I/4F2A1-3`.
  - Asteroid field: `AF <class><size>-<sector designation>-<n>`, where
    class is the letter from item 31 and size is floor(log10(radius in
    AU)) (a 0.001-1 ly field is 1-4). Example: `AF E3-4F2A1-02`.

## Asteroid field classes (item 31, built in v38)

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
