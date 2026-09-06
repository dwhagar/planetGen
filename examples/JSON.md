# System Specification Files

`../systemGen.py` can load a full (or partial) star system specification from a
JSON file with `--system-file` / `-f`:

```bash
python systemGen.py --system-file examples/solar_system.json
```

Every key is optional. Anything you leave out is generated normally (randomly,
subject to the usual physical rules). Anything given on the command line
*in addition to* `--system-file` overrides the corresponding value from the
file — the file sets a baseline, the command line wins.

See `solar_system.json` for a complete example that builds something
resembling our own solar system.

## Top-level keys

| Key | Type | Description |
|---|---|---|
| `star_type` | string or `null` | A specific spectral type + Yerkes class, e.g. `"G2V"`, `"M5V"`, `"B0IA"`. Format is `[OBAFGKM][0-9][Yerkes class]` (Yerkes classes: `0`, `IA+`, `IA`, `IAB`, `IB`, `II`, `III`, `IV`, `V`, `VI`, `VII`/`D`). If omitted, the star is generated randomly. |
| `name` | string or `null` | The star system's name. If omitted, a name is randomly generated. |
| `age` | `"young"`, `"old"`, or `null` | Biases the star's initial age toward the start or end of its lifespan. If omitted, age is random. |
| `num_orbits` | integer or `null` | The exact number of orbital slots (planets *and* asteroid belts, combined) to generate. If omitted, the count is estimated from the star's mass (see `max_planets` below) and randomized. |
| `slots` | array or `null` | Per-orbit specifications — see [Slots](#slots) below. |
| `habitable_world` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `asteroid_belt` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `large_star` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `moons` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `max_planets` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `intelligent_life` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `binary_system` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `planets` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `markdown` | `true`, `false`, or omitted | Output in Markdown instead of the default wikitext. |
| `flavor_chance_system` | float `0.0`-`1.0` | Overrides the odds of a system-level "sensor readings" flavor line appearing. |
| `flavor_chance_planet` | float `0.0`-`1.0` | Overrides the odds of a planet-level flavor line appearing. |
| `max_planet_flavor` | `true` or omitted | Raises the total flavor-text budget for the run to 99 and sets `flavor_chance_planet` to 1. |
| `output` | string or `null` | A file path to write the result to. Only used if `--output`/`-o` isn't also given on the command line. |

## Tri-state options

Most of the boolean-ish options above are **tri-state**: `true` forces the
feature to be present, `false` forces it to be absent, and omitting the key
(or, in JSON, setting it to `null`) leaves it up to chance. This mirrors the
command line's `+name` / `-name` flags (`+habitable_world` sets it `true`,
`-habitable_world` sets it `false`).

| Key | `true` forces... | `false` forces... |
|---|---|---|
| `habitable_world` | At least one habitable world exists. | No habitable world is generated. |
| `asteroid_belt` | At least one asteroid belt exists. | No asteroid belt is generated. |
| `large_star` | A large, massive star is generated (only matters when `star_type` is omitted). | Nothing beyond the normal default population — included for symmetry/explicitness. |
| `moons` | Every planet gets a chance at moons. | No planet in the system gets moons (unless a slot gives it an explicit `moons` count — see below). |
| `max_planets` | The system generates the maximum number of orbital objects its star can support. | The system generates the minimum (0, or whatever `habitable_world`/`asteroid_belt`/`planets` require). |
| `intelligent_life` | At least one planet reaches a technological civilization. Also forces `habitable_world` to `true`. | No planet reaches a technological civilization. Also forces `habitable_world` to `true` (a world can still be habitable without intelligent life). |
| `binary_system` | A binary (P-type, circumbinary) star system is generated instead of a single star. | A single star (the default). |
| `planets` | The system has at least one planet or asteroid belt (count won't be forced down to 0). | The system has no planets or asteroid belts at all — just the star. |

**Incompatible combinations** (the CLI will refuse these; take the same care
in a JSON file):
- `planets: false` with `moons: true`, `max_planets: true`, or `habitable_world: true`.
- `star_type` set together with `large_star: true`.
- `intelligent_life` (either value) together with `habitable_world: false`.
- `habitable_world: true` and `asteroid_belt: true` together with `large_star: false` (both objects need the extra room a large star provides).
- `num_orbits` together with `planets: false`.

## Slots

`slots` lets you dictate exactly what occupies specific orbital positions,
counted outward from the star starting at index 0. It's a JSON array; each
entry is either:

- `null` — generate that orbit normally (random, subject to the other options), or
- an object describing exactly what goes there:

  | Field | Required | Description |
  |---|---|---|
  | `type` | yes | `"planet"` or `"asteroid_belt"`. |
  | `planet_class` | no (`type: "planet"` only) | A specific planet class letter (see [Planet classes](#planet-classes)). If omitted, a class is chosen normally for whatever zone the orbit lands in. |
  | `moons` | no (`type: "planet"` only) | An exact number of moons to generate, `0` for none. If omitted, moon generation falls back to the `moons` option / random 50-50 chance. |

The array doesn't need to cover every orbit — it can be shorter than
`num_orbits` (or than however many orbits are otherwise generated), and any
entry can be `null`. Orbits past the end of the array, or with a `null`
entry, are generated normally.

Notes on how slots interact with the rest of generation:

- **Orbital distance isn't something you set directly.** Distances are
  computed sequentially, each orbit building outward from the last based on
  the star's mass and the previous object's size. If you request a
  `planet_class` that isn't valid at the distance an orbit would naturally
  land at (e.g. asking for `"M"`, an ecosphere-only class, in what would
  otherwise be a hot inner orbit), the generator nudges that orbit's distance
  into a zone that *does* support the requested class — the same way a
  forced habitable world snaps into the habitable zone. You don't need to
  pre-compute which index will land in which zone.
- **`moons` is best-effort.** A planet only has so much stable orbital room
  for satellites (governed by its Hill sphere). If you ask for more moons
  than fit, you'll get as many as fit and no more — no error.
- A slot with `planet_class` set to one of the habitable classes (see below)
  still counts toward satisfying a top-level `habitable_world: true`
  requirement.
- An invalid `planet_class` for *every* zone (a typo, or a class that exists
  in none of `h`/`e`/`c`) will raise a clear error at generation time rather
  than silently substituting something else.

## Planet classes

These are the class letters valid for `slots[].planet_class`. "Zone" says
where each class can naturally occur: **h**ot (inside the habitable zone's
inner edge), **e**cosphere (the habitable zone itself), **c**old (beyond the
habitable zone's outer edge). Classes marked "habitable" can host life and
count toward a `habitable_world: true` requirement; only `"M"` is
Earth-equivalent enough to be called out separately in the summary output as
a "class M world".

| Class | Description | Type | Zones | Habitable? |
|---|---|---|---|---|
| `A` | Small, barren, volcanic world | terrestrial | h | no |
| `B` | Small, molten world with a thin atmosphere | terrestrial | h | no |
| `C` | Dead world, no atmosphere | terrestrial | h, e, c | no |
| `D` | Small icy body | terrestrial | h, e, c | no |
| `E` | Molten core/crust world with a thin atmosphere | terrestrial | e | yes |
| `F` | Volcanic world with shallow seas and bacterial life | terrestrial | e | yes |
| `G` | Rocky, barren world with simple life | terrestrial | e | yes |
| `H` | Desert world, minimal water | terrestrial | e | yes |
| `I` | Ice giant with a tilted magnetic field | gas giant | c | no |
| `J` | Gas giant with a turbulent atmosphere and rings | gas giant | c | no |
| `K` | Adaptable world with a thin atmosphere | terrestrial | e | yes |
| `L` | Marginally habitable world with vegetation | terrestrial | e | yes |
| `M` | Terrestrial, Earth-like world | terrestrial | e | yes |
| `N` | Hot world with a dense, reducing atmosphere | terrestrial | e | no |
| `O` | Pelagic (ocean) world, >90% liquid water | terrestrial | e | yes |
| `P` | Cold, glaciated world | terrestrial | e | yes |
| `Q` | Eccentric orbit, extreme temperature swings | terrestrial | h, e, c | no |
| `R` | Ejected, geologically active world (rogue-like) | terrestrial | none¹ | no |
| `S` | Supergiant that shields the inner planets | gas giant | c | no |
| `T` | Gas dwarf with a thick atmosphere | gas giant | c | no |
| `U` | Ultragiant that could become a star | gas giant | c | no |
| `V` | Super-Earth with high gravity | terrestrial | e | yes |
| `W` | Tidally locked world, extreme temperature variation | terrestrial | h, e | yes |
| `X` | Stripped gas-giant core, no atmosphere | terrestrial | h | no |
| `Y` | "Demon" class world with a toxic atmosphere | terrestrial | h | no |

¹ `R` has no zone flags set (`h`/`e`/`c` are all `false`), so requesting it
explicitly in a slot will raise an error — it currently only appears from
free-form random generation via other code paths, not from a slot request.

Classes `Q`, `R`, `V`, `W`, `X`, and `Y` can't be generated as **moons**
(they're too large or physically implausible as satellites) — this only
matters for random moon generation, since a slot's own `planet_class` always
describes the planet itself, never its moons.

## Full example

```json
{
  "star_type": "G2V",
  "name": "Sol",
  "num_orbits": 9,
  "slots": [
    {"type": "planet", "planet_class": "B", "moons": 0},
    {"type": "planet", "planet_class": "N", "moons": 0},
    {"type": "planet", "planet_class": "M", "moons": 1},
    {"type": "planet", "planet_class": "C", "moons": 2},
    {"type": "asteroid_belt"},
    {"type": "planet", "planet_class": "J", "moons": 4},
    {"type": "planet", "planet_class": "J", "moons": 5},
    {"type": "planet", "planet_class": "I", "moons": 5},
    {"type": "planet", "planet_class": "I", "moons": 1}
  ]
}
```

This pins the star to a Sun-like `G2V`, lays out exactly 9 orbits, and
dictates each one: two small hot worlds, an Earth-like world with one moon, a
dead cold world with two moons, an asteroid belt, two gas giants, and two ice
giants — a rough analog of our own solar system. Anything not specified here
(the star's exact mass/radius/temperature within its type, each planet's
exact radius/mass/atmosphere details, orbital distances, names, ages, life
chemistry, flavor text, etc.) is still generated normally.

See [`EXAMPLES.md`](../README.md) for the full command-line reference.
