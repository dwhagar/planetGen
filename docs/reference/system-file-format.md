# System Specification Files

[`planetgen`](../../src/planetgen/cli/generate.py)'s `system` subcommand can load a full
(or partial) star system specification from a JSON file with
`--system-file` / `-f`:

```bash
planetgen system --system-file examples/systems/solar_system.json
```

Like any `planetgen system` run, this saves the system to the
configured database. Add `--output FILE` (or `--output -` for standard
output) to write the page to a file instead and leave the database alone.

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
| `num_orbits` | integer or `null` | The exact number of orbital slots (planets *and* asteroid belts, combined) to generate, at most 500. If omitted, the count is estimated from the star's mass (see `max_planets` below) and randomized. |
| `slots` | array or `null` | Per-orbit specifications — see [Slots](#slots) below. |
| `habitable_world` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `asteroid_belt` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `comets` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `large_star` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `moons` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `max_planets` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `intelligent_life` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `binary_system` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `wide_binary` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). Only meaningful together with `binary_system: true`. |
| `planets` | `true`, `false`, or omitted | See [Tri-state options](#tri-state-options). |
| `flavor_chance_system` | float `0.0`-`1.0` | Overrides the odds of a system-level "sensor readings" flavor line appearing. |
| `flavor_chance_planet` | float `0.0`-`1.0` | Overrides the odds of a planet-level flavor line appearing. |
| `max_planet_flavor` | `true` or omitted | Raises the total flavor-text budget for the run to 99 and sets `flavor_chance_planet` to 1. |

## Tri-state options

Most of the boolean-ish options above are **tri-state**: `true` forces the
feature to be present (except `habitable_world` and `intelligent_life`, where it
only raises the odds and never forces; the system may come out without and says so),
`false` forces it to be absent, and omitting the key
(or, in JSON, setting it to `null`) leaves it up to chance. This mirrors the
command line's `+name` / `-name` flags (`+habitable_world` sets it `true`,
`-habitable_world` sets it `false`).

| Key | `true` forces... | `false` forces... |
|---|---|---|
| `habitable_world` | Raises the odds of a habitable world; never forces one. | No habitable world is generated. |
| `asteroid_belt` | At least one asteroid belt exists. | No asteroid belt is generated. |
| `comets` | At least one star-bound comet exists. | No star-bound comet is generated. |
| `large_star` | A large, massive star is generated (only matters when `star_type` is omitted). | Nothing beyond the normal default population — included for symmetry/explicitness. |
| `moons` | Every planet gets a chance at moons. | No planet in the system gets moons (unless a slot gives it an explicit `moons` count — see below). |
| `max_planets` | The system generates the maximum number of orbital objects its star can support. | The system generates the minimum (0, or whatever `habitable_world`/`asteroid_belt`/`planets` require). |
| `intelligent_life` | Raises the odds of a planet with a technological civilization; never forces one. Also sets `habitable_world` to `true`. | No planet reaches a technological civilization. Also forces `habitable_world` to `true` (a world can still be habitable without intelligent life). |
| `binary_system` | A binary star system is generated instead of a single star (which of the two real configurations below is used is decided by `wide_binary`). | A single star (the default). |
| `wide_binary` | An S-type (wide) binary is generated: the two stars stay tens to thousands of AU apart, each keeping its own identity and independently-generated planets, with each star's maximum stable orbit limited by the companion's gravity. | A P-type (close/circumbinary) binary is generated instead: the two stars merge into one effective star for planet placement. If omitted, one of the two is chosen at random. Ignored entirely unless `binary_system` is also `true`. |
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
  | `planet_class` | no (`type: "planet"` only) | A specific planet class letter (see [Planet classes](#planet-classes)). If omitted, the planet is drawn normally and gets whatever class it comes out as. |
  | `moons` | no (`type: "planet"` only) | An exact number of moons to generate, `0` for none. If omitted, moon generation falls back to the `moons` option / random 50-50 chance. |

The array doesn't need to cover every orbit — it can be shorter than
`num_orbits` (or than however many orbits are otherwise generated), and any
entry can be `null`. Orbits past the end of the array, or with a `null`
entry, are generated normally.

Notes on how slots interact with the rest of generation:

- **Orbital distance isn't something you set directly.** Distances are
  computed sequentially, each orbit building outward from the last based on
  the star's mass and the previous object's size. A slot's `planet_class`
  never moves its orbit (see [Planet classes](#planet-classes)): the planet
  is drawn again at the orbit's own distance until it comes out as that class.
- **`moons` is best-effort.** A planet only has so much stable orbital room
  for satellites (governed by its Hill sphere). If you ask for more moons
  than fit, you'll get as many as fit and no more — no error.
- A slot with a `planet_class` counts toward a top-level `habitable_world: true`
  request only when the finished planet itself meets the request's conditions
  (a mild surface, air, and water at the surface); the class letter alone does not
  decide it.
- An unknown `planet_class` (a typo, or a letter no class uses) raises a clear
  error at generation time rather than silently substituting something else.
  A real class that the draws cannot reach at that place is listed as an
  unmet requirement in the generated system.

## Planet classes

These are the class letters valid for `slots[].planet_class`: every letter below (Q is the only reserved letter and is
not valid). A class is not drawn first. It is read off the finished body: the generator builds the whole planet (mass,
make, air, temperature, water), then gives it the class it fits best. Every class has its own conditions in
`src/planetgen/physics/data/planet_class_fit.json`; a condition just outside its edge lowers the match instead of ruling
the class out, and the class with the highest match wins (a tie goes to the class with more conditions, then to the
letter). No class relies on another having said no, so the order of the table does not matter. A slot that names a
class is generated again, up to 40 times, until the finished planet has that class.

A class that cannot come out at that place (a white dwarf's habitable zone holds no air) is listed in the system's
unmet requirements instead of being forced. An unknown letter raises an error at once. Only `"M"` is called out
separately in the summary output, as a "class M world". Whether life lives on a world is decided by its conditions and
history, never by its class letter.

The table gives the main conditions; the exact numbers and the other conditions (the temperature of the hottest and
coldest regions, wet bulb temperature, the phase of the liquid) are in `planet_class_fit.json`, and the class page of the
web interface lists them for each letter.

| Class | Name | Type | Main conditions |
|---|---|---|---|
| `T` | Gas dwarf | gas giant | Under 10 Earth masses. A rocky body with a hydrogen and helium envelope of 0.5 percent of its mass or more is judged by the giant classes. |
| `I` | Ice giant | gas giant | 10 to 35 Earth masses. |
| `J` | Gas giant | gas giant | 35 Earth masses or more. |
| `D` | Small airless world | terrestrial | No air (under 10 Pa), round, under 0.1 Earth masses. A body under a roundness ratio of 1 is an Asteroid and has no class. |
| `S` | Airless super-Earth | terrestrial | No air and 1.5 Earth masses or more. |
| `C` | Dead world | terrestrial | No air and 0.1 to 1.5 Earth masses. |
| `A` | Volcanic world | terrestrial | Sulfur dioxide 3 percent or more of the air. |
| `T` | Gas dwarf | terrestrial | Hydrogen 30 percent or more of the air and under 10 Earth masses. |
| `N` | Hot dense world | terrestrial | Air of 10 bar or more, a mean surface of 400 K or more, hydrogen under 30 percent. |
| `V` | High-gravity world | terrestrial | Air, 2 g or more, hydrogen under 30 percent, sulfur dioxide under 3 percent. |
| `L` | Marginal thin-air world | terrestrial | Air of 30 kPa to 3 bar, 273.15 to 320 K, carbon dioxide 1 kPa or more, oxygen 3 to 16 kPa, 0.5 to 2 g. |
| `K` | Thin-air world | terrestrial | Air under 3 kPa, 273.15 to 320 K, under 2 g, any oxygen. |
| `E` | Deep-frozen world | terrestrial | Air and a mean surface under 233.15 K (-40 C). |
| `P` | Glaciated world | terrestrial | The Earth-like world too cold: mean 240 to 270 K, the hottest region 270 K or more, poles under 245 K, oxygen 16 to 35.5 kPa, air 50 to 250 kPa, 0.75 to 2 g. |
| `W` | Frozen world | terrestrial | Air, a mean surface of 233.15 K to 273.15 K, under 2 g (class P wins where it fits). |
| `B` | Small hot world | terrestrial | Air, under 0.5 g, a mean surface above freezing. |
| `F` | Anoxic world | terrestrial | Oxygen under 3 kPa, carbon dioxide under 20 percent of the air, 273.15 to 400 K, 0.5 to 1.4 g. |
| `G` | Carbon dioxide world | terrestrial | Carbon dioxide 20 percent or more of the air, 273.15 to 400 K, 0.5 to 2 g. |
| `O` | Ocean world | terrestrial | M's too wet: water over 80 percent of the surface, mean 275 to 312 K, oxygen 16 to 35.5 kPa. |
| `H` | Desert world | terrestrial | M's too dry: under 20 percent water with groundwater, mean 275 to 315 K, oxygen 16 to 35.5 kPa. |
| `M` | Temperate world | terrestrial | Mean 275 to 305 K, liquid water over 20 to 80 percent of the surface, oxygen 16 to 35.5 kPa, carbon dioxide under 1 kPa, air 50 to 250 kPa, 0.75 to 2 g. |
| `R` | Hothouse world | terrestrial | M's too hot: mean 312 to 345 K, wet bulb 31 C or more. |
| `Z` | Hot thin-air world | terrestrial | 340 to 1500 K, air of 1 to 10 kPa, no ocean, 0.5 to 2 g. |
| `U` | Methane-ethane world | terrestrial | 91 to 130 K, air of 50 to 500 kPa with methane, methane-ethane liquid stable. |
| `Y` | Sulfuric acid world | terrestrial | 300 to 350 K, sulfur dioxide 3 percent or more of the air, sulfuric acid liquid stable. |
| `X` | Ammonia-water world | terrestrial | No air, a liquid cover over 20 percent, ammonia 10 percent or more of the water, 190 to 240 K. |

Moons are given a class by the same rules, so no letter is barred from a moon, and a slot's own `planet_class` always
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

See [`cli.md`](cli.md) for the full command-line reference,
or [`example-systems.md`](example-systems.md) for the rest of the example
files.
