# planetGen

**Version:** 5.1.0 &middot; [Changelog](CHANGELOG.md) &middot; [Repository](https://github.com/dwhagar/planetGen) &middot; License: [CC0 1.0 Universal](LICENSE.md)

A procedural planet and star system generator, designed for the Molten Aether FFRP game. The output is designed to be easily copied and pasted into the wiki.

## Features

*   **Star generation**: Stars are generated with physically-grounded mass, radius, temperature, and luminosity, either fully at random (weighted by realistic galactic spectral-class prevalence) or pinned to a specific spectral type and Yerkes luminosity class (e.g. `G2V`) via `--star-type`. `+large_star` biases generation toward hotter, more massive stars.
*   **Binary star systems**: `+binary_system` generates a P-type (circumbinary) binary pair — a primary and secondary star orbiting each other — represented as a single unified effective star (combined mass, luminosity, and habitable zone) for the purposes of planet placement, while still reporting each star's individual properties in the output.
*   **Planet and moon generation**: Planets are drawn from 25 distinct planet classes (terrestrial and gas giant), each with its own composition, atmosphere, and valid orbital zones (hot/ecosphere/cold). Planets can generate their own moons, with orbital placement, atmospheric conditions, and surface gravity calculated per body. Orbital spacing is validated against each object's Hill sphere to keep the system physically plausible.
*   **Asteroid belts**: Belts can appear between planets (or be forced/forbidden with `+asteroid_belt`/`-asteroid_belt`), each with a randomly generated density and mineral/gem composition.
*   **Explicit orbit/slot specification**: The number of orbital slots (planets and asteroid belts combined) can be pinned exactly with `--num-orbits`, and a `--system-file` JSON specification can dictate the exact contents of any specific orbital slot — whether it's a planet or an asteroid belt, the planet's class, and how many moons it has — leaving unspecified slots to normal random generation.
*   **Stellar age and lifespan modeling**: Star age and lifespan are derived from spectral and Yerkes class data, then adjusted (via `--age young|old`) or extended as needed so any habitable planets' life stages remain consistent with how long the star has existed and how much longer it has left.
*   **Life chemistry and evolutionary timelines**: Habitable planets are evaluated against the star's spectral class to determine which of several life chemistries (e.g. Chlorophyll a, Melanin, Retinal) are viable, each with its own evolutionary pace (fast/normal/slow). A speculative evolutionary timeline (abiogenesis through technological civilization) is generated for habitable worlds, influenced by `+intelligent_life` / `-intelligent_life`.
*   **Flavor text**: Randomly-selected descriptive "sensor readings" flavor text can be appended to systems and planets, with limits and chances controllable via `--flavor-chance-system`, `--flavor-chance-planet`, and `--max-planet-flavor`.
*   **Dual output formatting**: Every generated system can be rendered as either MediaWiki wikitext templates (default, ready to paste into the wiki) or Markdown (`--markdown`).
*   **Sector generation**: `sectorGen.py` generates a whole sector of independently-random star systems in one pass, with an optional guaranteed minimum number of habitable systems (`--min-habitable`) — see [Sector Generation](#sector-generation) below.

## Setup

This project uses the `nltk` library to generate phonetically pleasing names. To install the necessary dependencies, you will need to have `setuptools` installed.

First, ensure you have `setuptools` installed:

```bash
pip install setuptools
```

Then, run the setup script to install the project dependencies:

```bash
python setup.py install
```

This will install the `nltk` library and download the 'words' corpus, which is required for the name generation.

## Usage

To generate a new star system, run the `systemGen.py` script from the root of the project:

```bash
python systemGen.py [options]
```

To generate a whole sector of star systems at once, see [Sector Generation](#sector-generation) below.

### Options

Most generation options use a `+name`/`-name` tri-state syntax: `+name` forces that feature to be present, `-name` forces it to be absent, and leaving it off the command line leaves it up to chance.

*   `--version`: Prints the program's version, repository URL, and license, then exits immediately (also available on `sectorGen.py`).
*   `+habitable_world` / `-habitable_world`: Force / forbid the generation of a habitable world in the system.
*   `+asteroid_belt` / `-asteroid_belt`: Force / forbid the generation of an asteroid belt in the system.
*   `+large_star` / `-large_star`: Force / forbid the generation of a large star.
*   `+moons` / `-moons`: Force / forbid moons on the system's planets.
*   `+max_planets` / `-max_planets`: Force the system to the maximum, or the minimum, number of orbital objects it can support.
*   `+intelligent_life` / `-intelligent_life`: Ensure a planet with intelligent life is (or is not) generated. Either implies `+habitable_world`.
*   `+binary_system` / `-binary_system`: Force / forbid a binary star system (P-type).
*   `+planets` / `-planets`: Ensure the system has at least one planet or asteroid belt, or none at all (star only).
*   `--system-file`, `-f <path>`: Load a system generation specification from a JSON file (see below). Any of the options above, or the value options below, given on the command line override the corresponding value from the file.
*   `--num-orbits <int>`: Force an exact number of orbital slots (planets and asteroid belts combined) to be generated.
*   `--output`, `-o <file_path>`: Specifies a file path to write the output to.
*   `--markdown`, `-m`: Formats the output in Markdown instead of the default wikitext.
*   `--star-type <type>`: Force the generation of a specific star type (e.g., G2V).
*   `--name <name>`: Specifies a name for the star system, overriding the default random generation.
*   `--age <young|old>`: Specifies the age of the star system (young or old).
*   `--flavor-chance-system <float>`: Overrides the default system-level flavor text chance (0.0 to 1.0).
*   `--flavor-chance-planet <float>`: Overrides the default planet-level flavor text chance (0.0 to 1.0).
*   `--max-planet-flavor`: Sets the maximum flavor text total for planets to 99.

**Note on Incompatible Options:**

*   `-planets` cannot be combined with `+moons`, `+max_planets`, or `+habitable_world`.
*   `--star-type` cannot be combined with `+large_star`.
*   `+intelligent_life`/`-intelligent_life` cannot be combined with `-habitable_world`.
*   `+habitable_world` and `+asteroid_belt` together cannot be combined with `-large_star` (both objects require the room a large star provides).
*   `--num-orbits` cannot be combined with `-planets`.
*   `--flavor-chance-system` must be a float between 0.0 and 1.0.
*   `--flavor-chance-planet` must be a float between 0.0 and 1.0.

### System specification files

`--system-file`/`-f` loads a JSON file describing exactly what to generate. Any key may be omitted, in which case that aspect is generated normally. The tri-state options above become `true`/`false`/omitted JSON keys, e.g.:

```json
{
  "star_type": "G2V",
  "name": "Sol",
  "age": "young",
  "num_orbits": 5,
  "habitable_world": true,
  "asteroid_belt": true,
  "slots": [
    {"type": "planet", "planet_class": "M", "moons": 1},
    {"type": "asteroid_belt"},
    null,
    {"type": "planet", "planet_class": "J", "moons": 4},
    null
  ]
}
```

`slots` is an optional, per-orbit list: each entry is either `null` (generate that slot normally) or an object with a required `"type"` (`"planet"` or `"asteroid_belt"`) and, for planets, an optional `"planet_class"` (e.g. `"M"`) and/or `"moons"` (an exact moon count, `0` for none). The list doesn't need to cover every orbit — slots past the end of the list are generated normally too.

### Sector Generation

`sectorGen.py` generates a whole sector of independently-random star systems in one pass, reusing `systemGen.py`'s own generation logic for each one:

```bash
python sectorGen.py [options]
```

Most of `systemGen.py`'s options work here too, but apply *uniformly* to every system in the sector — `+asteroid_belt` guarantees a belt in every system, `--star-type G2V` makes every star in the sector a G2V, and so on. `--system-file`, `--num-orbits`, and systemGen.py's own per-system `--name` aren't offered here, since those describe one specific, hand-crafted system rather than a sector of varied ones; use `systemGen.py --system-file` directly for that.

Sector-specific options:

*   `--num-systems <int>`: How many star systems the sector contains. Defaults to 10.
*   `--name <name>`, `-n <name>`: Hard-sets the sector's own name, overriding the default random two-word name (e.g. `"Voranthis Kelmoor"`) generated the same phoneme-salad way as star/planet/moon names.
*   `--min-habitable <int>`: Guarantees at least this many systems in the sector have a habitable world, chosen randomly among them — without forcing *every* system to have one the way a uniform `+habitable_world` would. Extra systems can still turn out habitable by chance on top of this minimum. Cannot exceed `--num-systems`, and cannot be combined with a uniform `-habitable_world`.
*   `--db-path <path>`: Where the generated sector is saved in the SQLite database. Defaults to `db/planetgen.db` at the project root.

The output opens with a sector-wide summary and an index of every system's name and star type, followed by each system's full write-up in turn. Every run also saves the whole generated sector — every system, star, planet, moon, and asteroid belt, plus a rendered copy of the wiki page in both wikitext and Markdown — to the SQLite database described in [`db/README.md`](db/README.md), regardless of whether `--output` was given.

### Additional Information

This tool is designed as a personal tool for the Molten Aether FFRP game. The output is designed to be simply cut and paste from the program output into the wiki. See https://wiki.moltenaether.com for wiki and game information.

Using commands to force a habitable world and an asteroid belt will automatically force a large star to ensure there is room for both objects. The most common stars are small dwarf stars which make a smaller star system. Forcing a large star as well as the maximum number of planets will cause the generated system to be very large with an extremely high number of planets. Do not assume that just because it is generated here, it is accurate or possible, such large systems may require editing as some worlds may end up saying they are several thousand AU's from the central star.

## License

This project is licensed under the [CC0 1.0 Universal](LICENSE.md) license.